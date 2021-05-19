#include <random>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include "matrix.hpp"
#include "array.hpp"
#include "ligand.hpp"
#include "string.hpp"

ligand::ligand(const path& p, array<double, 3>& origin, const pka& pka, double ph)
	: xs{}
	, num_active_torsions(0)
{
	// Initialize necessary variables for constructing a ligand.
	lines.reserve(200); // A ligand typically consists of <= 200 lines.
	frames.reserve(30); // A ligand typically consists of <= 30 frames.
	frames.push_back(frame(0, 0, 1, 0, 0, 0)); // ROOT is also treated as a frame. The parent and rotorX of ROOT frame are dummy.
	heavy_atoms.reserve(100); // A ligand typically consists of <= 100 heavy atoms.
	hydrogens.reserve(50); // A ligand typically consists of <= 50 hydrogens.

	// Initialize helper variables for parsing.
	vector<vector<size_t>> bonds; //!< Covalent bonds.
	bonds.reserve(100); // A ligand typically consists of <= 100 heavy atoms.
	size_t current = 0; // Index of current frame, initialized to ROOT frame.
	frame* f = &frames.front(); // Pointer to the current frame.
	f->rotorYidx = 0; // Assume the rotorY of ROOT frame is the first atom.
	string line;

	// Start parsing.
	// To prepare pdbqt for ligand, please use vega. With prepare_ligand4.py, some necessary polar hydrogens are missing.
	// With OpenBabel, hydrogens are all treated polar.
	//
	// [VEGA ZZ Usage]
	// vega ligand.pdb -o ligand.pdbqt -f VINA -c Gasteiger -p VINA -l GEN -r APOLAR -j FLEX -w
	//   -o ligand.pdbqt    write to file with name ligand.pdbqt.
	//   -f VINA            output AudoDock Vina pdbqt format.
	//   -c Gasteiger       add Gasteiger-Marsili empirical atomic partial charges.
	//   -p VINA            assign atom types using the AutoDock Vina force field (based on AMBER) template.
	//   -l GEN             add hydrogens with generic organic molecule.
	//   -r APOLAR          remove non-polar hydrogens.
	//   -j FLEX            output as a flexible molecule with branches.
	//   -w                 remove water.
	for (ifstream ifs(p); safe_getline(ifs, line);)
	{
		const string record = line.substr(0, 6);
		if (record == "ATOM  " || record == "HETATM")
		{
			// Whenever an ATOM/HETATM line shows up, the current frame must be the last one.
			assert(current == frames.size() - 1);
			assert(f == &frames.back());

			// This line will be dumped to the output ligand file.
			lines.push_back(line);

			// Parse the line.
			atom a(line);

			// Harmonize a unsupported atom type to carbon.
			if (a.ad_unsupported()) a.ad = 2;

			if (a.is_hydrogen()) // Current atom is a hydrogen.
			{
				// For a polar hydrogen, the bonded hetero atom must be a hydrogen bond donor.
				if (a.is_polar_hydrogen())
				{
					for (size_t i = heavy_atoms.size(); i > f->habegin;)
					{
						atom& b = heavy_atoms[--i];
						if (!b.is_hetero()) continue; // Only a hetero atom can be a hydrogen bond donor.
						if (a.is_neighbor(b))
						{
							b.donorize();
							break;
						}
					}
				}

				// Save the hydrogen.
				hydrogens.push_back(a);
			}
			else // Current atom is a heavy atom.
			{
				// Test if the current atom should be ionized.
				if (pka.is_ionized(a.name, trim(line.substr(17, 3)), line[21], ph))
				{
					a.donorize();
				}

				// Find bonds between the current atom and the other atoms of the same frame.
				assert(bonds.size() == heavy_atoms.size());
				bonds.push_back(vector<size_t>());
				bonds.back().reserve(4); // An atom typically consists of <= 4 bonds.
				for (size_t i = heavy_atoms.size(); i > f->habegin;)
				{
					atom& b = heavy_atoms[--i];
					if (a.is_neighbor(b))
					{
						bonds[heavy_atoms.size()].push_back(i);
						bonds[i].push_back(heavy_atoms.size());

						// If carbon atom b is bonded to hetero atom a, b is no longer a hydrophobic atom.
						if (a.is_hetero() && !b.is_hetero())
						{
							b.dehydrophobicize();
						}
						// If carbon atom a is bonded to hetero atom b, a is no longer a hydrophobic atom.
						else if (!a.is_hetero() && b.is_hetero())
						{
							a.dehydrophobicize();
						}
					}
				}

				// Set rotorYidx if the serial number of current atom is rotorYsrn.
				if (current && (a.serial == f->rotorYsrn)) // current > 0, i.e. BRANCH frame.
				{
					f->rotorYidx = heavy_atoms.size();
				}

				// Save the heavy atom.
				heavy_atoms.push_back(a);
			}
		}
		else if (record == "BRANCH")
		{
			// This line will be dumped to the output ligand file.
			lines.push_back(line);

			// Parse "BRANCH   X   Y". X and Y are right-justified and 4 characters wide.
			const size_t rotorXsrn = stoul(line.substr( 6, 4));
			const size_t rotorYsrn = stoul(line.substr(10, 4));

			// Find the corresponding heavy atom with x as its atom serial number in the current frame.
			for (size_t i = f->habegin; true; ++i)
			{
				if (heavy_atoms[i].serial == rotorXsrn)
				{
					// Insert a new frame whose parent is the current frame.
					frames.push_back(frame(current, rotorXsrn, rotorYsrn, i, heavy_atoms.size(), hydrogens.size()));
					break;
				}
			}

			// Now the current frame is the newly inserted BRANCH frame.
			current = frames.size() - 1;

			// Update the pointer to the current frame.
			f = &frames[current];

			// The ending index of atoms of previous frame is the starting index of atoms of current frame.
			frames[current - 1].haend = f->habegin;
			frames[current - 1].hyend = f->hybegin;
		}
		else if (record == "ENDBRA")
		{
			// A frame may be empty, e.g. "BRANCH   4   9" is immediately followed by "ENDBRANCH   4   9".
			// This emptiness is likely to be caused by invalid input structure, especially when all the atoms are located in the same plane.
			if (f->habegin == heavy_atoms.size())
			{
				frames.pop_back();
				lines.pop_back();
			}
			else
			{
				// This line will be dumped to the output ligand file.
				lines.push_back(line);

				// If the current frame consists of rotor Y and a few hydrogens only, e.g. -OH and -NH2,
				// the torsion of this frame will have no effect on scoring and is thus redundant.
				if (current + 1 == frames.size() && f->habegin + 1 == heavy_atoms.size())
				{
					f->active = false;
				}
				else
				{
					++num_active_torsions;
				}

				// Set up bonds between rotorX and rotorY.
				bonds[f->rotorYidx].push_back(f->rotorXidx);
				bonds[f->rotorXidx].push_back(f->rotorYidx);

				// Dehydrophobicize rotorX and rotorY if necessary.
				atom& rotorY = heavy_atoms[f->rotorYidx];
				atom& rotorX = heavy_atoms[f->rotorXidx];
				if (rotorY.is_hetero() && !rotorX.is_hetero()) rotorX.dehydrophobicize();
				if (rotorX.is_hetero() && !rotorY.is_hetero()) rotorY.dehydrophobicize();

				// Calculate parent_rotorY_to_current_rotorY and parent_rotorX_to_current_rotorY.
				const frame& p = frames[f->parent];
				f->parent_rotorY_to_current_rotorY = rotorY.coord - heavy_atoms[p.rotorYidx].coord;
				f->parent_rotorX_to_current_rotorY = normalize(rotorY.coord - rotorX.coord);
			}

			// Now the parent of the following frame is the parent of current frame.
			current = f->parent;

			// Update the pointer to the current frame.
			f = &frames[current];
		}
		else if (record == "ROOT" || record == "ENDROO" || record == "TORSDO")
		{
			// This line will be dumped to the output ligand file.
			lines.push_back(line);
		}
	}
	assert(current == 0); // current should remain its original value if "BRANCH" and "ENDBRANCH" properly match each other.
	assert(f == &frames.front()); // The frame pointer should remain its original value if "BRANCH" and "ENDBRANCH" properly match each other.

	// Determine num_heavy_atoms and num_hydrogens.
	num_heavy_atoms = heavy_atoms.size();
	num_hydrogens = hydrogens.size();
	frames.back().haend = num_heavy_atoms;
	frames.back().hyend = num_hydrogens;

	// Determine num_frames, num_torsions, flexibility_penalty_factor.
	num_frames = frames.size();
	assert(num_frames >= 1);
	num_torsions = num_frames - 1;
	assert(num_torsions + 1 == num_frames);
	assert(num_torsions >= num_active_torsions);
	assert(num_heavy_atoms + num_hydrogens + (num_torsions << 1) + (num_torsions ? 3 : 0) == lines.size()); // ATOM/HETATM lines + BRANCH/ENDBRANCH lines + ROOT/ENDROOT/TORSDOF lines == lines.size()
	flexibility_penalty_factor = 1 / (1 + 0.05846 * (num_active_torsions + 0.5 * (num_torsions - num_active_torsions)));
	assert(flexibility_penalty_factor <= 1);

	// Detect the presence of XScore atom types.
	for (const auto& a : heavy_atoms)
	{
		xs[a.xs] = true;
	}

	// Save the coordinate of the first heavy atom to origin.
	origin = heavy_atoms[frames.front().rotorYidx].coord;

	// Update heavy_atoms[].coord and hydrogens[].coord relative to frame origin.
	for (const auto& f : frames)
	{
		const array<double, 3> origin = heavy_atoms[f.rotorYidx].coord;
		for (size_t i = f.habegin; i < f.haend; ++i)
		{
			heavy_atoms[i].coord -= origin;
		}
		for (size_t i = f.hybegin; i < f.hyend; ++i)
		{
			hydrogens[i].coord -= origin;
		}
	}

	// Find intra-ligand interacting pairs that are not 1-4.
	interacting_pairs.reserve(num_heavy_atoms * num_heavy_atoms);
	vector<size_t> neighbors;
	neighbors.reserve(10); // An atom typically consists of <= 10 neighbors.
	for (size_t k1 = 0; k1 < num_frames; ++k1)
	{
		const frame& f1 = frames[k1];
		for (size_t i = f1.habegin; i < f1.haend; ++i)
		{
			// Find neighbor atoms within 3 consecutive covalent bonds.
			for (const auto b1 : bonds[i])
			{
				if (find(neighbors.begin(), neighbors.end(), b1) == neighbors.end())
				{
					neighbors.push_back(b1);
				}
				for (const auto b2 : bonds[b1])
				{
					if (find(neighbors.begin(), neighbors.end(), b2) == neighbors.end())
					{
						neighbors.push_back(b2);
					}
					for (const auto b3 : bonds[b2])
					{
						if (find(neighbors.begin(), neighbors.end(), b3) == neighbors.end())
						{
							neighbors.push_back(b3);
						}
					}
				}
			}

			// Determine if interacting pairs can be possibly formed.
			for (size_t k2 = k1 + 1; k2 < num_frames; ++k2)
			{
				const frame& f2 = frames[k2];
				const frame& f3 = frames[f2.parent];
				for (size_t j = f2.habegin; j < f2.haend; ++j)
				{
					if ((k1 == f2.parent && (i == f2.rotorXidx || j == f2.rotorYidx)) // The former frame is the parent of the later and either atom is on the connector bond between.
						|| (k1 > 0 && f1.parent == f2.parent && i == f1.rotorYidx && j == f2.rotorYidx) // Both atoms are on a connector bond to the same parent frame.
						|| (f2.parent > 0 && k1 == f3.parent && i == f3.rotorXidx && j == f2.rotorYidx) // The former frame is the grandparent of the later and both atoms are on their connector bond.
						|| find(neighbors.cbegin(), neighbors.cend(), j) != neighbors.cend()) // The later atom is within 1-4 neighbors of the former.
					{
						continue;
					}
					interacting_pairs.push_back(interacting_pair(i, j, mp(heavy_atoms[i].xs, heavy_atoms[j].xs)));
				}
			}

			// Clear the current neighbor set for the next atom.
			neighbors.clear();
		}
	}
}

// This function does not require receptor::use_maps.
result ligand::complete_result_noconf(const array<double, 3>& origin, const scoring_function& sf, const receptor& rec, vector<bool>& mask) const
{
	// Initialize origin coordinate, which is rotorY.
	vector<array<double, 3>> orig(num_frames);

	// Initialize heavy atom and hydrogen coordinates.
	vector<array<double, 3>> heavy_atoms(num_heavy_atoms);
	vector<array<double, 3>> hydrogens(num_hydrogens);

	// Calculate the coor of both heavy atoms and hydrogens of ROOT frame.
	const frame& root = frames.front();
	orig.front() = origin;

	for (size_t i = root.habegin; i < root.haend; ++i)
	{
		heavy_atoms[i] = orig.front() + this->heavy_atoms[i].coord;
	}
	for (size_t i = root.hybegin; i < root.hyend; ++i)
	{
		hydrogens[i] = orig.front() + this->hydrogens[i].coord;
	}

	// Calculate the coor of both heavy atoms and hydrogens of BRANCH frames.
	for (size_t k = 1; k < num_frames; ++k)
	{
		const frame& f = frames[k];

		// Update origin.
		orig[k] = orig[f.parent] + f.parent_rotorY_to_current_rotorY;

		// Update coor.
		for (size_t i = f.habegin; i < f.haend; ++i)
		{
			heavy_atoms[i] = orig[k] + this->heavy_atoms[i].coord;
		}
		for (size_t i = f.hybegin; i < f.hyend; ++i)
		{
			hydrogens[i] = orig[k] + this->hydrogens[i].coord;
		}
	}

	double e = 0;
	assert(mask.size() == rec.residues.size());
	vector<array<double, 6>> e_residues(rec.residues.size());
	vector<double> e_heavy_atoms(num_heavy_atoms);

	// Calculate inter-ligand free energy.
	for (size_t i = 0; i < num_heavy_atoms; ++i)
	{
		const size_t xs = this->heavy_atoms[i].xs;
		const auto& coord = heavy_atoms[i];

		for (const auto& a : rec.atoms)
		{
			const double r2 = norm_sqr(coord - a.coord);
			if (r2 < scoring_function::cutoff_sqr)
			{
				const size_t nsr2 = static_cast<size_t>(sf.ns * r2);
				const double e0 = sf.e[mp(xs, a.xs)][nsr2]; // Read from template.
				
				// Aggregate the energy.
				e += e0; // Total energy.
				scoring_function::score(e_residues[a.residue].data(), xs, a.xs, r2); // Per residue energy component.
				e_residues[a.residue].back() += e0; // Per residue energy total.
				e_heavy_atoms[i] += e0; // Per ligand atom energy.
				mask[a.residue] = true; // Mark contributing residue.
			}
		}
	}

	// Weighten the score components.
	for (auto& e_residue : e_residues)
	{
		for (size_t i = 0; i < scoring_function::weights.size(); i++)
			e_residue[i] *= scoring_function::weights[i];
	}

	// Save inter-molecular free energy into f.
	double f = e;

	// Calculate intra-ligand free energy.
	for (const auto &p : interacting_pairs)
	{
		const double r2 = norm_sqr(heavy_atoms[p.i1] - heavy_atoms[p.i0]);
		if (r2 < scoring_function::cutoff_sqr)
		{
			const size_t nsr2 = static_cast<size_t>(sf.ns * r2);
			e += sf.e[p.p_offset][nsr2];
		}
	}

	return result(e, f, false, move(heavy_atoms), move(hydrogens), move(e_heavy_atoms), move(e_residues));
}

bool ligand::evaluate(const conformation& conf, const scoring_function& sf, const receptor& rec, const double e_upper_bound, double& e, double& f, change& g) const
{
	assert(rec.use_maps);
	if (!rec.within(conf.position))
		return false;

	// Initialize frame-wide conformational variables.
	vector<array<double, 3>> orig(num_frames); //!< Origin coordinate, which is rotorY.
	vector<array<double, 3>> axes(num_frames); //!< Vector pointing from rotor Y to rotor X.
	vector<array<double, 4>> oriq(num_frames); //!< Orientation in the form of quaternion.
	vector<array<double, 9>> orim(num_frames); //!< Orientation in the form of 3x3 matrix.
	vector<array<double, 3>> forc(num_frames); //!< Aggregated derivatives of heavy atoms.
	vector<array<double, 3>> torq(num_frames); //!< Torque of the force.

	// Initialize atom-wide conformational variables.
	vector<array<double, 3>> coor(num_heavy_atoms); //!< Heavy atom coordinates.
	vector<array<double, 3>> deri(num_heavy_atoms); //!< Heavy atom derivatives.

	// Apply position and orientation to ROOT frame.
	const frame& root = frames.front();
	orig.front() = conf.position;
	oriq.front() = conf.orientation;
	orim.front() = qtn4_to_mat3(conf.orientation);
	for (size_t i = root.habegin; i < root.haend; ++i)
	{
		coor[i] = orig.front() + orim.front() * heavy_atoms[i].coord;
		if (!rec.within(coor[i]))
			return false;
	}

	// Apply torsions to BRANCH frames.
	for (size_t k = 1, t = 0; k < num_frames; ++k)
	{
		const frame& f = frames[k];

		// Update origin.
		orig[k] = orig[f.parent] + orim[f.parent] * f.parent_rotorY_to_current_rotorY;
		if (!rec.within(orig[k]))
			return false;

		// If the current BRANCH frame does not have an active torsion, skip it.
		if (!f.active)
		{
			assert(f.habegin + 1 == f.haend);
			assert(f.habegin == f.rotorYidx);
			coor[f.rotorYidx] = orig[k];
			continue;
		}

		// Update orientation.
		assert(normalized(f.parent_rotorX_to_current_rotorY));
		axes[k] = orim[f.parent] * f.parent_rotorX_to_current_rotorY;
		assert(normalized(axes[k]));
		oriq[k] = vec4_to_qtn4(axes[k], conf.torsions[t++]) * oriq[f.parent];
		assert(normalized(oriq[k]));
		orim[k] = qtn4_to_mat3(oriq[k]);

		// Update coordinates.
		for (size_t i = f.habegin; i < f.haend; ++i)
		{
			coor[i] = orig[k] + orim[k] * heavy_atoms[i].coord;
			if (!rec.within(coor[i]))
				return false;
		}
	}

	// Check steric clash between atoms of different frames except for (rotorX, rotorY) pair.
	//for (size_t k1 = num_frames - 1; k1 > 0; --k1)
	//{
	//	const frame& f1 = frames[k1];
	//	for (size_t i1 = f1.habegin; i1 < f1.haend; ++i1)
	//	{
	//		for (size_t k2 = 0; k2 < k1; ++k2)
	//		{
	//			const frame& f2 = frames[k2];
	//			for (size_t i2 = f2.habegin; i2 < f2.haend; ++i2)
	//			{
	//				if ((distance_sqr(coor[i1], coor[i2]) < sqr(heavy_atoms[i1].covalent_radius() + heavy_atoms[i2].covalent_radius())) && (!((k2 == f1.parent) && (i1 == f1.rotorYidx) && (i2 == f1.rotorXidx))))
	//					return false;
	//			}
	//		}
	//	}
	//}

	e = 0;
	for (size_t i = 0; i < num_heavy_atoms; ++i)
	{
		const size_t xs = heavy_atoms[i].xs;

		// Find the index and fraction of the current coor.
		const auto index = rec.index(coor[i]);

		// Assert the validity of index.
		assert(index[0] + 1 < rec.num_probes[0]);
		assert(index[1] + 1 < rec.num_probes[1]);
		assert(index[2] + 1 < rec.num_probes[2]);

		// (x0, y0, z0) is the beginning corner of the partition.
		const size_t x0 = index[0];
		const size_t y0 = index[1];
		const size_t z0 = index[2];
		const double e000 = rec.e(xs, array<size_t, 3>{{x0    , y0    , z0    }});

		// The derivative of probe atoms can be precalculated at the cost of massive memory storage.
		const double e100 = rec.e(xs, array<size_t, 3>{{x0 + 1, y0    , z0    }});
		const double e010 = rec.e(xs, array<size_t, 3>{{x0    , y0 + 1, z0    }});
		const double e001 = rec.e(xs, array<size_t, 3>{{x0    , y0    , z0 + 1}});
		deri[i][0] = (e100 - e000) * rec.granularity_inverse;
		deri[i][1] = (e010 - e000) * rec.granularity_inverse;
		deri[i][2] = (e001 - e000) * rec.granularity_inverse;

		e += e000; // Aggregate the energy.
	}

	// Save inter-molecular free energy into f.
	f = e;

	// Calculate intra-ligand free energy.
	for (const auto &p : interacting_pairs)
	{
		const array<double, 3> r = coor[p.i1] - coor[p.i0];
		const double r2 = norm_sqr(r);
		if (r2 < scoring_function::cutoff_sqr)
		{
			const size_t nsr2 = static_cast<size_t>(sf.ns * r2);
			e += sf.e[p.p_offset][nsr2];
			const array<double, 3> d = sf.d[p.p_offset][nsr2] * r;
			deri[p.i0] -= d;
			deri[p.i1] += d;
		}
	}

	// If the free energy is no better than the upper bound, refuse this conformation.
	if (e >= e_upper_bound) return false;

	// Calculate and aggregate the force and torque of BRANCH frames to their parent frame.
	for (size_t k = num_frames - 1, t = 6 + num_active_torsions; k > 0; --k)
	{
		const frame& f = frames[k];

		for (size_t i = f.habegin; i < f.haend; ++i)
		{
			// The deri with respect to the position, orientation, and torsions
			// would be the negative total force acting on the ligand,
			// the negative total torque, and the negative torque projections, respectively,
			// where the projections refer to the torque applied to the branch moved by the torsion,
			// projected on its rotation axis.
			forc[k] += deri[i];
			torq[k] += cross(coor[i] - orig[k], deri[i]);
		}

		// Aggregate the force and torque of current frame to its parent frame.
		forc[f.parent] += forc[k];
		torq[f.parent] += torq[k] + cross(orig[k] - orig[f.parent], forc[k]);

		// If the current BRANCH frame does not have an active torsion, skip it.
		if (!f.active) continue;

		// Save the torsion.
		g[--t] = torq[k] * axes[k]; // dot product
	}

	// Calculate and aggregate the force and torque of ROOT frame.
	for (size_t i = root.habegin; i < root.haend; ++i)
	{
		forc.front() += deri[i];
		torq.front() += cross(coor[i] - orig.front(), deri[i]);
	}

	// Save the aggregated force and torque to g.
	g[0] = forc.front()[0];
	g[1] = forc.front()[1];
	g[2] = forc.front()[2];
	g[3] = torq.front()[0];
	g[4] = torq.front()[1];
	g[5] = torq.front()[2];

	return true;
}

result ligand::compose_result(const double e, const double f, const conformation& conf, bool from_docking) const
{
	vector<array<double, 3>> orig(num_frames);
	vector<array<double, 4>> oriq(num_frames);
	vector<array<double, 9>> orim(num_frames);
	vector<array<double, 3>> heavy_atoms(num_heavy_atoms);
	vector<array<double, 3>> hydrogens(num_hydrogens);

	orig.front() = conf.position;
	oriq.front() = conf.orientation;
	orim.front() = qtn4_to_mat3(conf.orientation);

	// Calculate the coor of both heavy atoms and hydrogens of ROOT frame.
	const frame& root = frames.front();
	for (size_t i = root.habegin; i < root.haend; ++i)
	{
		heavy_atoms[i] = orig.front() + orim.front() * this->heavy_atoms[i].coord;
	}
	for (size_t i = root.hybegin; i < root.hyend; ++i)
	{
		hydrogens[i]   = orig.front() + orim.front() * this->hydrogens[i].coord;
	}

	// Calculate the coor of both heavy atoms and hydrogens of BRANCH frames.
	for (size_t k = 1, t = 0; k < num_frames; ++k)
	{
		const frame& f = frames[k];

		// Update origin.
		orig[k] = orig[f.parent] + orim[f.parent] * f.parent_rotorY_to_current_rotorY;

		// Update orientation.
		oriq[k] = vec4_to_qtn4(orim[f.parent] * f.parent_rotorX_to_current_rotorY, f.active ? conf.torsions[t++] : 0) * oriq[f.parent];
		orim[k] = qtn4_to_mat3(oriq[k]);

		// Update coor.
		for (size_t i = f.habegin; i < f.haend; ++i)
		{
			heavy_atoms[i] = orig[k] + orim[k] * this->heavy_atoms[i].coord;
		}
		for (size_t i = f.hybegin; i < f.hyend; ++i)
		{
			hydrogens[i]   = orig[k] + orim[k] * this->hydrogens[i].coord;
		}
	}

	return result(e, f, from_docking, move(heavy_atoms), move(hydrogens));
}

double ligand::calculate_rf_score(const result& r, const receptor& rec, const forest& f) const
{
	array<double, tree::nv> x{};
	for (size_t i = 0; i < num_heavy_atoms; ++i)
	{
		const atom& la = heavy_atoms[i];
		for (const atom& ra : rec.atoms)
		{
			const auto ds = distance_sqr(r.heavy_atoms[i], ra.coord);
			if (ds >= 144) continue; // RF-Score cutoff 12A
			if (!la.rf_unsupported() && !ra.rf_unsupported())
			{
				++x[(la.rf << 2) + ra.rf];
			}
			if (ds >= 64) continue; // Vina score cutoff 8A
			if (!la.xs_unsupported() && !ra.xs_unsupported())
			{
				scoring_function::score(x.data() + 36, la.xs, ra.xs, ds);
			}
		}
	}
	x.back() = flexibility_penalty_factor;
	return f(x);
}

void ligand::write_models(const path& output_ligand_path, const vector<result>& results, const receptor& rec) const
{
	const size_t num_results = results.size();
	assert(num_results);

	// Dump binding conformations to the output ligand file.
	ofstream ofs(output_ligand_path);
	ofs.setf(ios::fixed, ios::floatfield);
	for (size_t k = 0; k < num_results; ++k)
	{
		const result& r = results[k];
		ofs << "MODEL     " << setw(4) << (k + 1) << endl << setprecision(2);
		if (!r.from_docking)
			ofs << "REMARK THIS MODEL IS IDENTICAL TO THE INPUT LIGAND AND NOT FROM DOCKING" << endl;
		ofs << "REMARK 921   NORMALIZED FREE ENERGY PREDICTED BY IDOCK:" << setw(8) << r.e_nd    << " KCAL/MOL" << endl
			<< "REMARK 922        TOTAL FREE ENERGY PREDICTED BY IDOCK:" << setw(8) << r.e       << " KCAL/MOL" << endl
			<< "REMARK 923 INTER-LIGAND FREE ENERGY PREDICTED BY IDOCK:" << setw(8) << r.f       << " KCAL/MOL" << endl
			<< "REMARK 924 INTRA-LIGAND FREE ENERGY PREDICTED BY IDOCK:" << setw(8) << r.e - r.f << " KCAL/MOL" << endl
			<< "REMARK 927      BINDING AFFINITY PREDICTED BY RF-SCORE:" << setw(8) << r.rf      << " PKD" << endl << setprecision(3);

		size_t heavy_atom = 0;
		size_t hydrogen = 0;
		for (const auto& line : lines)
		{
			if (line.size() >= 78) // This line starts with "ATOM" or "HETATM".
			{
				const bool is_hydrogen = line[77] == 'H' && (line.length() == 78 || line[78] == ' ' || line[78] == 'D');
				const double free_energy = is_hydrogen ? 0 : r.e_heavy_atoms[heavy_atom];
				const array<double, 3>& coordinate = is_hydrogen ? r.hydrogens[hydrogen++] : r.heavy_atoms[heavy_atom++];
				ofs << line.substr(0, 30)
					<< setw(8) << coordinate[0]
					<< setw(8) << coordinate[1]
					<< setw(8) << coordinate[2]
					<< line.substr(54, 16)
					<< setw(6) << free_energy
					<< setw(4) << left << line.substr(76) << right;
			}
			else // This line starts with "ROOT", "ENDROOT", "BRANCH", "ENDBRANCH", TORSDOF", which will not change during docking.
			{
				ofs << line;
			}
			ofs << '\n';
		}
		ofs << "ENDMDL\n";
		assert(heavy_atom == r.heavy_atoms.size());
		assert(hydrogen == r.hydrogens.size());
	}
}

//! Revisit a result and calculate inter-molecular free energy components to every ligand atom and from every receptor residue.
void ligand::calculate_by_comp(result& result, const scoring_function& sf, const receptor& rec, vector<bool>& mask) const
{
	auto& e_residues = result.e_residues;
	auto& e_heavy_atoms = result.e_heavy_atoms;
	assert(e_residues.empty());
	assert(e_heavy_atoms.empty());

	e_residues.resize(rec.residues.size());
	e_heavy_atoms.resize(num_heavy_atoms);

	for (size_t k = 0; k < num_heavy_atoms; ++k)
	{
		const auto& coord = result.heavy_atoms[k];
		const auto xs = heavy_atoms[k].xs;

		// Find the index of the current coord.
		const auto index = rec.index(coord);
		const auto coord_grid = rec.coord(index); // Aligned coordinate.

		// Iterate through all atom ids contributing energy to the grid.
		for (const auto& a : rec.atoms)
		{
			assert(!a.is_hydrogen());

			auto r2 = distance_sqr(a.coord, coord_grid);
			if (r2 < scoring_function::cutoff_sqr)
			{
				size_t r_offset = static_cast<size_t>(sf.ns * r2);
				assert(r_offset < sf.nr);

				size_t p = mp(a.xs, xs);
				assert(p < sf.np);

				assert(!sf.e[p].empty());
				const double e0 = sf.e[p][r_offset]; // Read from template.

				scoring_function::score(e_residues[a.residue].data(), a.xs, xs, r2);
				e_residues[a.residue].back() += e0; // Aggregate the energy for the residue the ligand atom locates in.
				e_heavy_atoms[k] += e0; // Aggregate the energy for the ligand atom.

				// Mark contributing residue.
				mask[a.residue] = true;
			}
		}
	}

	// Weighten the score components.
	for (auto& e_residue : e_residues)
	{
		for (size_t i = 0; i < scoring_function::weights.size(); i++)
			e_residue[i] *= scoring_function::weights[i];
	}
}

void ligand::monte_carlo(vector<result>& results, const size_t seed, const scoring_function& sf, const receptor& rec) const
{
	// Define constants.
	static const double pi = 3.1415926535897932; //!< Pi.
	static const size_t num_alphas = 5; //!< Number of alpha values for determining step size in BFGS
	const size_t num_mc_iterations = 100 * num_heavy_atoms; //!< The number of iterations correlates to the complexity of ligand.
	const size_t num_entities = 2 + num_active_torsions; // Number of entities to mutate.
	const size_t num_variables = 6 + num_active_torsions; // Number of variables to optimize.
	const double e_upper_bound = static_cast<double>(4 * num_heavy_atoms); // A conformation will be droped if its free energy is not better than e_upper_bound.
	const double required_square_error = static_cast<double>(1 * num_heavy_atoms); // Ligands with RMSD < 1.0 will be clustered into the same cluster.

	mt19937_64 rng(seed);
	uniform_real_distribution<double> u01(0, 1);
	uniform_real_distribution<double> u11(-1, 1);
	uniform_real_distribution<double> upi(-pi, pi);
	uniform_real_distribution<double> ub0(rec.corner0[0], rec.corner1[0]);
	uniform_real_distribution<double> ub1(rec.corner0[1], rec.corner1[1]);
	uniform_real_distribution<double> ub2(rec.corner0[2], rec.corner1[2]);
	uniform_int_distribution<size_t> uen(0, num_entities - 1);
	normal_distribution<double> n01(0, 1);

	// Generate an initial random conformation c0, and evaluate it.
	conformation c0(num_active_torsions);
	double e0, f0;
	change g0(num_active_torsions);
	bool valid_conformation = false;
	for (size_t i = 0; (i < 1000) && (!valid_conformation); ++i)
	{
		// Randomize conformation c0.
		c0.position = array<double, 3>{{ub0(rng), ub1(rng), ub2(rng)}};
		c0.orientation = normalize(array<double, 4>{{n01(rng), n01(rng), n01(rng), n01(rng)}});
		for (size_t i = 0; i < num_active_torsions; ++i)
		{
			c0.torsions[i] = upi(rng);
		}
		valid_conformation = evaluate(c0, sf, rec, e_upper_bound, e0, f0, g0);
	}
	if (!valid_conformation) return;
	double best_e = e0; // The best free energy so far.

	// Initialize necessary variables for BFGS.
	conformation c1(num_active_torsions), c2(num_active_torsions); // c2 = c1 + ap.
	double e1, f1, e2, f2;
	change g1(num_active_torsions), g2(num_active_torsions);
	change p(num_active_torsions); // Descent direction.
	double alpha, pg1, pg2; // pg1 = p * g1. pg2 = p * g2.
	size_t num_alpha_trials;

	// Initialize the inverse Hessian matrix to identity matrix.
	// An easier option that works fine in practice is to use a scalar multiple of the identity matrix,
	// where the scaling factor is chosen to be in the range of the eigenvalues of the true Hessian.
	// See N&R for a recipe to find this initializer.
	vector<double> h1(num_variables * (num_variables + 1) >> 1, 0); // Symmetric triangular matrix.
	for (size_t i = 0; i < num_variables; ++i)
		h1[mr(i, i)] = 1;

	// Initialize necessary variables for updating the Hessian matrix h.
	vector<double> h(h1);
	change y(num_active_torsions); // y = g2 - g1.
	change mhy(num_active_torsions); // mhy = -h * y.
	double yhy, yp, ryp, pco;

	for (size_t mc_i = 0; mc_i < num_mc_iterations; ++mc_i)
	{
		size_t mutation_entity;

		// Mutate c0 into c1, and evaluate c1.
		do
		{
			// Make a copy, so the previous conformation is retained.
			c1 = c0;

			// Determine an entity to mutate.
			mutation_entity = uen(rng);
			assert(mutation_entity < num_entities);
			if (mutation_entity < num_active_torsions) // Mutate an active torsion.
			{
				c1.torsions[mutation_entity] = upi(rng);
			}
			else if (mutation_entity == num_active_torsions) // Mutate position.
			{
				c1.position += array<double, 3>{{u11(rng), u11(rng), u11(rng)}};
			}
			else // Mutate orientation.
			{
				c1.orientation = vec3_to_qtn4(0.01 * array<double, 3>{{u11(rng), u11(rng), u11(rng)}}) * c1.orientation;
				assert(normalized(c1.orientation));
			}
		} while (!evaluate(c1, sf, rec, e_upper_bound, e1, f1, g1));

		// Initialize the Hessian matrix to identity.
		h = h1;

		// Given the mutated conformation c1, use BFGS to find a local minimum.
		// The conformation of the local minimum is saved to c2, and its derivative is saved to g2.
		// http://en.wikipedia.org/wiki/BFGS_method
		// http://en.wikipedia.org/wiki/Quasi-Newton_method
		// The loop breaks when an appropriate alpha cannot be found.
		while (true)
		{
			// Calculate p = -h*g, where p is for descent direction, h for Hessian, and g for gradient.
			for (size_t i = 0; i < num_variables; ++i)
			{
				double sum = 0;
				for (size_t j = 0; j < num_variables; ++j)
					sum += h[mp(i, j)] * g1[j];
				p[i] = -sum;
			}

			// Calculate pg = p*g = -h*g^2 < 0
			pg1 = 0;
			for (size_t i = 0; i < num_variables; ++i)
				pg1 += p[i] * g1[i];

			// Perform a line search to find an appropriate alpha.
			// Try different alpha values for num_alphas times.
			// alpha starts with 1, and shrinks to alpha_factor of itself iteration by iteration.
			alpha = 1.0;
			for (num_alpha_trials = 0; num_alpha_trials < num_alphas; ++num_alpha_trials)
			{
				// Obtain alpha from the precalculated alpha values.
				alpha *= 0.1;

				// Calculate c2 = c1 + ap.
				c2.position = c1.position + alpha * array<double, 3>{{p[0], p[1], p[2]}};
				assert(normalized(c1.orientation));
				c2.orientation = vec3_to_qtn4(alpha * array<double, 3>{{p[3], p[4], p[5]}}) * c1.orientation;
				assert(normalized(c2.orientation));
				for (size_t i = 0; i < num_active_torsions; ++i)
				{
					c2.torsions[i] = c1.torsions[i] + alpha * p[6 + i];
				}

				// Evaluate c2, subject to Wolfe conditions http://en.wikipedia.org/wiki/Wolfe_conditions
				// 1) Armijo rule ensures that the step length alpha decreases f sufficiently.
				// 2) The curvature condition ensures that the slope has been reduced sufficiently.
				if (evaluate(c2, sf, rec, e1 + 0.0001 * alpha * pg1, e2, f2, g2))
				{
					pg2 = 0;
					for (size_t i = 0; i < num_variables; ++i)
						pg2 += p[i] * g2[i];
					if (pg2 >= 0.9 * pg1)
						break; // An appropriate alpha is found.
				}
			}

			// If an appropriate alpha cannot be found, exit the BFGS loop.
			if (num_alpha_trials == num_alphas) break;

			// Update Hessian matrix h.
			for (size_t i = 0; i < num_variables; ++i) // Calculate y = g2 - g1.
				y[i] = g2[i] - g1[i];
			for (size_t i = 0; i < num_variables; ++i) // Calculate mhy = -h * y.
			{
				double sum = 0;
				for (size_t j = 0; j < num_variables; ++j)
					sum += h[mp(i, j)] * y[j];
				mhy[i] = -sum;
			}
			yhy = 0;
			for (size_t i = 0; i < num_variables; ++i) // Calculate yhy = -y * mhy = -y * (-hy).
				yhy -= y[i] * mhy[i];
			yp = 0;
			for (size_t i = 0; i < num_variables; ++i) // Calculate yp = y * p.
				yp += y[i] * p[i];
			ryp = 1 / yp;
			pco = ryp * (ryp * yhy + alpha);
			for (size_t i = 0; i < num_variables; ++i)
				for (size_t j = i; j < num_variables; ++j) // includes i
				{
					h[mr(i, j)] += ryp * (mhy[i] * p[j] + mhy[j] * p[i]) + pco * p[i] * p[j];
				}

			// Move to the next iteration.
			c1 = c2;
			e1 = e2;
			f1 = f2;
			g1 = g2;
		}

		// Accept c1 according to Metropolis criteria.
		const double delta = e0 - e1;
		if (delta > 0 || u01(rng) < exp(delta))
		{
			// best_e is the best energy of all the conformations in the container.
			// e1 will be saved if and only if it is even better than the best one.
			if (e1 < best_e || results.size() < results.capacity())
			{
				result::push(results, compose_result(e1, f1, c1, true), required_square_error);
				if (e1 < best_e) best_e = e0;
			}

			// Save c1 into c0.
			c0 = c1;
			e0 = e1;
		}
	}
}
