/*

   Copyright (c) 2011, The Chinese University of Hong Kong

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

*/

#include "fstream.hpp"
#include "parsing_error.hpp"
#include "ligand.hpp"

namespace idock
{
	ligand::ligand(const path& p) : num_active_torsions(0)
	{
		// Initialize necessary variables for constructing a ligand.
		lines.reserve(200); // A ligand typically consists of <= 200 lines.
		frames.reserve(30); // A ligand typically consists of <= 30 frames.
		frames.push_back(frame(0, 0, 0, 0)); // ROOT is also treated as a frame. The parent and rotorX of ROOT frame are dummy.
		heavy_atoms.reserve(100); // A ligand typically consists of <= 100 heavy atoms.
		hydrogens.reserve(50); // A ligand typically consists of <= 50 hydrogens.

		// Initialize helper variables for parsing.
		vector<size_t> numbers; ///< Atom serial numbers.
		numbers.reserve(100); // A ligand typically consists of <= 100 heavy atoms.
		size_t current = 0; // Index of current frame, initialized to ROOT frame.
		frame* f = &frames.front(); // Pointer to the current frame.
		f->rotorY = 0; // Assume the rotorY of ROOT frame is the first atom.
		size_t num_lines = 0; // Used to track line number for reporting parsing errors, if any.
		string line;
		line.reserve(79); // According to PDBQT specification, the last item AutoDock atom type locates at 1-based [78, 79].

		// Parse ROOT, ATOM/HETATM, ENDROOT, BRANCH, ENDBRANCH, TORSDOF.
		ifstream in(p); // Parsing starts. Open the file stream as late as possible.
		while (getline(in, line))
		{
			++num_lines;
			if (starts_with(line, "ATOM") || starts_with(line, "HETATM"))
			{
				// Whenever an ATOM/HETATM line shows up, the current frame must be the last one.
				BOOST_ASSERT(current == frames.size() - 1);
				BOOST_ASSERT(f == &frames.back());

				// This line will be dumped to the output ligand file.
				lines.push_back(line);

				// Parse and validate AutoDock4 atom type.
				const string ad_type_string = line.substr(77, isspace(line[78]) ? 1 : 2);
				const size_t ad = parse_ad_type_string(ad_type_string);
				if (ad == AD_TYPE_SIZE) throw parsing_error(p, num_lines, "Atom type " + ad_type_string + " is not supported by idock.");

				// Parse the Cartesian coordinate.
				const atom a(vec3(right_cast<fl>(line, 31, 38), right_cast<fl>(line, 39, 46), right_cast<fl>(line, 47, 54)), ad);

				// For a hydrogen, save it.
				if (a.is_hydrogen())
				{
					hydrogens.push_back(a);

					// For a polar hydrogen, the bonded hetero atom must be a hydrogen bond donor.
					if (ad == AD_TYPE_HD)
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
				}
				else // It is a heavy atom.
				{
					heavy_atoms.push_back(a);
					numbers.push_back(right_cast<size_t>(line, 7, 11));
				}
			}
			else if (starts_with(line, "BRANCH"))
			{
				// This line will be dumped to the output ligand file.
				lines.push_back(line);

				// Parse "BRANCH   X   Y". X and Y are right-justified and 4 characters wide.
				const size_t x = right_cast<size_t>(line, 7, 10);

				// Find the corresponding heavy atom with x as its atom serial number in the current frame.
				for (size_t i = f->habegin; true; ++i)
				{
					if (numbers[i] == x)
					{
						// Insert a new frame whose parent is the current frame.
						frames.push_back(frame(current, i, heavy_atoms.size(), hydrogens.size()));
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
			else if (starts_with(line, "ENDBRANCH"))
			{
				// This line will be dumped to the output ligand file.
				lines.push_back(line);

				// A frame may be empty, e.g. "BRANCH   4   9" is immediately followed by "ENDBRANCH   4   9".
				// This emptiness is likely to be caused by invalid input structure, especially when all the atoms are located in the same plane.
				if (f->habegin == heavy_atoms.size()) throw parsing_error(p, num_lines, "An empty BRANCH has been detected, indicating the input ligand structure is probably invalid.");

				// Parse "ENDBRANCH   X   Y". X and Y are right-justified and 4 characters wide.
				const size_t y = right_cast<size_t>(line, 14, 17);

				// Find the corresponding heavy atom with y as its atom serial number in the current frame.
				for (size_t i = f->habegin; true; ++i)
				{
					if (numbers[i] == y)
					{
						// Set rotorY of current frame.
						f->rotorY = i;
						break;
					}
				}

				// If the current frame consists of rotor Y and a few hydrogens only, e.g. -OH and -NH2,
				// the torsion of this frame will have no effect on scoring and is thus redundant.
				if ((current == frames.size() - 1) && (f->habegin + 1 == heavy_atoms.size()))
				{
					f->active = false;
				}
				else
				{
					++num_active_torsions;
				}

				// Now the parent of the following frame is the parent of current frame.
				current = frames[current].parent;

				// Update the pointer to the current frame.
				f = &frames[current];
			}
			else if (starts_with(line, "ROOT") || starts_with(line, "ENDROOT") || starts_with(line, "TORSDOF"))
			{
				// This line will be dumped to the output ligand file.
				lines.push_back(line);
			}
		}
		in.close(); // Parsing finishes. Close the file stream as soon as possible.
		BOOST_ASSERT(lines.size() <= num_lines); // Some lines like "REMARK", "WARNING", "TER" will not be dumped to the output ligand file.

		BOOST_ASSERT(current == 0); // current should remain its original value if "BRANCH" and "ENDBRANCH" properly match each other.
		BOOST_ASSERT(f == &frames.front()); // The frame pointer should remain its original value if "BRANCH" and "ENDBRANCH" properly match each other.		

		// Determine num_frames, num_torsions, flexibility_penalty_factor, and num_heavy_atoms_inverse.
		num_frames = frames.size();
		BOOST_ASSERT(num_frames >= 1);
		num_torsions = num_frames - 1;
		BOOST_ASSERT(num_torsions + 1 == num_frames);
		BOOST_ASSERT(num_torsions >= num_active_torsions);
		flexibility_penalty_factor = 1 / (1 + 0.05846 * (num_active_torsions + 0.5 * (num_torsions - num_active_torsions)));
		num_heavy_atoms = heavy_atoms.size();
		num_hydrogens = hydrogens.size();
		BOOST_ASSERT(num_heavy_atoms + num_hydrogens + (num_torsions << 1) + 3 <= num_lines); // ATOM/HETATM lines + BRANCH/ENDBRANCH lines + ROOT/ENDROOT/TORSDOF lines + REMARK lines (if any) == num_lines
		num_heavy_atoms_inverse = static_cast<fl>(1) / num_heavy_atoms;
		frames.back().haend = num_heavy_atoms;
		frames.back().hyend = num_hydrogens;

		// Dehydrophobicize carbons if necessary.
		for (size_t k = 0; k < num_frames; ++k)
		{
			frame& f = frames[k];
			for (size_t i = f.habegin; i < f.haend; ++i)
			{
				const atom& a = heavy_atoms[i];
				if (!a.is_hetero()) continue; // a is a hetero atom.

				for (size_t j = f.habegin; j < f.haend; ++j)
				{
					atom& b = heavy_atoms[j];
					if (b.is_hetero()) continue; // b is a carbon atom.

					// If the carbon atom b is bonded to the hetero atom a, it is no longer a hydrophobic atom.
					if (a.is_neighbor(b))
					{
						b.dehydrophobicize();
					}
				}
			}

			if (k) // k > 0
			{
				atom& rotorY = heavy_atoms[f.rotorY];
				atom& rotorX = heavy_atoms[f.rotorX];
				if ((rotorY.is_hetero()) && (!rotorX.is_hetero())) rotorX.dehydrophobicize();
				if ((rotorX.is_hetero()) && (!rotorY.is_hetero())) rotorY.dehydrophobicize();
			}
		}

		// Initialize parent_rotorY_to_current_rotorY and parent_rotorX_to_current_rotorY of BRANCH frames.
		for (size_t k = 1; k < num_frames; ++k)
		{
			frame& f = frames[k];
			const frame& p = frames[f.parent];
			f.parent_rotorY_to_current_rotorY =  heavy_atoms[f.rotorY].coordinate - heavy_atoms[p.rotorY].coordinate;
			f.parent_rotorX_to_current_rotorY = (heavy_atoms[f.rotorY].coordinate - heavy_atoms[f.rotorX].coordinate).normalize();
		}

		// Reserve enough capacity for bonds.
		vector<vector<size_t>> bonds(num_heavy_atoms);
		for (size_t i = 0; i < num_heavy_atoms; ++i)
		{
			bonds[i].reserve(4); // An atom typically consists of <= 4 bonds.
		}

		for (size_t k = 0; k < num_frames; ++k)
		{
			const frame& f = frames[k];
			for (size_t i = f.habegin; i < f.haend; ++i)
			{
				const atom& a1 = heavy_atoms[i];

				for (size_t j = i + 1; j < f.haend; ++j)
				{
					const atom& a2 = heavy_atoms[j];
					if (a1.is_neighbor(a2))
					{
						bonds[i].push_back(j);
						bonds[j].push_back(i);
					}
				}
			}
			if (k) // k > 0
			{
				bonds[f.rotorY].push_back(f.rotorX);
				bonds[f.rotorX].push_back(f.rotorY);
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
				const vector<size_t>& i0_bonds = bonds[i];
				const size_t num_i0_bonds = i0_bonds.size();
				for (size_t i0 = 0; i0 < num_i0_bonds; ++i0)
				{
					const size_t b1 = i0_bonds[i0];
					if (find(neighbors.begin(), neighbors.end(), b1) == neighbors.end())
					{
						neighbors.push_back(b1);
					}
					const vector<size_t>& i1_bonds = bonds[b1];
					const size_t num_i1_bonds = i1_bonds.size();
					for (size_t i1 = 0; i1 < num_i1_bonds; ++i1)
					{
						const size_t b2 = i1_bonds[i1];
						if (find(neighbors.begin(), neighbors.end(), b2) == neighbors.end())
						{
							neighbors.push_back(b2);
						}
						const vector<size_t>& i2_bonds = bonds[b2];
						const size_t num_i2_bonds = i2_bonds.size();
						for (size_t i2 = 0; i2 < num_i2_bonds; ++i2)
						{
							const size_t b3 = i2_bonds[i2];
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
					for (size_t j = f2.habegin; j < f2.haend; ++j)
					{
						if (((k1 == f2.parent) && ((j == f2.rotorY) || (i == f2.rotorX))) || (find(neighbors.begin(), neighbors.end(), j) != neighbors.end())) continue;
						const size_t type_pair_index = triangular_matrix_permissive_index(heavy_atoms[i].xs, heavy_atoms[j].xs);
						interacting_pairs.push_back(interacting_pair(i, j, type_pair_index));
					}
				}

				// Clear the current neighbor set for the next atom.
				neighbors.clear();
			}
		}

		// Update heavy_atoms[].coordinate and hydrogens[] to relative coordinates.
		for (size_t k = 0; k < num_frames; ++k)
		{
			frame& f = frames[k];
			const vec3 origin = heavy_atoms[f.rotorY].coordinate;
			for (size_t i = f.habegin; i < f.haend; ++i)
			{
				heavy_atoms[i].coordinate -= origin;
			}
			for (size_t i = f.hybegin; i < f.hyend; ++i)
			{
				hydrogens[i].coordinate -= origin;
			}
		}
	}

	vector<size_t> ligand::get_atom_types() const
	{
		vector<size_t> atom_types;
		atom_types.reserve(10); // A ligand typically consists of <= 10 XScore atom types.
		for (size_t i = 0; i < num_heavy_atoms; ++i)
		{
			const size_t t = heavy_atoms[i].xs;
			if (find(atom_types.begin(), atom_types.end(), t) == atom_types.end()) atom_types.push_back(t);
		}
		return atom_types;
	}

	bool ligand::evaluate(const conformation& conf, const scoring_function& sf, const box& b, const vector<array3d<fl>>& grid_maps, const fl e_upper_bound, fl& e, fl& f, change& g) const
	{
		if (!b.within(conf.position))
			return false;
		
		// Define frame-wide conformational variables.
		vector<vec3> origin; ///< Origin coordinate, which is rotorY.
		vector<qt>   orientation_q; ///< Orientation in the form of quaternion.
		vector<mat3> orientation_m; ///< Orientation in the form of 3x3 matrix.
		vector<vec3> axis; ///< Vector pointing from rotor Y to rotor X.
		vector<vec3> force; ///< Aggregated derivatives of heavy atoms.
		vector<vec3> torque; /// Torque of the force.
		origin.resize(num_frames);
		orientation_q.resize(num_frames);
		orientation_m.resize(num_frames);
		axis.resize(num_frames);
		force.resize(num_frames);
		torque.resize(num_frames);

		// Define atom-wide conformational variables.
		vector<vec3> coordinates; ///< Heavy atom coordinates.
		vector<vec3> derivatives; ///< Heavy atom derivatives.
		vector<fl> energies; ///< Heavy atom free energies.
		coordinates.resize(num_heavy_atoms);
		derivatives.resize(num_heavy_atoms);
		energies.resize(num_heavy_atoms);


		// Apply position and orientation to ROOT frame.
		const frame& root = frames.front();
		origin.front() = conf.position;
		orientation_q.front() = conf.orientation;
		orientation_m.front() = quaternion_to_matrix(conf.orientation);
		for (size_t i = root.habegin; i < root.haend; ++i)
		{
			coordinates[i] = origin.front() + orientation_m.front() * heavy_atoms[i].coordinate;
			if (!b.within(coordinates[i]))
				return false;
		}

		// Apply torsions to BRANCH frames.
		for (size_t k = 1, t = 0; k < num_frames; ++k)
		{
			const frame& f = frames[k];

			// Update origin.
			origin[k] = origin[f.parent] + orientation_m[f.parent] * f.parent_rotorY_to_current_rotorY;
			if (!b.within(origin[k]))
				return false;

			// If the current BRANCH frame does not have an active torsion, skip it.
			if (!f.active)
			{
				BOOST_ASSERT(f.habegin + 1 == f.haend);
				BOOST_ASSERT(f.habegin == f.rotorY);
				coordinates[f.rotorY] = origin[k];
				continue;
			}

			// Update orientation.
			BOOST_ASSERT(f.parent_rotorX_to_current_rotorY.normalized());
			axis[k] = orientation_m[f.parent] * f.parent_rotorX_to_current_rotorY;
			BOOST_ASSERT(axis[k].normalized());
			orientation_q[k] = axis_angle_to_quaternion(axis[k], conf.torsions[t++]) * orientation_q[f.parent];
			BOOST_ASSERT(quaternion_is_normalized(orientation_q[k]));
			orientation_m[k] = quaternion_to_matrix(orientation_q[k]);

			// Update coordinates.
			for (size_t i = f.habegin; i < f.haend; ++i)
			{
				coordinates[i] = origin[k] + orientation_m[k] * heavy_atoms[i].coordinate;
				if (!b.within(coordinates[i]))
					return false;
			}
		}

		// Check steric clash between atoms of different frames except for (rotorX, rotorY) pair.
		//for (size_t k1 = num_frames - 1; k1 > 0; --k1)
		//{
		//	const frame& f1 = frames[k1];
		//	const size_t num_heavy_atoms1 = f1.heavy_atoms.size();
		//	for (size_t i1 = 0; i1 < num_heavy_atoms1; ++i1)
		//	{
		//		for (size_t k2 = 0; k2 < k1; ++k2)
		//		{
		//			const frame& f2 = frames[k2];
		//			const size_t num_heavy_atoms2 = f2.heavy_atoms.size();
		//			for (size_t i2 = 0; i2 < num_heavy_atoms2; ++i2)
		//			{
		//				if ((distance_sqr(f1.coordinates[i1], f2.coordinates[i2]) < sqr(f1.heavy_atoms[i1].covalent_radius() + f2.heavy_atoms[i2].covalent_radius())) && (!((k2 == f1.parent) && (i1 == 0) && (i2 == f1.rotorX))))
		//					return false;
		//			}
		//		}
		//	}
		//}

		e = 0;
		for (size_t i = 0; i < num_heavy_atoms; ++i)
		{
			// Retrieve the grid map in need.
			const array3d<fl>& grid_map = grid_maps[heavy_atoms[i].xs];
			BOOST_ASSERT(grid_map.initialized());

			// Find the index and fraction of the current coordinates.
			const array<size_t, 3> index = b.grid_index(coordinates[i]);

			// Assert the validity of index.
			BOOST_ASSERT(index[0] < b.num_grids[0]);				
			BOOST_ASSERT(index[1] < b.num_grids[1]);
			BOOST_ASSERT(index[2] < b.num_grids[2]);

			// (x0, y0, z0) is the beginning corner of the partition.
			const size_t x0 = index[0];
			const size_t y0 = index[1];
			const size_t z0 = index[2];
			const fl e000 = grid_map(x0, y0, z0);
			energies[i] = e000;

			// The derivative of probe atoms can be precalculated at the cost of massive memory storage.
			const fl e100 = grid_map(x0 + 1, y0,     z0    );
			const fl e010 = grid_map(x0,     y0 + 1, z0    );
			const fl e001 = grid_map(x0,     y0,     z0 + 1);
			derivatives[i][0] = (e100 - e000) * b.grid_granularity_inverse;
			derivatives[i][1] = (e010 - e000) * b.grid_granularity_inverse;
			derivatives[i][2] = (e001 - e000) * b.grid_granularity_inverse;

			e += e000; // Aggregate the energy.
		}

		// Save inter-molecular free energy into f.
		f = e;

		// Calculate intra-ligand free energy.
		const size_t num_interacting_pairs = interacting_pairs.size();
		for (size_t i = 0; i < num_interacting_pairs; ++i)
		{
			const interacting_pair& p = interacting_pairs[i];
			const vec3 r = coordinates[p.i2] - coordinates[p.i1];
			const fl r2 = r.norm_sqr();
			if (r2 < scoring_function::Cutoff_Sqr)
			{
				const scoring_function_element element = sf.evaluate(p.type_pair_index, r2);
				e += element.e;
				const vec3 derivative = element.dor * r;
				derivatives[p.i1] -= derivative;
				derivatives[p.i2] += derivative;
			}
		}

		// If the free energy is no better than the upper bound, refuse this conformation.
		if (e >= e_upper_bound) return false;

		// Initialize force and torque. TODO: try assign().
		for (size_t k = 0; k < num_frames; ++k)
		{
			force[k]  = zero3; // Initialize force to zero.
			torque[k] = zero3; // Initialize torque to zero.
		}

		// Calculate and aggregate the force and torque of BRANCH frames to their parent frame.
		for (size_t k = num_frames - 1, t = num_active_torsions; k > 0; --k)
		{
			const frame&  f = frames[k];

			for (size_t i = f.habegin; i < f.haend; ++i)
			{
				// The derivatives with respect to the position, orientation, and torsions
				// would be the negative total force acting on the ligand,
				// the negative total torque, and the negative torque projections, respectively,
				// where the projections refer to the torque applied to the branch moved by the torsion,
				// projected on its rotation axis.
				force[k]  += derivatives[i];
				torque[k] += cross_product(coordinates[i] - origin[k], derivatives[i]);
			}

			// Aggregate the force and torque of current frame to its parent frame.
			force[f.parent]  += force[k];
			torque[f.parent] += torque[k] + cross_product(origin[k] - origin[f.parent], force[k]);

			// If the current BRANCH frame does not have an active torsion, skip it.
			if (!f.active) continue;

			// Save the aggregated torque to torsion.
			g.torsions[--t] = torque[k] * axis[k]; // dot product
		}

		// Calculate and aggregate the force and torque of ROOT frame.
		for (size_t i = root.habegin; i < root.haend; ++i)
		{
			force.front()  += derivatives[i];
			torque.front() += cross_product(coordinates[i] - origin.front(), derivatives[i]);
		}

		// Save the aggregated force and torque to g.position and g.orientation.
		g.position    = force.front();
		g.orientation = torque.front();

		return true;
	}

	result ligand::compose_result(const fl e, const fl f, const conformation& conf) const
	{
		vector<vec3> origins(num_frames);
		vector<qt>   orientations_q(num_frames);
		vector<mat3> orientations_m(num_frames);
		vector<vec3> heavy_atoms(num_heavy_atoms);
		vector<vec3> hydrogens(num_hydrogens);

		origins.front() = conf.position;
		orientations_q.front() = conf.orientation;
		orientations_m.front() = quaternion_to_matrix(conf.orientation);

		// Calculate the coordinates of both heavy atoms and hydrogens of ROOT frame.
		const frame& root = frames.front();
		for (size_t i = root.habegin; i < root.haend; ++i)
		{
			heavy_atoms[i] = origins.front() + orientations_m.front() * this->heavy_atoms[i].coordinate;
		}
		for (size_t i = root.hybegin; i < root.hyend; ++i)
		{
			hydrogens[i]   = origins.front() + orientations_m.front() * this->hydrogens[i].coordinate;
		}

		// Calculate the coordinates of both heavy atoms and hydrogens of BRANCH frames.
		for (size_t k = 1, t = 0; k < num_frames; ++k)
		{
			const frame& f = frames[k];

			// Update origin.
			origins[k] = origins[f.parent] + orientations_m[f.parent] * f.parent_rotorY_to_current_rotorY;

			// Update orientation.
			orientations_q[k] = axis_angle_to_quaternion(orientations_m[f.parent] * f.parent_rotorX_to_current_rotorY, f.active ? conf.torsions[t++] : 0) * orientations_q[f.parent];
			orientations_m[k] = quaternion_to_matrix(orientations_q[k]);

			// Update coordinates.
			for (size_t i = f.habegin; i < f.haend; ++i)
			{
				heavy_atoms[i] = origins[k] + orientations_m[k] * this->heavy_atoms[i].coordinate;
			}
			for (size_t i = f.hybegin; i < f.hyend; ++i)
			{
				hydrogens[i]   = origins[k] + orientations_m[k] * this->hydrogens[i].coordinate;
			}
		}

		return result(e, f, static_cast<vector<vec3>&&>(heavy_atoms), static_cast<vector<vec3>&&>(hydrogens));
	}

	void ligand::write_models(const path& output_ligand, const ptr_vector<result>& results, const size_t num_conformations)
	{
		BOOST_ASSERT(num_conformations > 0);
		BOOST_ASSERT(num_conformations <= results.size());

		const size_t num_lines = lines.size();

		// Dump binding conformations to the output ligand file.
		using namespace std;
		ofstream out(output_ligand); // Dumping starts. Open the file stream as late as possible.
		out.setf(ios::fixed, ios::floatfield);
		out << setprecision(3);
		for (size_t i = 0; i < num_conformations; ++i)
		{
			const result& r = results[i];
			out << "MODEL     " << setw(4) << (i + 1) << '\n'
				<< "REMARK     FREE ENERGY PREDICTED BY IDOCK:" << setw(8) << r.e << " KCAL/MOL\n";
			for (size_t j = 0, heavy_atom = 0, hydrogen = 0; j < num_lines; ++j)
			{
				const string& line = lines[j];
				if (line.size() >= 79) // This line starts with "ATOM" or "HETATM"
				{
					const vec3& coordinate = line[77] == 'H' ? r.hydrogens[hydrogen++] : r.heavy_atoms[heavy_atom++];
					out << line.substr(0, 30)
						<< setw(8) << coordinate[0]
						<< setw(8) << coordinate[1]
						<< setw(8) << coordinate[2]
						<< line.substr(54);
				}
				else // This line starts with "ROOT", "ENDROOT", "BRANCH", "ENDBRANCH", TORSDOF", which will not change during docking.
				{
					out << line;
				}
				out << '\n';
			}
			out << "ENDMDL\n";
		}
		out.close(); // Dumping finishes. Close the file stream as soon as possible.
	}
}
