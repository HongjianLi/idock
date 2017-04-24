#pragma once
#ifndef IDOCK_LIGAND_HPP
#define IDOCK_LIGAND_HPP

#include <boost/filesystem/fstream.hpp>
#include "scoring_function.hpp"
#include "random_forest.hpp"
#include "atom.hpp"
#include "receptor.hpp"
using namespace boost::filesystem;

//! Represents a ROOT or a BRANCH in PDBQT structure.
class frame
{
public:
	size_t parent; //!< Frame array index pointing to the parent of current frame. For ROOT frame, this field is not used.
	size_t rotorXsrn; //!< Serial atom number of the parent frame atom which forms a rotatable bond with the rotorY atom of current frame.
	size_t rotorYsrn; //!< Serial atom number of the current frame atom which forms a rotatable bond with the rotorX atom of parent frame.
	size_t rotorXidx; //!< Index pointing to the parent frame atom which forms a rotatable bond with the rotorY atom of current frame.
	size_t rotorYidx; //!< Index pointing to the current frame atom which forms a rotatable bond with the rotorX atom of parent frame.
	size_t childYidx; //!< The exclusive ending index to the heavy atoms of the current frame.
	bool active; //!< Indicates if the current frame is active.
	array<float, 3> yy; //!< Vector pointing from the origin of parent frame to the origin of current frame.
	array<float, 3> xy; //!< Normalized vector pointing from rotor X of parent frame to rotor Y of current frame.
	vector<size_t> branches; //!< Indexes to child branches.

	//! Constructs an active frame, and relates it to its parent frame.
	explicit frame(const size_t parent, const size_t rotorXsrn, const size_t rotorYsrn, const size_t rotorXidx, const size_t rotorYidx) : parent(parent), rotorXsrn(rotorXsrn), rotorYsrn(rotorYsrn), rotorXidx(rotorXidx), rotorYidx(rotorYidx), active(true) {}

	//! Outputs a BRANCH line in PDBQT format.
	void output(boost::filesystem::ofstream& ofs) const;
};

//! Represents a ligand.
class ligand
{
public:
	path filename; //!< Filename of the input ligand.
	vector<frame> frames; //!< ROOT and BRANCH frames.
	vector<atom> atoms; //!< Heavy atoms. Coordinates are relative to frame origin, which is the first atom by default. Hydrogens are saved under heavy atoms.
	array<bool, scoring_function::n> xs; //!< Presence of XScore atom types.
	size_t nv; //!< Number of variables to optimize, which equals 6 plus the number of active frames.
	size_t nf; //!< Number of frames, both active and inactive.
	size_t na; //!< Number of heavy atoms.
	size_t np; //!< Number of non 1-4 interacting pairs.
	vector<float> affinities; //!< Binding affinities of predicted conformations.

	//! Constructs a ligand by parsing a ligand file in PDBQT format.
	explicit ligand(const path& p);

	//! Encodes the current ligand into an array of integers.
	void encode(int* const p) const;

	//! Writes conformations in PDBQT format to file.
	void write(const float* const ex, const path& output_folder_path, const size_t max_conformations, const size_t num_tasks, const receptor& rec, const forest& f, const scoring_function& sf);

	//! Gets the number of elements of the current ligand.
	size_t get_lig_elems() const;

	//! Gets the number of elements of a solution.
	size_t get_sln_elems() const;

	//! Gets the number of elements of a conformation.
	size_t get_cnf_elems() const;

	void load_from_path(const path& p);
private:
	//! Represents a pair of interacting atoms that are separated by 3 consecutive covalent bonds.
	class interacting_pair
	{
	public:
		size_t i0; //!< Index of atom 0.
		size_t i1; //!< Index of atom 1.
		size_t p_offset; //!< Type pair index to the scoring function. It can be precalculated to save the creating time of grid maps.

		//! Constructs a pair of non 1-4 interacting atoms.
		interacting_pair(const size_t i0, const size_t i1, const size_t p_offset) : i0(i0), i1(i1), p_offset(p_offset) {}
	};

	vector<interacting_pair> interacting_pairs; //!< Non 1-4 interacting pairs.
};

#endif
