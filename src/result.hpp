#pragma once
#ifndef IDOCK_RESULT_HPP
#define IDOCK_RESULT_HPP

#include <vector>

//! Represents a result found by BFGS local optimization for later clustering.
class result
{
public:
	double e; //!< Free energy.
	double f; //!< Inter-molecular free energy.
	double e_nd; //!< Normalized free energy.
	double rf; //!< RF-Score binding affinity.
	bool from_docking; //!< Indicates if the result is created from docking.
	vector<array<double, 3>> heavy_atoms; //!< Heavy atom coordinates.
	vector<array<double, 3>> hydrogens; //!< Hydrogen atom coordinates.
	vector<double> e_heavy_atoms; //!< Per heavy atom inter-molecular free energy.
	vector<array<double, 6>> e_residues; //!< Per residue total of inter-molecular free energy, including weighted 5 term score components plus the weighted overall score.

	//! Constructs a result from free energy e, force f, whether created from docking, heavy atom coordinates, hydrogen atom coordinates, per heavy atom energy and per residue energy.
	explicit result(const double e, const double f, bool from_docking, vector<array<double, 3>>&& heavy_atoms, vector<array<double, 3>>&& hydrogens, vector<double>&& e_heavy_atoms, vector<array<double, 6>>&& e_residues)
		: e(e)
		, f(f)
		, e_nd()
		, rf()
		, from_docking(from_docking)
		, heavy_atoms(move(heavy_atoms))
		, hydrogens(move(hydrogens))
		, e_heavy_atoms(move(e_heavy_atoms))
		, e_residues(move(e_residues))
	{
	}

	//! Constructs a result from free energy e, force f, whether created from docking, heavy atom coordinates and hydrogen atom coordinates.
	explicit result(const double e, const double f, bool from_docking, vector<array<double, 3>>&& heavy_atoms, vector<array<double, 3>>&& hydrogens)
		: result(e, f, from_docking, move(heavy_atoms), move(hydrogens), {}, {})
	{
	}

	result(const result&) = default;
	result(result&&) = default;
	result& operator=(const result&) = default;
	result& operator=(result&&) = default;

	//! For sorting vector<result>.
	bool operator<(const result& r) const
	{
		return e < r.e;
	}

	//! Clusters a result into an existing result set with a minimum RMSD requirement.
	static void push(vector<result>& results, result&& r, const double required_square_error);
};

#endif
