#pragma once
#ifndef IDOCK_SCORING_FUNCTION_HPP
#define IDOCK_SCORING_FUNCTION_HPP

#include <vector>
#include <array>
using namespace std;

//! Represents a scoring function.
class scoring_function
{
public:
	static const size_t n = 15; //!< Number of XScore atom types.
	static const size_t np = n*(n+1)>>1; //!< Number of XScore atom type pairs.
	static const size_t ns = 1024; //!< Number of samples in a unit distance.
	static const size_t cutoff = 8; //!< Atom type pair distance cutoff.
	static const size_t nr = ns*cutoff*cutoff+1; //!< Number of samples within the entire cutoff.
	static const double cutoff_sqr; //!< Cutoff square.
	static const array<double, 5> weights; //!< Weight constants for 5 terms.

	//! Constructs an empty scoring function.
	explicit scoring_function();

	//! Returns the score between two atoms of XScore atom types t0 and t1 with distance r.
	static double score(const size_t t0, const size_t t1, const double r);

	//! Accumulates the unweighted score between two atoms of XScore atom types t0 and t1 with square distance r2.
	static void score(double* const v, const size_t t0, const size_t t1, const double r2);

	//! Precalculates the scoring function values of sample points for the type combination of t0 and t1.
	void precalculate(const size_t t0, const size_t t1);

	//! Clears precalculated values.
	void clear();

	vector<vector<double>> e; //!< Scoring function values.
	vector<vector<double>> d; //!< Scoring function derivatives divided by distance.

private:
	static const array<double, n> vdw; //!< Van der Waals distances for XScore atom types.
	vector<double> rs; //!< Distance samples.
};

#endif
