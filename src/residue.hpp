#pragma once
#ifndef IDOCK_RESIDUE_HPP
#define IDOCK_RESIDUE_HPP

#include <vector>
#include <string>
#include <set>
using namespace std;

//! Represents a residue formed by a few atoms.
class residue
{
public:
	//! Constructs a residue from an ATOM/HETATM line in PDBQT format.
	explicit residue(const string& line);

	string name; //!< Residue name.
	char chain; //!< Chain identifier.
	int seq; //!< Residue sequence number. Could be negative.

	//! Determine if the residue is an amino acid.
	bool is_amino_acid();

private:
	static const set<string> amino_acids;
};

#endif
