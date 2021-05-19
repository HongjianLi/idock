#pragma once
#ifndef IDOCK_PKA_HPP
#define IDOCK_PKA_HPP

#include <string>
#include <vector>
#include <filesystem>
using namespace std;
using namespace std::filesystem;

//! Represents a pKa line of an atom.
class pka_line
{
public:
	string residue_name;
	string atom_name;
	char chain_id;
	double pka_value;

	//! Creates from a summary line in pka file.
	pka_line(const string& line);
};

//! Represents pKa values of a molecule.
class pka
{
private:
	vector<pka_line> pka_lines;

public:
	//! Creates a dummy object where no pKa lines are defined.
	pka() {}

	//! Creates from a pka file.
	pka(const path& p);

	//! Returns a pKa value for the given query criteria.
	double get_pka(const string& atom_name, const string& residue_name, char chain_id) const;

	//! Returns whether the query atom should be inoized in the given pH environment.
	bool is_ionized(const string& atom_name, const string& residue_name, char chain_id, double ph) const;
};
#endif