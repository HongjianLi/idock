#include "pka.hpp"
#include <fstream>
#include "string.hpp"

//! Creates from a summary line in pka file.
pka_line::pka_line(const string& line)
	: residue_name(trim(line.substr(3, 3)))
	, atom_name(trim(line.substr(7, 3)))
	, chain_id(line[11])
	, pka_value(stod(line.substr(13, 8)))
{
	// line sample: "   2GM  N6 A    10.00      10.00                N33   "
	//               012345678901234567890123456789012345678901234567890123
}

//! Creates from a pka file.
pka::pka(const path& p)
{
	string line;
	bool found = false;
	for (ifstream ifs(p); safe_getline(ifs, line);)
	{
		if (!found)
		{
			if (starts_with(line, "SUMMARY OF THIS PREDICTION"))
				found = true;
		}
		else
		{
			if (starts_with(line, "--"))
				break;
			if (line.size() == 54)
			{
				pka_lines.push_back(pka_line(line));
			}
		}
	}
}

//! Returns a pKa value for the given query criteria.
double pka::get_pka(const string& atom_name, const string& residue_name, char chain_id) const
{
	for (auto& line : pka_lines)
	{
		if (line.atom_name == atom_name && line.residue_name == residue_name && line.chain_id == chain_id)
			return line.pka_value;
	}
	return 0;
}

//! Returns whether the query atom should be inoized in the given pH environment.
bool pka::is_ionized(const string& atom_name, const string& residue_name, char chain_id, double ph) const
{
	return ph < get_pka(atom_name, residue_name, chain_id);
}
