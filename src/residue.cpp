#include "residue.hpp"
#include "string.hpp"

const set<string> residue::amino_acids
{{
	"ALA", // A, Alanine
	"ARG", // R, Arginine
	"ASN", // N, Asparagine
	"ASP", // D, Aspartic acid
	"CYS", // C, Cysteine
	"GLU", // E, Glutamic acid
	"GLN", // Q, Glutamine
	"GLY", // G, Glycine
	"HIS", // H, Histidine
	"ILE", // I, Isoleucine
	"LEU", // L, Leucine
	"LYS", // K, Lysine
	"MET", // M, Methionine
	"PHE", // F, Phenylalanine
	"PRO", // P, Proline
	"SER", // S, Serine
	"THR", // T, Threonine
	"TRP", // W, Tryptophan
	"TYR", // Y, Tyrosine
	"VAL", // V, Valine
}};


//! Constructs a residue from an ATOM/HETATM line in PDBQT format.
residue::residue(const string& line)
	: name(trim(line.substr(17, 3)))
	, chain(line[21])
	, seq(stoi(line.substr(22, 4)))
{
}

//! Determine if the residue is an amino acid.
bool residue::is_amino_acid()
{
	return amino_acids.count(name);
}
