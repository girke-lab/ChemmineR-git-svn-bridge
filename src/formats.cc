#include "desc.h"
#include <string>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <sstream>

#define MAX_LINE_LENGTH 100000

#ifndef HAS_OPENBABEL
static char elements[112][3] = {
	    "R",
	    "H",
	    "He",
	    "Li",
	    "Be",
	    "B",
	    "C",
	    "N",
	    "O",
	    "F",
	    "Ne",
	    "Na",
	    "Mg",
	    "Al",
	    "Si",
	    "P",
	    "S",
	    "Cl",
	    "Ar",
	    "K",
	    "Ca",
	    "Sc",
	    "Ti",
	    "V",
	    "Cr",
	    "Mn",
	    "Fe",
	    "Co",
	    "Ni",
	    "Cu",
	    "Zn",
	    "Ga",
	    "Ge",
	    "As",
	    "Se",
	    "Br",
	    "Kr",
	    "Rb",
	    "Sr",
	    "Y",
	    "Zr",
	    "Nb",
	    "Mo",
	    "Tc",
	    "Ru",
	    "Rh",
	    "Pd",
	    "Ag",
	    "Cd",
	    "In",
	    "Sn",
	    "Sb",
	    "Te",
	    "I",
	    "Xe",
	    "Cs",
	    "Ba",
	    "La",
	    "Ce",
	    "Pr",
	    "Nd",
	    "Pm",
	    "Sm",
	    "Eu",
	    "Gd",
	    "Tb",
	    "Dy",
	    "Ho",
	    "Er",
	    "Tm",
	    "Yb",
	    "Lu",
	    "Hf",
	    "Ta",
	    "W",
	    "Re",
	    "Os",
	    "Ir",
	    "Pt",
	    "Au",
	    "Hg",
	    "Tl",
	    "Pb",
	    "Bi",
	    "Po",
	    "At",
	    "Rn",
	    "Fr",
	    "Ra",
	    "Ac",
	    "Th",
	    "Pa",
	    "U",
	    "Np",
	    "Pu",
	    "Am",
	    "Cm",
	    "Bk",
	    "Cf",
	    "Es",
	    "Fm",
	    "Md",
	    "No",
	    "Lr",
	    "Rf",
	    "Db",
	    "Sg",
	    "Bh",
	    "Hs",
	    "Mt",
	    "Ds",
	    "Rg"	
};


void parse_line_4(const char* buf, int& n_atoms, int& n_bonds)
{
	char num[4];
	strncpy(num, buf, 3);
	num[3] = '\0';
	n_atoms = atoi(num);
	DEBUG_VAR(n_atoms);
	strncpy(num, buf + 3, 3);
	num[3] = '\0';
	n_bonds = atoi(num);
	DEBUG_VAR(n_bonds);
}

int parse_atoms(const char* buf, Molecule& mol, int line)
{
	char ele[4] = {0x0, 0x0, 0x0, 0x0};
	int i;
	for (i = 31; i < 34; i ++) {
		if (buf[i] != ' ') break;
	}
	int j = 0;
	while (i < 34) {
		if (buf[i] != ' ') ele[j++] = buf[i];
		i ++;
	}
	for (i = 0; i < 112; i ++) {
		if (strcmp(ele, elements[i]) == 0) {
			Atom a(line - 4, i);
			mol.add_atom(a);
			DEBUG_VAR(i);
			break;
		}
	}
	if (i == 112) {
		std::cerr << "Cannot understand atom type " << ele << " on line " << line << std::endl;
		return 0;
	}
	return 1;
}

int parse_bonds(const char* buf, Molecule& mol, int line)
{
	int left, right, bond_type;
	char num[4];
	strncpy(num, buf, 3);
	num[3] = '\0';
	left = atoi(num);
	DEBUG_VAR(left);
	strncpy(num, buf + 3, 3);
	num[3] = '\0';
	right = atoi(num);
	DEBUG_VAR(right);
	strncpy(num, buf + 6, 3);
	num[3] = '\0';
	bond_type = atoi(num);
	DEBUG_VAR(bond_type);
	
	if (!(left != 0 && right != 0 && bond_type != 0)) 
		throw "invalid bond line";
	
	Atom* a1 = mol.GetAtom(left);
	Atom* a2 = mol.GetAtom(right);
	if (a1 == NULL) {
		std::cerr << "Bond definition contains unknown atom : " << left << " on line " << line << std::endl;
		return 0;
	}
	if (a2 == NULL) {
		std::cerr << "Bond definition contains unknown atom : " << right << " on line " << line << std::endl;
		return 0;
	}
	return mol.add_bond(*a1, *a2, bond_type);
}

void parse_sdf(std::istream & ifs, Molecule ** mol)
{
	int line_cntr = 0;
	if (ifs.good()) {
		char buf[MAX_LINE_LENGTH+2];
		int n_atoms = 0;
		int n_bonds = 0;
		while (true) {
			ifs.getline(buf, MAX_LINE_LENGTH+2);
			line_cntr ++;
			if (ifs.fail() or strlen(buf) > MAX_LINE_LENGTH) {
				if (strlen(buf) > MAX_LINE_LENGTH) 
					std::cerr << "SDF not well-formatted : line exceeds "<<MAX_LINE_LENGTH<<" characters" << " len=" << strlen(buf) << " last=" << buf[strlen(buf) - 1] << std::endl;
				else
					std::cerr << "SDF not well-formatted : error when reading line " << line_cntr << std::endl;
				delete (*mol);
				*mol = NULL;
				break;				
			}
			if (line_cntr <= 3) continue;
			if (line_cntr == 4) {
				parse_line_4(buf, n_atoms, n_bonds);
				if (n_atoms == 0 || n_bonds == 0) {
					std::cerr << "SDF not well-formatted : failed when reading number of atoms and number of bonds on line " << line_cntr << std::endl;
					std::cerr << " line is: " << buf << std::endl;
					delete (*mol);
					*mol = NULL;
					break;				
				}
			}
			else if (line_cntr <= n_atoms + 4) {
				if (parse_atoms(buf, **mol, line_cntr) == 0) {
					std::cerr << "SDF not well-formatted: bond contains unknown atoms on line " << line_cntr << std::endl;
					std::cerr << " line is: " << buf << std::endl;
					delete (*mol);
					*mol = NULL;
					break;				
				}
			}
			else if (line_cntr <= n_atoms + 4 + n_bonds) {
				if (parse_bonds(buf, **mol, line_cntr) == 0) {
					std::cerr << "SDF not well-formatted: bond contains unknown atoms on line " << line_cntr << std::endl;
					std::cerr << " line is: " << buf << std::endl;
					delete (*mol);
					*mol = NULL;
					break;				
				}
			}
			else break;
		}
	}
}

#endif

Molecule * 
new_mol_from_sdfile(const char* sdfile)
{
	Molecule *mol = new Molecule;
#ifdef HAS_OPENBABEL
	mol->Clear();
	OpenBabel::OBConversion conv;
    conv.SetInAndOutFormats("SDF", "SDF");
	conv.ReadFile(mol, std::string(sdfile));
#else
	std::ifstream ifs(sdfile, std::ifstream::in);
	try {
		parse_sdf(ifs, &mol);
	} catch (char const* err) {
		std::cerr << err << std::endl;
		delete mol;
		return NULL;
	}
	ifs.close();
#endif
	return mol;
}

Molecule*
new_mol_from_sdf(const char* sdf)
{
	Molecule *mol = new Molecule;
#ifdef HAS_OPENBABEL
	mol->Clear();
	OpenBabel::OBConversion conv;
    conv.SetInAndOutFormats("SDF", "SDF");
	conv.ReadString(mol, std::string(sdf));
#else
	std::string _sdf(sdf);
	std::istringstream ifs(_sdf);
	try {
		parse_sdf(ifs, &mol);
	} catch (char const* err) {
		std::cerr << err << std::endl;
		delete mol;
		return NULL;
	}
#endif
	return mol;	
}


Molecule *
new_mol_from_smiles(const char* smiles)
{
#ifdef HAS_OPENBABEL
	Molecule *mol = new Molecule;
	mol->Clear();
	OpenBabel::OBConversion conv;
	conv.SetInAndOutFormats("SMILES", "SMILES");
	conv.ReadString(mol, std::string(smiles));
	return mol;
#else
	std::cerr << "`descriptor' tool is not compiled with OpenBabel and cannot read SMILES format. Return NULL molecule." << std::endl;
	return NULL;
#endif
}

int sdf_iter(std::fstream& ifs, std::string& sdf, int& line_cntr)
{
	sdf.clear();
	char line[MAX_LINE_LENGTH+2]; line[MAX_LINE_LENGTH+1] = '\0';
	char buf_4[5]; buf_4[4] = '\0';

	ifs.getline(line, MAX_LINE_LENGTH+2);
	line_cntr ++;

	while (ifs.good() or ifs.fail()) {
		if (strlen(line) > MAX_LINE_LENGTH) {
			std::cerr << "Line exceeds "<<MAX_LINE_LENGTH<<" characters when reading line "
				<< line_cntr << std::endl;
			sdf.clear();
			return 0;
		}
		if (ifs.fail()) return 0;
		sdf += line;
		sdf += '\n';
		strncpy(buf_4, line, 4);
		if (strcmp(buf_4, "$$$$") == 0) {
			return 1;
		}

		ifs.getline(line, MAX_LINE_LENGTH+2);
		line_cntr ++;
	}
	return 0;
}
