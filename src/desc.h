#ifndef __DESC_H
#define __DESC_H
#ifdef HAS_OPENBABEL
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/obiter.h>
typedef OpenBabel::OBAtom Atom; 
typedef OpenBabel::OBBond Bond;
typedef OpenBabel::OBMol Molecule;
#define _FOR_BONDS_OF_ATOM(b,p) for( OpenBabel::OBAtomBondIter b(p); b; ++b )
#else
#include "molecule.h"
typedef SimpleAtom Atom;
typedef SimpleBond Bond;
typedef SimpleMolecule Molecule;
#define _FOR_BONDS_OF_ATOM(b,p) Bond *b; for ( std::vector<SimpleBond*>::iterator ii = p.get_bonds_iter(); ((b=*ii) != NULL) && (ii != p.get_bonds_iter_end()); ++ii )
#endif
#include <iostream>
#include <set>
#include <vector>
#include "debug.h"

unsigned int calc_desc(Molecule & mol, std::multiset<unsigned int> & descriptors);
unsigned int calc_desc(Molecule & mol, std::vector<unsigned int> & descriptors);

/* formats */
Molecule * new_mol_from_sdf(const char* sdf); /* remember to delete the memory after use */
Molecule * new_mol_from_sdfile(const char* sdfile); /* remember to delete the memory after use */
Molecule * new_mol_from_smiles(const char* smile); /* remember to delete the memory after use */
int sdf_iter(std::fstream& ifs, std::string& sdf, int& line_cntr);

/* similarity */
double similarity(Molecule* mol1, Molecule* mol2);
double similarity(std::multiset<unsigned int> & descriptors, Molecule* mol2);
double similarity(std::multiset<unsigned int> & descriptors, std::multiset<unsigned int> & descriptors2);
double similarity(std::vector<unsigned int> & descs1, std::vector<unsigned int> & descs2, int sorted=1);


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


#endif

