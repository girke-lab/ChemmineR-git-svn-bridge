#ifndef __DESC_H
#define __DESC_H
#ifdef HAS_OPENBABEL
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/obiter.h>
#include <openbabel/oberror.h>
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
#endif

