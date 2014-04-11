#include "script.h"
#include "desc.h"
#include <iostream>
#include <algorithm>
#include <assert.h>
#include <string.h>


Descriptors::Descriptors()
{
}

Descriptors::~Descriptors()
{
	//std::cout << "descriptors object released" << std::endl;
}

int Descriptors::parse_sdf(const char* sdf)
{
	Molecule* mol = new_mol_from_sdf(sdf);
	descs.clear();
	if (mol == NULL) return 0;
	int ret = calc_desc(*mol, descs);
	delete mol;
	return ret;
}

int Descriptors::parse_sdfile(const char* sdfile)
{
	Molecule* mol = new_mol_from_sdfile(sdfile);
	descs.clear();
	if (mol == NULL) return 0;
	int ret = calc_desc(*mol, descs);
	delete mol;
	return ret;
}

int Descriptors::parse_smiles(const char* smiles)
{
	Molecule* mol = new_mol_from_smiles(smiles);
	descs.clear();
	if (mol == NULL) return 0;
	int ret = calc_desc(*mol, descs);
	delete mol;
	return ret;
}

unsigned int Descriptors::get_descriptor(unsigned int i)
{
	if (i >= descs.size()) return 0;
	return descs[i];
}

unsigned int Descriptors::get_len()
{
	return descs.size();
}

double similarity(Descriptors* d1, Descriptors* d2)
{
	if (d1 == NULL || d2 == NULL) {
		std::cerr << "one or both input compounds are invalid" << std::endl;
		return 0;
	}
	
	return similarity(d1->descs, d2->descs, 1);
}


