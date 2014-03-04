#include "molecule.h"
#include "assert.h"
#include <iostream>

SimpleAtom::SimpleAtom()
{
	aid = 0;
	number = 0;
}

SimpleAtom::SimpleAtom(unsigned short int _aid, unsigned int _number)
{
	aid = _aid;
	number = _number;
}

unsigned short int SimpleAtom::get_id()
{
	if (aid == 0)
		throw "atom sequence id cannot be 0";
	return aid;
}


unsigned int SimpleAtom::get_number()
{
	if (number == 0)
		throw "atom 'R' is not allowed";
	return number;
}

unsigned int SimpleAtom::GetAtomicNum()
{
	return get_number();
}

void SimpleAtom::register_bond(SimpleBond* b)
{
	if (aid == 0)
		throw "atom sequence id cannot be 0";
	if (number == 0)
		throw "atom 'R' is not allowed";
	if (b->get_nbr_atom(this) == NULL)
		throw "invalid bond";
	bonds.push_back(b);
}

std::vector<SimpleBond*>::iterator SimpleAtom::get_bonds_iter()
{
	return bonds.begin();
}

std::vector<SimpleBond*>::iterator SimpleAtom::get_bonds_iter_end()
{
	return bonds.end();
}

SimpleAtom::~SimpleAtom()
{
}

// ----------------------

SimpleBond::SimpleBond()
{
	atom1 = NULL;
	atom2 = NULL;
	bond_order = 0;
}

SimpleBond::SimpleBond(SimpleAtom* a1, SimpleAtom* a2, unsigned int order)
{
	atom1 = a1;
	atom2 = a2;
	bond_order = order;
}

unsigned int SimpleBond::get_bond_order()
{
	if (bond_order == 0)
		throw "bond order cannot be 0";
	return bond_order;
}

unsigned int SimpleBond::GetBondOrder()
{
	return get_bond_order();
}

SimpleAtom* SimpleBond::get_nbr_atom(SimpleAtom* a)
{
	if (bond_order == 0)
		throw "bond order cannot be 0";
	if (atom1->get_id() == a->get_id())
		return atom2;
	if (atom2->get_id() == a->get_id())
		return atom1;
	return NULL;
}

SimpleAtom* SimpleBond::GetNbrAtom(SimpleAtom* atom)
{
	return get_nbr_atom((SimpleAtom*) atom);
}

SimpleBond::~SimpleBond()
{
}

// --------------------------

SimpleMolecule::SimpleMolecule()
{
}

unsigned int SimpleMolecule::get_num_atoms()
{
	return atoms.size();
}

unsigned int SimpleMolecule::NumAtoms()
{
	return get_num_atoms();
}

void SimpleMolecule::add_atom(SimpleAtom& a)
{
	atoms[a.get_id()] = a;
}

int SimpleMolecule::add_bond(SimpleAtom& a1, SimpleAtom& a2, unsigned int bond_type)
{
	// check existence
	if (atoms.find(a1.get_id()) == atoms.end()) return 0;
	if (atoms.find(a2.get_id()) == atoms.end()) return 0;
	SimpleAtom& _a1 = atoms[a1.get_id()];
	SimpleAtom& _a2 = atoms[a2.get_id()];
	assert(_a1.get_number() == a1.get_number());
	assert(_a2.get_number() == a2.get_number());
	
	SimpleBond b(&_a1, &_a2, bond_type);
	int key;
	if (_a1.get_id() < _a2.get_id())
		key = (_a1.get_id() << 16) + _a2.get_id();
	else
		key = (_a2.get_id() << 16) + _a1.get_id();
	bonds[key] = b;
	
	_a1.register_bond(&(bonds[key]));
	_a2.register_bond(&(bonds[key]));
	return 1;
}

SimpleAtom* SimpleMolecule::get_atom(unsigned short int aid)
{
	if (atoms.find(aid) != atoms.end())
		return &(atoms[aid]);
	else
		return NULL;
}

SimpleAtom* SimpleMolecule::GetAtom(unsigned int id)
{
	return get_atom((unsigned short int) id);
}

SimpleBond* SimpleMolecule::get_bond(unsigned short int a1, unsigned short int a2)
{
	unsigned int key;
	if (a1 < a2)
		key = (a1 << 16) + a2;
	else
		key =  (a2 << 16) + a1;
	if (bonds.find(key) != bonds.end())
		return &(bonds[key]);
	else
		return NULL;
}

SimpleBond* SimpleMolecule::GetBond(unsigned int a_id1, unsigned int a_id2)
{
	return get_bond((unsigned short int) a_id1, (unsigned short int) a_id2);
}

SimpleMolecule::~SimpleMolecule()
{
}

