#ifndef MOLECULE_H_
#define MOLECULE_H_
#include <vector>
#include <map>
class SimpleBond;
class SimpleAtom
{
private:
	short unsigned int aid;
	unsigned int number;
	std::vector<SimpleBond*> bonds;
public:
	SimpleAtom();
	SimpleAtom(unsigned short int, unsigned int);
	unsigned short int get_id();
	unsigned int get_number();
	void register_bond(SimpleBond* b);
	std::vector<SimpleBond*>::iterator get_bonds_iter();
	std::vector<SimpleBond*>::iterator get_bonds_iter_end();
	virtual ~SimpleAtom();
	/* overriding methods */
	unsigned int GetAtomicNum();
};

class SimpleBond
{
private:
	SimpleAtom* atom1;
	SimpleAtom* atom2;
	unsigned int bond_order;
public:
	unsigned int get_bond_order();
	SimpleAtom* get_nbr_atom(SimpleAtom* atom);
	SimpleBond();
	SimpleBond(SimpleAtom*, SimpleAtom*, unsigned int);
	virtual ~SimpleBond();
	/* overriding methods */
	SimpleAtom* GetNbrAtom(SimpleAtom* atom);
	unsigned int GetBondOrder();
};

class SimpleMolecule
{
private:
	std::map<unsigned short int, SimpleAtom> atoms;
	std::map<unsigned int, SimpleBond> bonds;
public:
	unsigned int get_num_atoms();
	SimpleAtom* get_atom(unsigned short int id);
	SimpleBond* get_bond(unsigned short int a_id1, unsigned short int a_id2);
	SimpleMolecule();
	void add_atom(SimpleAtom&);
	int add_bond(SimpleAtom&, SimpleAtom&, unsigned int);
	virtual ~SimpleMolecule();
	/* overriding methods */
	unsigned int NumAtoms();
	SimpleAtom* GetAtom(unsigned int id);
	SimpleBond* GetBond(unsigned int a_id1, unsigned int a_id2);
};

#endif /*MOLECULE_H_*/
