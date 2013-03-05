#include <iostream>
#include <fstream>
#include <set>
#include <algorithm>
#include <iterator>
#include "desc.h"

static void check_bonds(Atom * atom, char & degree, char & pi_electrons)
{
  degree = 0;
  pi_electrons = 0;
  _FOR_BONDS_OF_ATOM(b, (*atom)) {
    if (b->GetNbrAtom(atom)->GetAtomicNum() != 1) {
      char bond_order = b->GetBondOrder();
			if (bond_order >=4)
				std::cerr << "pi bond will be ignored" << std::endl;
			else
				pi_electrons += (bond_order - 1);
      degree ++;
    } else {
      //std::cerr << "ERROR unexpected hydrogen. Ignored.\n";
      continue;
    }
  }
}

static unsigned int calc_atom_desc(char atom, char degree, char pi_electrons)
{
  unsigned int atom_desc = ((unsigned int) atom << 6) & 0x00001FC0;
  if (degree > 7) degree = 7;
  if (pi_electrons > 7) pi_electrons = 7;

  atom_desc = atom_desc | (((unsigned int) degree << 3) & 0x00000038);
  atom_desc = atom_desc | ((unsigned int) pi_electrons & 0x00000007);

  return atom_desc;
}

static int allpairshort (int **ct, int n, int ***ret)
{
  int **scrtch;
  scrtch = new int *[n];
  for (int i = 0; i < n; i++) {
    scrtch[i] = new int[n];
  }
  // init
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (ct[i][j] == 0) {
        scrtch[i][j] = 256;
      } else
        scrtch[i][j] = ct[i][j];
    }
  }

  // use Floyd–Warshall algorithm
  int i, j, k;
  for (k = 0; k < n; k++) {
    for (i = 0; i < n; i++) {
			if (i == k) continue;
      for (j = 0; j < i; j++) {
        if (scrtch[i][j] > scrtch[i][k] + scrtch[k][j]) {
          scrtch[i][j] = scrtch[i][k] + scrtch[k][j];
					scrtch[j][i] = scrtch[i][j];
        } 
			}
    }
  }

  *ret = scrtch;

  return 1;
}

static unsigned int _core(Molecule & mol, std::vector<unsigned int> & descriptors)
{
  // an array of pointers to OBAtom objects
  Atom ** atoms;
  int num_atoms = (int) mol.NumAtoms ();
  atoms = new Atom *[num_atoms];
  for (int c = 1; c <= num_atoms; c++) {
    atoms[c - 1] = mol.GetAtom (c);
  }


  int **connection_tbl = new int *[num_atoms];
  for (int i = 0; i < num_atoms; i ++)
    connection_tbl[i] = new int[num_atoms];
  // read the bonds to a connection table
  for (int i = 0; i < num_atoms; i++) {
    for (int j = i; j < num_atoms; j++) {
      if (mol.GetBond (i + 1, j + 1) == NULL)
    	connection_tbl[j][i] = connection_tbl[i][j] = 0;
      else
        connection_tbl[j][i] = connection_tbl[i][j] = 1;
    }
  }
  // calculate all-pair shortest distances
  int **shortest;
  allpairshort(connection_tbl, num_atoms, &shortest);

  char atom, degree, pi_electrons;

  for (int i = 0; i < num_atoms; i++) {
    for (int j = i + 1; j < num_atoms; j++) {
      int dist = shortest[i][j];
      if (dist >= 128) continue;

      atom = (char) atoms[i]->GetAtomicNum();
			if (atom == 1) continue;
      check_bonds(atoms[i], degree, pi_electrons);
      unsigned int atom_desc1 = calc_atom_desc(atom, degree, pi_electrons);

      atom = (char) atoms[j]->GetAtomicNum ();
			if (atom == 1) continue;
      check_bonds(atoms[j], degree, pi_electrons);
      unsigned int atom_desc2 = calc_atom_desc(atom, degree, pi_electrons);

      unsigned int aa_dist = shortest[i][j] & 0x3F;

      unsigned int desc = 0;
      /** note: if both atoms's atomic numbers are beyond 64, then descriptors will lose the highest bit
       * because the smaller atom only has 12 instead of 13 bits */
      if (atom_desc1 < atom_desc2) 
        desc = desc | (atom_desc1 << 20) | (aa_dist << 13) | atom_desc2;
      else
        desc = desc | (atom_desc2 << 20) | (aa_dist << 13) | atom_desc1;
      descriptors.push_back(desc);
    }
  }

  for (int i = 0; i < num_atoms; i++) {
    delete[]connection_tbl[i];
    delete[]shortest[i];
  }
  delete[]connection_tbl;
  delete[]shortest;
  delete[]atoms;
  return 1;
}

unsigned int calc_desc(Molecule & mol, std::vector<unsigned int> & descriptors)
{
	int ret = _core(mol, descriptors);
	sort(descriptors.begin(), descriptors.end());
	return ret;
}
unsigned int calc_desc(Molecule & mol, std::multiset<unsigned int> & descriptors)
{
	std::vector<unsigned int> descs;
	int ret = calc_desc(mol, descs);
	std::copy(descs.begin(), descs.end(), std::inserter(descriptors, descriptors.begin()));
	return ret;
}
