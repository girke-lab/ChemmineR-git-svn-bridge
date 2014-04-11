#include "desc.h"
#include <set>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <vector>
#include <string.h>

/*
double similarity(Molecule* mol1, Molecule* mol2)
{
	 if (mol1 == NULL || mol2 == NULL) {
		  std::cerr << "Similarity calculation function receives an empty molecule." << std::endl;
		  return -1;
	  }
	  
	  std::multiset<unsigned int> descs1;
	  calc_desc(*mol1, descs1);
	  
	  return similarity(descs1, mol2);
}
double similarity(std::multiset<unsigned int> & descs1, Molecule* mol2)
{
	 if (mol2 == NULL) {
		  std::cerr << "Similarity calculation function receives an empty molecule." << std::endl;
		  return -1;
	  }
	  
	  std::multiset<unsigned int> descs2;
	  calc_desc(*mol2, descs2);
	  return similarity(descs1, descs2);
}

double similarity(std::multiset<unsigned int> & descs1, std::multiset<unsigned int> & descs2)
{
	  std::multiset<unsigned int> intersect, union_set;
	  std::set_intersection(descs1.begin(), descs1.end(), descs2.begin(), descs2.end(), std::inserter(intersect, intersect.begin()));
	  std::set_union(descs1.begin(), descs1.end(), descs2.begin(), descs2.end(), std::inserter(union_set, union_set.begin()));
	  double sim = intersect.size() * 1.0 / union_set.size();
	  return sim;
}
*/

double similarity(std::vector<unsigned int> & descs1, std::vector<unsigned int> & descs2, int sorted)
{
	unsigned int i = 0; unsigned int j = 0;
	unsigned int set_union = 0; unsigned int set_intersection = 0;
	if (! sorted) {
		sort(descs1.begin(), descs1.end());
		sort(descs2.begin(), descs2.end());
	}
	while (i < descs1.size() && j < descs2.size()) {
		if (descs1[i] == descs2[j]) {
			set_intersection ++;
			i ++;
			j ++;
		} else if (descs1[i] < descs2[j]) {
			i ++;
		} else {
			j ++;
		}
		set_union ++;
	}
	if (i < descs1.size()) set_union += (descs1.size() - i);
	else set_union += (descs2.size() - j);
	return set_intersection * 1.0 / set_union;
}
