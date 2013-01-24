/*
 * Patrick-Jarvis Clustering
 * Input: a neighbor file, each line is of this format:
 *        entry_id,neighbor_1,neighbor_2,... neighbor_p
 * Output: to STDOUT
 *        cluster_id : entry_id1, entry_id 2, ....
 *        cluster_id : entry_id, ...
 */
#include "debug.h"
#include "DisjointSets.h"
#include <map>
#include <vector>
#include <fstream>
#include <iostream>
#include <cstring>
#include <cassert>
#include <algorithm>

#ifdef NO_MAIN
#include <R.h>
#include <Rinternals.h>

extern "C" {
	SEXP jarvis_patrick(SEXP neighbors,SEXP minNbrs,
		SEXP fast,SEXP bothDirections);
}
#endif

void print_clusters(DisjointSets &s);
void print_clusters(DisjointSets &s,int N);

#define LINE_BUF_SIZE 10240

//This is not thread safe
std::map<std::string, int> name_to_id;
std::vector<std::string> names;
std::vector<std::vector<int> > nbr_list;

/* read in a neighbor file */
void static prepare_neighbors(const char* nbr_file, int skip, int p)
{
	// open file
	std::fstream f(nbr_file, std::ios::in);
	if (not f.good()) {
		std::cerr << "**Error in opening file " << nbr_file << std::endl;
		exit(1);
	}

	// read line by line
	char line_buf[LINE_BUF_SIZE];
	for (int i = 0; i < skip; i ++) {
		f.getline(line_buf, LINE_BUF_SIZE);
	}
	int line_cntr = 0;
	names.clear();
	name_to_id.clear();
	while (f.good()) {
		f.getline(line_buf, LINE_BUF_SIZE);
		// parse the first before comma
		char* raw = strtok(line_buf, ",");
		if (raw == NULL) continue;
		DEBUG_VAR(line_cntr);
		DEBUG_VAR(raw);
		names.push_back(std::string(raw));
		name_to_id[std::string(raw)] = line_cntr;
		line_cntr ++;
	}
	
	assert(f.eof());
	f.clear();
	f.seekg(0, std::ios::beg);
	for (int i = 0; i < skip; i ++) {
		f.getline(line_buf, LINE_BUF_SIZE);
	}

	nbr_list.clear();
	line_cntr = 0;
	while (f.good()) {
		f.getline(line_buf, LINE_BUF_SIZE);
		if (f.eof()) break;
		// parse the first before comma
		char* raw = strtok(line_buf, ",");
		if (raw == NULL) break;
		// now parse the rest
		std::vector<int> nbrs;
		while (true) {
			char* raw = strtok(NULL, ",");
			if (raw == NULL) break;
			DEBUG_VAR(raw);
			int nbr = name_to_id[std::string(raw)];
			if (nbr != line_cntr) nbrs.push_back(nbr);
			DEBUG_VAR(nbr);
		}
		// apply p
		int to_erase = nbrs.size() - p;
		if (to_erase > 0)
			nbrs.erase(nbrs.begin(), nbrs.begin() + to_erase);
		// sort nbrs
		std::sort(nbrs.begin(), nbrs.end());
		nbr_list.push_back(nbrs);
		line_cntr ++;
	}
}

void print_neighbors()
{
	for (int i = 0; i < nbr_list.size(); i ++) {
		std::cout << i;
		for (int j = 0; j < nbr_list[i].size(); j ++) {
			std::cout << " " << nbr_list[i][j];
		}
		std::cout << std::endl;
	}
}

/* size of neighbor set intersection. Expect a sorted list */
int nbr_intersect(std::vector<int>& nbrs1, std::vector<int>& nbrs2)
{
	int i = 0, j = 0, intrsct = 0;
	while (i < nbrs1.size() and j < nbrs2.size()) {
		if (nbrs1[i] == nbrs2[j]) {
			intrsct ++; 
			i ++; j ++;
		} else if (nbrs1[i] > nbrs2[j])
			j ++;
		else
			i ++;
	}
	return intrsct;
}
int contains(int x, std::vector<int> &list)
{
	for(int i=0; i < list.size(); i++)
		if(x==list[i])
			return 1;
	return 0;
}

void checkPair(DisjointSets &s,int i, int j,int m)
{
	//printf("%d:%d\n",i,j);
	if (s.FindSet(i) == s.FindSet(j)) return;

	//printf("different sets\n");
	// check condition 2
	if (nbr_intersect(nbr_list[i], nbr_list[j]) < m)
		return;

	//printf("met intersection req\n");
	// merging clusters
	//printf("merged %d and %d\n",s.FindSet(i), s.FindSet(j));
	s.Union(s.FindSet(i), s.FindSet(j));
}

DisjointSets clusterAllPairs(int n,int m)
{
	DisjointSets s;
	s.AddElements(n);
	for (int i = 0; i < n; i ++) {
		for (int j = i+1; j < n; j ++) { 
			checkPair(s,i,j,m);
		}
		//print_clusters(s,n);
	}
	return s;
}

DisjointSets cluster(int n,int m,int bothDirections)
{
	DisjointSets s;
	s.AddElements(n);
	for (int i = 0; i < n; i ++) {
		for (int j = 0; j < nbr_list[i].size(); j ++) {
			if(!bothDirections || contains(i,nbr_list[j]))
				checkPair(s,i,nbr_list[i][j],m);
		}
	//	print_clusters(s,n);
	}
	return s;
}
DisjointSets cluster(int m)
{
	return cluster(names.size(),m,0);
}

#ifdef NO_MAIN
SEXP jarvis_patrick(SEXP neighbors,SEXP minNbrs,
		SEXP fast,SEXP bothDirections)
{
	// neightbors is NxKx2  last 2 are (id,distance)
	
	//load nbr_list with data from neighbors
   SEXP dims= getAttrib(neighbors,R_DimSymbol);
   int N = INTEGER(dims)[0]; // num compounds
   int K = INTEGER(dims)[1]; // num neighbors given

	//Rprintf("N:%d, K:%d,m:%d\n",N,K,*INTEGER(minNbrs));

	nbr_list.clear();
	for(unsigned i=0; i<N; i++) //rows
	{
		std::vector<int> nbrs;
		for(int j=0; j<K; j++)  //cols
		{  // R arrays are column major 
			int n=INTEGER(neighbors)[j*N+i];
			if( n!= NA_INTEGER && n != -1)
				nbrs.push_back(INTEGER(neighbors)[j*N+i]-1);
		}
		std::sort(nbrs.begin(), nbrs.end());
		nbr_list.push_back(nbrs);
	}
	//Rprintf("loaded nbr_list\n");

	// do actual clustering
	DisjointSets s= *INTEGER(fast)? 
		cluster(N,*INTEGER(minNbrs),*INTEGER(bothDirections)):
		clusterAllPairs(N,*INTEGER(minNbrs));

	//Rprintf("done clustering\n");
	//print_clusters(s,N);

	//pull result out of s
	SEXP result;
   PROTECT(result = allocVector(INTSXP,N));
	//Rprintf("allocated vector\n");


	for(unsigned i=0; i<N; i++){
		//Rprintf("cluster(%u)=%d\n",i,s.FindSet(i));
		INTEGER(result)[i]=s.FindSet(i)+1;
	}
	//Rprintf("done\n");

	UNPROTECT(1);

	return result;
}
#endif

void print_clusters(DisjointSets &s)
{
	print_clusters(s,names.size());
}
void print_clusters(DisjointSets &s,int N)
{
	for (int i = 0; i < N; i ++) 
		std::cout << i<<","<<s.FindSet(i) <<"  ";
	std::cout << std::endl;	
}

void print_cluster_stat(DisjointSets &s,int print_pair = 0)
{
	std::map<int, std::vector<int> > lg_cls;
	int *st = new int[names.size()];
	bzero(st, sizeof(int) * names.size());
	assert(st[0] == 0 and st[1] == 0);
	for (int i = 0; i < names.size(); i ++) {
		int set = s.FindSet(i);
		if (st[set] == 0) {
			st[set] = i + 1;
		} else if (st[set] > 0) {
			std::vector<int> _t;
			_t.push_back(st[set] - 1);
			_t.push_back(i);
			st[set] = -1;
			lg_cls[set] = _t;
		} else {
			assert(st[set] == -1);
			lg_cls[set].push_back(i);
		}
	}

	if (print_pair == 0) {
		for (std::map<int, std::vector<int> >::iterator i = lg_cls.begin();
				i != lg_cls.end();
				i ++) {
			std::cout << "cluster " << (*i).first << ":";
			for (int j = 0; j < (*i).second.size(); j ++)
				std::cout << " " << names[(*i).second[j]];
			std::cout << std::endl;
		}
	} else {
		for (std::map<int, std::vector<int> >::iterator i = lg_cls.begin();
				i != lg_cls.end();
				i ++) {
			for (int j = 0; j < (*i).second.size(); j ++) {
				for (int k = j + 1; k < (*i).second.size(); k ++) {
					long _i = (*i).second[j] * names.size() + (*i).second[k];
					long _j = (*i).second[k] * names.size() + (*i).second[j];
					if (_i > _j) 
						std::cout << " " << _j;
					else
						std::cout << " " << _i;
				}
			}
		}
		std::cout << std::endl;
	}
}

/*
 *
 * Modes:
 * 1: print non-singleton clusters and their members
 * 2: for each point, print its cluster
 * other: print pairs, encoded in a single long int
 */
#ifndef NO_MAIN
int main(int argc, char* argv[])
{
	if (argc < 5) {
		std::cerr << "Usage: " << argv[0] << " input.nbr skip_line p m [mode]" << std::endl;
		return 1;
	}
	int mode = 0;
	if (argc == 6)
		mode = atoi(argv[5]);
	prepare_neighbors(argv[1], atoi(argv[2]), atoi(argv[3]));
	//print_neighbors();
	DisjointSet s =cluster(atoi(argv[4]));
	if (mode == 2)
		print_clusters(s);
	else if (mode == 1)
		print_cluster_stat(s);
	else
		print_cluster_stat(s,1);
	return 0;
}
#endif
