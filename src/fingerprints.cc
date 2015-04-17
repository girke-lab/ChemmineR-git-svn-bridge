#include <stdlib.h>
#include <math.h>
#include <vector>
#include <algorithm>
/*
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
*/

#include <Rcpp.h>


//extern "C"{

using namespace Rcpp;
using namespace std;

struct IndexedValue{
	int index;
	long long int value;
	int dupCount;
	IndexedValue(){}
	IndexedValue(int i, int v): index(i), value(v),dupCount(0) {}
};


bool byValue (IndexedValue *a, IndexedValue *b){
	return a->value < b->value;
}
bool byIndex (IndexedValue *a, IndexedValue *b){
	return a->value < b->value;
}
RcppExport SEXP uniquifyAtomPairs(SEXP atomPairsS){
	
	NumericVector atomPairs(atomPairsS);
	vector<IndexedValue*> aps(atomPairs.size());

	for(int i=0; i < aps.size(); i++)
		aps[i] = new IndexedValue(i,atomPairs[i]);
	
	stable_sort(aps.begin(),aps.end(),byValue);

	long long int lastValue = -1;
	int dupCount=0;
	for(int i=0; i < aps.size(); i++){

		if(lastValue == aps[i]->value) {//found dup
			dupCount++;
		}else{
			dupCount=0;
		}
		aps[i]->dupCount = dupCount;

		lastValue = aps[i]->value;

	}

	for(int i=0; i < aps.size(); i++){
		atomPairs(aps[i]->index) = (aps[i]->value << 7) + aps[i]->dupCount;
		delete aps[i];
	}

	return atomPairs;
}


//}
