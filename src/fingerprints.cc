
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

#define  TANIMOTO 0
#define  EUCLIDIAN 1
#define  TVERSKY 2
#define  DICE 3

#define M00 0 //neither bit set
#define M01 1 //bit only set in target
#define M10 2 //bit only set in query
#define M11 3 //both bits set

int* features(NumericVector &query, NumericMatrix &targets, int targetRow);


RcppExport SEXP similarity(SEXP queryS, SEXP targetsS,SEXP typeS,SEXP addoneS,SEXP alphaS,SEXP betaS){

	NumericVector query(queryS);

//	printf("query:  ");
//	for(int k=0;k < query.size() && k < 30;k++)
//		printf("%d,",(int)query(k));
//	printf("\n");


	NumericMatrix targets(targetsS);
	NumericVector result(targets.nrow());
	NumericVector t(typeS), t1(addoneS),t2(alphaS),t3(betaS);
	int type = t(0);
	int addone = t1(0);
	int alpha= t2(0);
	int beta = t3(0);
		

	for(int i=0; i < targets.nrow(); i++){
		int *counts = features(query,targets,i);
		
		switch(type){
			case TANIMOTO:
					//if(method=="Tanimoto" | method=="tanimoto") method <- function(a,b,c,d) (c+addone)/(a+b+c+addone)
				result(i) = (double)(counts[M11] + addone)	 / (double)(counts[M10]+counts[M01]+counts[M11]+addone);
				//printf("result(%d): %f\n",i,result(i));
				break;
			case EUCLIDIAN:
				//if(method=="Euclidean" | method=="euclidean") method <- function(a,b,c,d) sqrt((c+d+addone)/(a+b+c+d+addone))
				result(i) = sqrt( (double)(counts[M11]+counts[M00]+addone) / (double)(counts[M00]+counts[M01]+counts[M10]+counts[M11]+addone));
				break;
			case TVERSKY:
				//if(method=="Tversky" | method=="tversky") method <- function(a,b,c,d) (c+addone)/(alpha*a + beta*b+c+addone)
				result(i) = (double)(counts[M11] + addone)	 / (double)(alpha*counts[M10]+beta*counts[M01]+counts[M11]+addone);
				break;
			case DICE:
				//if(method=="Dice" | method=="dice") method <- function(a,b,c,d) (2*c+addone)/(a+c+b+c+addone)
				result(i) = (double)(2*counts[M11] + addone)	 / (double)(counts[M10]+counts[M01]+2*counts[M11]+addone);
				break;
			default:
				Rf_error("unknown similarity type");
		}

	delete [] counts;
	}
	
	return result;
}
int* features(NumericVector &query, NumericMatrix &targets, int targetRow){
	//int m00=0, m01=0, m10=0, m11=0;
	
//	printf("target: ");
//	for(int k=0;k < targets.ncol() && k < 30; k++)
//		printf("%d,",(int)targets(targetRow,k));
//	printf("\n");
	
	int i,j;
	int *counts= new int[4];
	for(int k=0; k < 4;k++) counts[k]=0;

	int lookup[2][2] = {{M00,M01}, 
						     {M10,M11}}; 
	if(query.size() != targets.ncol())
		Rf_error("query size must match the target size");

	int querySize = query.size();
	int targetSize = targets.ncol();
	for(i=0, j=0; i < querySize && j < targetSize; i++, j++ ){
		//printf("(i,j)=(%d,%d)\n",i,j);
		int q=query[i];
		int t = targets(targetRow,j);
		//printf("query: %d, target: %d, lookup: %d\n",q,t,lookup[q][t]);
		if( (q!=0 && q!=1) || (t!=0 && t!=1))
			Rf_error("non binary digits found");
		counts[lookup[q][t]]++;
	}
//	printf("counts: ");
//	for(int k=0;k < 4; k++)
//		printf("%d,",counts[k]);
//	printf("\n");
	

	return counts;
}
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
