#include <Rcpp.h>
#include <math.h>
#include <bigmemory/MatrixAccessor.hpp>

// Written by Tyler Backman April 2015

using namespace Rcpp;

#define  TANIMOTO 0
#define  EUCLIDIAN 1
#define  TVERSKY 2
#define  DICE 3

RcppExport SEXP bigMatrixSimilarity(SEXP pQueryVector, SEXP pTargetMatrix, 
    SEXP typeS, SEXP addoneS, SEXP alphaS, SEXP betaS){
BEGIN_RCPP
    // tell Rcpp what class to use for function parameters
    NumericVector queryVector(pQueryVector);
    XPtr<BigMatrix> targetMatrix(pTargetMatrix);
    MatrixAccessor<char> mat(*targetMatrix); // use char to save memory over int
    IntegerVector addoneVector(addoneS);
    NumericVector typeV(typeS), alphaV(alphaS),betaV(betaS);
    int addone = addoneVector[0]; // add to numerator and denominator to avoid divide by zero
    int type = typeV[0];
    double alpha = alphaV[0];
    double beta = betaV[0];
    
    // create R vector to store results
    NumericVector similarities(targetMatrix->nrow());

    int sharedBits;
    int unsharedBits;
    int neither;
    int queryOnly;
    int targetOnly;

    // loop over targets (rows) and compute simlarity
    for(size_t i=0; i < targetMatrix->nrow(); ++i){
        sharedBits = 0;
        neither = 0;
        queryOnly = 0;
        targetOnly = 0;

        // loop over each bit (columns) and count shared vs unshared 
        for(size_t j=0; j < targetMatrix->ncol(); ++j){
            if(queryVector[j] && mat[j][i])
                sharedBits++;
            else if(queryVector[j] && ! mat[j][i])
                queryOnly++;
            else if(! queryVector[j] && mat[j][i])
                targetOnly++;
            else
                neither++;
        }

        unsharedBits = queryOnly + targetOnly;

        switch(type){
            case TANIMOTO:
                similarities[i] = (double)(sharedBits + addone) / (double)(sharedBits + unsharedBits + addone);
                break;
            case EUCLIDIAN:
                similarities[i] = sqrt( (double)(sharedBits + neither + addone) / (double)(neither + 
                    unsharedBits + sharedBits + addone)); 
                break;
            case TVERSKY:
                similarities[i] = (double)(sharedBits + addone) / (double)((alpha * queryOnly) +
                    (beta * targetOnly) + sharedBits + addone);
                break;
            case DICE:
                similarities[i] = (double)((2 * sharedBits) + addone) /
                    (double)(unsharedBits + (2 * sharedBits) + addone);
                break;
        }
    }

    return similarities;
END_RCPP
}
