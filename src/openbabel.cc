


#include <iostream>
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/oberror.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <string>

using namespace std;
using namespace OpenBabel;

extern "C" {
   SEXP smile2sdf_file(SEXP smileFile, SEXP sdfFile);
   SEXP smile2sdf_string(SEXP smileFile);
}

SEXP smile2sdf_file(SEXP smileFile, SEXP sdfFile)
{

   ifstream ifs(CHAR(STRING_ELT(smileFile,0)));
   if(!ifs){
      error("cannot open smile file for input: %s",CHAR(STRING_ELT(smileFile,0)));
      return R_NilValue;
   }

   ofstream ofs(CHAR(STRING_ELT(sdfFile,0)));
   if(!ofs){
      error("cannot open sdf file for output: %s",CHAR(STRING_ELT(sdfFile,0)));
      return R_NilValue;
   }

   OpenBabel::OBConversion conv(&ifs,&ofs);

   if(!conv.SetInAndOutFormats("SMI","SDF")) {
      error("conversion from smile to sdf not available");
      return R_NilValue;
   }

   conv.AddOption("gen2D",OBConversion::GENOPTIONS);
	
   conv.Convert();

   return R_NilValue;

}

SEXP smile2sdf_string(SEXP smileFile)
{

   ifstream ifs(CHAR(STRING_ELT(smileFile,0)));
   if(!ifs){
      error("cannot open smile file for input: %s",CHAR(STRING_ELT(smileFile,0)));
      return R_NilValue;
   }


	ostringstream ofs;

   OpenBabel::OBConversion conv(&ifs,&ofs);

   if(!conv.SetInAndOutFormats("SMI","SDF")) {
      error("conversion from smile to sdf not available");
      return R_NilValue;
   }

   conv.AddOption("gen2D",OBConversion::GENOPTIONS);
	
   conv.Convert();

   /* start at returning data in a list of vectors
    * problem is we may not be able to fit entire file into memory

   istringstream smileStream(ofs.str());
   string line;
   vector<string>
   while(getline(smileStream,line))
   {
   }
   */

   SEXP sdfString;
   PROTECT(sdfString = NEW_STRING(1));
   SET_STRING_ELT(sdfString,0,mkChar(ofs.str().c_str()));
   UNPROTECT(1);


   return sdfString;

}
