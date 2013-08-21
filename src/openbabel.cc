


#include <iostream>
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/oberror.h>
#include <R.h>
#include <Rinternals.h>

using namespace std;
using namespace OpenBabel;

extern "C" {
   SEXP smile2sdf_file(SEXP smileFile, SEXP sdfFile);
}

SEXP smile2sdf_file(SEXP smileFile, SEXP sdfFile)
{

   obErrorLog.SetOutputLevel(obDebug);
	OpenBabel::OBConversion conv;
   conv.SetInAndOutFormats("SDF", "SDF");
	//OpenBabel::OBMol *mol = new OpenBabel::OBMol;
	OpenBabel::OBMol mol ;
	mol.Clear();

	bool ret = conv.ReadFile(&mol, "test.sdf");
	std::cout<<"ret: "<<ret<<endl;
   std::cout<<"NumAtoms: "<<mol.NumAtoms()<<std::endl;
   std::cout<<"weight: "<<mol.GetMolWt()<<std::endl;





/*

   ifstream ifs(CHAR(STRING_ELT(smileFile,0)));
   if(!ifs){
      error("cannot open smile file for input: %s",CHAR(STRING_ELT(smileFile,0)));
      return;
   }

   ofstream ofs(CHAR(STRING_ELT(sdfFile,0)));
   if(!ofs){
      error("cannot open sdf file for output: %s",CHAR(STRING_ELT(sdfFile,0)));
      return;
   }

   OpenBabel::OBConversion conv(&ifs,&ofs);
	OpenBabel::OBConversion conv;
   if(!conv.SetInAndOutFormats("SMILE","SDF")) {
      error("conversion from smile to sdf not available");
      return;
   }
	
   conv.Convert();
   */

   return R_NilValue;

}
