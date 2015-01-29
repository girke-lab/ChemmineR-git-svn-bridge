
#include <Rcpp.h>

using namespace Rcpp;

#include <R.h>
#include <boost/algorithm/string.hpp>


RcppExport SEXP cstrsplit( SEXP l ){

  std::vector<std::string> strs;
  const char *line = CHAR(STRING_ELT(l,0));
  boost::split(strs, line, boost::is_any_of("\t "),boost::token_compress_on);

  CharacterVector output(strs.begin(),strs.end());
  return output;

}
