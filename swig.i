%{
#include "script.h"
#include <vector>
%}
class Descriptors{
public:
   Descriptors();
   int parse_sdf(const char* sdf);
   int parse_sdfile(const char* sdfile);
   int parse_smiles(const char* smile);
   unsigned int get_descriptor(unsigned int i);
   unsigned int get_len();
};

double similarity(Descriptors* d1, Descriptors* d2);

