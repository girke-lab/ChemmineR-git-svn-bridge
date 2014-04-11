/* Simple Objects that provide an easier interface for SWIG wrapper generation */

#ifndef SCRIPT_H_
#define SCRIPT_H_
#include <vector>
#include <fstream>
class Descriptors;


class Descriptors {
private:
	std::vector<unsigned int> descs;
public:
	Descriptors();
	int parse_sdf(const char* sdf);
	int parse_sdfile(const char* sdfile);
	int parse_smiles(const char* smile);
	unsigned int get_descriptor(unsigned int i);
	unsigned int get_len();
	virtual ~Descriptors();
	friend double similarity(Descriptors*, Descriptors*);
};

double similarity(Descriptors* d1, Descriptors* d2);


#endif /*SCRIPT_H_*/
