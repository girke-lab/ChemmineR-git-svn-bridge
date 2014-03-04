#ifndef DEBUG_H_
#define DEBUG_H_
#ifdef DEBUG
#define LOC std::cout << "debug:" << __FILE__ << ":" << __LINE__ << " ";
#define DEBUG_PRINT(fmt,...) LOC printf(fmt,__VA_ARGS__);
#define DEBUG_MSG(text) LOC std::cout << text << std::endl;
#define DEBUG_VAR(text) LOC std::cout << (#text) << "=" << text << std::endl;
#else
#define DEBUG_PRINT(fmt,...)
#define DEBUG_MSG(text)
#define DEBUG_VAR(text)
#endif
#include <iostream>

#endif /*DEBUG_H_*/
