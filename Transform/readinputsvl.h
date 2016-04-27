
#ifndef __READ_INPUT_SVL_H_
#define __READ_INPUT_SVL_H_

#include "svArray.h"

#define  SVL_FORMAT  100 // default
#define  RBF_INFOR_FORMAT 110
using namespace __svl_lib;

void readinputsvl(char *infile,
                  svVector3Array *vec3profile, 
                  svVector3Array *vec3colorProfile, 
		  int *seed_num);

#endif //__READ_INPUT_SVL_H_
