#ifndef __SV_ARRAY_H
#define __SV_ARRAY_H

#include "MArray.h"
#include "svVectorMatrix.h"
#include "svType.h"
//#include "svParticle.h"
namespace __svl_lib {
typedef MArray<svInt> svIntArray;
typedef MArray<svScalar> svScalarArray;

typedef svVector3 svVector3;
typedef MArray<svVector3> svVector3Array;
typedef svVector3Array* svVector3ArrayPtr;

typedef MArray<svVector4> svVector4Array;
typedef svVector4Array* svVector4ArrayPtr;

typedef svVector3* svVector3Ptr;
typedef MArray<svVector3Ptr> svVector3PtrArray;
}
#endif // __SV_ARRAY_H
