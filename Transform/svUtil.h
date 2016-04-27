
#ifndef __SV_UTIL_H
#define __SV_UTIL_H

//#include "svVectorField.h"
#include "svType.h"
#include "svArray.h"

namespace __svl_lib {
#define SWAP(T, a, b)   { T t; t=a; a=b; b=t;  }

#define SV_PI	  (3.14159265358979323846)
#define SV_2_PI   (3.14159265358979323846 * 2.0)
#define SV_SMALL  (1e-6)

#define svToDegree(x)             ((x)*(180.)/SV_PI)
#define svToRadian(x)	           ((x)*SV_PI/(180.))
#define svMax(a,b)                (((a) > (b)) ? (a) : (b))
#define svMin(a,b)                (((a) < (b)) ? (a) : (b))
#define svClamp(v, min_v, max_v)  { v = svMax(svMin(v, max_v), min_v); }

	
template <class T> inline
  void svSwap(T & a, T & b)
{  T t = a; a = b; b = t;  }

svVector3 svGetPerpendicularVector(const svVector3 & v);
svVector3 svGetRotatePoint(const svVector3& start, const svVector3& org, 
                           const svVector3& axis, svScalar rotate_degree);
svVector3 svGetNorm(const svVector3& v1, const svVector3& v0,
                    const svVector3& v2);
svVector3 svAverage(const svVector3& v1, const svVector3& v2);

}


#endif // __SV_UTIL_H
