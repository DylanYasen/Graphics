
#include "svUtil.h"

namespace __svl_lib {

svVector3 svGetPerpendicularVector(const svVector3 & v)
{
    float x,y,z;
    v.getValue(x,y,z);
    
    svVector3 proj(x,y,0);
    
    return normalize(proj.cross(v));
}

svVector3 svGetRotatePoint(const svVector3& start, const svVector3& org, 
			const svVector3& axis, svScalar rotate_degree)
{
// CMSC: to be implemented 
   return svVector3(0,0,0);
}

// get the normal vector
//      v1-------------v0
//                    /
//                   /
//                  / v2
//
svVector3 svGetNorm(const svVector3& v1, const svVector3& v0,
                    const svVector3& v2)
{
  return normalize(cross((v1-v0), (v2-v0)));
}

svVector3 svAverage(const svVector3& v1, const svVector3& v2)
{
  return normalize(v1+v2);
}

}
