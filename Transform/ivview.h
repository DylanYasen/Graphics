

#ifndef __IV_VIEW_H
#define __IV_VIEW_H

#include "svType.h"

enum PROJECTION_TYPE {
        PERSPECTIVE,
        PARALLEL
};

namespace __svl_lib {
class view3d{
public:
  view3d();
  virtual ~view3d();

  virtual void init(); // init
  void setEyePos(svScalar x, svScalar y, svScalar z) // eye pos
  { eye[0]=x; eye[1]=y; eye[2]=z; }
  void getEyePos(svScalar *x, svScalar *y, svScalar *z) const
  { *x = eye[0]; *y = eye[1]; *z = eye[2]; }

  void setCoi(svScalar x, svScalar y, svScalar z) // eye pos
  { coi[0]=x; coi[1]=y; coi[2]=z; }
  void getCoi(svScalar *x, svScalar *y, svScalar *z) const
  { *x = coi[0]; *y = coi[1]; *z = coi[2]; }

//private:
  PROJECTION_TYPE type;// projection type, either PERSPECTIVE or PARALLEL
  svScalar view_angle;  // angle between line of sight and
                       //   top of view pyramid
  svScalar yon;         // distance from eye to far plane
  svScalar hither;      // distance from eye to near plane
  svScalar image_plane; // distance from eye to image plane
  svScalar head_tilt;    // (CCW) angle of head tilt in degrees
  svScalar aspect_ratio; // width/height of view pyramid

  svScalar eye[3];      // eye position
  svScalar coi[3];      // center of interest

  svScalar left;                // for parallel
  svScalar right;
  svScalar bottom;
  svScalar top;
};
}
#endif

