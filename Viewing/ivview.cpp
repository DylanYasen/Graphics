
#include "ivview.h"

namespace __svl_lib {
view3d::view3d()
{
  init();
}

view3d::~view3d()
{
}

void view3d::init()
{
  //type = PARALLEL;
  type = PERSPECTIVE;
  yon  = 1000;
  hither = 0.001;
  //yon  = 200;
  //hither = 1;
  image_plane = hither;
  view_angle = 45;
  head_tilt  = 0;
  aspect_ratio = 1.0;
}

}
