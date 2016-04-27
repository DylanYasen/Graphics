
/**----------------------------------------------
 * Author: Jian Chen
 * Nov 2013
 *
 * utilities  routines
 */
#ifndef __SVL_MGL_H
#define __SVL_MGL_H

#ifdef __APPLE__
  #include <GLUT/glut.h>
#else
  #include <GL/glut.h>
#endif

namespace __svl_lib {
void WriteToPNM(int w, int h, char *fn);
void drawAxis(GLfloat x, GLfloat y, GLfloat z, GLfloat l);
}

#endif //__MGL_H
