

/**----------------------------------------------
 * Author: Jian Chen
 * Nov 2013
 */

#include "MGL.h"

#include <stdio.h>
#include <fstream>

namespace __svl_lib {

void WriteToPNM(int w, int h, char *fn)
{
  char *image_buf = new char [w*h*3];
  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  glReadPixels(0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, image_buf);

  FILE* fp;
  if (NULL != (fp = fopen(fn, "wb"))){
    // Write the 'header' information
    fprintf(fp, "P6 %d %d 255\n", w, h);
    for (int i=h-1; i >= 0; --i)
    {
       // write binary data
       fwrite(image_buf+3*i*w, sizeof(unsigned char), 3*w, fp); 
    } 
    fclose(fp);
  }
  delete [] image_buf;
}

void drawAxis(GLfloat x, GLfloat y, GLfloat z, GLfloat l)
{
  glBegin(GL_LINES);
    glColor3f(1, 0, 0);
    glVertex3f(x,y,z);	glVertex3f(x+l, y, z);
    glColor3f(0, 1, 0);
    glVertex3f(x,y,z);	glVertex3f(x, y+l, z);
    glColor3f(0, 0, 1);
    glVertex3f(x,y,z);	glVertex3f(x, y, z+l);
  glEnd();
}

}
