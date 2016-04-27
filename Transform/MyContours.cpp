#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <assert.h>
#include <GLUT/glut.h>

#include "svUtil.h"
#include "ivTrackball.h"
#include "ivview.h"
#include "MGL.h"

#include "svVectorField.h"

#define IMAGE_WIDTH  512
#define IMAGE_HEIGHT 512

#define CONE_BASE_R 0.2f
#define CONE_HEIGHT 0.1f
#define CYL_BASE_R 0.03f
#define CYL_HEIGHT 0.65f

using namespace __svl_lib;

void reshape(int w, int h);
void display(void);
void key(unsigned char key, int x, int y);


view3d view_info;

GLuint image_width;
GLuint image_height;

GLint nx, ny, nz;
Trackball trackball;

svVectorField  *myVectorField;

void InitLight()
{
  //setting of lighting
  GLfloat mat_diffuse[] = { 0.8, .0, 0.0, .4};
  GLfloat mat_specular[] = { 1, 1, 1, 1 };
  GLfloat mat_shininess[] = { 20.0 };
  //GLfloat light_position[] = { 24, 24, 60, 0.0 };
  GLfloat light_position[] = { 50,50,50, 0.0 };
  GLfloat white_light[] = { 1.0, 1.0, 1.0, 1.0 };
  glClearColor(0.0, 0.0, 0.0, 0.0);
  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat_diffuse);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
  glLightfv(GL_LIGHT0, GL_POSITION, light_position);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, white_light);
  glLightfv(GL_LIGHT0, GL_SPECULAR, white_light);
}

//***************************
// GLUT callback functions
//****************************
void reshape(int w, int h)
{

        image_width  = w;
        image_height = h;

        glViewport (0, 0, w, h);

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();

        if (view_info.type == PARALLEL){
                glOrtho(view_info.left, view_info.right, view_info.bottom, view_info.top,
                        view_info.hither, view_info.yon);

                glMatrixMode(GL_MODELVIEW);
                glLoadIdentity();

        }
        else{ // perspective
                float aspect = view_info.aspect_ratio * float(w)/float(h);
                float GLfloat_fov;
                if ( aspect < 1 ){
                  // fovy is a misnomer.. we GLfloatly mean the fov applies to the
                  // smaller dimension
                  float fovx, fovy;
                  fovx = fovy = view_info.view_angle;
                  GLfloat_fov = svToDegree(2*atan(tan(svToRadian(fovx/2))/aspect));
                }
                else{
                        GLfloat_fov = view_info.view_angle;
                }

                gluPerspective(GLfloat_fov,
                                           aspect,
                                           view_info.hither,
                                           view_info.yon);

                glMatrixMode(GL_MODELVIEW);
                glLoadIdentity();
                gluLookAt(view_info.eye[0],view_info.eye[1],view_info.eye[2],
                                view_info.coi[0],view_info.coi[1],view_info.coi[2],
                                0, 1, 0);
        }
        InitLight();

}

void display(void)
{
        // set new model view and projection matrix
        glClearColor(0.3, 0.3, 0.3, 1);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glEnable(GL_DEPTH_TEST);

        glEnable(GL_COLOR_MATERIAL);
        glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);

        GLfloat m[16];
        trackball.getMatrix().getValue(m);

        glPushMatrix();
        glMultMatrixf(m);

        glEnable(GL_LIGHTING);
        glEnable(GL_LIGHT0);
        glEnable(GL_TEXTURE_2D);
    
    // CMSC:  rendering your geometry here...

    myVectorField->Render();
    
    
    drawAxis(0., 0., 0., 1.);
    
    
    
    glDisable(GL_TEXTURE_2D);
    glDisable(GL_BLEND);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);
    glDisable(GL_LIGHT0);

    glPopMatrix();
    glFlush();
    glutSwapBuffers();
}

void key(unsigned char key, int x, int y)
{
  switch (key) {
        case '\033':
        case 'q':
        case 'Q':
                exit(0);
                break;
  }
}

void mouse(int button, int state, int x, int y)
{
        long s=0x00000000;
    s |= (button == GLUT_LEFT_BUTTON)   ? ((state == GLUT_DOWN) ? Trackball::LBUTTON_DOWN : Trackball::LBUTTON_UP) : 0;
    s |= (button == GLUT_MIDDLE_BUTTON) ? ((state == GLUT_DOWN) ? Trackball::MBUTTON_DOWN : Trackball::MBUTTON_UP) : 0;
    s |= (button == GLUT_RIGHT_BUTTON)  ? ((state == GLUT_DOWN) ? Trackball::RBUTTON_DOWN : Trackball::RBUTTON_UP) : 0;

    int key_state = glutGetModifiers();
    s |= (key_state & GLUT_ACTIVE_CTRL)  ? Trackball::CTRL_DOWN  : 0;
    s |= (key_state & GLUT_ACTIVE_ALT)   ? Trackball::ALT_DOWN   : 0;
    s |= (key_state & GLUT_ACTIVE_SHIFT) ? Trackball::SHIFT_DOWN : 0;

        if (s & Trackball::BUTTON_DOWN){
        trackball.mouseDown(s, x, y);
        }

        if (s & Trackball::BUTTON_UP){
        trackball.mouseUp(s, x, y);
        }
}

void motion(int x, int y)
{
        trackball.mouseMotion(x, y);
        glutPostRedisplay();
}

//******************************************
// initialization code for GL and NV_EXT
//******************************************

void init(char *infname)
{
// CMSC: init your geometry if any

  svVector3 center = svVector3(0, 0, 0);
  center.getValue(view_info.coi);

  // CMSC: please replace the numbers with the real geometry boundary 
  GLfloat x,y,z;
  myVectorField = new svVectorField(infname);
  myVectorField->Generate();
  myVectorField->GetPhysicalDimension(&x,&y, &z);

  view_info.eye[0] = x/2.0;
  view_info.eye[1] = y/2.0;
  view_info.eye[2] = z*2.0;

  trackball.setEye(view_info.eye);
  trackball.setFocus(center);
  trackball.setWindowSize(IMAGE_WIDTH, IMAGE_HEIGHT);
}

//**********************
// program entry
//**********************
int main(int argc, char** argv)
{
        glutInit(&argc, argv);
        glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);

        image_width  = IMAGE_WIDTH;
        image_height = IMAGE_HEIGHT;

        glutInitWindowSize(image_width, image_height);
        glutCreateWindow("sweeping vectors");
        glutInitWindowPosition(200, 200);

    
        char* fileName= "/Users/Yadikaer/Projects/OpenGL/Project1/Project1/data/4Contours.svl";
        init(fileName);
        //init(argv[1]);

        glutDisplayFunc(display);
        glutReshapeFunc(reshape);
        glutKeyboardFunc(key);
        glutMouseFunc(mouse);
        glutMotionFunc(motion);

        glutMainLoop();
        return 0;
}


