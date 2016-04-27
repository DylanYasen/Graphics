#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <assert.h>
#include <GLUT/GLUT.h>

#include "svUtil.h"
#include "ivTrackball.h"
#include "ivview.h"
#include "jcUtil/MGL.h"

#include "svVectorField.h"

#define IMAGE_WIDTH  512
#define IMAGE_HEIGHT 512

using namespace __svl_lib;
using namespace std;

void reshape(int w, int h);
void display(void);
void key(unsigned char key, int x, int y);

#define checkImageWidth 64
#define checkImageHeight 64

#define SPHERE_SLICE 18

static GLubyte checkImage[checkImageWidth][checkImageWidth][4];
static GLuint texName;

view3d view_info;

GLuint image_width;
GLuint image_height;

GLfloat nx, ny, nz;
Trackball trackball;


// ========= My Vars ========= //
GLfloat r,l,t,b,n,f;
svVector4 points[360]; // shpere outline
svVector3 viewPoint(5,4,3);
svVector3 viewDir(-5,-4,-3);
svVector3 viewUp(0,1,0);
svVector3 projNormal(5,4,3);
GLfloat projDistance = 5;
GLfloat viewWidth = 2.5f;
GLfloat viewHeight = 2.5f;
svVector3 sphereCenter(0,0,0);
GLfloat sphereRadius = 1;
GLfloat imagePlaneOffsetVal = -4;
svVector3 imagePlaneOffset;
GLfloat nearPlaneOffset = 3.5f;

void AlignWithAxis(const svVector3 axis, GLfloat &rx,GLfloat &ry,GLfloat &rz,GLfloat &angle){
    svVector3 dirVec = axis;
    svVector3 normDirVec = normalize(dirVec);
    float dirX,dirY,dirZ;
    dirVec.getValue(dirX, dirY, dirZ);
    
    // calculate rotation axis and angle
    svVector3 up (0,0,1);
    svVector3 axisTo = normDirVec;
    svVector3 rotationAxis = up.cross(axisTo);
    rotationAxis.getValue(rx, ry, rz);
    angle = acos(up.dot(axisTo));
    angle *= (180/3.1415926);
}

void DrawFrustum(){
    glColor3f(1,1,1);
    
    GLfloat length = viewPoint.length() + imagePlaneOffset.length();

    // bottom left
    glBegin(GL_LINES);
    glVertex3f(l,b,0); glVertex3f(l,b,length);
    glVertex3f(r,b,0); glVertex3f(r,b,length);
    glVertex3f(l,t,0); glVertex3f(l,t,length);
    glVertex3f(r,t,0); glVertex3f(r,t,length);
    glEnd();
}

void DrawViewPlane(){
    
    glColor3f(0.5f, 0.5f, 0.5f);
    
    float rx,ry,rz;
    float angle;
    AlignWithAxis(viewDir,rx, ry, rz, angle);
    
    svVector3 viewDirNorm = normalize(viewDir);
    svVector3 nearPoint = viewPoint + viewDirNorm * nearPlaneOffset;
    svVector3 farPoint =  nearPoint + viewDirNorm * projDistance;

    // near plane
    glPushMatrix();
    glTranslated(nearPoint[0],nearPoint[1],nearPoint[2]);
    glRotatef(angle, rx, ry, rz);
    glBegin(GL_QUADS);                      // Draw A Quad
    glVertex3f(l, t, 0.0f);              // Top Left
    glVertex3f(r, t, 0.0f);              // Top Right
    glVertex3f(r, b, 0.0f);              // Bottom Right
    glVertex3f(l, b, 0.0f);              // Bottom Left
    glEnd();
    glPopMatrix();
    
    glColor3f(0,0,0);
    
    glPushMatrix();
    glTranslated(farPoint[0],farPoint[1],farPoint[2]);
    glRotatef(angle, rx, ry, rz);
    
    glBegin(GL_QUADS);                      // Draw A Quad
    glVertex3f(l, t, 0.0f);              // Top Left
    glVertex3f(r, t, 0.0f);              // Top Right
    glVertex3f(r, b, 0.0f);              // Bottom Right
    glVertex3f(l, b, 0.0f);              // Bottom Left
    glEnd();
    
    glPopMatrix();
}

void DrawCamera(){
    
    float rx,ry,rz;
    float angle;
    AlignWithAxis(viewDir,rx, ry, rz, angle);
    
    glPushMatrix();
    glTranslatef(viewPoint[0],viewPoint[1],viewPoint[2]);
    glRotatef(angle, rx, ry, rz);
    DrawFrustum();

    // =========== drawing camera ============ //
    // draw camera axis
    glPushMatrix();
    drawAxis(0, 0, 0, 2);
    glPopMatrix();
    
    // draw cam
    glColor3f(0,0,1);
    glutSolidCube(0.5f);
    
    // draw cam body
    glPushMatrix();
    glColor3f(0,0,0);
    glTranslated(0, 0, -0.5);
    glutSolidCube(1);
    glPopMatrix();
    // ============================================= //
    
    glPopMatrix();
}

bool sortByY(const svVector4 &lhs, const svVector4 &rhs) { return lhs[1] < rhs[1]; }


// CMSC: replace the texture image generate with the calculated pixel
// values.
void makeCheckImage(void)
{
   int i, j;
    
   for (i = 0; i < checkImageHeight; i++) {
      for (j = 0; j < checkImageWidth; j++) {
         //c = ((((i&0x8)==0)^((j&0x8)==0)))*255;
         checkImage[i][j][0] = (GLubyte) 255;
         checkImage[i][j][1] = (GLubyte) 255;
         checkImage[i][j][2] = (GLubyte) 255;
         checkImage[i][j][3] = (GLubyte) 255;
      }
   }
    
    for (int i = 0; i < 360; i++) {
        
        float fx,fy,fz;
        points[i].getValue(fx, fy, fz, fz);
        
        int x = (int)fx;
        int y = (int)fy;
        int z = (int)fz;
        
        int nextY = 0;
        int nextX = 0;
        if(i < 360){
            nextX = (int)points[i+1][0];
            nextY = (int)points[i+1][1];
        }
        
        // same line
        if(nextY == y){
            for (int j = x; j <= nextX; j++) {
                checkImage[j][y][0] = (GLubyte) 255;
                checkImage[j][y][1] = (GLubyte) 0;
                checkImage[j][y][2] = (GLubyte) 0;
            }
        }
    }
}

void ScanFill(){
    
}

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
        glEnable(GL_TEXTURE_2D);

        GLfloat m[16];
        trackball.getMatrix().getValue(m);

        glPushMatrix();
        glMultMatrixf(m);
    
        glEnable(GL_TEXTURE_2D);
    
    // draw texture quad
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
    glBindTexture(GL_TEXTURE_2D, texName);
    
    glPushMatrix();
    GLfloat rx,ry,rz;
    GLfloat angle;
    AlignWithAxis(viewDir, rx, ry, rz, angle);

    glTranslated(imagePlaneOffset[0],imagePlaneOffset[1],imagePlaneOffset[2]);
    glRotated(angle, rx, ry, rz);
    
    glBegin(GL_QUADS);
    glTexCoord2f(0.0, 0.0); glVertex3f(-1.0, -1.0, 0.0);
    glTexCoord2f(0.0, 1.0); glVertex3f(-1.0, 1.0, 0.0);
    glTexCoord2f(1.0, 1.0); glVertex3f(0.0, 1.0, 0.0);
    glTexCoord2f(1.0, 0.0); glVertex3f(0.0, -1.0, 0.0);
    glEnd();
    
    glPopMatrix();
    
    glFlush();
    glDisable(GL_TEXTURE_2D);

    glDisable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);
    glDisable(GL_LIGHT0);
    
    // draw shpere
    glColor3f(1, 0, 0);
    glutSolidSphere(sphereRadius, SPHERE_SLICE, SPHERE_SLICE);

    DrawViewPlane();
    
    DrawCamera();
    
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

void init()
{
    // CMSC: init your geometry if any
    nx = 64;
    ny = 64;
    
    l = -1.25f;
    r = 1.25f;
    n = 0;
    f = -5;
    t = 1.25f;
    b = -1.25f;
    
    imagePlaneOffset = svVector3(0,0,0) - normalize(viewDir) * imagePlaneOffsetVal;
    
    svMatrix4 mvp(svVector4(nx/2, 0, 0, (nx-1)/2),
                  svVector4(0, ny/2, 0, (ny-1)/2),
                  svVector4(0, 0, 1, 0),
                  svVector4(0, 0, 0, 1)
                  );
    
//    svMatrix4 mor(svVector4(2/(r-l), 0, 0, 0),
//                  svVector4(0, 2/(t-b), 0, 0),
//                  svVector4(0, 0, 2/(n-f), 0),
//                  svVector4(-(r-l)/(r+l),-(t+b)/(t-b),-(n+f)/(f-n),1)
//                  );
    
    svMatrix4 mor(svVector4(1/r, 0, 0, 0),
                  svVector4(0, 1/t, 0, 0),
                  svVector4(0, 0, 2/(n-f), -(n+f)/(f-n)),
                  svVector4(0, 0, 0, 1)
                  );
 
    svMatrix4 M = mvp * mor;
    
    svVector4 topPoint = svVector4(0,1,0,1);
    GLfloat rotateAngle = 360 / 360 * (3.1415926 / 180);
  
    for (int i = 0; i < 360; i++) {
        
        GLfloat deg = rotateAngle * i;
        
        svMatrix4 rotM(
                       svVector4(cos(deg), sin(deg),0,0),
                       svVector4(-sin(deg),cos(deg),0,0),
                       svVector4(0,0,1,0),
                       svVector4(0,0,0,1)
                       );
        
        svVector4 point = rotM * topPoint;
        svVector4 pointProject = M * point ;
        
        points[i] = pointProject;
    }

    
    // get points ready for scan line algorithm
    std::sort(std::begin(points),std::end(points), sortByY);
    
    
  svVector3 center = svVector3(0, 0, 0);
  center.getValue(view_info.coi);

  // CMSC: need to make your own image
  makeCheckImage();

  glShadeModel(GL_FLAT);
  glEnable(GL_DEPTH_TEST);

  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

#ifdef GL_VERSION_1_1
   glGenTextures(1, &texName);
   glBindTexture(GL_TEXTURE_2D, texName);
#endif

   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
#ifdef GL_VERSION_1_1
   glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, checkImageWidth, checkImageHeight,
                0, GL_RGBA, GL_UNSIGNED_BYTE, checkImage);
#else
   glTexImage2D(GL_TEXTURE_2D, 0, 4, checkImageWidth, checkImageHeight,
                0, GL_RGBA, GL_UNSIGNED_BYTE, checkImage);
#endif

  GLfloat x=2, y=2, z=2;
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
        glutCreateWindow("My excellent viewing project");
        glutInitWindowPosition(200, 200);

        init();

        glutDisplayFunc(display);
        glutReshapeFunc(reshape);
        glutKeyboardFunc(key);
        glutMouseFunc(mouse);
        glutMotionFunc(motion);

        glutMainLoop();
        return 0;
}


