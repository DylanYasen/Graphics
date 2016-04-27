#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <assert.h>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#define PI 3.1415926

GLint mouse_button = -1;
GLint mouse_state = -1;
void Find_Nearest(int, int);
int ni, nj;
int width;
GLdouble wx[4], wy[4], wz[4];

int SELECTED = GL_FALSE;

GLfloat CtrlPoints[10][2] = {
        {26, 33},{70, 51},{100, 83},{113, 137},
        {154, 56},{211, 56},{248, 135},
        {273, 106},{309, 98},{347, 99}
    };

GLfloat BezMatrix[4][4] = {
    { -1.0,  3.0, -3.0, 1.0},
    {  3.0, -6.0,  3.0, 0.0},
    { -3.0,  3.0,  0.0, 0.0},
    {  1.0,  0.0,  0.0, 0.0}};

GLfloat BP[10][2];

int in_box = -1;
int off_x, off_y;

// circle radius
GLfloat r = 25;

void init()
{
    int i, j, k;
    
    for(i=0; i<10; i++)
        for(j=0; j<2; j++) {
            BP[i][j] = 0.0;
            for(k=0; k < 4; k++) {
                BP[i][j] += BezMatrix[i%4][k] * CtrlPoints[k][j];
        }
    }
}

float* CalculateBezierPoint(float t, float startPoint[], float controlPoint1[], float controlPoint2[], float endPoint[]){
    float u = 1 - t;
    float tt = t*t;
    float uu = u*u;
    float uuu = uu * u;
    float ttt = tt * t;
    
    float *p = new float[2];
    p[0] = startPoint[0] * uuu;
    p[1] = startPoint[1] * uuu;
    
    p[0] += 3 * uu * t * controlPoint1[0];
    p[1] += 3 * uu * t * controlPoint1[1];
    p[0] += 3 * u * tt * controlPoint2[0];
    p[1] += 3 * u * tt * controlPoint2[1];
    p[0] += ttt * endPoint[0];
    p[1] += ttt * endPoint[1];
    
    return p;
}

void draw_Ctrl_Points()
{
    int i;
    glColor3f(1.0, 1.0, 1.0);
    glBegin(GL_POINTS);
    glPointSize(4.0);
    for(i=0; i<10; i++)
        glVertex2f(CtrlPoints[i][0], CtrlPoints[i][1]);
    glEnd();
    
    glBegin(GL_LINE_STRIP);
    for(i=0; i<10; i++)
        glVertex2f(CtrlPoints[i][0], CtrlPoints[i][1]);
    glEnd();
    
    for(i=0; i<10; i++) {
        glBegin(GL_LINE_LOOP);
        
        if(i== in_box)
            glColor3f(0.0,  1.0, 1.0);
        else
            glColor3f(1.0, 0.0, 0.0);
        
        glVertex2f(CtrlPoints[i][0] -5., CtrlPoints[i][1] - 5.);
        glVertex2f(CtrlPoints[i][0] -5., CtrlPoints[i][1] + 5);
        glVertex2f(CtrlPoints[i][0] + 5.,  CtrlPoints[i][1] + 5.);
        glVertex2f(CtrlPoints[i][0] +5. ,  CtrlPoints[i][1] -5.);
        
        glEnd();
    }
    
    // draw the highlighted one
    if(SELECTED)
    {
        glPointSize(5.0);
        glColor3f(1, 0, 0);
        glBegin(GL_POINTS);
        glVertex2f(CtrlPoints[in_box][0], CtrlPoints[ni][1]);
        glEnd();
    }
}


void DrawFilledCircle(GLfloat x,GLfloat y,GLfloat r){
    int i;
    int triangleAmount = 36;
    
    GLfloat twicePi = 2.0f * PI;
    
    glBegin(GL_TRIANGLE_FAN);
    glVertex2f(x, y); // center of circle
    for(i = 0; i <= triangleAmount;i++) {
        glVertex2f(
                   x + (r * cos(i *  twicePi / triangleAmount)),
                   y + (r * sin(i * twicePi / triangleAmount))
                   );
    }
    glEnd();
}

void DrawCircle(GLfloat x, GLfloat y, GLfloat r){
    int i;
    int lineAmount = 100;
    
    GLfloat twicePi = 2.0f * PI;
    
    glBegin(GL_LINE_LOOP);
    for(i = 0; i <= lineAmount;i++) {
        glVertex2f(
                   x + (r * cos(i *  twicePi / lineAmount)),
                   y + (r* sin(i * twicePi / lineAmount))
                   );
    }
    glEnd();
}


void draw_curve(float p1[],float p2[],float p3[],float p4[])
{
    glBegin(GL_POINTS);
    glColor3f(1.0, 1.0, 0.0);
    for(float u=0.0; u<=1.0; u+=0.01) {
        float *p = CalculateBezierPoint(u, p1, p2, p3, p4);
        glVertex2f(p[0], p[1]);
        delete p;
    }
    
    glEnd();
}


float moveStep = 100; // pixels per frame
float currentPercent = 0;
float lastUpdateTime = 0;
int path = 0;
float yOffset = 200;
float lastY = 0;
void display(void)
{
    // calculate delta time
    float timeSinceStart = glutGet(GLUT_ELAPSED_TIME);
    float dt = (timeSinceStart - lastUpdateTime) / 1000; // millisec to sec
    lastUpdateTime = timeSinceStart;
    
    glClear (GL_COLOR_BUFFER_BIT);

    draw_Ctrl_Points();
    draw_curve(CtrlPoints[0],CtrlPoints[1],CtrlPoints[2],CtrlPoints[3]);
    draw_curve(CtrlPoints[3],CtrlPoints[4],CtrlPoints[5],CtrlPoints[6]);
    draw_curve(CtrlPoints[6],CtrlPoints[7],CtrlPoints[8],CtrlPoints[9]);
    
    if(currentPercent <= 100){
        
        // get current pos on curve
        float *position = CalculateBezierPoint(currentPercent / 100, CtrlPoints[path], CtrlPoints[path+1], CtrlPoints[path+2], CtrlPoints[path+3]);
        
        // render ball at position
        DrawCircle(position[0],yOffset + position[1],r);
        
        // going upwards with less speed
        if (position[1] < lastY) {
            currentPercent += (moveStep/3) * dt ;
        }
        // drop it faster
        else{
            // update progress
            currentPercent += moveStep*dt;
        }
        
        // keep record of the y
        lastY = position[1];
        
        // reset progress when finish
        if (currentPercent > 100) {
            currentPercent = 0;
            
            // update current path
            path+=3;
            if(path>6)
                path = 0;
            
        }
        
        delete position;
    }
    
    glFlush();
    glutPostRedisplay();
}

void reshape(int w, int h)
{
    width=w;
    glViewport(0, 0, (GLsizei) w, (GLsizei) h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0.0, 500.0, 500.0, 0.0, -10.0, 10.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

static unsigned int m_off_x, m_off_y;
void move_ctrlpoint(int button, int x, int y)
{
    int i;
    //m_off_x = x - (int)(CtrlPoints[i][0]);
    //m_off_y = y - (int)(CtrlPoints[i][1]);
    printf("off_x=%d, off_y=%d\n", off_x, off_y);
    switch (button) {
        case GLUT_LEFT_BUTTON:
            if(in_box != -1) {
                CtrlPoints[in_box][0] = (float)(x-off_x);
                CtrlPoints[in_box][1] = (float)(y-off_y);
                init();
            };
        default:
            break;
    }
}


void mouse(int button, int state, int x, int y)
{
    int i;
    mouse_state = state;
    mouse_button = button;
    
	   //move_ctrlpoint(button, x, y);
    switch (button) {
        case GLUT_LEFT_BUTTON:
            if (mouse_state == GLUT_UP) {
                if(in_box != -1) {
                    CtrlPoints[in_box][0] = (float)(x-off_x);
                    CtrlPoints[in_box][1] = (float)(y-off_y);
                    init();
                }
                in_box = -1;
            }
            if (mouse_state == GLUT_DOWN) {
                printf("x=%d, y=%d\n", x, y);
                in_box = -1;
                for(i=0; i<10; i++) {
                    if( (x > (int)(CtrlPoints[i][0] - 15.))
                       && (y > (int)(CtrlPoints[i][1] -  15.))
                       && (x < (int)(CtrlPoints[i][0] +  15.))
                       && (y < (int)(CtrlPoints[i][1] +  15.)))
                    {
                        in_box = i;
                        
                        off_x = x - (int)(CtrlPoints[i][0]);
                        off_y = y - (int)(CtrlPoints[i][1]);
                        printf("inside_box %d!\n", in_box);
                        printf("%d %d %d %d \n",  (int)(CtrlPoints[i][0] -5.),
                               (int)(CtrlPoints[i][1] - 5.),
                               (int)(CtrlPoints[i][0] + 5.),
                               (int)(CtrlPoints[i][1] + 5.));
                        break;
                    }
                }
            }
            break;
        default:
            break;
    }
    /*
     */
    
    glFlush();
    display();
    
}


void Find_Nearest(int x, int y)
{
    int i, j;
    GLint viewport[4];
    GLdouble mvmatrix[16], projmatrix[16];
    GLdouble td, dd;
    
    glGetIntegerv(GL_VIEWPORT, viewport);
    glGetDoublev(GL_MODELVIEW_MATRIX, mvmatrix);
    glGetDoublev(GL_PROJECTION_MATRIX, projmatrix);
    
    for(i=0; i<10; i++)
    {
        gluProject((GLdouble)CtrlPoints[i][0],
                   (GLdouble)CtrlPoints[i][1],
                   0,
                   mvmatrix, projmatrix, viewport,
                   &(wx[i]), &(wy[i]), &(wz[i]));
        wy[i]=(GLdouble)width-wy[i];
        printf("wx=%lf, wy=%lf, wz%lf\n", wx[i], wy[i], wz[i]);
    };
    
    printf("\n");
    printf("x=%d, y=%d\n", x, y);
    dd=9e+9;
    
    ni=0;
    for(i=0; i<10; i++)
    {
        td=((GLdouble)x -wx[i]) * ((GLdouble)x-wx[i]) +
	       ((GLdouble)y -wy[i]) * ((GLdouble)y-wy[i]);
        if(td<dd) {
            dd=td;
            ni=i;
            printf("ni=%d\n", ni);
        };
    }
    if(dd<50.) {
        SELECTED=GL_TRUE;
    }
}

void motion(int x, int y)
{
    printf("here\n");
    //if((mouse_button == GLUT_LEFT_BUTTON) && (mouse_state==GLUT_DOWN))
    {
        if(!SELECTED)
        {
            Find_Nearest(x, y);
        }
        if(SELECTED) {
            move_ctrlpoint(mouse_button, x, y);
            if(in_box != -1) {
                CtrlPoints[in_box][0] = (float)(x-off_x);
                CtrlPoints[in_box][1] = (float)(y-off_y);
                init();
            }
        }
    }
    //else
    //{
    //}; // end if(mouse_button)
}

void keyboard(unsigned char key, int x, int y)
{
    switch (key) {
        case 27:
            exit(0);
            break;
    }
}

int main(int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB);
    glutInitWindowSize (500, 500);
    glutInitWindowPosition (100, 100);
    glutCreateWindow (argv[0]);
    init();
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    
    glutKeyboardFunc (keyboard);
    glutMainLoop();
    return 0;
}