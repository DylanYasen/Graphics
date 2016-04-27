
//********************************************
// CMSC: To be implemented

#ifndef __SV_VECTOR_FIELD_H
#define __SV_VECTOR_FIELD_H

#include "readinputsvl.h"
#include "MArray.h"

using namespace __svl_lib;

class svVectorField
{
public:
    svVectorField(svChar *inf) {
        
        readinputsvl(inf, position, orientation, &seedNum);
        
        // init physical space size
        x = 8;
        y = 8;
        z = 8;
        
        printf("%d\n",seedNum);
        printf("%d\n",position->size());
        
        
        // diverging color
        color0 = svVector3(0.230f, 0.299f, 0.754f);
        color1 = svVector3(0.865f,0.865f,0.865f);
        color2 = svVector3(0.706f, 0.016f, 0.150f);
        color10Diff = normalize(color1 - color0);
        color12Diff = normalize(color1 - color2);
    
        svVector3 test = color1 - color12Diff * 0.5f;
    
        printf("%f %f %f \n",color1[0],color1[1],color1[2]);
        printf("%f %f %f \n",color2[0],color2[1],color2[2]);
        printf("%f %f %f \n",test[0],test[1],test[2]);
    }
    
	virtual ~svVectorField(){};

    void Generate(){
     
    }
    
    void Render(){
        
        for (int i = 0; i < seedNum; i++) {
            for (int j = 0; j < position[i].size(); j++) {
              
                // get positions
                float posX,posY,posZ;
                position[i][j].getValue(posX,posY,posZ);
                
                // get direction vector
                svVector3 dirVec = orientation[i][j] - svVector3(0,0,0);
                svVector3 normDirVec = normalize(dirVec);
                float dirX,dirY,dirZ;
                float dirVecMag = dirVec.length();
                dirVec.getValue(dirX, dirY, dirZ);
                
                // calculate rotation axis and angle
                svVector3 up(0,0,1);
                svVector3 axisTo = dirVec;
                float rx,ry,rz;
                svVector3 rotationAxis = up.cross(axisTo);
                rotationAxis.getValue(rx, ry, rz);
                float angle = acos(up.dot(axisTo));
                angle *= (180/3.1415);
                
                glPushMatrix();
                glTranslatef(posX, posY, posZ);
                glRotatef(angle, rx, ry, rz);
                DrawArrowByOrient(angle);
                glPopMatrix();
            }   
        }
    }
    
    void DrawArrow(){
        
        // draw cylinder
        glColor3f(0, 0, 1);
        GLUquadricObj *quadratic;
        quadratic = gluNewQuadric();
        gluCylinder(quadratic,0.03f, 0.03f, 0.65f, 32, 32);
    
        // draw cone
        glPushMatrix();
        glColor3f(0, 0, 1);
        glTranslatef(0, 0, 0.65f);
        glutSolidCone(0.1f, 0.35f, 32, 32);
        glPopMatrix();
    }
    
    void DrawArrowByOrient(float deg){
        
        float c = cosf(deg);
        
        svVector3 color;
        
        // lower palette
        if(c < 0){
            color = color1 + color10Diff * c;
        }
        // upper palette
        else{
            color = color1 - color12Diff * c;
        
        }
        
        // draw cylinder
        glColor3f(color[0], color[1], color[2]);
        GLUquadricObj *quadratic;
        quadratic = gluNewQuadric();
        gluCylinder(quadratic,0.03f, 0.03f, 0.65f, 16, 16);
        
        // draw cone
        glPushMatrix();
        glColor3f(color[0], color[1], color[2]);
        glTranslatef(0, 0, 0.65f);
        glutSolidCone(0.1f, 0.35f, 16, 16);
        glPopMatrix();
    }
    
	void GetPhysicalDimension(GLfloat *xx, GLfloat *yy, GLfloat *zz)
	{
	  *xx=x; *yy=y; *zz=z;
	}

  private:
    GLfloat x;
    GLfloat y;
    GLfloat z;
    
    svVector3Array position[10];
    svVector3Array orientation[10];
    int seedNum;
    
    svVector3 color0;
    svVector3 color1;
    svVector3 color10Diff;
    svVector3 color2;
    svVector3 color12Diff;
};


#endif //__SV_VECTOR_FIELD_H
