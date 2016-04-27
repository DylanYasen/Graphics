
//********************************************
// CMSC: To be implemented

#ifndef __SV_VECTOR_FIELD_H
#define __SV_VECTOR_FIELD_H

using namespace __svl_lib;

class svVectorField
{
  public: 
    svVectorField(svChar *inf) {};
    virtual ~svVectorField(){};

    void Generate(){};
    void Render(){};
	void GetPhysicalDimension(GLfloat *xx, GLfloat *yy, GLfloat *zz)
	{
	  *xx=x; *yy=y; *zz=z;
	}

  private:
    GLfloat x, y, z; 
	// physical dimension of the vector field

	// CMSC: add your own private member here

};


#endif //__SV_VECTOR_FIELD_H
