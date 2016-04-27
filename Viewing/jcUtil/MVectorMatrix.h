
/**----------------------------------------------
 * Author: Jian Chen
 * Nov 2013
 */

#ifndef __SVL_M_VECTOR_MATRIX_H
#define __SVL_M_VECTOR_MATRIX_H

#include <iostream>
#include <math.h>
#include <assert.h>

typedef float Scalar;

#define  SCALAR_EPSILON Scalar(0.000001)

////
// vector base class
////


namespace __svl_lib {
template <int N, class T> 
class VBase{				
 protected:
  T _v[N];

 public:
  VBase(const T & t = T()) 
    { for(int i = 0; i < N; i++) _v[i] = t; }

  VBase(const T * t)
    { for(int i = 0; i < N; i++) _v[i] = t[i]; }

  VBase(const VBase<N-1, T> & v, const T & t = T(1)) 
    { for(int i = 0; i < N-1; i++) _v[i]=v[i]; _v[N-1] = t; }

  VBase(const VBase<N, T> & v) 
    { for(int i = 0; i < N; i++) _v[i] = v[i]; }

  const T * getValue( T * t = NULL ) const
    { if (t!=NULL) for(int i = 0; i < N; i++) t[i] = _v[i]; return _v; }

  void setValue( const T * t )
    { for(int i = 0; i < N; i++) _v[i] = t[i]; }

  int size() const 
    { return N; }
  
  void zero()
    { for(int i = 0; i < N; i++) _v[i] = T(0); }
  
  void negate()
    { for(int i = 0; i < N; i++) _v[i] = -_v[i]; }

  bool agree( const VBase<N,T> & v) const
    {
      for(int i = 0; i < N; i++) 
        if (T(fabs(_v[i]-v[i]))>SCALAR_EPSILON) return false;
      return true;
    }

  T length() const
    {
      T r = T(0);
      for(int i = 0; i < N; i++) r += _v[i]*_v[i]; 
      return T(sqrt(r));
    }	
  
  T length2() const
    {
      T r = T(0);
      for(int i = 0; i < N; i++) r += _v[i]*_v[i]; 
      return r;
    }	

  T dot( const VBase<N,T> & v ) const
    {
      T r = T(0);
      for(int i = 0; i < N; i++) r += _v[i]*v[i];
      return r;
    }

  void normalize() 
    { 
      T l;
      
      l = length();
      if (l > SCALAR_EPSILON){
	for(int i = 0; i < N; i++) 
	  _v[i] /= l;
      }
    }

  void print() const
    {
      //cout << "[ ";
      //for (int i=0; i < N; i++) cout << _v[i] << " ";
      //cout << "]" << endl;
    }

  // operator functions

  T operator () ( int i ) const
    { return _v[i]; }

  T & operator [] ( int i )
    { return _v[i]; }
  
  const T & operator [] ( int i ) const
    { return _v[i]; }

  VBase<N,T> & operator += ( const VBase<N,T> & v )
    { for(int i = 0; i < N; i++) _v[i] += v[i]; return (*this); }

  VBase<N,T> & operator -= ( const VBase<N,T> & v )
    { for(int i = 0; i < N; i++) _v[i] -= v[i]; return (*this); }

  VBase<N,T> & operator *= ( T t )
    { for(int i = 0; i < N; i++) _v[i] *= t; return (*this); }

  VBase<N,T> & operator *= ( const VBase<N,T> & v )
    { for(int i = 0; i < N; i++) _v[i] *= v[i]; return (*this);}
 
  VBase<N,T> & operator /= ( T t )
    {
      assert(t); 
      for(int i = 0; i < N; i++) _v[i] /= t;
      return (*this);
    }
  
  VBase<N,T> operator + () const
    { return VBase<N,T> (*this); }
    
  VBase<N,T> operator - () const
    { VBase<N,T> rv(*this); rv.negate(); return rv; }
  
  VBase<N,T> operator + ( const VBase<N,T> &u ) const
    { VBase<N,T> rv(*this); return rv += u; }
  
  VBase<N,T> operator - ( const VBase<N,T> &u ) const
    { VBase<N,T> rv(*this); return rv -= u; }
  
  VBase<N,T> operator * ( T t ) const
    { VBase<N,T> rv(*this); return rv *= t; }

  VBase<N,T> operator / ( T t ) const
    { VBase<N,T> rv(*this); return rv /= t; }
};
	
// VBase<N,T> friend operators
	
template <int N, class T> inline
  VBase<N,T> operator + ( const VBase<N,T> & v, const VBase<N,T> & u )
{
  VBase<N,T> rv(v);
  return rv += u;
}

template <int N, class T> inline
  VBase<N,T> operator - ( const VBase<N,T> & v, const VBase<N,T> & u )
{
  VBase<N,T> rv(v);
  return rv -= u;
}

template <int N, class T> inline
  VBase<N,T> operator * ( const VBase<N,T> & v, T t )
{
  VBase<N,T> rv(v);
  return rv *= t;
}

template <int N, class T> inline
  VBase<N,T> operator * ( T t, const VBase<N,T> & v )
{ 
  return v*t; 
}

template <int N, class T> inline
  VBase<N,T> operator * ( const VBase<N,T> & v, const VBase<N,T> & u )
{
  VBase<N,T> rv(v);
  return rv *= u;
}

template <int N, class T> inline
  VBase<N,T> operator / ( const VBase<N,T> & v, T t )
{ 
  VBase<N,T> rv(v); 
  return rv /= t; 
}

template <int N, class T> inline
  bool operator == ( const VBase<N,T> & v, const VBase<N,T> & u )
{
  for(int i = 0; i < N; i++){
    if(v[i] != u[i])  return false;
  }
  return true;
}

template <int N, class T> inline
  bool operator != ( const VBase<N,T> & v, const VBase<N,T> & u )
{  
  return !(v == u); 
}

template <int N, class T> inline
  T dot ( const VBase<N,T> & v, const VBase<N,T> & u )
{
  return v.dot(u);
}

template <int N, class T> inline
  VBase<N,T> normalize ( const VBase<N,T> & v )
{
  VBase<N,T> rv(v);
  rv.normalize();
  return rv;
}

////
// Vector3 class
////

class Vector3 : public VBase<3,Scalar>
{
 public:
  Vector3(const Scalar & t = Scalar()) : VBase<3,Scalar>(t) 
    {}
  Vector3(const Scalar * t ) : VBase<3,Scalar>(t) 
    {}
  Vector3(const VBase<2,Scalar> & v, Scalar t) : VBase<3,Scalar>(v, t)
    {}
  Vector3(const VBase<3,Scalar> & v) : VBase<3,Scalar>(v) 
    {}
  Vector3(Scalar x, Scalar y, Scalar z)  
    { _v[0] = x; _v[1] = y; _v[2] = z; }

  
  void getValue(Scalar & x, Scalar & y, Scalar & z) const
    { x = _v[0]; y = _v[1]; z = _v[2]; }
  
  const Scalar * getValue( Scalar * t = NULL ) const
    { if (t!=NULL) for(int i = 0; i < 3; i++) t[i] = _v[i]; return _v; }
  
  void setValue(const Scalar & x, const Scalar & y, const Scalar & z)
    { _v[0] = x; _v[1] = y; _v[2] = z; }

  Vector3 cross(const Vector3 & rhs ) const
    {
      Vector3 rv;
      rv[0] = _v[1]*rhs[2] - _v[2]*rhs[1];
      rv[1] = _v[2]*rhs[0] - _v[0]*rhs[2];
      rv[2] = _v[0]*rhs[1] - _v[1]*rhs[0];	
      return rv;
    }
};

inline Vector3 cross(const Vector3 & l, const Vector3 & r)
{
  return l.cross(r);
}

////
// Vector4 class
////

class Vector4 : public VBase<4,Scalar>
{
 public:
  Vector4(const Scalar & t = Scalar()) : VBase<4,Scalar>(t) 
    {}
  Vector4(const Scalar * t ) : VBase<4,Scalar>(t) 
    {}
  Vector4(const VBase<3,Scalar> & v, Scalar t) : VBase<4,Scalar>(v, t)
    {}
  Vector4(const VBase<4,Scalar> & v) : VBase<4,Scalar>(v) 
    {}
  Vector4(Scalar x, Scalar y, Scalar z, Scalar w)  
    { _v[0] = x; _v[1] = y; _v[2] = z; _v[3] = w; }
  
  void getValue(Scalar & x, Scalar & y, Scalar & z, Scalar & w) const
    { x = _v[0]; y = _v[1]; z = _v[2];  w = _v[3]; }
  
  const Scalar * getValue( Scalar * t = NULL ) const
    { if (t!=NULL) for(int i = 0; i < 4; i++) t[i] = _v[i]; return _v; }
  
  void setValue(const Scalar & x, const Scalar & y, const Scalar & z, 
		const Scalar & w)
    { _v[0] = x; _v[1] = y; _v[2] = z; _v[3] = w; }

  void homogenize()
    {  
       assert(_v[3] != 0.0); 
       _v[0] /= _v[3]; _v[1] /= _v[3]; _v[2] /= _v[3]; _v[3] = 1.0; 
    }
};

inline Vector3 homogenize(const Vector4 & v)
{
    assert(v[3] != 0.0); 
	return Vector3(v[0]/v[3], v[1]/v[3], v[2]/v[3]);
}

////
// MxN matrix base class
////

template <int M, int N, class T> 
class MBase{
 protected:
  VBase<N, T> _m[M];
  
 public:
  MBase(const T & t = T())
    { for(int i = 0; i < M; i++) _m[i] = VBase<N,T>(t); }

  MBase(const T * t)
    { setValue(t); }

  MBase(const MBase<M, N, T> & n)
    { for(int i = 0; i < M; i++) _m[i] = n[i]; }

  MBase(const VBase<N, T> * p)
    { for(int i = 0; i < M; i++) _m[i] = p[i]; }

  void getValue( T * rt ) const
    // return 1D array conforming with OpenGL matrix' element order
    {
      for(int i = 0; i < M; i++)
        for(int j = 0; j < N; j++)
          rt[j*M+i] = _m[i][j];
    }

  void setValue( const T * t)
    {
      for(int i = 0; i < M; i++)
        for(int j = 0; j < N; j++)
          _m[i][j] = t[j*M+i];
    }

  VBase<N,T> getRowVector( int i , T * t = NULL ) const
    { _m[i].getValue(t);  return _m[i];  }

  void setRowVector( int i, T * t )
    { _m[i].setValue(t); }

  void setRowVector( int i, const VBase<N,T> & v )
    { _m[i] = v; }

  VBase<M,T> getColVector( int j , T * t = NULL ) const
    { 
      VBase<M,T> rv;
      for( int i = 0; i < M; i++ ) rv[i] = _m[i][j];
      rv.getValue(t);
      return rv; 
    }

  void setColVector( int j, T * t )
    { for( int i = 0; i < M; i++) _m[i][j] = t[i]; }

  void setColVector( int j, const VBase<M,T> & v )
    { for( int i = 0; i < M; i++) _m[i][j] = v[i]; }
  
  int row () const
    { return M; }

  int col () const
    { return N; }

  void zero()
    { for(int i = 0; i < M; i++)  _m[i].zero(); }

  void negate()
    { for(int i = 0; i < M; i++)  _m[i].negate(); }

  bool agree(const MBase<M,N,T> & n)
    {
      for(int i = 0; i < M; i++) if (!_m[i].agree(n[i]))  return false;
      return true;
    }

  void print() 
    { for(int i = 0; i < M; i++) _m[i].print(); }

  // operator functions

  T operator () ( int i, int j ) const
    { return _m[i][j]; }

  VBase<N,T> operator () ( int i ) const
    // return row vector
    { return _m[i]; }

  VBase<N,T> & operator [] ( int i )
    // return row vector
    { return _m[i]; }

  const VBase<N,T> & operator [] ( int i ) const
    // return row vector
    { return _m[i]; }

  MBase<M,N,T> & operator = ( const MBase<M,N,T> & n )
    { for(int i = 0; i < M; i++)  _m[i] = n[i]; return (*this); }

  MBase<M,N,T> & operator += ( const MBase<M,N,T> & n )
    { for(int i = 0; i < M; i++)  _m[i] += n[i]; return (*this); }

  MBase<M,N,T> & operator -= ( const MBase<M,N,T> & n )
    { for(int i = 0; i < M; i++)  _m[i] -= n[i]; return (*this); }

  MBase<M,N,T> & operator *= ( T t )
    { for(int i = 0; i < M; i++)  _m[i] *= t; return (*this); }

  MBase<M,N,T> & operator /= ( T t )
    {
      assert(t != 0.0);
      for(int i = 0; i < M; i++) _m[i] /= t;
      return (*this);
    }

  MBase<M,N,T> operator + () const
    { return MBase<M,N,T> (*this); }

  MBase<M,N,T> operator - () const
    { MBase<M,N,T> rm(*this); rm.negate(); return rm; } 
  
  MBase<M,N,T> operator + ( const MBase<M,N,T> & n ) const
    { MBase<M,N,T> rm(*this); return rm += n; } 

  MBase<M,N,T> operator - ( const MBase<M,N,T> & n ) const
    { MBase<M,N,T> rm(*this); return rm -= n; } 

  MBase<M,N,T> operator * ( T t ) const
    { MBase<M,N,T> rm(*this); return rm *= t; }

  MBase<M,N,T> operator / ( T t ) const
    { MBase<M,N,T> rm(*this); return rm /= t; }
 
  VBase<M,T> operator * ( const VBase<N,T> & v ) const
    { 
      VBase<M,T> rv; 
      for (int i=0; i<M; i++) rv[i] = _m[i].dot(v);  
      return rv; 
    } 
};

// MBase<M,N,T> friend operators
	
template <int M, int N, class T> inline
  MBase<M,N,T> operator + ( const MBase<M,N,T> & m, const MBase<M,N,T> & n )
{ MBase<M,N,T> rm(m); return rm += n; }

template <int M, int N, class T> inline
  MBase<M,N,T> operator - ( const MBase<M,N,T> & m, const MBase<M,N,T> & n )
{ MBase<M,N,T> rm(m); return rm -= n; }

template <int M, int N, class T> inline
  MBase<M,N,T> operator * ( const MBase<M,N,T> & m, T t )
{ MBase<M,N,T> rm(m); return rm *= t; }

template <int M, int N, class T> inline
  MBase<M,N,T> operator * ( T t, const MBase<M,N,T> & m )
{ return m*t; }

template <int M, int N, class T> inline
  MBase<M,N,T> operator / ( const MBase<M,N,T> & m, T t )
{ MBase<M,N,T> rm(m); return rm /= t; }

template <int M, int N, class T> inline
  bool operator == ( const MBase<M,N,T> & m, const MBase<M,N,T> & n )
{
  for(int i = 0; i < M; i++){
    if(m[i] != n[i])  return false;
  }
  return true;
}

template <int M, int N, class T> inline
  bool operator != ( const MBase<M,N,T> & m, const MBase<M,N,T> & n )
{ return !(m == n); }

// transpose function

template <int M, int N, class T> inline
  MBase<N,M,T> transpose ( const MBase<M,N,T> & m )
{
  MBase<N,M,T> rm;
  for(int i = 0; i < M; i++) rm.setColVector(m.getRowVector());
  return rm;
}


////
// MxM matrix class
////
template <int M, class T> 
class Matrix: public MBase<M,M,T>
{
 public:
  Matrix() 
    { identity(); }
  Matrix(const T t) : MBase<M,M,T> (t)
    {}
  Matrix(const T * t) : MBase<M,M,T> (t) 
    {}
  Matrix(const Matrix<M, T> & n) : MBase<M,M,T> (n) 
    {}
  Matrix(const VBase<M, T> * p) : MBase<M,M,T> (p) 
    {}
  
  void identity()
  {
      for(int i = 0; i < M; i++)
      {
          this->_m[i].zero();
          this->_m[i][i] = T(1);
      }
  }
  
  void transpose()
    { 
      for(int i = 0; i < M; i++)
        for(int j = i; j < M; j++)
          { T t = this->_m[i][j]; this->_m[i][j] = this->_m[j][i]; this->_m[j][i] = t; }
    }

  Matrix<M,T> & operator = ( const MBase<M,M,T> & n )
    { for(int i = 0; i < M; i++)  this->_m[i] = n[i]; return (*this); }

  Matrix<M,T> & operator = ( const Matrix<M,T> & n )
    { for(int i = 0; i < M; i++)  this->_m[i] = n[i]; return (*this); }
  
  Matrix<M,T> & operator += ( const Matrix<M,T> & n )
    { for(int i = 0; i < M; i++)  this->_m[i] += n[i]; return (*this); }

  Matrix<M,T> & operator -= ( const Matrix<M,T> & n )
    { for(int i = 0; i < M; i++)  this->_m[i] -= n[i]; return (*this); }

  Matrix<M,T> & operator *= ( T t )
    { for(int i = 0; i < M; i++)  this->_m[i] *= t; return (*this); }

  Matrix<M,T> & operator /= ( T t )
    {
      if(t == 0) return (*this);
      for(int i = 0; i < M; i++) this->_m[i] /= t;
      return (*this);
    }

  Matrix<M,T> & operator *= ( const Matrix<M,T> & n )
    { 
      Matrix<M,T> l(*this);
      for(int i = 0; i < M; i++)  
        for(int j = 0; j < M; j++)
          this->_m[i][j] = l[i].dot(n.getColVector(j));
      return (*this); 
    }

  Matrix<M,T> operator + () const
    { return Matrix<M,T> (*this); }

  Matrix<M,T> operator - () const
    { Matrix<M,T> rm(*this); rm.negate(); return rm; }

  Matrix<M,T> operator + ( const Matrix<M,T> & n ) const
    { Matrix<M,T> rm(*this); return rm += n; }
  
  Matrix<M,T> operator - ( const Matrix<M,T> & n ) const
    { Matrix<M,T> rm(*this); return rm -= n; }

  Matrix<M,T> operator * ( T t ) const
    { Matrix<M,T> rm(*this); return rm *= t; } const

  Matrix<M,T> operator / ( T t ) const
    { Matrix<M,T> rm(*this); return rm /= t; }

  VBase<M,T> operator * ( const VBase<M,T> & v ) const
    { 
      VBase<M,T> rv; 
      for(int i = 0; i < M; i++) rv[i] = this->_m[i].dot(v);
      return rv; 
    }

  Matrix<M,T> operator * ( const Matrix<M,T> & n ) const
    { Matrix<M,T> rm(*this);  return rm *= n; }
};

// Matrix<M,T> friend operator

template <int M, class T> inline
  Matrix<M,T> operator + ( const Matrix<M,T> & m, const Matrix<M,T> & n )
{ Matrix<M,T> rm(*m); return rm += n; }

template <int M, class T> inline
  Matrix<M,T> operator - ( const Matrix<M,T> & m, const Matrix<M,T> & n )
{ Matrix<M,T> rm(*m); return rm -= n; }

template <int M, class T> inline
  Matrix<M,T> operator * ( const Matrix<M,T> & m, T t )
{ Matrix<M,T> rm(*m); return rm *= t; }

template <int M, class T> inline
  Matrix<M,T> operator * ( T t, const Matrix<M,T> & m )
{ return m*t; }

template <int M, class T> inline
  Matrix<M,T> operator / ( const Matrix<M,T> & m, T t )
{ Matrix<M,T> rm(*m); return rm /= t; }

//
// inverse function 
//
template <int M, class T>
  Matrix<M,T> inverse( const Matrix<M,T> & m )
{
  register short i,j,k;
  double work_mat[M][M+M];
  double *r[M];
  Matrix<M,T> result(0.0);
		
  //
  // setup a work matrix of Mx2M, the left part (MxM) is input matrix, 
  // the MxM right part as an identity matrix
  //  3x6 workmat:
  //                [ m0 m3 m6   1 0 0 ] 
  //                [ m1 m4 m7   0 1 0 ]
  //                [ m2 m5 m8   0 0 1 ]
  //
  for(i=0; i<M; i++) {
    for(j=0; j<M; j++) {
      work_mat[i][j] = m[i][j];
      work_mat[i][j+M] = 0.0;
    }
    work_mat[i][i+M] = 1.0;
  }

  // setup pointer for each row
  for(i=0; i<M; i++){
    r[i] = &work_mat[i][0];
  }

  // find out the element with the largest absolute value for each row
  double row_max[M];
  for(i=0; i<M; i++) {
    row_max[i] = fabs(r[i][0]);
    for(j=1; j<M; j++)
      if (fabs(r[i][j]) > row_max[i]) row_max[i] = fabs(r[i][j]);

    // singular matrix if containing all-0s row
    if(row_max[i] == 0.0) return result; 
  }
		
  //
  // transform the left part into an upper triangle matrix using Gaussian
  // elimination.
  //     [ m0 m3 m6   1 0 0 ]         [ m0 m3 m6   r0 r3 r6 ]
  //     [ m1 m4 m7   0 1 0 ]   ==>   [ 0  m4 m7   r1 r4 r7 ]
  //     [ m2 m5 m8   0 0 1 ]         [ 0  0  m8   r2 r5 r8 ]
  //           
  for(i=0; i<M; i++) {

    // find target row for top
    int    target = i;
    double col_row_max = fabs(r[i][i]/row_max[i]);
    for(j=i+1; j<M; j++)
      if (fabs(r[j][i]/row_max[j]) > col_row_max){
        target = j; 
	col_row_max = fabs(r[j][i]/row_max[j]); 
      }

    // swap the row with the largest col_row_max to the current work row
    if (target != i) {
      double *tmp = r[i];
      r[i] = r[target];
      r[target] = tmp;

      double tmp2 = row_max[i];
      row_max[i] = row_max[target];
      row_max[target] = tmp2;
    }
				
    // execute gaussian elimination row by row
    for(j = i+1; j < M; j++) {
      double factor = r[j][i]/r[i][i];
      r[j][i] = 0.0;
      for(k = i+1; k < M<<1; k++)
	r[j][k] -= r[i][k] * factor;
    }
  }

  // singular if the rank is less than M
  if(r[M-1][M-1] == 0.0) return result; 
		
  //
  // transfrom the left upper triangle into identity, the right part is
  // the inverse what we need
  //    [ m0 m3 m6   r0 r3 r6 ]         [ 1 0 0   r0 r3 r6 ]
  //    [ 0  m4 m7   r1 r4 r7 ]   ==>   [ 0 1 0   r1 r4 r7 ]
  //    [ 0  0  m8   r2 r5 r8 ]         [ 0 0 1   r2 r5 r8 ]
  //
  for(i = M-1; i > 0; --i) {
    for(j = i-1; j >= 0; --j) {
      double factor = r[j][i]/r[i][i];
      for(k = j+1; k < M<<1; k++)
        r[j][k] -= factor * r[i][k];
    }
  }
		
  // set the resulting inverse matrix
  for(i=0; i<M; i++)
    for(j=0; j<M; j++)
      result[i][j] = T(r[i][j+M] / r[i][i]);
			
  return result;
}

//
// safe inverse, verfying the resulting inverse matrix
//
template <int M, class T>
  bool safe_inverse(Matrix<M,T> & m)
{
  Matrix<M,T> inv;
  inv = inverse(m);
  
  Matrix<M,T> id;
  Matrix<M,T> l, r;

  l = inv * m;
  r = m * inv;
 
  if (l.agree(id) && r.agree(id)) {
    m = inv;
    return true;
  }
  else{
    return false;
  }
}

////
// Matrix3 class
////

class Matrix3 : public Matrix<3,Scalar> {
 public:
  Matrix3() : Matrix<3,Scalar> () 
    {};
  Matrix3(const Scalar t) : Matrix<3,Scalar> (t)
    {};
  Matrix3(const Scalar * t) : Matrix<3,Scalar> (t)
    {}
  Matrix3(const Matrix<3, Scalar> & n) : Matrix<3,Scalar> (n)
    {}
  Matrix3(const VBase<3,Scalar> * p) : Matrix<3,Scalar> (p)
    {}
  Matrix3(const VBase<3,Scalar> &v0, const VBase<3,Scalar> &v1,
          const VBase<3,Scalar> &v2)
    {  _m[0] = v0; _m[1] = v1; _m[2] = v2;  }

  Matrix3(const Scalar t0, const Scalar t1, const Scalar t2,
	  const Scalar t3, const Scalar t4, const Scalar t5,
	  const Scalar t6, const Scalar t7, const Scalar t8)
    { 
       _m[0][0] = t0; _m[0][1] = t3; _m[0][2] = t6; 
       _m[1][0] = t1; _m[1][1] = t4; _m[1][2] = t7; 
       _m[2][0] = t2; _m[2][1] = t5; _m[2][2] = t8; 
    }
  double det()
  {
   return (_m[0][0]* (_m[1][1]*_m[2][2]-_m[2][1]*_m[1][2]) - 
           _m[0][1]* (_m[1][0]*_m[2][2]-_m[2][0]*_m[1][2]) +
           _m[0][2]* (_m[1][0]*_m[2][1]-_m[2][0]*_m[1][1]));
	  };
  Matrix3  inv()
  {
   return inverse(Matrix<3,Scalar>(_m));
  }
};

////
// Matrix4 class
////

class Matrix4 : public Matrix<4,Scalar> {
 public:
  Matrix4() : Matrix<4,Scalar> ()
    {};
  Matrix4(const Scalar t) : Matrix<4,Scalar> (t)
    {};
  Matrix4(const Scalar * t) : Matrix<4,Scalar> (t)
    {}
  Matrix4(const Matrix<4, Scalar> & n) : Matrix<4,Scalar> (n)
    {}
  Matrix4(const VBase<4,Scalar> * p) : Matrix<4,Scalar> (p)
    {}
  Matrix4(const VBase<4,Scalar> &v0, const VBase<4,Scalar> &v1, 
          const VBase<4,Scalar> &v2, const VBase<4,Scalar> &v3) 
    {  _m[0] = v0; _m[1] = v1; _m[2] = v2; _m[3] = v3;  }

  Matrix4(const Scalar t0, const Scalar t1, const Scalar t2, const Scalar t3,
          const Scalar t4, const Scalar t5, const Scalar t6, const Scalar t7,
          const Scalar t8, const Scalar t9, const Scalar t10, const Scalar t11,
          const Scalar t12, const Scalar t13, const Scalar t14,const Scalar t15)
    { 
       _m[0][0] = t0; _m[0][1] = t4; _m[0][2] = t8; _m[0][3] = t12;
       _m[1][0] = t1; _m[1][1] = t5; _m[1][2] = t9; _m[1][3] = t13;
       _m[2][0] = t2; _m[2][1] = t6; _m[2][2] = t10; _m[2][3] = t14;
       _m[3][0] = t3; _m[3][1] = t7; _m[3][2] = t11; _m[3][3] = t15;  
    }
};

typedef Vector3 Vector3f;
typedef Vector4 Vector4f;
typedef Matrix3 Matrix3f;
typedef Matrix4 Matrix4f;
}

#endif // __M_VECTOR_MATRIX_H
