
/**----------------------------------------------
 * Author: Jian Chen
 * Nov 2013
 */

#ifndef __SVL_MARRAY_H 
#define __SVL_MARRAY_H


#include <stdlib.h>
#include <stdio.h>

namespace __svl_lib {

template<class T> class MArray {
 private:
    T* objs;               // the actual array pointer
    int _size;             // the actual size of array , so, _size <= nalloc
    int nalloc;            // the number of cells have been allocated 
    int default_grow_size; // the growing size at a time 

 public:
    MArray(const MArray&);
    MArray(int size=0, int default_grow_size=500, int asize=-1);
    virtual ~MArray();

    MArray<T>& operator=(const MArray&);

    // Accesses the nth element of the array
    inline T& operator[](int n) {
        if (n<0 || n>=_size) {
	        printf("array access out of range: n = %d\n", n); 
	        exit(1); 
        }
        return objs[n];
    }

    // Returns the size of the array
    inline int size() const{ return _size;}
    inline int real_size() const{return nalloc; }

    // Returns the distance between two arrays
    /*
    inline T& dis(const MArray a) 
    {
      int asize = a.size(), bsize=b.size();
      float dis=0;
      if(asize != _size) 
      {
        printf("array sizes doesn't agree: %d, %d\n", asize, _size);
	exit(1);
      }
      for(int i=0; i<_size; i++)
      {
	dis += (objs[i]-a[i]) * (objs[i]-a[i]);
      }
      return (sqrt(dis));
    }
    */

    // Make the array larger by count elements
    void grow(int count, int grow_size=100);

    // Add one element to the array.  equivalent to:
    //  grow(1)
    //  array[array.size()-1]=data;
    void add(const T&);

    // Insert one element to the ith position in the array
    void insert(const T& obj, int idx);

    // Remove one element from the array.  This is very inefficient
    // if you remove anything besides the last element.
    void remove(int);

    // Remove the last element
    void remove();

    // Remove all elements in the array.  The array is not freed,
    // and the number of allocated elements remains the same.
    void remove_all();

    // Free the whole array
    void free();
};

//
// copy constructor
//
template<class T>
MArray<T>::MArray(const MArray<T>& a) 
{
  _size=a._size;
  nalloc=_size;
  objs=new T[_size];
  for(int i=0;i<_size;i++)objs[i]=a.objs[i];
  default_grow_size=a.default_grow_size;
}

template<class T>
MArray<T>::MArray(int size, int gs, int asize)
{
  if (size < 0) {
    printf(" can't declare a negative size array.\n"); 
    abort(); 
  }

  default_grow_size=gs;

  if(size){
    if(asize<size){
      objs=new T[size];
      _size=size;
      nalloc=_size;
    } else {
      objs=new T[asize];
      _size=size;
      nalloc=asize;
    }
  } else {
    if(asize==-1){
      objs=0;
      _size=0;
      nalloc=0;
    } else {
      objs=new T[asize];
      _size=0;
      nalloc=asize;
    }
  }

  nalloc=_size;
}	

template<class T>
MArray<T>::~MArray()
{
  if(objs) delete [] objs;
}

template<class T>
void MArray<T>::grow(int count, int grow_size)
{
  int newsize=_size+count;
  if(newsize>nalloc){
    // Reallocate...
    int gs1=newsize>>2;
    int gs=gs1>grow_size?gs1:grow_size;
    int newalloc=newsize+gs;
    T* newobjs=new T[newalloc];
    if(objs){
      for(int i=0;i<_size;i++){
	newobjs[i]=objs[i];
      }
      delete[] objs;
    }
    objs=newobjs;
    nalloc=newalloc;
  }
  _size=newsize;
}

template<class T>
void MArray<T>::add(const T& obj)
{
  grow(1, default_grow_size);
  objs[_size-1]=obj;
}

template<class T>
void MArray<T>::insert(const T& obj, int idx)
{
  if (idx>_size || idx<0){
    printf("The index is out of the array range. index = %d\n", idx);
    abort();
  }
  else if (idx==_size){
    add(obj);
  }
  else {
    grow(1, default_grow_size);
    for (int i=_size-1; i>idx; i--) objs[i] = objs[i-1];
    objs[idx] = obj;
  }
}

template<class T>
void MArray<T>::remove(int idx)
{
  _size--;
  for(int i=idx;i<_size;i++)objs[i]=objs[i+1];
}

template<class T>
void MArray<T>::remove()
{
  if (_size>0) _size --;
}

template<class T>
void MArray<T>::remove_all()
{
  _size=0;
}

template<class T>
void MArray<T>::free()
{
  if(objs) delete [] objs;
  objs = NULL;
  _size = 0;
  nalloc = 0;
}

//
// assignment operation for array
//
template<class T>
MArray<T>& MArray<T>::operator=(const MArray<T>& copy)
{
  if (objs) delete [] objs;
  _size=copy._size;
  nalloc=_size;
  objs=new T[_size];
  for(int i=0;i<_size;i++) objs[i]=copy.objs[i];
  default_grow_size=copy.default_grow_size;
  return(*this);
}

template class MArray<int>; 
template class MArray<float>;
template class MArray<double>;
typedef MArray<int> IntArray;
typedef MArray<float> FloatArray;
typedef MArray<double> DoubleArray;

}

#endif // __MARRAY_H
