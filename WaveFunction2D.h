#ifndef WaveFunction2D_h
#define WaveFunction2D_h

#include <iostream>
#include <complex>
#include <cmath>
#include <new>
#include "Constants.h"

using namespace std;
typedef double mydouble;
typedef complex<mydouble> mycomplex;

template <typename T>
class Vector2D{
public:
   T x;
   T y;

   Vector2D(){;}

   Vector2D(T X, T Y){
      x=X;
      y=Y;
   }

   void set(T X, T Y){
      x=X;
      y=Y;
   }
   
   ~Vector2D(){
   }
};

class WaveFunction2D{
private:
   mycomplex* psi;
   int size;
   int size1;
   mydouble width;
public:
   WaveFunction2D(){;}

   WaveFunction2D(int Size, mydouble Width){
      size=Size;
      size1=size+1;
      width=Width;
      psi = new mycomplex [size1*size1];
   }

   ~WaveFunction2D(){
      delete [] psi;
   }

   mycomplex& operator()(int i, int j){
      return psi[i+size1*j];
   }

   mycomplex* ArrayPointer(void){
      return psi;
   }

   void Clear(void);

   void SetGaussian(Vector2D<mydouble>& X, Vector2D<mydouble>& P, mydouble vX);

   mydouble Norm(void);

   void ShowWaveFunction(void);
};

#endif
