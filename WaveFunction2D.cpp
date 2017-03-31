#include "WaveFunction2D.h"

void WaveFunction2D::Clear(void){
   for(int j=0;j<size1*size1;j++){
      psi[j]=0;
   }
}

void WaveFunction2D::SetGaussian(Vector2D<mydouble>& X, Vector2D<mydouble>& P, mydouble vX){
   mycomplex I(0.,1.);
   mydouble x, y, dx, dy;
   mycomplex psi_x, psi_y;
   mydouble vP;
   
   vP=1./vX; // vP = sqrt(2*<p^2>), vX = sqrt(2*<x^2>). vP*vX = 1

   for(int j=0;j<size1;j++){
      y=j*width/size;
      dy = y - X.y;
      psi_y = exp(-.5*vP*vP*dy*dy+I*P.y*dy);
      for(int i=0;i<size1;i++){
         x=i*width/size;
         dx = x - X.x;
         psi_x = exp(-.5*vP*vP*dx*dx+I*P.x*dx);
         psi[i+size1*j] = vP*I_SQRTPI*psi_y*psi_x;
      }
   }
}

mydouble WaveFunction2D::Norm(void){
   mydouble norm=0;
   mydouble c;
   int Ntot=size1*size1;
   for(int i=0;i<Ntot;i++){
      c=abs(psi[i]);
      norm+=c*c;
   }
   return norm*width*width/size/size;
}

void WaveFunction2D::ShowWaveFunction(void){
   mydouble x, y, p;
   int i, j;
   cout.setf(ios::scientific);

   for(i=0;i<size1;i++){
      x=i*width/size;
      for(j=0;j<size1;j++){
         y=j*width/size;
         p=abs(psi[i+size1*j]);
         cout << x << " " << y << " " << p*p << endl;
      }
      cout << endl;
   }
}
