#ifndef QuantumDyamics2D_h
#define QuantumDyamics2D_h

#include "WaveFunction2D.h"
#include <new>
#include <string>
#include <fstream>

//////////////////////
// 
// hbar = e = c = m = 1
//
class QuantumDynamics2D{
private:
   mydouble width; //  box width
   int size; // size = width/dx
   int size1; // size1 = size+1
   mydouble dt; // time discretization 
   mycomplex** A; // vector potential (A0,A1,A2)=(phi,Ax,Ay)
   mycomplex* V; // scalar potential
   mydouble q; // charge (-1 for electron)
   WaveFunction2D psi;
   WaveFunction2D psi_tmp;

   mycomplex* U_p;  // U_p = (1+I*H(t)*dt/2)
   mycomplex* U_m;  // U_m = (1-I*H(t)*dt/2)

   string ExportFileName;
public:
   QuantumDynamics2D(void){;};
   QuantumDynamics2D(int Size, mydouble Width);

   ~QuantumDynamics2D(){
      for(int i=0;i<3;i++){
         delete [] A[i];
      }
      delete [] A;
      delete [] V;

      delete [] U_p;
      delete [] U_m;
   }
   
   void SetExportFileName(string str){
      ExportFileName = str;
   }

   void SetDt(mydouble DT){
      dt=DT;
   }

   void SetWidth(mydouble w){
      width=w;
   }

   mycomplex& VectorPotential(int n, int i, int j){
      return (A[n%3][i+size1*j]);
   }

   mycomplex& ScalarPotential(int i, int j){
      return (V[i+size1*j]);
   }

   void SetCharge(mydouble Q){
      q=Q;
   }

   void SetGaussian(Vector2D<mydouble>& X, Vector2D<mydouble>& P, mydouble vX);

   //void CreateHamiltonian(void);

   void CreateEvolutionMatrices(void);

   void Evolve(void);

   void SetFixedBoundary(void);

   void ShowWaveFunction(void);
};


#endif
