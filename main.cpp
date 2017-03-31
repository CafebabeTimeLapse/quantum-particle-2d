#include "WaveFunction2D.h"
#include "QuantumDynamics2D.h"
#include <sstream>

string int2str(int i){
   stringstream ss;
   ss << i;
   return ss.str();
}

int main(int argc, char** argv){
   int      N     = 30;
   mydouble width = 10.;
   mydouble dt    = 0.001;
   string DataDirectory="./data/";
   string FileName;
   
   Vector2D<mydouble> x0(6.,5.);
   Vector2D<mydouble> p0(5.,0.);
   mydouble vX=1.5;

   ////////////////////////////////////////
   ////////////////////////////////////////
   // construct simulator class
   //   width: length of the box
   //   N:     # of meshes in {x,y}-axis 
   //   dt:    mesh size in t-axis 
   ////////////////////////////////////////
   QuantumDynamics2D q(N,width);
   q.SetDt(dt);

   ////////////////////////////////////////
   ////////////////////////////////////////
   // set Gaussian wavepacket,
   //   x0: center of the packet
   //   vX: width of the packet
   //   p0: momentum
   //
   ////////////////////////////////////////
   q.SetGaussian(x0,p0,vX);

   ////////////////////////////////////////
   ////////////////////////////////////////
   // set magnetic fields,
   //   Ax(x,y) = -B*y/2,
   //   Ay(x,y) =  B*x/2
   //
   ////////////////////////////////////////
   mydouble x,y,E,B;
   B=0.;
   E=0.;
   for(int i=0;i<N;i++){
      x=width/N*i;
      for(int j=0;j<N;j++){
         y=width/N*j;
         q.VectorPotential(0,i,j) =  E*x;
         q.VectorPotential(1,i,j) = -B*y/2;
         q.VectorPotential(2,i,j) =  B*x/2;
      }
   }

   ////////////////////////////////////////
   ////////////////////////////////////////
   // set scalar potential
   //
   ////////////////////////////////////////
   mydouble V0=0;
   for(int i=5*N/10;i<6*N/10;i++){
      x=width/N*i;
      for(int j=0;j<N;j++){
         q.ScalarPotential(i,j) = V0;
      }
   }

   ////////////////////////////////////////
   ////////////////////////////////////////
   // create hamiltonian
   //
   ////////////////////////////////////////
   q.CreateEvolutionMatrices();


   ////////////////////////////////////////
   ////////////////////////////////////////
   // time evolution 
   //
   ////////////////////////////////////////
   for(int i=0;i<100;i++){
      cout << i << endl;
      q.SetExportFileName(DataDirectory+int2str(i)+".dat");
      q.ShowWaveFunction();
      for(int j=0;j<10;j++){
         q.Evolve();
         q.SetFixedBoundary();
      }
   }

   return 0;
}
