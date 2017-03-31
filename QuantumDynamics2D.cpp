#include "QuantumDynamics2D.h"

extern "C"{
   void zgemv(const char& trans, const int& m, const int& n, 
     const mycomplex* alpha, const mycomplex* A, const int& ldA, const mycomplex* x, const int& incx, 
     const mycomplex* beta , const mycomplex* y, const int& incy);

   void zgesv(const int& N, const int& NRHS,
     const mycomplex* A, const int& LDA, const int* IPIV,
     const mycomplex* B, const int& LDB, int& INFO );
}

////////////////////////////////////////////////
////////////////////////////////////////////////
// Calculate w = M*v 
//   M  : n*n matrix
//  v,w : n*1 vector
//
void MatrixVectorProduct(int n, mycomplex* M, mycomplex* v, mycomplex* w){
   int One=1;
   mycomplex alpha(1.,0);
   mycomplex beta(0.,0);
   zgemv('N', n, n, &alpha, M, n, v, One, &beta, w, One);
}

////////////////////////////////////////////////
////////////////////////////////////////////////
// Solve M*v=w for v
//   M  : n*n matrix
//  v,w : n*1 vector
//
void SolveLinearEquations(int n, mycomplex* M, mycomplex* v, mycomplex* w){
   int info, *ipiv;
   ipiv = new int [n];
   for(int i=0;i<n;i++){
      v[i] = w[i];
   }
   zgesv( n, 1, M, n, ipiv, v, n, info );

   delete [] ipiv;
}

void QuantumDynamics2D::SetGaussian(Vector2D<mydouble>& X, Vector2D<mydouble>& P, mydouble vX){
   psi.SetGaussian(X,P,vX);
}

QuantumDynamics2D::QuantumDynamics2D(int Size, mydouble Width){
   width = Width;
   size = Size;
   size1 = size+1;
   dt=0.1;
   ExportFileName="default.dat";

   ///////////////////////////////////////////
   ///////////////////////////////////////////
   // Wave function variables
   //
   new(&psi) WaveFunction2D(size,width);
   new(&psi_tmp) WaveFunction2D(size,width);
  
   ///////////////////////////////////////////
   ///////////////////////////////////////////
   // potential variables
   //
   A = new mycomplex* [3];
   for(int i=0;i<3;i++){
      A[i] = new mycomplex [size1*size1];
   }
   V = new mycomplex [size1*size1];

   ///////////////////////////////////////////
   ///////////////////////////////////////////
   // charge 
   //
   q=-1.;

   ///////////////////////////////////////////
   ///////////////////////////////////////////
   // time evolution operator 
   //
   U_p = new mycomplex [size1*size1*size1*size1];
   U_m = new mycomplex [size1*size1*size1*size1];
}

////////////////////////////////////////
////////////////////////////////////////
// Fix the boundary of psi: 
//   psi(0,i) = psi(size,i) = psi(i,0) = psi(i,size) = 0
//   for i = 0,1,...,size
//
////////////////////////////////////////
void QuantumDynamics2D::SetFixedBoundary(){
   for(int i=0;i<size1;i++){
      psi(i,0)    = 0.;
      psi(i,size) = 0.;
      psi(0,i)    = 0.;
      psi(size,i) = 0.;
   }
}

////////////////////////////////////////
////////////////////////////////////////
// create Caylay operator
//   U_p^(-1) * U_m,
//
// where
//   U_p = 1+I*H(t)*dt/2,
//   U_m = 1-I*H(t)*dt/2,
//
// and the Hamiltonian H(t) given by
//   H_p2   = P^2/(2m)            // obi-matrix
//   H_Ap   = -2q(A*p)/(2m)       // obi-matrix
//   H_A2   = (qA)^2/(2m)         // diagonal
//   H_divA = iq(hbar)(divA)/(2m) // diagonal
//   H_phi  = q(phi)              // diagonal
//   H_V    = V                   // diagonal
//
void QuantumDynamics2D::CreateEvolutionMatrices(){
   int dummy,s2,diag;
   s2=size1*size1;
   mydouble idx=size/width; // idx = 1/(dx)
   mydouble idx2=size*size/width/width; // idx2 = 1/(dx)^2
   mycomplex mesh(0.,0.5*dt); // mesh = I*dt/2
   mycomplex I(0.,1.);
   mycomplex cdummy;
   mycomplex Ax, Ay; // A{x,y} = differential of A{x,y} w.r.t {x,y}
   for(int i=0;i<size1;i++){
      for(int j=0;j<size1;j++){
         dummy=i+size1*j;
         diag = dummy+s2*dummy;

         ////////////////////////////////////////
         ////////////////////////////////////////
         // 1 (s2*s2 identity matrix)
         //
         U_p[diag] = 1.;
         U_m[diag] = 1.;
   
   
         ////////////////////////////////////////
         ////////////////////////////////////////
         // H_p2
         //
         cdummy = mesh*2.0*idx2;
         U_p[diag] += cdummy;
         U_m[diag] -= cdummy;
         if(i!=0){
            cdummy = -mesh*0.5*idx2; 
            U_p[dummy+s2*(dummy-1) ] += cdummy;
            U_m[dummy+s2*(dummy-1) ] -= cdummy;
         }
         if(i!=size1-1){
            cdummy = -mesh*0.5*idx2; 
            U_p[dummy+s2*(dummy+1) ] += cdummy;
            U_m[dummy+s2*(dummy+1) ] -= cdummy;
         }
         if(j!=0){
            cdummy = -mesh*0.5*idx2; 
            U_p[dummy+s2*(dummy-size1) ] += cdummy;
            U_m[dummy+s2*(dummy-size1) ] -= cdummy;
         }
         if(j!=size1-1){
            cdummy = -mesh*0.5*idx2; 
            U_p[dummy+s2*(dummy+size1) ] += cdummy;
            U_m[dummy+s2*(dummy+size1) ] -= cdummy;
         }

         ////////////////////////////////////////
         ////////////////////////////////////////
         // H_A2
         //
         cdummy = mesh*0.5*q*q*( A[1][dummy]*A[1][dummy] + A[2][dummy]*A[2][dummy]);
         U_p[diag]+=cdummy;
         U_m[diag]-=cdummy;
         
         // differential of Ax
         if(i==0){
            Ax = (-1.5*A[1][dummy]+2.0*A[1][dummy+1]-0.5*A[1][dummy+2])*idx;
         }else if(i==size1-1){
            Ax = (-1.5*A[1][dummy]+2.0*A[1][dummy-1]-0.5*A[1][dummy-2])*idx;
         }else{
            Ax = 0.5*(A[1][dummy+1]-A[1][dummy-1])*idx;
         } 
         // differential of Ay
         if(j==0){
            Ay = (-1.5*A[2][dummy]+2.0*A[2][dummy+size1]-0.5*A[2][dummy+2*size1])*idx;
         }else if(j==size1-1){
            Ay = (-1.5*A[2][dummy]+2.0*A[2][dummy-size1]-0.5*A[2][dummy-2*size1])*idx;
         }else{
            Ay = 0.5*(A[2][dummy+size1]-A[2][dummy-size1])*idx;
         } 
         
         ////////////////////////////////////////
         ////////////////////////////////////////
         // H_divA
         //
         cdummy = mesh*0.5*I*q*(Ax+Ay);
         U_p[diag] += cdummy; 
         U_m[diag] -= cdummy; 
        
         ////////////////////////////////////////
         ////////////////////////////////////////
         // H_phi
         //
         cdummy = mesh*q*A[0][dummy];
         U_p[diag] += cdummy; 
         U_m[diag] -= cdummy;

         ////////////////////////////////////////
         ////////////////////////////////////////
         // H_V
         //
         cdummy = mesh*V[dummy];
         U_p[diag] += cdummy; 
         U_m[diag] -= cdummy;


         ////////////////////////////////////////
         ////////////////////////////////////////
         // H_Ap
         //
         if(i!=size1-1){
            cdummy = mesh*I*q*A[1][dummy]*0.5*idx;
            U_p[dummy+s2*(dummy+1)] += cdummy;
            U_m[dummy+s2*(dummy+1)] -= cdummy;
         }
         if(i!=0){
            cdummy = -mesh*I*q*A[1][dummy]*0.5*idx;
            U_p[dummy+s2*(dummy-1)] += cdummy; 
            U_m[dummy+s2*(dummy-1)] -= cdummy; 
         }
         if(j!=size1-1){
            cdummy = mesh*I*q*A[2][dummy]*0.5*idx;
            U_p[dummy+s2*(dummy+size1)] += cdummy;
            U_m[dummy+s2*(dummy+size1)] -= cdummy;
         }
         if(j!=0){
            cdummy = -mesh*I*q*A[2][dummy]*0.5*idx;
            U_p[dummy+s2*(dummy-size1)] += cdummy;
            U_m[dummy+s2*(dummy-size1)] -= cdummy;
         }

      }
   }
}

////////////////////////////////////////
////////////////////////////////////////
// time evolution 
//   psi(t+dt) = exp( -I*H(t)*dt ) * psi(t) 
//            ~=   ( 1+I*H(t)*dt/2 )^{-1}
//               * ( 1-I*H(t)*dt/2 )
//               * psi(t)
////////////////////////////////////////
void QuantumDynamics2D::Evolve(){
   int s2 = size1*size1;
   MatrixVectorProduct(s2, U_m, psi.ArrayPointer(), psi_tmp.ArrayPointer());
   //MatrixVectorProduct(s2, U_m, psi_tmp_diag.ArrayPointer(), psi.ArrayPointer());
   SolveLinearEquations(s2, U_p, psi.ArrayPointer(), psi_tmp.ArrayPointer());
   double snorm = sqrt(psi.Norm());
   for(int i=0;i<size1;i++){
      for(int j=0;j<size1;j++){
         psi(i,j) /=snorm;
      }
   }
}


void QuantumDynamics2D::ShowWaveFunction(){
   mydouble x, y, p;
   int i, j;
   ofstream ofs;
   ofs.open(ExportFileName.c_str());
   ofs.setf(ios::scientific);
   for(i=0;i<size1;i++){
      x=i*width/size;
      for(j=0;j<size1;j++){
         y=j*width/size;
         p=abs(psi(i,j));
         ofs << x << " " << y << " " << p*p << endl;
      }
      ofs << endl;
   }

   ofs.close();
}
