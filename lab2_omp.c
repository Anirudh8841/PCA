#include <malloc.h>
#include <omp.h>
#include "func.h"
#include <stdio.h>
#include <iostream>

// /*
// 	*****************************************************
// 		TODO -- You must implement this function
// 	*****************************************************
// */
void SVD(int M, int N, float* D, float** U, float** SIGMA, float** V_T)
{

   double* D_Tras = new double[N*M];

   double* r = new double[N*N];

   trans_mat_p_f(M,N,D,D_Tras);

   mult_mat_p_f(M,N,D_Tras,D,r);

   delete[] D_Tras;

   double* E_i = new double[N*N](); 
   
	double* E_f = new double[N*N]();

   // double* q_fin = new double[N*N]();
   double* q_real = new double[N*N]();

   double* D_o = new double[N*N];

   double* D_i = new double[N*N]();  

   for(int i =0;i<N;i++){
    	E_i[i*N+i] = 1;
   }

   int boo=1;
   double change =1e6;
   double* D_check = new double[N*N](); 
   // double change[1] ={1e6};
   copy_mat(N,N,D_o,r);

   while(change>1e-6){
   // while(boo<100){

      double q_fin[N*N]; //= new double[N*N]();

      #pragma omp parallel
      {

      #pragma omp for schedule(runtime)
      for(int i =0;i<N;i++){

      for(int j =0;j<N;j++){
         if(i==j){
            q_fin[i*N+i] = 1;
         }
         else{
            q_fin[i*N+j] = 0;
         }
    	   
      }
      }
      }

      qr_decomp_givens_p(N,N,r,q_fin,boo); 

      mult_mat_p_spld(N,E_i,q_fin,E_f,r,D_i);
 
      change = con_val_p(N,D_o,D_i,r,E_i,E_f);  // both copy but parallel
      
      cout<< "change: "<<boo<<" " << change<<endl;

      boo++;
      // delete[] q_fin;
    }

   //  delete[] q_fin;
    delete[] E_i;
    delete[] D_o;
    delete[] D_i;
    
   double* v_fin = new double[N*N]; 
   double* sig_inv = new double[N*N]();
   double* sig_1 = new double[N*N]();
 
   sort_diag(M,N,r,sig_inv,sig_1,v_fin,E_f);

   delete[] r;
   delete[] E_f;

    double* u_dum = new double[M*N]; 
    double* u_fin = new double[M*N]; 

    mult_mat_gen_f(M,N,N,D,v_fin,u_dum);
    mult_mat_gen(M,N,N,u_dum,sig_inv,u_fin);
    delete[] u_dum;
    delete[] sig_inv;



      for (int j = 0; j < N; ++j)
      {
         SIGMA[0][j] =  sig_1[j*N+j];   /* code */
             
      }
   
   delete[] sig_1;
   
  
   
   // making it for Dt 

    
    double* v_t_dum = new double[N*M];
  

     trans_mat_p(M,N,u_fin,v_t_dum);
     delete[] u_fin;
   

     for (int i = 0; i < N; ++i)
     {
     	for (int j = 0; j < M; ++j)
     	{
     		V_T[0][i*M+j] = v_t_dum[i*M+j];
     		/* code */
     	}
     	/* code */
     }
    delete[] v_t_dum;
    copy_mat_f(N,N,*U,v_fin);
     
    delete[] v_fin;
 
  
  
   



}

// /*
// 	*****************************************************
// 		TODO -- You must implement this function
// 	*****************************************************
// */
void PCA(int retention, int M, int N, float* D, float* U, float* SIGMA, float** D_HAT, int *K)
{
    double* e_val = new double[N];
    double eigen_sum =0;

    for (int i = 0; i < N; ++i)
    {
    	e_val[i] = SIGMA[i]*SIGMA[i];
    	eigen_sum += e_val[i];
    	
    }
    
    if(eigen_sum==0){
    	printf("no eigen values");
    }

    for (int i = 0; i < N; ++i)
    {
    	e_val[i] = e_val[i] /eigen_sum;
    	/* code */
    }
     
   
    eigen_sum = 0;
    for (int i = 0; i < N ; ++i)
    {
       eigen_sum += e_val[i];
       if(100*eigen_sum >= retention){
       
       	K[0] = i+1;
       	// printf("val_in is k: %d",K[0]);
       	break;
       }
    	/* code */
    }
     delete[] e_val;
    double* w = new double[N*K[0]];

    for (int i = 0; i < N; ++i)
    {
    	for (int j = 0; j < K[0]; ++j)
    	{
    		w[i*K[0]+j] = U[i*M+j];
    		/* code */
    	}
    	/* code */
    }
    *D_HAT = new float[M*K[0]];
    mult_mat_gen_f2(M,N,K[0],D,w,*D_HAT);
    delete[] w; 



}
