#ifndef FUNC_H
#define FUNC_H

#include <stdio.h>
#include <utility>
using namespace std;

void trans_mat_p(int M,int N,double* mat,double* res);
void trans_mat_p_f(int M,int N,float* mat,double* res);

void mult_mat_p(int M,int N,double* mat1,double* mat2,double* res);

void mult_mat_p_spld(int N,double* mat1,double* mat2,double* res,double* r,double* D_i);
void mult_mat_p_f(int M,int N,double* mat1,float* mat2,double* res);

double len_mat(int size,int ind,double* mat);

void mult_mat_gen(int M1,int N1,int N2,double* mat1,double* mat2,double* res);
void print_mat(int M,int N,double* mat);
void print_mat_f(int M,int N,float* mat);

void mult_mat_gen_f(int M1,int N1,int N2,float* mat1,double* mat2,double* res);
void mult_mat_gen_f2(int M1,int N1,int N2,float* mat1,double* mat2,float* res);
void mult_mat_gen_f_f_d(int M1,int N1,int N2,float* mat1,float* mat2,double* res);
void mult_mat_gen_d_f_d(int M1,int N1,int N2,double* mat1,float* mat2,double* res);
void mult_p_mat_gen_f(int M1,int N1,int N2,float* mat1,double* mat2,double* res);
double copy_mat_change_copy(int M,int N,double* mat1,double* mat2,double* d_o);

void given_ab_pdf(double a,double b,double* gt,double* p);
void given_inv(double p,double* gt);

void copy_mat(int M,int N,double* mat1,double* mat2);
void copy_mat_f(int M,int N,float* mat1,double* mat2);
double con_copy_mat_change_copy(int N,double* D_o,double* D_i,double* r,double* E_i,double* E_f);
void qr_decomp_givens(int M,int N,double*a,double* q_fin,int boo);
void qr_decomp_givens_p(int M,int N,double* a,double* q_fin,int boo);
bool compare_fun(pair<double,int> P1,pair<double,int> P2);
double con_val(int M,int N,double* D_o,double* D_i);

double con_val_p(int N,double* D_o,double* D_i,double* r,double* E_i,double* E_f);
void sort_diag(int s,int s_v,double D[],double* sig_inv,double* sig,double* v_fin,double *E_f);
void opt_mult_mat(int M,int N,double* gt,double* A,int i,int j);

#endif



// func.cpp


#include "func.h"
#include <omp.h>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <utility>
#include<bits/stdc++.h> 

using namespace std;

void trans_mat_p_f(int M,int N,float* mat,double *res){ 

   #pragma omp parallel for schedule(runtime) 
    for(int i=0;i<N;i++){
   	 for(int j=0;j<M;j++){     
   	 	// res[i][j] = mat[j][i];
   		res[i*M + j] = mat[j*N +i]; 
     }
	}   
}
void trans_mat_p(int M,int N,double* mat,double *res){ 

   #pragma omp parallel for schedule(runtime) 
    for(int i=0;i<N;i++){
   	 for(int j=0;j<M;j++){
   	 	// res[i][j] = mat[j][i];
   		res[i*M + j] = mat[j*N +i]; 
     }
	}   
}
void mult_mat_gen_f2(int M1,int N1,int N2,float* mat1,double* mat2,float* res){
  
   double mat2_t[N1*N2];
   trans_mat_p(N1,N2,mat2,mat2_t);

 #pragma omp parallel for schedule(runtime) 
 // M1*N1 N1*N2
 for (int i1 = 0; i1 < M1; i1++)
 {
  for(int i=0;i<N2;i++){
    double sum = 0.0;
    for(int j=0;j<N1;j++){

       sum +=  mat1[i1*N1+j]*mat2_t[i*N1+j];
    }
    res[i1*N2+i] = sum;
    }
  }
}

void mult_mat_gen_f(int M1,int N1,int N2,float* mat1,double* mat2,double* res){
  
   double mat2_t[N1*N2];
   trans_mat_p(N1,N2,mat2,mat2_t);

  #pragma omp parallel for schedule(runtime) 
 // M1*N1 N1*N2
  for (int i1 = 0; i1 < M1; i1++)
  {
  for(int i=0;i<N2;i++){
    double sum = 0.0;
    for(int j=0;j<N1;j++){
       sum +=  mat1[i1*N1+j]*mat2_t[i*N1+j];
    }
    res[i1*N2+i] = sum;  
    }
  }
}
void mult_mat_gen_f_f_d(int M1,int N1,int N2,float* mat1,float* mat2,double* res){
  
   double mat2_t[N1*N2];
   trans_mat_p_f(N1,N2,mat2,mat2_t);

  #pragma omp parallel for schedule(runtime) 

 // M1*N1 N1*N2
 for (int i1 = 0; i1 < M1; i1++)
 {
    for(int i=0;i<N2;i++){
    double sum = 0.0;
    for(int j=0;j<N1;j++){
       sum +=  mat1[i1*N1+j]*mat2_t[i*N1+j];
    }
    res[i1*N2+i] = sum;
   }
  }
}
void mult_mat_gen_d_f_d(int M1,int N1,int N2,double* mat1,float* mat2,double* res){
  
   double mat2_t[N1*N2];
   trans_mat_p_f(N1,N2,mat2,mat2_t);

  #pragma omp parallel for schedule(runtime) 

 // M1*N1 N1*N2
 for (int i1 = 0; i1 < M1; i1++)
 {
  for(int i=0;i<N2;i++){
    double sum = 0.0;
    for(int j=0;j<N1;j++){
      //  sum +=  mat1[i1*N1+j]*mat2[j*N2+i];
       sum +=  mat1[i1*N1+j]*mat2_t[i*N1+j];

    }
    res[i1*N2+i] = sum; 
    }
  }
}
void mult_mat_gen(int M1,int N1,int N2,double* mat1,double* mat2,double* res){
  
 // M1*N1 N1*N2
 double mat2_t[N1*N2];
 trans_mat_p(N1,N2,mat2,mat2_t);
 #pragma omp parallel for schedule(runtime)
 for (int i1 = 0; i1 < M1; i1++)
 {
   for(int i=0;i<N2;i++)
   {
    double sum = 0.0;
    for(int j=0;j<N1;j++){
      //  sum +=  mat1[i1*N1+j]*mat2[j*N2+i];
      sum +=  mat1[i1*N1+j]*mat2_t[i*N1+j];
    }
    res[i1*N2+i] = sum;
    }
  }
}
void mult_mat_p_f(int M,int N,double* mat1,float* mat2,double* res){
 
  #pragma omp parallel for schedule(runtime)   // almost same with / without collapse
 // N*M   M*N  
	// for(int j=0;j<N;j++){
	  for(int i=0;i<N;i++){
	    for(int j=0;j <= i;j++){                     // major diff  interchanged i j 
	  	double sum =0;
    	for(int k=0;k<M;k++){
   	 	   // sum += mat1[i][k] * mat2[k][j];
   	 	sum += mat1[i*M +k] * mat1[j*M+k];   	 	     
   		// res[i*M + j] = mat[j*N +i]; 
    }
     // res[i][j] = sum;
     if(i==j){
        res[i*N+j] = sum;
     }
     else{
        res[i*N+j] = sum;
        res[j*N+i] = sum;
     }
   
	}

    }
}

void mult_mat_p_spld(int N,double* mat1,double* mat2,double* res,double* r,double* D_i){
 
  #pragma omp parallel   // almost same with / without collapse
 // N*M   M*N  
	// for(int j=0;j<N;j++){
    {
    #pragma omp for schedule(runtime)
	  // for(int i=0;i<N;i++){
    	for(int j=0;j<N;j++){     
         for(int i=0;i<N;i++){                // major diff  interchanged i j 
		  double sum =0;
      double sum2 =0;
   	 	for(int k=0;k<N;k++){
   	 	   // sum += mat1[i][k] * mat2[k][j];
   	 	   sum += mat1[i*N +k] * mat2[j*N+k];   
         sum2 += r[i*N +k] * mat2[j*N+k];
   		// res[i*M + j] = mat[j*N +i]; 
    	}
     // res[i][j] = sum;
    	res[i*N+j] = sum;
      D_i[i*N+j] = sum2;
	  }
  }
  }
}
void mult_mat_p(int M,int N,double* mat1,double* mat2,double* res){
 
   double mat2_t[M*N];
   trans_mat_p(M,N,mat2,mat2_t);

  #pragma omp parallel for collapse(2) schedule(runtime)   // almost same with / without collapse
 // N*M   M*N  
	for(int j=0;j<N;j++){
	  for(int i=0;i<N;i++){
	// for(int j=0;j<N;j++){                     // major diff  interchanged i j 
		  double sum =0;
   	 	for(int k=0;k<M;k++){
   	 	   // sum += mat1[i][k] * mat2[k][j];
   	 	   sum += mat1[i*M +k] * mat2_t[j*M+k];   	 	     
   		// res[i*M + j] = mat[j*N +i]; 
    	}
     // res[i][j] = sum;
    	res[i*N+j] = sum;
	  }
  }
}
double len_mat(int size,int ind,double* mat){

   double sum_sq=0;
    for (int j = 0; j <  size; j++)
		{
			// sum_sq += (a[ind][j]*a[ind][j]);
			sum_sq += (mat[ind+j]*mat[ind+j]);  // not ind*N since while passing converted
		}
    double ans = sqrt(sum_sq);
	return ans;
}

double len_mat_new(int size,int ind,double* mat){

  double sum_sq=0;
  for (int j = 0; j <  size; j++)
	{
			// sum_sq += (a[j][ind]*a[j][ind]);
		sum_sq += (mat[ind+j*size]*mat[ind+j*size]);  // not ind*N since while passing converted
	}

  double ans = sqrt(sum_sq);
	return ans;
}


void given_ab_pdf(double a,double b,double* gt,double *p){

  double r=0;
  double sig=0;
  double gamma = 0;

  if(fabs(b)<=1e-10){
    sig = 0;
    gamma=1;
    p[0] = 0;
  }
  else if(fabs(a)<=1e-10){
    sig =1;
    gamma = 0;
    p[0] =1;
  }
  else{
    if(fabs(b)>fabs(a)){
      r = -1*a/b;
      sig = (double)1/sqrt(1+r*r);  
      gamma = sig*r;
      p[0] = 2/gamma;
    }
    else{
      r= -1*b/a;
      gamma = (double)1/sqrt(1+r*r);
      sig = (gamma*r);
      p[0]=sig/2;
    }
  }
  gt[0] = gamma;
  // gt[1]=sig;
  gt[1] = -1*sig;
  gt[2] = sig;
  gt[3] = gamma;

}
void given_inv(double p,double* gt){
  double gamma;
  double sig;

  if(fabs(p)<=1e-10){
    gamma =1;
    sig =0;
  }
  else if(fabs(p-1)<1e-10){
    gamma =0;
    sig =1;
  }
  else if(fabs(p)>2){
    gamma = 2/p;
    sig = sqrt(1-(gamma*gamma));
  }
  else{
    sig = 2*p;
    gamma = sqrt(1-(sig*sig));
  }

  gt[0] = gamma;
  gt[1] = sig;
  gt[2] = -1*sig; //changed for serial part
  gt[3] = gamma;
}
void print_mat_f(int M,int N,float* mat){
  cout<< "printing mat "<< M << " * " << N << endl;
  for (int i = 0; i < M; ++i)
  {
    for (int j = 0; j < N; ++j)
    {
      cout<<" "<< mat[i*N+j];
    }
    cout<< " "<<endl;
  }
}
void print_mat(int M,int N,double* mat){
  cout<< "printing mat "<< M << " * " << N << endl;
  for (int i = 0; i < M; ++i)
  {
    for (int j = 0; j < N; ++j)
    {
      cout<<" "<< mat[i*N+j];
    }
    cout<< " "<<endl;
  }
}
void copy_mat_f(int M,int N,float* mat1,double* mat2){
  #pragma omp parallel for schedule(runtime)

  for (int i = 0; i < M*N; ++i)
  {
     mat1[i] = mat2[i] ;
  }
}
double copy_mat_change_copy(int M,int N,double* mat1,double* mat2,double* d_o){
  
  double ans =0;
  for (int i = 0; i < M; ++i)
  {
    for (int j = 0; j < N; ++j)
    {
      int k = i*N+j;
      mat1[k] = mat2[k];
      double diff = fabs(mat2[k]-d_o[k]);
      d_o[k] = mat1[k];
      
      if(ans<=diff){
        ans = diff;
      }
    }
  }
  return ans;
}
void copy_mat(int M,int N,double* mat1,double* mat2){

  #pragma omp parallel for schedule(runtime)
  for (int i = 0; i < M*N; ++i)
  {
     mat1[i] = mat2[i] ;
  }
}
void opt_mult_mat(int M,int col,double* gt,double*A,int m,int n){
      
      // double* A_dum1 = new double[col];
      // double* A_dum2 = new double[col];
      double A_dum1[col];
      double A_dum2[col];
       for (int i = n; i < col; ++i)
       {
         A_dum1[i] = gt[0]*A[(m-1)*col+i] + gt[1]* A[m*col+i];
         A_dum2[i] = gt[2]* A[(m-1)*col+i]+gt[3]*A[m*col+i] ;
       }
        for (int i = n; i < col; ++i)
       {
         A[(m-1)*col+i] = A_dum1[i];
         A[(m)*col+i] = A_dum2[i];
       }   
  


}

void qr_decomp_givens(int M,int N,double* a,double* q_fin,int boo){

  double* gt = new double[4];
  double* p = new double[1];
 
  for (int i = 0; i < N;i++)
  {
    for(int j = M-1; (j>i);j--){

          given_ab_pdf(a[i+(j-1)*N],a[i+j*N],gt,p);

          opt_mult_mat(M,N,gt,a,j,i);  
         
          a[j*N+i] = p[0];

      }
  }
  // exit(0);
  for(int i=0;i<M;i++){
    for(int j=0;j<M;j++)
    {
      if(i==j){
         q_fin[i*M+i] = 1;
      }
      else
      {
        q_fin[i*M+j] =0;
      }  
    }
  }
  for (int j = M-2; j >= 0 ; j--)
  {
    for (int i = j+1; i < M; ++i)
    {
      given_inv(a[i*M+j],gt);
      // gt[1]=-gt[1];
      opt_mult_mat(M,M,gt,q_fin,i,0);
      a[i*M+j] =0;
    }
  }
  delete p;
  delete gt;
      
}

void qr_decomp_givens_p(int M,int N,double* a,double* q_fin,int boo){

  int iteration[16]={0};
  int steps[16]={0};
  // double* q_fin_dum= new double[M*M]();
  // for(int i = 0; i < M; i++)
  // {
  //   q_fin_dum[i*M+i] = 1;
  //   /* code */
  // }

  #pragma omp parallel shared(iteration,steps)
  {
 
  int num_threads =  omp_get_num_threads();
  int tid = omp_get_thread_num();
  int tid_prereq = ((tid-1+num_threads)%num_threads);
  int notmaster = (tid_prereq < tid);
  int master =!notmaster;

   for (int i = tid; i < N-1;i=i+num_threads)
   {
        //  double gt[2];

    double gt[4];
    double p[1];
    for(int j = M-1; (j>=i+1);j--){
      int same_iter = (iteration[tid]==iteration[tid_prereq]);
      int not_same_iter = !(same_iter);
      // int not_same_iter = (iteration[tid_prereq]<iteration[tid]);
      int dependent = (steps[tid_prereq]< (steps[tid] + 2)); 

      while((notmaster&&same_iter&&dependent) || (master && not_same_iter &&dependent))
      { 
         same_iter = (iteration[tid]==iteration[tid_prereq]);
         not_same_iter = !(same_iter);
        //  not_same_iter = (iteration[tid_prereq]<iteration[tid]);
         dependent = (steps[tid_prereq]< (steps[tid]+2)); 
      }

      given_ab_pdf(a[i+(j-1)*N],a[i+j*N],gt,p);
   
      opt_mult_mat(M,N,gt,a,j,i);  
      opt_mult_mat(M,M,gt,q_fin,j,0); // for qt 

      // opt_mult_mat(M,M,gt,q_fin_dum,j,0);  // for q
             a[j*N+i] = 0;
      //  a[j*N+i] = p[0];

      steps[tid]++;
      
      }
      steps[tid]=0;
      iteration[tid]++;
    }
  }

  // trans_mat_p(M,M,q_fin_dum,q_fin);
  // delete q_fin_dum;

  //   double gt[4];
  //  for(int i=0;i<M;i++){
  //   for(int j=0;j<M;j++)
  //   {
  //     if(i==j){
  //        q_fin[i*M+i] = 1;
  //     }
  //     else
  //     {
  //       q_fin[i*M+j] =0;
  //     }  
  //   }
  // }
  // for (int j = M-2; j >= 0 ; j--)
  // {
  //   for (int i = j+1; i < M; ++i)
  //   {
  //     given_inv(a[i*M+j],gt);
  //     opt_mult_mat(M,M,gt,q_fin,i,0);
  //     a[i*M+j] =0;
  //   }
  // }

}

double con_val(int M,int N,double* D_o,double* D_i)
{
  double ans=0;
  for (int i = 0; i < M; i++)
  {
    for (int j = 0; j < N; j++)
    {
      double diff = fabs(D_o[i*N+j] - D_i[i*N+j]);
      if(ans <=diff){
        ans = diff;
      }/* code */
    }
    /* code */
  }
  return ans;
}
double con_copy_mat_change_copy(int N,double* D_o,double* D_i,double* r,double* E_i,double* E_f)
{
  double ans =0;
  for (int i = 0; i < N; ++i)
  {
    for (int j = 0; j < N; ++j)
    {
      int k= i*N+j;
     r[k] = D_i[k];
      double diff = fabs(D_o[k]-D_i[k]);
      D_o[k] = r[k];
     	E_i[k] = E_f[k];
      
      if(ans<=diff){
        ans = diff;
      }
    }
  }
 return ans;
}
double con_val_p(int N,double* D_o,double* D_i,double* r,double* E_i,double* E_f)
{
  double ans=-1e9;

  #pragma omp parallel
  {
  double ansl=-1e9;
  #pragma omp for nowait
  for (int i = 0; i < N; ++i)
  {   
    for (int j = 0; j < N; ++j)
    {
      int k=i*N+j;
      r[k] = D_i[k];
      double diff = fabs(D_o[k] - D_i[k]);
      D_o[k] = r[k];
      E_i[k] = E_f[k];
      if(ansl <= diff){
        ansl = diff;
      }
    }
    /* code */
  }
   #pragma omp critical   
   if(ans <= ansl){
      ans = ansl;
    }
  } 
 return ans;
  // return ans;
}

bool compare_fun(pair <double,int> P1, pair <double,int> P2){

  return P1.first>P2.first;
}
void sort_diag(int M,int N,double D[],double* sig_inv,double* sig,double* v_fin,double *E_f){

    pair <double,int> PAIR[N];
    
    #pragma omp parallel for 
    for (int i = 0; i < N; ++i)
     {
        D[i*N+i] = sqrt(fabs(D[i*N+i]));
    
        PAIR[i] = make_pair(D[i*N+i],i);

       /* code */
     } 

    sort(PAIR,PAIR+N,compare_fun);
    for (int i = 0; i < N; ++i)
    {

        sig[i*N+i] = PAIR[i].first;
        if(PAIR[i].first!=0){
        sig_inv[i*N+i] = 1/PAIR[i].first;
    }

    }
     for (int i = 0; i < N; ++i)
     {
      for (int j = 0; j < N; ++j)
      {       
        v_fin[j*N+i] = E_f[PAIR[i].second+N*j]; 

      }
              
     }

}


// end

// lab2_omp.c

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

// end 