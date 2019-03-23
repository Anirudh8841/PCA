#include "lab2_io.h"
#include "lab2_omp.h"

#include <stdlib.h>
#include <omp.h>
#include <iostream>
#include <stdio.h> 
#include "func.h" 
  
using namespace std; 

int main(int argc, char const *argv[])
{
	if (argc < 2){
		printf("\nLess Arguments\n");
		return 0;
	}

	if (argc > 2){
		printf("\nTOO many Arguments\n");
		return 0;
	}

	//---------------------------------------------------------------------
	int M;			//no of rows (samples) in input matrix D (input)
	int N;
    float* D;		//1D array of M x N matrix to be reduced (input)
	float* U;

	double start_time1,start_time, end_time,end_time1;
	double computation_time,computation_time1;

	omp_set_num_threads(4); 
	read_matrix (argv[1], &M, &N, &D);
	double D_Tras_p[M*N];
  	// D_Tras_p = (double*) malloc(sizeof(double) * M*N);
	double D_t_d[M*N];
	double D_t_d1[M*N];
  	// D_t_d = (double*) malloc(sizeof(double) * M*N);
    trans_mat_p_f(M,N,D,D_Tras_p);
    start_time1 = omp_get_wtime();

// mult_mat(M,N,D_Tras_p,D,res);
// mult_mat_gen(M,N,M,D,D_Tras_p,res);
//  double change2 = con_val_p(M,N,D,D_Tras_p);
	mult_mat_p_f(M,N,D_Tras_p,D,D_t_d);
	end_time1 = omp_get_wtime();

	computation_time1 = ( (end_time1 - start_time1));
//    cout<< " printing Dt1"<<endl;
// 	print_mat(N,M,D_t_d);
	
		// cout<< " for serial change"<< change2 <<endl;
	cout<< " for serial "<< computation_time1 <<endl;
	start_time = omp_get_wtime();
	// double change = con_val(M,N,D,D_Tras_p);
	
	mult_mat_gen_d_f_d(N,M,N,D_Tras_p,D,D_t_d1);

	end_time = omp_get_wtime();
	computation_time = ( (end_time - start_time));
	//  cout<< " printing Dt2"<<endl;
	// print_mat(N,M,D_t_d1);
	 double change  = con_val(N,M,D_t_d1,D_t_d);
	cout<< " change "<< change <<endl;
	cout<< " for parallel "<< computation_time <<endl;

	return 0;
}

// #include <malloc.h>
// #include <omp.h>
// #include "func.h"
// #include <stdio.h>
// #include <iostream>

// // /*
// // 	*****************************************************
// // 		TODO -- You must implement this function
// // 	*****************************************************
// // */
// void SVD(int M, int N, float* D, float** U, float** SIGMA, float** V_T)
// {

//    double* D_Tras = new double[N*M];

//    double* r = new double[N*N];

//    trans_mat_p_f(M,N,D,D_Tras);

//    mult_mat_p_f(M,N,D_Tras,D,r);

//    delete[] D_Tras;

//    double* E_i = new double[N*N](); 
   
// 	double* E_f = new double[N*N]();

//    double* q_fin = new double[N*N]();

//    double* D_o = new double[N*N];

//    double* D_i = new double[N*N]();  

//    for(int i =0;i<N;i++){
//     	E_i[i*N+i] = 1;
//    }

//    int boo=1;
//    double change =1e6;
//    double* D_check = new double[N*N](); 
//    // double change[1] ={1e6};
//    copy_mat(N,N,D_o,r);

//    while(change>1e-6){
//    // while(boo<2000){

//       // copy_mat(N,N,D_o,r);

//       // qr_decomp_givens(N,N,r,q_fin,boo);

//       qr_decomp_givens_p(N,N,r,q_fin,boo);

//    // cout<<"printing q "<<endl;
//    // print_mat(N,N,q_fin);

//    // cout<<"printing r "<<endl;
//    // print_mat(N,N,r);

//    //    mult_mat_p(N,N,q_fin,r,D_check);
//    //    // // //  cout<<"printing D_check "<<endl;
//    //    // // //  print_mat(N,N,D_check);
      
//    //    double change2 = con_val(N,N,D_check,D_o);
//    //    cout<< "change_QR : " << change2<<endl;
   
//       mult_mat_p_spl(N,r,q_fin,D_i);
//       // change = con_val_spl(N,D_o,D_i);
//       // change[0] = con_val(N,N,D_o,D_i);
//       // copy_mat(N,N,r,D_i);
//       //   printf("change : %d %lf\n",boo,change);
//       change =  copy_mat_change_copy(N,N,r,D_i,D_o);
//       // con_val_p(N,D_o,D_i,r,change,E_i,q_fin,E_f);
//       // con_copy_mat_change_copy(N, D_o, D_i, r, change, E_i, q_fin, E_f);

//       cout<< "change : "<< boo <<"  "<<change<<endl;
//       mult_mat_p(N,N,E_i,q_fin,E_f);
//       copy_mat(N,N,E_i,E_f);

//          //   cout<<"printing D_i "<<endl;
//          //   print_mat(N,N,D_i);
    
//        boo++;
//       //  exit(0);
//     }
   
// //  exit(0);
//     delete[] q_fin;
//     delete[] E_i;
//     delete[] D_o;
//     delete[] D_i;
//     // exit(0);
//       // printf("here eigen val: %d \n",boo);
     
//       // print_mat(N,N,r);
//       // print_mat()  // e
//       // printf("here U\n");

//      // print_mat(N,N,E_f);    // v
//      // float* sig_nn = new float[N*]();
//    double* v_fin = new double[N*N]; 
//    double* sig_inv = new double[N*N]();
//    double* sig_1 = new double[N*N]();
//       // eigenvec
//       //  printf("here eigen val: \n");
     
//       // print_mat(N,N,r);

//     sort_diag(M,N,r,sig_inv,sig_1,v_fin,E_f);

//    delete[] r;
//    delete[] E_f;

//     double* u_dum = new double[M*N]; 
//     double* u_fin = new double[M*N]; 

//     mult_mat_gen_f(M,N,N,D,v_fin,u_dum);
//     mult_mat_gen(M,N,N,u_dum,sig_inv,u_fin);
//     delete[] u_dum;
//     delete[] sig_inv;

     

//   // printf("SIGMA in: \n");
//   // print_mat(N,N,sig_1);

//     for (int i = 0; i < N; ++i)
//     {
//       /* code */
//       for (int j = 0; j < M; ++j)
//       {
//               if(j<N){
//                  SIGMA[0][i*M+j] =  sig_1[i*N+j];   /* code */
//               }
//               else{
//                  SIGMA[0][i*M+j] =0;
//               }
             
//       }
//     }
//    delete[] sig_1;
   
  
   
//    // making it for Dt 

    
//     double* v_t_dum = new double[N*M];
  

//      trans_mat_p(M,N,u_fin,v_t_dum);
//      delete[] u_fin;
     
//      // // copy_mat(M,N,)

//      for (int i = 0; i < N; ++i)
//      {
//      	for (int j = 0; j < M; ++j)
//      	{
//      		V_T[0][i*M+j] = v_t_dum[i*M+j];
//      		/* code */
//      	}
//      	/* code */
//      }
//     delete[] v_t_dum;
//     copy_mat_f(N,N,*U,v_fin);
     
//     delete[] v_fin;
 
  
  
   



// }

// // /*
// // 	*****************************************************
// // 		TODO -- You must implement this function
// // 	*****************************************************
// // */
// void PCA(int retention, int M, int N, float* D, float* U, float* SIGMA, float** D_HAT, int *K)
// {
//     double* e_val = new double[N];
//     double eigen_sum =0;
//     // print_mat(N,M,SIGMA);
//     for (int i = 0; i < N; ++i)
//     {
//     	// printf("here loop i : %d %f\n",i,SIGMA[i*M+i]);
//     	e_val[i] = SIGMA[i*M+i]*SIGMA[i*M+i];
//     	eigen_sum += e_val[i];
//     	/* code */
//     }

//     // printf("came0\n");

    
//     if(eigen_sum==0){
//     	printf("no eigen values");
//     }

//     for (int i = 0; i < N; ++i)
//     {
//     	e_val[i] = e_val[i] /eigen_sum;
//     	/* code */
//     }
     
   
//     eigen_sum = 0;
//     // printf("came1 \n");
//     for (int i = 0; i < N ; ++i)
//     {
//        eigen_sum += e_val[i];
//        if(100*eigen_sum >= retention){
       
//        	K[0] = i+1;
//        	// printf("val_in is k: %d",K[0]);
//        	break;
//        }
//     	/* code */
//     }
//      delete[] e_val;
//     // printf("came2\n");
//     double* w = new double[N*K[0]];

//     for (int i = 0; i < N; ++i)
//     {
//     	for (int j = 0; j < K[0]; ++j)
//     	{
//     		w[i*K[0]+j] = U[i*M+j];
//     		/* code */
//     	}
//     	/* code */
//     }
//     // printf("came3\n");
//     // print_mat(M,N,D);
//     // printf("came4\n");
//     // print_mat(N,K[0],w);
//     *D_HAT = new float[M*K[0]];
//     mult_mat_gen_f2(M,N,K[0],D,w,*D_HAT);
//     delete[] w; 



// }
