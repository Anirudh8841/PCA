#include "lab2_io.h"
#include "lab2_omp.h"

#include <stdlib.h>
#include <omp.h>

#include <stdio.h> 
#include "func.h" 
#include <iostream>



int main(int argc, char const *argv[])
{
	if (argc < 3){
		printf("\nLess Arguments\n");
		return 0;
	}

	if (argc > 3){
		printf("\nTOO many Arguments\n");
		return 0;
	}

	//---------------------------------------------------------------------
	int M;			//no of rows (samples) in input matrix D (input)
	int N;			//no of columns (features) in input matrix D (input)
	float* D;		//1D array of M x N matrix to be reduced (input)
	float* U;		//1D array of N x N matrix U (to be computed by SVD)
	float* SIGMA;	//1D array of N x M diagonal matrix SIGMA (to be computed by SVD)
	float* V_T;		//1D array of M x M matrix V_T (to be computed by SVD)
	int K;			//no of coulmns (features) in reduced matrix D_HAT (to be computed by PCA)
	float *D_HAT;	//1D array of M x K reduced matrix (to be computed by PCA)
	int retention;	//percentage of information to be retained by PCA (command line input)
	//---------------------------------------------------------------------

	retention = atoi(argv[2]);	//retention = 90 means 90% of information should be retained

	float start_time, end_time;
	double computation_time;

	/*
		-- Pre-defined function --
		reads matrix and its dimentions from input file and creats array D
	    #elements in D is M * N
        format - 
        --------------------------------------------------------------------------------------
        | D[0][0] | D[0][1] | ... | D[0][N-1] | D[1][0] | ... | D[1][N-1] | ... | D[M-1][N-1] |
        --------------------------------------------------------------------------------------
	*/
	read_matrix (argv[1], &M, &N, &D);

  
 


	U = (float*) malloc(sizeof(float) * N*N);
	SIGMA = (float*) malloc(sizeof(float) * N);
	V_T = (float*) malloc(sizeof(float) * M*M);

	start_time = omp_get_wtime();
	
	// /*
	// 	*****************************************************
	// 		TODO -- You must implement these two function
	// 	*****************************************************
	// */
	// omp_set_num_threads(4);
	SVD(M, N, D, &U, &SIGMA, &V_T);

	PCA(retention, M, N, D, U, SIGMA, &D_HAT, &K);

	end_time = omp_get_wtime();
	computation_time = ((double) (end_time - start_time));
	
	/*
		--Pre-defined functions --
		checks for correctness of results computed by SVD and PCA
		and outputs the results
	*/
	write_result(M, N, D, U, SIGMA, V_T, K, D_HAT, computation_time);


// by me for extra check 

printf("k val is %d\n",K);

  double* D_Tras_p;
  D_Tras_p = (double*) malloc(sizeof(double) * M*N);

  trans_mat_p_f(M,N,D,D_Tras_p);

  
  double* dup_ch = new double[N*M];
  double* dup_dt = new double[N*M];
	float* dup_sig = new float[N*M]();
	// for(int i = 0; i < N; i++)
	// {
		for(int j = 0; j < N; j++)
		{

			dup_sig[j*M+j] = SIGMA[j]; 
			/* code */
		}
		
		/* code */
	// }
	

  mult_mat_gen_f_f_d(N,N,M,U,dup_sig,dup_ch);

  mult_mat_gen_d_f_d(N,M,M,dup_ch,V_T,dup_dt);

  delete[] dup_ch;


  //  printf("DT: \n");
  // print_mat(N,M,SIGMA);

  double change = 0;
  change = con_val(N,M,dup_dt,D_Tras_p);


  //  printf("DT in main : \n");
  // print_mat(N,M,D_Tras_p);


  delete[] dup_dt;
  delete[] D_Tras_p;
  // printf("change val is %lf\n",change);
  // printf("time taken is %lf\n",computation_time);

//  printf(" U \n");
//   print_mat_f(N,N,U);
//   printf("SIGMA: \n");
//   for (int i = 0; i < N; ++i)
//   {
//   	cout<<SIGMA[i*M+i]<<endl;
//   	/* code */
//   }


  // print_mat(N,M,SIGMA);
 
cout<< " val k is "<<K<<endl;
cout<< "time taken is  "<< computation_time<<endl;
cout<< "change val is "<< change<<endl;
  





	return 0;
}
