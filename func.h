#ifndef FUNC_H
#define FUNC_H

#include <stdio.h>
#include <utility>
using namespace std;

// int mat_eq(int size,double* mat1,double* mat2 );
void trans_mat_p(int M,int N,double* mat,double* res);
void trans_mat_p_f(int M,int N,float* mat,double* res);

// void trans_mat(int M,int N,double* mat,double* res);
// void mult_mat(int M,int N,double* mat1,double* mat2,double* res);
void mult_mat_p(int M,int N,double* mat1,double* mat2,double* res);
// void mult_mat_p_spld(int M,int N,double* mat1,double* mat2,double* res);
void mult_mat_p_spld(int N,double* mat1,double* mat2,double* res,double* r,double* D_i);
void mult_mat_p_f(int M,int N,double* mat1,float* mat2,double* res);
// void mult_mat_p_spl(int N,double* r,double* q_fin,double* D_i);
double len_mat(int size,int ind,double* mat);
// double mult_mat_ind(int N,double* mat1,double* mat2,int s1,int s2);
void mult_mat_gen(int M1,int N1,int N2,double* mat1,double* mat2,double* res);
void print_mat(int M,int N,double* mat);
void print_mat_f(int M,int N,float* mat);

void mult_mat_gen_f(int M1,int N1,int N2,float* mat1,double* mat2,double* res);
void mult_mat_gen_f2(int M1,int N1,int N2,float* mat1,double* mat2,float* res);
void mult_mat_gen_f_f_d(int M1,int N1,int N2,float* mat1,float* mat2,double* res);
void mult_mat_gen_d_f_d(int M1,int N1,int N2,double* mat1,float* mat2,double* res);
void mult_p_mat_gen_f(int M1,int N1,int N2,float* mat1,double* mat2,double* res);
double copy_mat_change_copy(int M,int N,double* mat1,double* mat2,double* d_o);
// void qr_decomp(int M,int N,double* a,double* q,double* q_dum,double* r);
// void given_ab(double a,double b,double* gt,double* p);
void given_ab_pdf(double a,double b,double* gt,double* p);
void given_inv(double p,double* gt);
// void add_to_i(int m,int s,double mat[],double* gt);
void copy_mat(int M,int N,double* mat1,double* mat2);
void copy_mat_f(int M,int N,float* mat1,double* mat2);
double con_copy_mat_change_copy(int N,double* D_o,double* D_i,double* r,double* E_i,double* E_f);
void qr_decomp_givens(int M,int N,double*a,double* q_fin,int boo);
void qr_decomp_givens_p(int M,int N,double* a,double* q_fin,int boo);
bool compare_fun(pair<double,int> P1,pair<double,int> P2);
double con_val(int M,int N,double* D_o,double* D_i);
// double con_val_spl(int N,double* D_o,double* D_i);
double con_val_p(int N,double* D_o,double* D_i,double* r,double* E_i,double* E_f);
// void opt_my_mult(int M,double* q_fin,double* gt,int m);
void sort_diag(int s,int s_v,double D[],double* sig_inv,double* sig,double* v_fin,double *E_f);
void opt_mult_mat(int M,int N,double* gt,double* A,int i,int j);
// void opt_gt_mult(int M,double* q_fin,double* gt,int m);


#endif