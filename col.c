#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <fftw3.h>

#define lambda 0.0
#define C_lambda 1.0
#define this_rho 1e-26
#define this_K 0.6
#define pi M_PI


#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

double calculate_T(double u, double y) {
  if (u*y == 0)
    return C_lambda*64*M_PI*pow(2.0*u,lambda);
  else
    return C_lambda*64*M_PI*pow(2.0*u,lambda)*sin(u*y)/(u*y);
}
double Krook_Wu_dist_1D(double v, double K, double beta) {
  return 0.5*exp(-0.5*v*v/(K*beta*beta))/(sqrt(2.0*pi*K*beta*beta))*((3.*K-1.)/K + (1.-K)*v*v/(K*K*beta*beta));
}
double Krook_Wu_dist_3D(double vx, double vy, double vz, double K, double beta) {
  return 0.5*exp(-0.5*(vx*vx+vy*vy+vz*vz)/(K*beta*beta))/((2.0*pi*K*beta*beta)*sqrt(2.0*pi*K*beta*beta))*((5.*K-3.)/K + (1.-K)*(vx*vx+vy*vy+vz*vz)/(K*K*beta*beta));
}
double calculate_B(double u) {
  return 4.*M_PI*C_lambda*pow(u,lambda);
}
void real_fftshift(double *a, int n) {
  // return;
  double *b;
  b = malloc(sizeof(double)*(n)*(n)*(n));
  for (int i = 0; i < n; i++) {
    int ii = (i+n/2)%(n);
    for (int j = 0; j < n; j++) {
      int jj = (j+n/2)%(n);
      for (int k = 0; k < n; k++) {
	int kk = (i+n/2)%(n);
	b[i*(n)*(n) + j*(n) + k] =  a[ii*(n)*(n) + jj*(n) + kk];
      }
    }
  }
  memcpy(a, b, sizeof(double)*(n)*(n)*(n));
  free(b);
}
void complex_fftshift(_Complex double *a, int n) {
  // return;
  _Complex double *b;
  b = malloc(sizeof(_Complex double)*(n)*(n)*(n));
  for (int i = 0; i < n; i++) {
    int ii = (i+n/2)%(n);
    for (int j = 0; j < n; j++) {
      int jj = (j+n/2)%(n);
      for (int k = 0; k < n; k++) {
	int kk = (i+n/2)%(n);
	b[i*(n)*(n) + j*(n) + k] =  a[ii*(n)*(n) + jj*(n) + kk];
      }
    }
  }
  memcpy(a, b, sizeof(_Complex double)*(n)*(n)*(n));
  free(b);
}
void create_grids(int* Cn_grid, double** Cv_grid, double** Cu_grid, double** Cy_grid, double *V, int n, double h_v, double h_y) {
  Cn_grid = malloc(sizeof(int)*(n));
  Cv_grid = malloc(sizeof(double*)*3);
  Cu_grid = malloc(sizeof(double*)*3);
  Cy_grid = malloc(sizeof(double*)*3);
  for (int i = 0; i < n; i++) {
    Cn_grid[i] = -n/2+i;
  }
  for (int i = 0; i < 3; i++) {		
    Cv_grid[i] = malloc(sizeof(double)*(n));
    Cu_grid[i] = malloc(sizeof(double)*(n));
    Cy_grid[i] = malloc(sizeof(double)*(n));
    for (int j = 0; j < n; j++) {
      Cv_grid[i][j] = V[i] + Cn_grid[j]*h_v;
      Cu_grid[i][j] = Cn_grid[j]*h_v;
      Cy_grid[i][j] = Cn_grid[j]*h_y;
    }
  }
}

double *init_dist_fn(double** Cv_grid, double* V, int n, double rho, double K, double sigma, double h_v3) {
  double *dist_fn; //3D velocity distribution function.
  dist_fn = malloc(sizeof(double)*(n)*(n)*(n));
  double sum = 0.0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      for (int k = 0; k < n; k++) {
	dist_fn[i*(n)*(n) + j*(n) + k] =
	  Krook_Wu_dist_3D(Cv_grid[0][i]-V[0],Cv_grid[1][j]-V[1],Cv_grid[2][k]-V[2],this_K, sigma);
	sum += dist_fn[i*(n)*(n) + j*(n) + k]*h_v3;
      }
    }
  }
  double check_sum = 0.0;
  printf("sum = %lf\n",sum);
   for (int i =
	  0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      for (int k = 0; k < n; k++) {
	dist_fn[i*(n)*(n) + j*(n) + k] *= rho/sum;
	check_sum += dist_fn[i*(n)*(n) + j*(n) + k]*h_v3;
      }
    }
  }
  return dist_fn;
}
int main()  {
  int n = 16;
  double L = 50.0;
  double h_v = 2.*L/n;
  double h_v3 = h_v*h_v*h_v;
  double h_y = M_PI/L;
  double h_y3 = h_y*h_y*h_y;
  double sigma = 10.0;
  double normal_coeff = 1./(2*M_PI*sigma*sigma*sqrt(2*M_PI*sigma*sigma));
  double fft_coeff = 1./(8.*L*L*L);
  double V[3] = {10.0, 0.0, 0.0};


  int *Cn_grid;
  double **Cv_grid;
  double **Cu_grid;
  double **Cy_grid;

  (void) create_grids(Cn_grid, Cv_grid, Cu_grid, Cy_grid, V, n, h_v, h_y);
  double *dist_fn = init_dist_fn(Cv_grid,V, n, this_rho, this_K, sigma, h_v3);
  return 0;
}
