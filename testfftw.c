#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <fftw3.h>

void real_fftshift(double *a, int n) {
  double *b;
  b = malloc(sizeof(double)*(n));
  for (int i = 0; i < n; i++) {
    int ii = (i+n/2)%(n);
    b[i] =  a[ii];
  }
  memcpy(a, b, sizeof(double)*(n));
  free(b);
}
void complex_fftshift(_Complex double *a, int n) {
  // return;
  _Complex double *b;
  b = malloc(sizeof(_Complex double)*(n));
  for (int i = 0; i < n; i++) {
    int ii = (i+n/2)%(n);
    b[i] =  a[ii];
  }
  memcpy(a, b, sizeof(_Complex double)*(n));
  free(b);
}


int main() {
  fftw_complex *a;
  fftw_complex *x;
  fftw_complex *b,*c;
  double abs_x = 5.0;
  int N = 129;
  x = fftw_malloc(sizeof(fftw_complex)*N);
  a = fftw_malloc(sizeof(fftw_complex)*N);
  b = fftw_malloc(sizeof(fftw_complex)*N);
  c = fftw_malloc(sizeof(fftw_complex)*N);
  fftw_plan plan = fftw_plan_dft_1d(N,a,b,FFTW_BACKWARD,FFTW_ESTIMATE);
  fftw_plan plan2 = fftw_plan_dft_1d(N,b,c,FFTW_FORWARD,FFTW_ESTIMATE);
  for(int i=0; i<N; i++) {
    x[i] = -1.*abs_x+2.*abs_x/N*i;
    a[i] = exp(-1.*x[i]*x[i]);
  }
  for(int i=0; i<N; i++) {
    printf("%g\t",creal(a[i]));
  }
  complex_fftshift(a,N);
  printf("\n");
  for(int i=0; i<N; i++) {
    printf("%g\t",creal(a[i]));
  }
  fftw_execute(plan);
  complex_fftshift(a,N);
  complex_fftshift(b,N);

  // backward
  complex_fftshift(b,N);
  fftw_execute(plan2);
  complex_fftshift(c,N);
  printf("\n");
  for(int i=0; i<N; i++) {
    printf("%g\t",creal(a[i])-creal(c[i])/N);
  }
  return 0;
  
}
