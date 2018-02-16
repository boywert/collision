#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <fftw3.h>

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

double calculate_T(double u, double y) {
  double C = 1.0;
  if (u*y == 0)
    return C*64*M_PI*u;
  else
    return C*64*M_PI*u*sin(u*y)/(u*y);
}

void real_fftshift(double *a, int n) {
  // return;
  double *b;
  b = malloc(sizeof(double)*(n+1)*(n+1)*(n+1));
  for (int i = 0; i < n+1; i++) {
    int ii = (i+n/2)%(n+1);
    //if(ii < 0) ii+=(n+1);
    for (int j = 0; j < n+1; j++) {
      int jj = (j+n/2)%(n+1);
      //if(jj < 0) jj+=(n+1);
      for (int k = 0; k < n+1; k++) {
	int kk = (i+n/2)%(n+1);
	//if(kk < 0) kk+=(n+1);
	b[i*(n+1)*(n+1) + j*(n+1) + k] =  a[ii*(n+1)*(n+1) + jj*(n+1) + kk];
      }
      // memcpy(&b[i*(n+1)*(n+1) + j*(n+1)], &a[i*(n+1)*(n+1) + (j+n/2)%(n+1)*(n+1)], sizeof(double)*(n+1));
    }
    // memcpy(&b[i*(n+1)*(n+1)], &a[(i+n/2)%(n+1)*(n+1)*(n+1)], sizeof(double)*(n+1)*(n+1));
  }
  memcpy(a, b, sizeof(double)*(n+1)*(n+1)*(n+1));
  free(b);
}
void complex_fftshift(_Complex double *a, int n) {
  // return;
  _Complex double *b;
  b = malloc(sizeof(_Complex double)*(n+1)*(n+1)*(n+1));
  for (int i = 0; i < n+1; i++) {
    int ii = (i+n/2)%(n+1);
    //if(ii < 0) ii+=(n+1);
    for (int j = 0; j < n+1; j++) {
      int jj = (j+n/2)%(n+1);
      //if(jj < 0) jj+=(n+1);
      for (int k = 0; k < n+1; k++) {
	int kk = (i+n/2)%(n+1);
	//if(kk < 0) kk+=(n+1);
	b[i*(n+1)*(n+1) + j*(n+1) + k] =  a[ii*(n+1)*(n+1) + jj*(n+1) + kk];
      }
      // memcpy(&b[i*(n+1)*(n+1) + j*(n+1)], &a[i*(n+1)*(n+1) + (j+n/2)%(n+1)*(n+1)], sizeof(double)*(n+1));
    }
    // memcpy(&b[i*(n+1)*(n+1)], &a[(i+n/2)%(n+1)*(n+1)*(n+1)], sizeof(double)*(n+1)*(n+1));
  }
  memcpy(a, b, sizeof(_Complex double)*(n+1)*(n+1)*(n+1));
  free(b);
}
int main()  {
  int n = 100;
  double L = 100.0;
  double h_v = 2.*L/n;
  double h_v3 = h_v*h_v*h_v;
  double h_y = M_PI/L;
  double h_y3 = h_y*h_y*h_y;
  double sigma = 10.0;
  double normal_coeff = 1./(2*M_PI*sigma*sigma*sqrt(2*M_PI*sigma*sigma));
  double fft_coeff = 1./(8.*L*L*L);
  double V[3] = {0.0, 0.0, 0.0};

	
  int *Cn_grid;
  Cn_grid = malloc(sizeof(int)*(n+1));
  for (int i = 0; i < n+1; ++i) {
    Cn_grid[i] = -n/2+i;
    // printf("%d\n", Cn_grid[i]);
  }
  double **Cv_grid;
  double **Cu_grid;
  double **Cy_grid;
  Cv_grid = malloc(sizeof(double*)*3);
  Cu_grid = malloc(sizeof(double*)*3);
  Cy_grid = malloc(sizeof(double*)*3);
  for (int i = 0; i < 3; i++) {		
    Cv_grid[i] = malloc(sizeof(double)*(n+1));
    Cu_grid[i] = malloc(sizeof(double)*(n+1));
    Cy_grid[i] = malloc(sizeof(double)*(n+1));
    for (int j = 0; j < n+1; j++) {
      Cv_grid[i][j] = V[i] + Cn_grid[j]*h_v;
      Cu_grid[i][j] = Cn_grid[j]*h_v;
      Cy_grid[i][j] = Cn_grid[j]*h_y;
      // printf("%lf \t %lf \t %lf\n", Cv_grid[i][j],Cu_grid[i][j],Cy_grid[i][j]);
    }
  }
  double *dist_fn; //3D velocity distribution function.
  dist_fn = malloc(sizeof(double)*(n+1)*(n+1)*(n+1));
  double sum = 0.0;
  for (int i = 0; i < n+1; i++) {
    for (int j = 0; j < n+1; j++) {
      for (int k = 0; k < n+1; k++) {
	// printf("%d %d %d %d\n",i,j,k, i*(n+1)*(n+1) + j*(n+1) + k);

	dist_fn[i*(n+1)*(n+1) + j*(n+1) + k] =
	  normal_coeff
	  *exp(-0.5/(sigma*sigma)*
	       ((Cv_grid[0][i]-V[0])*(Cv_grid[0][i]-V[0])
		+(Cv_grid[1][j]-V[1])*(Cv_grid[1][j]-V[1])
		+(Cv_grid[2][k]-V[2])*(Cv_grid[2][k]-V[2])
		)
	       );
	/* dist_fn[i*(n+1)*(n+1) + j*(n+1) + k] =  */
	/*   exp(-0.5/(sigma*sigma)* */
	/*        ((Cv_grid[0][i]-V[0])*(Cv_grid[0][i]-V[0]) */
	/* 	+(Cv_grid[1][j]-V[1])*(Cv_grid[1][j]-V[1]) */
	/* 	+(Cv_grid[2][k]-V[2])*(Cv_grid[2][k]-V[2]) */
	/* 	) */
	/*        ) */
	/*    + exp(-0.5/(sigma*sigma)* */
	/*        ((Cv_grid[0][i]+V[0])*(Cv_grid[0][i]+V[0]) */
	/* 	+(Cv_grid[1][j]-V[1])*(Cv_grid[1][j]-V[1]) */
	/* 	+(Cv_grid[2][k]-V[2])*(Cv_grid[2][k]-V[2]) */
	/* 	) */
	/*        ); */
	sum += dist_fn[i*(n+1)*(n+1) + j*(n+1) + k]*h_v3;

      }
    }
  }
  // printf("sum = %lf\n", sum);
  // exit(0);

  // Set up g3 for each y in C_y

  fftw_complex *g3;
  g3 = fftw_malloc(sizeof(fftw_complex)*(n+1)*(n+1)*(n+1));
  fftw_complex *Q;
  Q = fftw_malloc(sizeof(fftw_complex)*(n+1)*(n+1)*(n+1));
  fftw_plan plan2 = fftw_plan_dft_3d(n+1, n+1, n+1, g3, Q, FFTW_BACKWARD,FFTW_ESTIMATE);
  for(int i = 0; i < (n+1)*(n+1)*(n+1); i++) {
    g3[i] = 0.0;
  }
  // for each u_k in C_u


  for (int i = 0; i < n+1; i++) {
    for (int j = 0; j < n+1; j++) {
      for (int k = 0; k < n+1; k++) {
	// super nested loop
	double abs_u = sqrt(Cu_grid[0][i]*Cu_grid[0][i]+Cu_grid[1][j]*Cu_grid[1][j]+Cu_grid[2][k]*Cu_grid[2][k]);
	double *g1;
	g1 = malloc(sizeof(double)*(n+1)*(n+1)*(n+1));
	fftw_complex *g2;
	g2 = fftw_malloc(sizeof(_Complex double)*(n+1)*(n+1)*(n+1));
	fftw_plan plan = fftw_plan_dft_r2c_3d(n+1, n+1, n+1, g1, g2, FFTW_ESTIMATE);
	// for each z_l in C_v
	for (int ii = 0; ii < n+1; ii++) {
	  for (int jj = 0; jj < n+1; jj++) {
	    for (int kk = 0; kk < n+1; kk++) {
	      int index_i = MIN(MAX(Cn_grid[ii]-Cn_grid[i],Cn_grid[0]),Cn_grid[n])+Cn_grid[n];
	      int index_j = MIN(MAX(Cn_grid[jj]-Cn_grid[j],Cn_grid[0]),Cn_grid[n])+Cn_grid[n];
	      int index_k = MIN(MAX(Cn_grid[kk]-Cn_grid[k],Cn_grid[0]),Cn_grid[n])+Cn_grid[n];
	      // printf("index1 %d %d %d\n", index_i, index_j, index_k);
	      // if ((index_i > n) || (index_j > n) || (index_k > n)
	      // 	 || (index_i < 0) || (index_j < 0) || (index_k < 0))   
	      // {
	      // 	printf("index not correct\n");
	      // 	exit(1);
	      // }
	      g1[ii*(n+1)*(n+1) + jj*(n+1) + kk] = dist_fn[index_i*(n+1)*(n+1) + index_j*(n+1) + index_k];
	      index_i = MIN(MAX(Cn_grid[ii]+Cn_grid[i],Cn_grid[0]),Cn_grid[n])+Cn_grid[n];
	      index_j = MIN(MAX(Cn_grid[jj]+Cn_grid[j],Cn_grid[0]),Cn_grid[n])+Cn_grid[n];
	      index_k = MIN(MAX(Cn_grid[kk]+Cn_grid[k],Cn_grid[0]),Cn_grid[n])+Cn_grid[n];
	      // if ((index_i > n) || (index_j > n) || (index_k > n)
	      // 	 || (index_i < 0) || (index_j < 0) || (index_k < 0))   
	      // {
	      // 	printf("index not correct\n");
	      // 	exit(1);
	      // }
	      // printf("index2 %d %d %d\n", index_i, index_j, index_k);
	      g1[ii*(n+1)*(n+1) + jj*(n+1) + kk] *= dist_fn[index_i*(n+1)*(n+1) + index_j*(n+1) + index_k];
	      // printf("%lf\n", g1[ii*(n+1)*(n+1) + jj*(n+1) + kk]);
	    }
	    // printf("\n");
	  }
	}
	// perform forward (-1) FFTW
	real_fftshift(g1,n);
	fftw_execute(plan);
	fftw_destroy_plan(plan);
	free(g1);
	complex_fftshift(g2,n);

	// for all element in C_y
	for (int ii = 0; ii < n+1; ii++) {
	  for (int jj = 0; jj < n+1; jj++) {
	    for (int kk = 0; kk < n+1; kk++) {
	      // printf("%lf\t", creal(g2[ii*(n+1)*(n+1) + jj*(n+1) + kk]));
	      double theta = Cy_grid[0][ii]*V[0] + Cy_grid[1][jj]*V[1] + Cy_grid[2][kk]*V[2];
	      _Complex double itheta = cos(theta) + I*sin(theta);
	      // printf("itheta = %g+%gi, theta = %g\n",creal(itheta),cimag(itheta),theta);
	      /* double realpart =  */
	      /* 	g2[ii*(n+1)*(n+1) + jj*(n+1) + kk][0]*sin(theta) */
	      /* 	-g2[ii*(n+1)*(n+1) + jj*(n+1) + kk][1]*cos(theta); */
	      /* double imgpart =  */
	      /* 	g2[ii*(n+1)*(n+1) + jj*(n+1) + kk][0]*cos(theta) */
	      /* 	+g2[ii*(n+1)*(n+1) + jj*(n+1) + kk][1]*sin(theta); */
	      double abs_y = sqrt(Cy_grid[0][ii]*Cy_grid[0][ii]+Cy_grid[1][jj]*Cy_grid[1][jj]+Cy_grid[2][kk]*Cy_grid[2][kk]);
	      // printf("abs_y = %g\n",abs_y);
	      double T = calculate_T(abs_u,abs_y);
	      // use g3 to collect each value for u_k
	      g3[ii*(n+1)*(n+1) + jj*(n+1) + kk] += fft_coeff*g2[ii*(n+1)*(n+1) + jj*(n+1) + kk]*itheta*T*h_v3;
	 
	    }
	    // printf("\n");
	  }
	  // printf("\n");
	}
	  		
	fftw_free(g2);	
      }
    }
  }

  // exit(0);
  for (int i = 0; i < n+1; i++) {
    for (int j = 0; j < n+1; j++) {
      for (int k = 0; k < n+1; k++) {
	//printf("%lf\t", creal(g3[i*(n+1)*(n+1) + j*(n+1) + k])); //, g3[i*(n+1)*(n+1) + j*(n+1) + k][1]);			
      }
      //printf("\n");
    }
    //printf("\n\n");
  }
  //printf("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
  complex_fftshift(g3,n);
  for (int i = 0; i < n+1; i++) {
    for (int j = 0; j < n+1; j++) {
      for (int k = 0; k < n+1; k++){
	//printf("%lf\t",crealf(g3[i*(n+1)*(n+1) + j*(n+1) + k])); //, g3[i*(n+1)*(n+1) + j*(n+1) + k][1]);			
      }
      //printf("\n");
    }
    //printf("\n\n");
  }
  //exit(0);
  fftw_execute(plan2);
  complex_fftshift(Q,n);
  // fftw_free(g3);
  /* for (int i = 0; i < 1; i++) { */
  /*   for (int j = 0; j < n+1; j++) { */
  /*     for (int k = 0; k < n+1; k++) { */
  /* 	 printf("%lf\t", creal(g3[i*(n+1)*(n+1) + j*(n+1) + k])); //, g3[i*(n+1)*(n+1) + j*(n+1) + k][1]);			 */
  /*     } */
  /*     printf("\n"); */
  /*   } */
  /*   printf("\n\n"); */
  /* } */
  fftw_destroy_plan(plan2);
  // end for each u_k
  for (int i = 0; i < n+1; i++) {
    for (int j = 0; j < n+1; j++) {
      for (int k = 0; k < n+1; k++) {
	// printf("%lf\n", Q[i*(n+1)*(n+1) + j*(n+1) + k][0]);
	double theta = Cu_grid[0][i]*V[0] + Cu_grid[1][j]*V[1] + Cu_grid[2][k]*V[2];
	_Complex double itheta = cos(theta) + I*sin(theta);
	// printf("%lf\n",creal(itheta)*h_y*h_y*h_y);
	Q[i*(n+1)*(n+1) + j*(n+1) + k] *= itheta*h_y*h_y*h_y; //((n+1)*(n+1)*(n+1));
      }
    }	
  }
  for (int i = 0; i < n+1; i++) {
    for (int j = 0; j < n+1; j++) {
      for (int k = 0; k < n+1; k++) {
	//printf("%lf\t", creal(Q[i*(n+1)*(n+1) + j*(n+1) + k]));						
      }
      //printf("\n");
    }
    //printf("\n\n");
  }


  double *g;
  g = malloc(sizeof(double)*(n+1)*(n+1)*(n+1));
  for (int i = 0; i < n+1; i++) {
    for (int j = 0; j < n+1; j++) {
      for (int k = 0; k < n+1; k++) {
	g[i*(n+1)*(n+1) + j*(n+1) + k]	 = 0.0;	
	for (int ii = 0; ii < n+1; ii++) {
	  for (int jj = 0; jj < n+1; jj++) {
	    for (int kk = 0; kk < n+1; kk++) {
	      g[i*(n+1)*(n+1) + j*(n+1) + k] += 
		dist_fn[ii*(n+1)*(n+1) + jj*(n+1) + kk]
		*h_v*sqrt((double) (abs(ii - i)*abs(ii - i)
				    + abs(jj - j)*abs(jj - j)
				    + abs(kk - k)*abs(kk - k)));
	    }
	  }
	}
	g[i*(n+1)*(n+1) + j*(n+1) + k] *= h_v3*dist_fn[i*(n+1)*(n+1) + j*(n+1) + k];
      }
			
    }
  }
  #define Myr 1e6*36000*24*365.25
  for (int i = 0; i < n+1; i++) {
    for (int j = 0; j < n+1; j++) {
      for (int k = 0; k < n+1; k++) {
	Q[i*(n+1)*(n+1) + j*(n+1) + k] -= g[i*(n+1)*(n+1) + j*(n+1) + k];
	printf("%lf\t",Myr*creal(Q[i*(n+1)*(n+1) + j*(n+1) + k]));//dist_fn[i*(n+1)*(n+1) + j*(n+1) + k]/dist_fn[i*(n+1)*(n+1) + j*(n+1) + k]);
      }
      printf("\n");
    }
    printf("\n\n");
  }


  return 0;
}
