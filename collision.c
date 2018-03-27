#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <fftw3.h>

#define lambda 0.0
#define C_lambda (187.0/M_PI/4.0) //cm^2/g km/s
#define rho_Msun_per_pc3 0.02   
#define this_K 0.6
#define pi M_PI
#define Myr (1000000.0*3600*24*365.25)
#define UnitVelocity_in_cm_per_s 100000.0 // km/s
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define RIGHTBOUND(x,y) (((x) < (y)) ? (x) : (0))
//#define RIGHTBOUND(x,y) (MIN(x,y))
#define third 0.3333333333333333333333333333333333333333333333333333333333333
#define cube(x) (x*x*x)
inline  double simpson_third(int i) {
  return third*((i%2)*2.0+2.0);
}
double calculate_T(double u, double y) {
  double crosssection = C_lambda* UnitVelocity_in_cm_per_s; // cm^3/g/sec
  if (u*y == 0)
    return crosssection*32*M_PI*pow(2.0*u,lambda);
  else
    return crosssection*32*M_PI*pow(2.0*u,lambda)*sin(u*y)/(u*y);
}
double maxwell_dist_3D(double vx, double vy, double vz, double beta) {
  return exp(-0.5*(vx*vx+vy*vy+vz*vz)/(beta*beta))/(sqrt(2.0*pi*beta*beta));
}
double Krook_Wu_dist_1D(double v, double K, double beta) {
  return 0.5*exp(-0.5*v*v/(K*beta*beta))/(sqrt(2.0*pi*K*beta*beta))*((3.*K-1.)/K + (1.-K)*v*v/(K*K*beta*beta));
}
double Krook_Wu_dist_3D(double vx, double vy, double vz, double K, double beta) {
  return 0.5*exp(-0.5*(vx*vx+vy*vy+vz*vz)/(K*beta*beta))/((2.0*pi*K*beta*beta)*sqrt(2.0*pi*K*beta*beta))*((5.*K-3.)/K + (1.-K)*(vx*vx+vy*vy+vz*vz)/(K*K*beta*beta));
}
double calculate_B(double u) {
  double crosssection = C_lambda* UnitVelocity_in_cm_per_s; // cm^3/g/sec
  return 4.*M_PI*crosssection*pow(u,lambda);
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
	int kk = (k+n/2)%(n);
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
	int kk = (k+n/2)%(n);
	b[i*(n)*(n) + j*(n) + k] =  a[ii*(n)*(n) + jj*(n) + kk];
	//printf("%d %d %d -> %d %d %d\t",i,j,k,ii,jj,kk);
      }
      //printf("\n");
    }
    //printf("\n\n");
  }
  memcpy(a, b, sizeof(_Complex double)*(n)*(n)*(n));
  free(b);
}

int main()  {
  int n = 16;
  double L = 800.0;
  L *=  UnitVelocity_in_cm_per_s;
  double sigma = 187.0;
  sigma *=  UnitVelocity_in_cm_per_s;
  double h_v = 2.*L/n;
  double h_v3 = h_v*h_v*h_v;
  double h_y = M_PI/L;
  double h_y3 = h_y*h_y*h_y;
  double normal_coeff = 1./(2*M_PI*sigma*sigma*sqrt(2*M_PI*sigma*sigma));
  double fft_coeff = 1./(8.*M_PI*M_PI*M_PI);
  double V[3] = {0.0, 0.0, 0.0};
  double t = 0.0;
  double rho = rho_Msun_per_pc3*1.989e33/((3.1e18)*(3.1e18)*(3.1e18)); //g/cm^3
  double tau = 1./(4*M_PI*rho*C_lambda*UnitVelocity_in_cm_per_s);
  double t_max = tau*12.0;
  double dt = tau*0.01;
  double *simpson_third_w;

  printf("tau = %g Myr\n",tau/Myr);
  printf("Myr = %g, dt = %g\n, t = %g",Myr,dt/Myr,t);
  int *Cn_grid;
  Cn_grid = malloc(sizeof(int)*(n));
  simpson_third_w = malloc(sizeof(double)*(n));
  for (int i = 0; i < n; i++) {
    Cn_grid[i] = -n/2+i;
    simpson_third_w[i] = simpson_third(i);
  }
  double **Cv_grid;
  double **Cu_grid;
  double **Cy_grid;

  Cv_grid = malloc(sizeof(double*)*3);
  Cu_grid = malloc(sizeof(double*)*3);
  Cy_grid = malloc(sizeof(double*)*3);
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
  for (int j = 0; j < n; j++) {		
    printf("%g %g %g\n",Cv_grid[0][j],Cv_grid[1][j],Cv_grid[2][j]);
  }

  double *dist_fn; //3D velocity distribution function.
  dist_fn = fftw_malloc(sizeof(double)*(n)*(n)*(n));
  double sum = 0.0, sum_vx =0.0, sum_vy = 0.0,sum_vz = 0.0;
  for (int i = 0; i < n; i++) {
    double w_i,w_j,w_k;
    w_i = simpson_third_w[i];
    for (int j = 0; j < n; j++) {
      w_j = simpson_third_w[j];
      for (int k = 0; k < n; k++) {
	w_k = simpson_third_w[k];
	dist_fn[i*(n)*(n) + j*(n) + k] =
	  //maxwell_dist_3D(Cv_grid[0][i]-V[0],Cv_grid[1][j]-V[1],Cv_grid[2][k]-V[2], sigma);
	  Krook_Wu_dist_3D(Cv_grid[0][i]-V[0],Cv_grid[1][j]-V[1],Cv_grid[2][k]-V[2],this_K, sigma);
	sum += dist_fn[i*(n)*(n) + j*(n) + k]*h_v3*w_i*w_j*w_k;
	sum_vx += dist_fn[i*(n)*(n) + j*(n) + k]*Cu_grid[0][i]*h_v3*w_i*w_j*w_k;
	sum_vy += dist_fn[i*(n)*(n) + j*(n) + k]*Cu_grid[1][j]*h_v3*w_i*w_j*w_k;
	sum_vz += dist_fn[i*(n)*(n) + j*(n) + k]*Cu_grid[2][k]*h_v3*w_i*w_j*w_k;
	// printf("Vx = %g, Vy = %g, Vz = %g\n",sum_vx/sum,sum_vy/sum,sum_vz/sum);
      }
    }
  }
  /* printf("Vx = %g, Vy = %g, Vz = %g\n",sum_vx/sum,sum_vy/sum,sum_vz/sum); */
  /* exit(1); */
  double check_sum = 0.0;
  printf("sum = %g\n",sum);
  for (int i =0; i < n; i++) {
    double w_i,w_j,w_k;
    w_i = simpson_third_w[i];
    for (int j = 0; j < n; j++) {
      w_j = simpson_third_w[j];
      for (int k = 0; k < n; k++) {
	w_k = simpson_third_w[k];
	dist_fn[i*(n)*(n) + j*(n) + k] *= rho/sum;
	check_sum += dist_fn[i*(n)*(n) + j*(n) + k]*h_v3*w_i*w_j*w_k;
      }
    }
  }

  printf("check sum = %g\n", check_sum);
  // Set up g3 for each y in C_y
  FILE *fp_time, *fp_x;
  fp_time = fopen("time.out","w");
  fp_x = fopen("vx.out","w");
  while(t < t_max) {
    fftw_complex *g3;
    g3 = fftw_malloc(sizeof(fftw_complex)*(n)*(n)*(n));
    fftw_complex *Q;
    Q = fftw_malloc(sizeof(fftw_complex)*(n)*(n)*(n));
    fftw_plan plan2 = fftw_plan_dft_3d(n, n, n, g3, Q, FFTW_BACKWARD,FFTW_ESTIMATE);
    for(int i = 0; i < (n)*(n)*(n); i++) {
      g3[i] = 0.0;
    }
    // for each u_k in C_u


    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
	for (int k = 0; k < n; k++) {
	  
	  //if((k==0) || (k==n)) wk = 0.5;
	  double abs_u = sqrt(Cu_grid[0][i%n]*Cu_grid[0][i%n]+Cu_grid[1][j%n]*Cu_grid[1][j%n]+Cu_grid[2][k%n]*Cu_grid[2][k%n]);
	  fftw_complex *g1;
	  g1 = fftw_malloc(sizeof(_Complex double)*(n)*(n)*(n));
	  fftw_complex *g2;
	  g2 = fftw_malloc(sizeof(_Complex double)*(n)*(n)*(n));
	  fftw_plan plan = fftw_plan_dft_3d(n, n, n, g1, g2, FFTW_FORWARD,FFTW_ESTIMATE);
	  for (int ii = 0; ii < n; ii++) {
	    for (int jj = 0; jj < n; jj++) {
	      for (int kk = 0; kk < n; kk++) {
		int index_i = RIGHTBOUND(MAX(Cn_grid[ii]-Cn_grid[i%n]+n/2,0),n-1);
		int index_j = RIGHTBOUND(MAX(Cn_grid[jj]-Cn_grid[j%n]+n/2,0),n-1);
		int index_k = RIGHTBOUND(MAX(Cn_grid[kk]-Cn_grid[k%n]+n/2,0),n-1);
		//printf("%d %d %d\n",index_i,index_j,index_k);
		g1[ii*(n)*(n) + jj*(n) + kk] = dist_fn[index_i*(n)*(n) + index_j*(n) + index_k];
		index_i = RIGHTBOUND(MAX(Cn_grid[ii]+Cn_grid[i%n]+n/2,0),n-1);
		index_j = RIGHTBOUND(MAX(Cn_grid[jj]+Cn_grid[j%n]+n/2,0),n-1);
		index_k = RIGHTBOUND(MAX(Cn_grid[kk]+Cn_grid[k%n]+n/2,0),n-1);
		g1[ii*(n)*(n) + jj*(n) + kk] *= dist_fn[index_i*(n)*(n) + index_j*(n) + index_k]*simpson_third_w[ii]*simpson_third_w[jj]*simpson_third_w[kk];
	      }
	    }
	  }

	  // perform forward (-1) FFTW
	  complex_fftshift(g1,n);
	  fftw_execute(plan);
	
	  complex_fftshift(g2,n);
	  fftw_free(g1);
	  
	  // for all element in C_y
	  for (int ii = 0; ii < n; ii++) {
	    for (int jj = 0; jj < n; jj++) {
	      for (int kk = 0; kk < n; kk++) {
		double theta = Cy_grid[0][ii]*V[0] + Cy_grid[1][jj]*V[1] + Cy_grid[2][kk]*V[2];
		_Complex double itheta = cos(theta) + I*sin(theta);
		double abs_y = sqrt(Cy_grid[0][ii]*Cy_grid[0][ii]+Cy_grid[1][jj]*Cy_grid[1][jj]+Cy_grid[2][kk]*Cy_grid[2][kk]);
		double T = calculate_T(abs_u,abs_y);
		// use g3 to collect each value for u_k
		g3[ii*(n)*(n) + jj*(n) + kk] += fft_coeff*h_v3*h_v3*g2[ii*(n)*(n) + jj*(n) + kk]*itheta*T;//*simpson_third_w[i]*simpson_third_w[j]*simpson_third_w[k];

	      }
	    }

	  }
	  fftw_destroy_plan(plan);		
	  fftw_free(g2);	
	}
      }
    }
    

    complex_fftshift(g3,n);
    fftw_execute(plan2);
    fftw_destroy_plan(plan2);
    complex_fftshift(Q,n);
    // end for each u_k
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
	for (int k = 0; k < n; k++) {
    
	  double theta = Cu_grid[0][i]*V[0] + Cu_grid[1][j]*V[1] + Cu_grid[2][k]*V[2];
	  _Complex double itheta = cos(theta) + I*sin(theta);
	  
	  Q[i*(n)*(n) + j*(n) + k] *= itheta*h_y3;
	  //if(creal(Q[i*(n)*(n) + j*(n) + k]) < 0)
	    // printf ("Q < 0; %g\n",creal(Q[i*(n)*(n) + j*(n) + k])*dt/dist_fn[i*(n)*(n) + j*(n) + k]);
	  //printf ("%g\t",creal(Q[i*(n)*(n) + j*(n) + k]));
	}
	//printf("\n");
      }
      //printf("\n");
    }

    double *g;
    g = malloc(sizeof(double)*(n)*(n)*(n));
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
	for (int k = 0; k < n; k++) {
	  g[i*(n)*(n) + j*(n) + k] = 0.0;	
	  for (int ii = 0; ii < n; ii++) {
	    for (int jj = 0; jj < n; jj++) {
	      for (int kk = 0; kk < n; kk++) {
		g[i*(n)*(n) + j*(n) + k] +=
		  //simpson_third_w[ii]*simpson_third_w[jj]*simpson_third_w[kk]
		  dist_fn[ii*(n)*(n) + jj*(n) + kk]
		  *calculate_B(h_v*sqrt((double)((ii - i)*(ii - i)
						 + (jj - j)*(jj - j)
						  + (kk - k)*(kk - k))));
	      }
	    }
	  }
	  g[i*(n)*(n) + j*(n) + k] *= h_v3*dist_fn[i*(n)*(n) + j*(n) + k];
	}
			
      }
    }
    double sum_Qp = 0.0;
    double sum_Qm = 0.0;
    _Complex double *temp_dist_fn1,*temp_dist_fn2,*temp_dist_fn3;
    temp_dist_fn1 = fftw_malloc(sizeof(_Complex double)*n*n*n);
    temp_dist_fn2 = fftw_malloc(sizeof(_Complex double)*n*n*n);
    temp_dist_fn3 = fftw_malloc(sizeof(_Complex double)*n*n*n);
    fftw_plan plan3 = fftw_plan_dft_3d(n, n, n, temp_dist_fn1, temp_dist_fn2, FFTW_BACKWARD,FFTW_ESTIMATE);
    fftw_plan plan4 = fftw_plan_dft_3d(n, n, n, temp_dist_fn2, temp_dist_fn3, FFTW_FORWARD,FFTW_ESTIMATE);
    //double max_t = 10*Myr;

    double check_Q = 0.0;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
	for (int k = 0; k < n; k++) {
	  sum_Qp += creal(Q[i*(n)*(n) + j*(n) + k])*h_v3;
	  sum_Qm += g[i*(n)*(n) + j*(n) + k]*h_v3;
	  Q[i*(n)*(n) + j*(n) + k] -= g[i*(n)*(n) + j*(n) + k];
	  //temp_dist_fn1[i*(n)*(n) + j*(n) + k] = (dist_fn[i*(n)*(n) + j*(n) + k] + creal(Q[i*(n)*(n) + j*(n) + k])*dt);//*simpson_third_w[i]*simpson_third_w[j]*simpson_third_w[k];
	  check_Q += Q[i*(n)*(n) + j*(n) + k]*h_v3;
	  //max_t = MIN(max_t,(0.1*dist_fn[i*(n)*(n) + j*(n) + k]/fabs(creal(Q[i*(n)*(n) + j*(n) + k]))));
	
	}
      }
    }
    //dt = max_t;
    printf("dt = %g Myr\t check_Q = %g Q+ = %g Q- = %g\n",dt/Myr,check_Q,sum_Qp,sum_Qm);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
	for (int k = 0; k < n; k++) {
	  temp_dist_fn1[i*(n)*(n) + j*(n) + k] = (dist_fn[i*(n)*(n) + j*(n) + k] + creal(Q[i*(n)*(n) + j*(n) + k])*dt);//*simpson_third_w[i]*simpson_third_w[j]*simpson_third_w[k];
	  if (creal(temp_dist_fn1[i*(n)*(n) + j*(n) + k]) < 0) {
	    printf("Stop: distribution function < 0; dist = %g Q = %g, new_dist = %g pos = (%d,%d,%d)\n", dist_fn[i*(n)*(n) + j*(n) + k],dt*creal(Q[i*(n)*(n) + j*(n) + k]),creal(temp_dist_fn1[i*(n)*(n) + j*(n) + k]),i,j,k);
	    //exit(0);
	  }
	}
      }
    }
    double check_rho1 = 0.0;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
	for (int k = 0; k < n; k++) {
	  check_rho1 += creal(temp_dist_fn1[i*(n)*(n) + j*(n) + k])*h_v3;
	}
      }
    }
    printf("First Check:\tTime: %g rho = %g, check_rho = %g\n",t/Myr,rho, check_rho1);
    // Conservation correction
    complex_fftshift(temp_dist_fn1,n);
    fftw_execute(plan3);
    fftw_destroy_plan(plan3);
    // origin
    printf("old = %g, new = %g\n",creal(temp_dist_fn2[0]), rho/h_v3);
    temp_dist_fn2[0] = rho/h_v3;
    
    // 1st axis
    for (int i = 0; i < n; i++) {
     
    }
    // 2nd axis
    for (int j = 0; j < n; j++) {
     
    }
    // 3rd axis
    for (int k = 0; k < n; k++) {
      
    }

    fftw_execute(plan4);
    fftw_destroy_plan(plan4);
    
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
	for (int k = 0; k < n; k++) {
	  dist_fn[i*(n)*(n) + j*(n) + k] = creal(temp_dist_fn3[i*(n)*(n) + j*(n) + k])/(n*n*n);
	}
      }
    }
    real_fftshift(dist_fn,n);
    fftw_free(temp_dist_fn1);
    fftw_free(temp_dist_fn2);
    fftw_free(temp_dist_fn3);
    double check_rho = 0.0;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
	for (int k = 0; k < n; k++) {
	  check_rho += dist_fn[i*(n)*(n) + j*(n) + k]*h_v3;
	}
      }
    }
    printf("After Correction:\tTime: %g Myr\t  rho = %g, check_rho = %g\n",t/Myr,rho, check_rho);
    t += dt;
    fftw_free(Q);
    fftw_free(g3);
    free(g);
    fprintf(fp_time,"%g\n",t/tau);
    for (int i = 0; i < n; i++) {
      double gg = 0.0;
      double gc = 0.0;
      for (int j = 0; j < n; j++) {
	for (int k = 0; k < n; k++) {
	  gg += dist_fn[i*(n)*(n) + j*(n) + k]*h_v*h_v;
	  gc += Q[i*(n)*(n) + j*(n) + k]*h_v*h_v*dt;
	}
      }
      fprintf(fp_x,"%g\t",gg);
    }
    fprintf(fp_x,"\n");
    fflush(fp_x);
    fflush(fp_time);
  }
  fclose(fp_time);
  fclose(fp_x);
  return 0;
}
