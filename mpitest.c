#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <time.h>
#include <unistd.h>

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define Nhist 16
double P(double x) {
  return x;
}
int binarySearch(double *arr, int l, int r, double x) {
  printf("%d %d\n",l,r);
  if (r > l) {
        int mid = l + (r - l)/2;
 
        // If the element is present at the middle 
        // itself
        if (arr[mid] == x)  
            return mid;
 
        // If element is smaller than mid, then 
        // it can only be present in left subarray
        if (arr[mid] > x) 
            return binarySearch(arr, l, mid-1, x);
 
        // Else the element can only be present
        // in right subarray
        return binarySearch(arr, mid+1, r, x);
   }
 
   // We reach here when element is not 
   // present in array
  printf("x = %g, find %g %g\n",x,arr[l],arr[l+1]);
  return l;
}
double interpolate(double *rf, double *rx, double x, int Np) {
  int i = binarySearch(rx,0,Np-1,x);
  return rf[i]+(rf[i+1]-rf[i])/(rx[i+1]-rx[i])*(x-rx[i]);
}
// A utility function to swap to integers
void swap (int *a, int *b)
{
    int temp = *a;
    *a = *b;
    *b = temp;
}
 
// A utility function to print an array
void printArray (int arr[], int n)
{
    for (int i = 0; i < n; i++)
        printf("%d ", arr[i]);
    printf("\n");
}
 
// A function to generate a random permutation of arr[]
void randomize ( int arr[], int n )
{
    // Use a different seed value so that we don't get same
    // result each time we run this program
    srand ( time(NULL) );
 
    // Start from the last element and swap one by one. We don't
    // need to run for the first element that's why i > 0
    for (int i = n-1; i > 0; i--)
    {
        // Pick a random index from 0 to i
        int j = rand() % (i+1);
 
        // Swap arr[i] with the element at random index
        swap(&arr[i], &arr[j]);
    }
}

double *init_vel(int N, int Np, double *rf, double *rx, double xmin, double xmax) {
  int MaxSize = 10000;
  double *vel;
  int largeset = 0;
  vel = (double*)malloc(sizeof(double)*N*3);
  int accepted = 0;
  double ymax = 0.0;
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus2);
  gsl_rng_set(rng, time(NULL)); // Seed with time

  for(int i=0; i < Np; i++){
    ymax = MAX(ymax,rf[i]);
  }
  ymax *= 1.1;
  {
    gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, Np);
    gsl_spline_init (spline, rx, rf, Np);

    while(accepted < N) {
      double x = xmax*gsl_rng_uniform(rng)+xmin;
      double y = ymax*gsl_rng_uniform(rng);
      //printf("\tx = %g, y = %g\n",x,y);
      double v = gsl_spline_eval (spline, x, acc);
      //printf("\tx = %g, v = %g\n",x,v);
      if(y < v) {
	double phi = 2.0*M_PI*gsl_rng_uniform(rng);
	double costheta = -1.0 + 2.0*gsl_rng_uniform(rng);
	double theta = acos(costheta);
	//printf("theta = %g, phi = %g, s(theta) = %g, s(phi) = %g\n",theta,phi,sin(theta),sin(phi));
	vel[accepted*3] = x*cos(theta);
	vel[accepted*3+1] = x*sin(theta)*sin(phi);
	vel[accepted*3+2] = x*sin(theta)*cos(phi);
	accepted++;
      }
      //printf("accepted = %d / %d\n",accepted, N);
    }
    //exit(9);
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);
  }
  gsl_rng_free (rng);

  return vel;
}
void write_histogram(FILE *fp, double *vel, int N, double time) {
  long int_hist[Nhist];
  for(int i=0; i < Nhist; i++)
    int_hist[i] = 0;
  double vmin = -400.0*1e5;
  double vmax = 400.0*1e5;
  double dv = (vmax-vmin)/Nhist;
  for(int i=0; i<N; i++) {
    if((vel[3*i] > vmin ) && (vel[3*i] < vmax)) {
      int index = (int)((vel[3*i]-vmin)/dv);
      int_hist[index]++;
    }
  }
  fprintf(fp,"%g\t",time);
  for(int k=0; k<Nhist; k++)
    fprintf(fp,"%g\t",(float)int_hist[k]/N/dv);
  fprintf(fp,"\n");
  fflush(fp);
}
double KrookWu_f_speed(double v, double K, double beta) {
  return 2.0*M_PI*v*v*exp(-0.5*v*v/(K*beta*beta))/((2.0*M_PI*K*beta*beta)*sqrt(2.0*M_PI*K*beta*beta))*((5.*K-3.)/K + (1.-K)*v*v/(K*K*beta*beta));
}

int main(int argc, char** argv) {
  MPI_Init(NULL, NULL);
  
  // Get the number of processes
  int MPI_nRank;
  MPI_Comm_size(MPI_COMM_WORLD, &MPI_nRank);
  
  // Get the rank of the process
  int MPI_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &MPI_rank);

  int N = 32;
  N = N*N*N;
  time_t p_t;
  struct tm* tm;
  char Time[23],fp_out_name[256];
  time(&p_t);
  tm = localtime(&p_t);
  strftime(Time, sizeof Time, "%Y%m%d_%I%M%S", tm);
  sprintf(fp_out_name,"output_%d_%s",N,Time);

  int Np = 10000;
  double sigma = 187.0*1e5; //cm/s
  double xmin = 0.0;
  double xmax = sigma*10;
  double cspm = 1.0*sigma;
  double rho_Msun_per_pc3 = 0.02;
  double rho = rho_Msun_per_pc3*1.989e33/((3.1e18)*(3.1e18)*(3.1e18)); //g/cm^3
  double Mp = rho/N;
  double kappa = cspm/4.0/M_PI;
  double *x, *f;
  double tau = 1/(rho*4.0*M_PI*kappa);
#define Myr (3600*24*365.25*1e6)
  printf("tau = %g Myr\n",tau/Myr);
 
  double dt = 0.01*tau;
  x = (double *)malloc(sizeof(double)*Np);
  f = (double *)malloc(sizeof(double)*Np);
  for(int i = 0; i < Np; i++) { 
    x[i] = i*(xmax-xmin)/(Np-1) + xmin; 
    f[i] = KrookWu_f_speed(x[i],0.6,sigma/sqrt(3.0));
    //printf("%g\t",f[i]);
  }
  //printf("\n");
  double *vel;
  vel = init_vel(N, Np, f, x, xmin, xmax);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&vel, 3*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus2);
  gsl_rng_set(rng, time(NULL)); // Seed with time
  double t = 6.0*log(2.5)*tau;

  int s_index[N];
  for(int i = 0; i< N; i++)
    s_index[i] = i;
  FILE *fp_out;
  if(MPI_rank == 0) {
    fp_out = fopen (fp_out_name,"w");
    fprintf(fp_out, "%g\t%d\n", tau, Nhist );
    write_histogram(fp_out, vel, N, t/tau);
  }
  while (t < 12*tau) {
    //randomize (s_index, N);
    for(int ii=0; ii<N; ii++) {
      int i = s_index[ii];
      for(int jj=ii+1; jj<N; jj++) {
	int j = s_index[jj];
	double v1[3],v2[3];
	v1[0] = vel[i*3];
	v1[1] = vel[i*3+1];
	v1[2] = vel[i*3+2];
	v2[0] = vel[j*3];
	v2[1] = vel[j*3+1];
	v2[2] = vel[j*3+2];      
	double rel_v = (v1[0]-v2[0])*(v1[0]-v2[0]);
	rel_v += (v1[1]-v2[1])*(v1[1]-v2[1]);
	rel_v += (v1[2]-v2[2])*(v1[2]-v2[2]);
	rel_v = sqrt(rel_v);
	double crosssection = kappa*4.*M_PI;
	double P = rho/N*crosssection*dt;
	double finger = gsl_rng_uniform(rng);
	//printf("P = %g, finger = %g\n",P, finger);
	if(finger < P) {
	  double n[3];
	  double phi = 2.0*M_PI*gsl_rng_uniform(rng);
	  double costheta = -1.0 + 2.0*gsl_rng_uniform(rng);
	  double theta = acos(costheta);
	  n[0] = cos(theta);
	  n[1] = sin(theta)*sin(phi);
	  n[2] = sin(theta)*cos(phi);
	  double dotproduct = (v1[0]-v2[0])*n[0] + (v1[1]-v2[1])*n[1] + (v1[2]-v2[2])*n[2];
	  //printf("i = %d, dot = %g\n",i,dotproduct);
	  vel[i*3] = v1[0]-dotproduct*n[0];
	  vel[i*3+1] = v1[1]-dotproduct*n[1];
	  vel[i*3+2] = v1[2]-dotproduct*n[2];
	  vel[j*3] = v2[0]+dotproduct*n[0];
	  vel[j*3+1] = v2[1]+dotproduct*n[1];
	  vel[j*3+2] = v2[2]+dotproduct*n[2];
	}
      }
    }
    t+=dt;
    if(MPI_rank == 0)
      write_histogram(fp_out, vel, N, t/tau);
 

  }
  MPI_Finalize();
  return(0);
}
