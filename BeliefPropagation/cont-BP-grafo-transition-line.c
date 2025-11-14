#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_hyperg.h>


#define FNORM   (2.3283064365e-10)
#define RANDOM  ((ira[ip++] = ira[ip1++] + ira[ip2++]) ^ ira[ip3++])
#define FRANDOM (FNORM * RANDOM) //it generates random numbers between 0 and 1
#define pm1 ((FRANDOM > 0.5) ? 1 : -1)
#define sign(x) ((x) > 0 ? 1 : -1)
#define max(a,b) ((a) > (b) ? (a) : (b))


#define C 3  // Set C to the connectivity (number of neighbors per node)
#define NITER 10000
#define Ninit 100
#define NPRINT_OBS 1
#define PI 3.14159265358979323846
#define Infinite_N_MAX 2000
//#define EPS (4e-10)//(4e-324)

double LAMBDA= 1e-6;
double T;
double BETA;
double MU;
double BL;
double step;
double delta;
int N=256;

/* global variables */
unsigned myrand, ira[256];
unsigned char ip, ip1, ip2, ip3;
double DAMP, damp2;
double mu_in, mu_fin, delta_mu;
int N_tot;

void error(char *string) {
  fprintf(stderr, "ERROR: %s\n", string);
  exit(EXIT_FAILURE);
}

#include "functions-cont-BP-grafo-transition-line.h"



unsigned randForInit(void) {
  unsigned long long y;
  
  y = myrand * 16807LL;
  myrand = (y & 0x7fffffff) + (y >> 31);
  if (myrand & 0x80000000) {
    myrand = (myrand & 0x7fffffff) + 1;
  }
  return myrand;
}

void initRandom(void) {
  int i;
  
  ip = 128;    
  ip1 = ip - 24;    
  ip2 = ip - 55;    
  ip3 = ip - 61;
  
  for (i = ip3; i < ip; i++) {
    ira[i] = randForInit();
  }
}




int main(int argc, char *argv[])
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();
  Tree_type tree;
  Conf_type conf;
  Par_type param;
  Obs_type obs;
  int numbofiter=NITER;
 
 if (argc != 5&&argc!=4) {
    fprintf(stderr, "usage: %s MU T N_tot [seed]\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  
  T = atof(argv[2]);
  N_tot = atoi(argv[3]);

  if (argc == 5) {
    myrand = (unsigned)atoi(argv[4]);
    if (myrand == 2147483647)
      error("seed must be less than 2147483647");
  } else {
    FILE *devran = fopen("/dev/random","r");
    fread(&myrand, 4, 1, devran);
    fclose(devran);
  }

int seed_put= myrand;

  initRandom();
  allocate_tree(&tree);
  
  allocate_obs(&obs);
  allocate_conf(&conf, &tree);
  allocate_par(&param);
  create_grid(&conf);
  compute_simpson_weights(&conf);

  MakeRandomRegularGraph(&tree);
  
    if (tree.signal)
    {
      fprintf(stderr,"GRAFO SBAGLIATO\n");
      return 1;
    }

  double step_MU=0.001;
  double small_step_T=0.0005;
  double large_step_T=0.001;
  double step_T;

  FILE *file_MU_c = fopen("MU_crit.dat", "w"); //open file in writing mode */

  if (file_MU_c == NULL) { 
      printf("Impossibile aprire il file.\n"); 
      return 1; 
  } 


while (T<0.0341){

  if (T<0.032)
  {
    step_T=large_step_T;
  }else{
    step_T=small_step_T;
  }
  
  printf("T= %g\n",T);
  BETA= 1./T;
  BL=BETA*LAMBDA;
  MU = atof(argv[1]);
  compute_field(&param,&conf);

while (MU<0.5){
    printf("# N= %d    T=%g (Beta= %g)    delta= %g    Mu= %g    N_tot= %d   seed= %i\n", N, 1./BETA, BETA, delta, MU, N_tot, seed_put);   

    compute_exp_interaction(&param, &conf);
    computeZ_hat_ij_0(&conf);
  

    printf("0 %g %g \n", obs.mean_AV, obs.variance_AV);
    

    int consecutive_iterations = 0;
    int consecutive_threshold = 10;
    numbofiter=NITER;

    int iter;
    for (iter = 1; iter < NITER + 1; iter++) {
        double diff_means = 0.;
        double diff_variances = 0.;
        BP_iteration(&param, &conf, &tree);
        compute_eta_hat_MOMENTS(&obs, &conf);
        printf("%i %g %g %g %g \n",iter, obs.mean_AV, obs.variance_AV, obs.mean[3],obs.variance[3]);

        for (int i = 0; i < N*C; i++)
        {
          if (fabs((obs.mean[i]-obs.previous_mean[i])/obs.previous_mean[i])>diff_means)
          {
            diff_means=fabs((obs.mean[i]-obs.previous_mean[i])/obs.previous_mean[i]);          
          }      
          if (fabs((obs.variance[i]-obs.previous_variance[i])/obs.previous_variance[i])>diff_variances)
          {
            diff_variances=fabs((obs.variance[i]-obs.previous_variance[i])/obs.previous_variance[i]);
          }    
        }
        
        if (diff_means < 0.000001 && diff_variances <0.000001) {
            consecutive_iterations++;
            if (consecutive_iterations > consecutive_threshold) {
                numbofiter=iter;
                break;
            }
        }
        else {
              consecutive_iterations = 0; 
        }

        for (int i = 0; i < N*C; i++)
        {        
            obs.previous_mean[i]=obs.mean[i];
            obs.previous_variance[i]=obs.variance[i];
        }

    }

    if (numbofiter==NITER)
    {
      fprintf(file_MU_c,"%g\t%g\n",T,MU);
      fflush(file_MU_c); 
      break;
    }

    MU+=step_MU;
  }
T+=step_T;
}

  free_all_memory(&param, &conf, &tree, &obs);

    end = clock();
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("#CPU Time: %f seconds\n", cpu_time_used);

    return 0;

}

