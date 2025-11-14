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
#define NPRINT_OBS 200
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
int N;

/* global variables */
unsigned myrand, ira[256];
unsigned char ip, ip1, ip2, ip3;
double DAMP;
int N_tot;

void error(char *string) {
  fprintf(stderr, "ERROR: %s\n", string);
  exit(EXIT_FAILURE);
}

typedef struct parameters{ 

  double *field;
  double *field_gauss_part;
  double **exp_interaction;
  double **integrand;
  double *integrand_zi;
  double *integrand_gaussian_part;

}Par_type;


typedef struct tree{

  int **ivic;
  int **ineighpos;
  int *n_neigh;
  int *unchosen;
  int signal;
  double minAlpha;
  
}Tree_type;


typedef struct configuration{

  double ** eta_hat;
  double ** eta_i;
  double ** gauss_i;
  double ** eta_hat_new;
  double ** NN_eta_hat;
  double * grid;
  double * simps_weight;
  double * NNetahatij;
  
}Conf_type;

typedef struct observables{

  double *mean;
  double *previous_mean;
  double mean_AV;
  
  double *variance;
  double *previous_variance;
  double variance_AV;


  double *true_mean;
  double true_mean_AV;
  
  double *true_variance;
  double true_variance_AV;

  double *gauss_part_mean;
  double gauss_part_mean_AV;
  
  double *gauss_part_variance;
  double gauss_part_variance_AV;

  double *gauss_part_skewness;
  double gauss_part_skewness_AV;

  double * av_eta_i;
  double *av_gauss_i;

  
}Obs_type;


void allocate_tree(Tree_type * tree){
  
  int i;
  tree->ivic = (int **)calloc(N,sizeof(int *));
  for(i=0;i<N;i++) tree->ivic[i] = (int *)calloc(C,sizeof(int));

  tree->ineighpos = (int **)calloc(N,sizeof(int *));
  for(i=0;i<N;i++) tree->ineighpos[i] = (int *)calloc(C,sizeof(int));

  tree->n_neigh= (int *)calloc(N,sizeof(int));     
  tree->unchosen= (int *)calloc(N*C,sizeof(int));  

  tree->signal=0;   

  return;
  
}


void allocate_obs(Obs_type * obs){

  int l,t;
  obs->mean = (double *)calloc(N*C,sizeof(double));
  obs->previous_mean = (double *)calloc(N*C,sizeof(double));

  for ( l = 0; l < N*C; l++)
  {
    obs->previous_mean[l]=1.;
  }
  
  obs->variance = (double *)calloc(N*C,sizeof(double));
  obs->previous_variance = (double *)calloc(N*C,sizeof(double));

  for ( t = 0; t < N*C; t++)
  {
    obs->previous_variance[t]=1.;
  }


  obs->mean_AV=0;

  obs->variance_AV=0;


  obs->true_mean = (double *)calloc(N,sizeof(double));
  obs->true_variance = (double *)calloc(N,sizeof(double));

  obs->av_eta_i = (double *)calloc(N_tot,sizeof(double));

  obs->gauss_part_mean = (double *)calloc(N,sizeof(double));
  obs->gauss_part_variance = (double *)calloc(N,sizeof(double));
  obs->gauss_part_skewness = (double *)calloc(N,sizeof(double));

  obs->av_gauss_i = (double *)calloc(N_tot,sizeof(double));
  
  return;
  
}

void allocate_conf(Conf_type * conf, Tree_type * tree){
 
  int i;
  
  conf->eta_hat = (double **)calloc(N*C,sizeof(double *));
  for(i=0;i<N*C;i++) conf->eta_hat[i] = (double *)calloc(N_tot,sizeof(double));

  conf->eta_i = (double **)calloc(N,sizeof(double *));
  for(i=0;i<N;i++) conf->eta_i[i] = (double *)calloc(N_tot,sizeof(double));

  conf->gauss_i = (double **)calloc(N,sizeof(double *));
  for(i=0;i<N;i++) conf->gauss_i[i] = (double *)calloc(N_tot,sizeof(double));

  conf->eta_hat_new = (double **)calloc(N*C,sizeof(double *));
  for(i=0;i<N*C;i++) conf->eta_hat_new[i] = (double *)calloc(N_tot,sizeof(double));

  conf->NN_eta_hat = (double **)calloc(N*C,sizeof(double *));
  for(i=0;i<N*C;i++) conf->NN_eta_hat[i] = (double *)calloc(N_tot,sizeof(double));

  conf->NNetahatij= (double *)calloc(N_tot,sizeof(double));

    
  return;

}

void allocate_par(Par_type * param){

  int i;  
  param->field = (double *)calloc(N_tot,sizeof(double));
  param->field_gauss_part = (double *)calloc(N_tot,sizeof(double));

  param->exp_interaction = (double **)malloc(N_tot*sizeof(double*));
  for (i = 0; i < N_tot; i++)
  {
    param->exp_interaction[i] = (double *)calloc(N_tot,sizeof(double));
  }   
  param->integrand = (double **)malloc((N_tot)*sizeof(double*));
  for (i = 0; i < N_tot; i++)
  {
    param->integrand[i] = (double *)calloc(N_tot,sizeof(double));    
  }

  param->integrand_zi = (double *)calloc(N_tot,sizeof(double)); 
  param->integrand_gaussian_part = (double *)calloc(N_tot,sizeof(double)); 

  return;

}


double gaussian_pdf(double x, double mean) {
    double coefficient = 1.0 / (sqrt(2.0 * PI * 25));
    double exponent = -((x - mean) * (x - mean)) / (2.0 * 25);
    return coefficient * exp(exponent);
}

void compute_exp_interaction(Par_type *param, Conf_type *conf){

int i,j;
    for (i = 0; i < N_tot; i++)
    {
      for (j = 0; j < N_tot; j++)
      {
        param->exp_interaction[i][j] = exp(- BETA * MU * conf->grid[i] * conf->grid[j]);
      }
    } 

return;

}

double parabolic_cylinder_D(double v, double z) {
    const double z2 = 0.5 * z * z;            

    const double M1 = gsl_sf_hyperg_1F1(-0.5 * v, 0.5, z2);   
    const double M2 = gsl_sf_hyperg_1F1((1.0 - v) * 0.5, 1.5, z2);

    const double sqrt_pi = sqrt(PI);

    const double pref = pow(2.0, 0.5 * v) * exp(-0.25 * z * z);
    const double term1 = (sqrt_pi / gsl_sf_gamma((1.0 - v) * 0.5)) * M1;
    const double term2 = (sqrt(2.0) * sqrt_pi * z / gsl_sf_gamma(-0.5 * v)) * M2;

    return pref * (term1 - term2);
}

double analytic_first_step(double beta, double lambda, double alpha_ij, double n_j) {
    const double bl = beta * lambda;
    const double z = sqrt(beta) * (alpha_ij * n_j - 1.0);

    const double pref_beta = pow(beta, -0.5 * bl);
    const double gamma_bl = gsl_sf_gamma(bl);
    const double ez = exp(0.25 * z * z);

    const double D = parabolic_cylinder_D(-bl, z);

    return pref_beta * gamma_bl * ez * D;
}


double computeZ_hat_ij_0(Conf_type *conf)
{
    int i;
    double z = 0.;
    
    for (i = 0; i < N_tot; i++)
    {
        conf->NN_eta_hat[0][i]=analytic_first_step(BETA, LAMBDA, MU, conf->grid[i]);
        z += conf->NN_eta_hat[0][i];
    }

    for (i = 0; i < N_tot; i++)
    {
        conf->eta_hat[0][i]= conf->NN_eta_hat[0][i]/z;
    }

    for (int j = 1; j < N*C; j++)
    {
      for (i = 0; i < N_tot; i++)
      {
        conf->NN_eta_hat[j][i]=conf->NN_eta_hat[0][i];
        conf->eta_hat[j][i]=conf->eta_hat[0][i];
      }
      
    }
    
    return 0;
}


void compute_eta_hat_MOMENTS(Obs_type *obs, Conf_type *conf){

int j;
double sum_mean=0.;
double sum_var=0.;

for (int i = 0; i < N*C; i++)
{
  obs->mean[i]=0.;
  obs->variance[i]=0.;
}


for (int i = 0; i < N*C; i++)
{
  for(j = 1; j < N_tot; j++){
    obs->mean[i]+=conf->simps_weight[j] * conf->grid[j] * conf->eta_hat[i][j];
  }

  sum_mean += obs->mean[i];

  for(j = 1; j < N_tot; j++){
    obs->variance[i]+= conf->simps_weight[j] * pow(conf->grid[j]-obs->mean[i],2) * conf->eta_hat[i][j];
  }

  sum_var += obs->variance[i];

}

obs->mean_AV=sum_mean / (N*C);
obs->variance_AV=sum_var / (N*C);
    
return;
  
}


void compute_true_MOMENTS(Obs_type *obs, Conf_type *conf){

int n_i;
double sum_mean=0.;
double sum_var=0.;

for (int i = 0; i < N; i++)
{
  obs->true_mean[i]=0.;
  obs->true_variance[i]=0.;
}


for (int i = 0; i < N; i++)
{
  for(n_i = 1; n_i < N_tot; n_i++){
    obs->true_mean[i]+=conf->simps_weight[n_i] * conf->grid[n_i] * conf->eta_i[i][n_i];
  }

  sum_mean += obs->true_mean[i];

  for(n_i = 1; n_i < N_tot; n_i++){
    obs->true_variance[i]+= conf->simps_weight[n_i] * pow(conf->grid[n_i] - obs->true_mean[i],2) * 
                            conf->eta_i[i][n_i];
  }

  sum_var += obs->true_variance[i];

}

obs->true_mean_AV=sum_mean / N;
obs->true_variance_AV=sum_var / N;
    
return;
  
}


void compute_gauss_part_MOMENTS(Obs_type *obs, Conf_type *conf){

int n_i;
double sum_mean=0.;
double sum_var=0.;
double sum_skew=0.;

for (int i = 0; i < N; i++)
{
  obs->gauss_part_mean[i]=0.;
  obs->gauss_part_variance[i]=0.;
  obs->gauss_part_skewness[i]=0.;
}


for (int i = 0; i < N; i++)
{
  for(n_i = 1; n_i < N_tot; n_i++){
    obs->gauss_part_mean[i]+=conf->simps_weight[n_i] * conf->grid[n_i] * conf->gauss_i[i][n_i];
  }

  sum_mean += obs->gauss_part_mean[i];

  for(n_i = 1; n_i < N_tot; n_i++){
    obs->gauss_part_variance[i]+= conf->simps_weight[n_i] * pow(conf->grid[n_i] - obs->gauss_part_mean[i],2) * 
                                  conf->gauss_i[i][n_i];
  }

  sum_var += obs->gauss_part_variance[i];

  for(n_i = 1; n_i < N_tot; n_i++){
    obs->gauss_part_skewness[i]+= conf->simps_weight[n_i] *  pow(conf->grid[n_i] - obs->gauss_part_mean[i], 3) * 
                                  conf->gauss_i[i][n_i];
  }

  sum_skew += obs->gauss_part_skewness[i];

}

obs->gauss_part_mean_AV=sum_mean / N;
obs->gauss_part_variance_AV=sum_var / N;
obs->gauss_part_skewness_AV=sum_skew / N;
    
return;
  
}

void free_all_memory(Par_type *param, Conf_type *conf,Tree_type * tree, Obs_type *obs){

    int i;

     for (i = 0; i < N; i++)
    {
        free(tree->ivic[i]);
    }
    free(tree->ivic);
    

    for (i = 0; i < C*N; i++)
    {
        free(conf->eta_hat[i]);
        free(conf->eta_hat_new[i]);
        free(conf->NN_eta_hat[i]);
    }
    free(conf->eta_hat);
    free(conf->eta_hat_new);
    free(conf->NN_eta_hat);

    
    free(obs->mean);
    free(obs->variance);
    free(obs->previous_mean);
    free(obs->previous_variance);
    free(obs->av_eta_i);
    free(obs->gauss_part_mean);
    free(obs->gauss_part_variance);
    free(obs->gauss_part_skewness);
    free(obs->av_gauss_i);
    free(tree->ineighpos);
    free(tree->n_neigh);
    free(tree->unchosen);
    free(conf->eta_i);
    free(conf->gauss_i);
    free(conf->grid);
    free(conf->simps_weight);
    free(conf->NNetahatij);
    free(param->field);
    free(param->field_gauss_part);
    free(param->integrand_zi);
    free(param->integrand_gaussian_part);
    for (i = 0; i < N_tot; i++) {
        free(param->integrand[i]);
    }
    free(param->integrand);

for (i = 0; i < N_tot; i++) {
        free(param->exp_interaction[i]);
}
free(param->exp_interaction);


    return;
}

void create_grid(Conf_type *conf){

    conf->grid = (double *)calloc(N_tot,sizeof(double));

    double start = 0.0;
    double end = 2.0;

    step = (end - start) / N_tot;
    delta= step/50.0;

    int i;
    conf->grid[0] = start; 
    conf->grid[1] = delta;
    for (i = 2; i < N_tot; i++) {
        conf->grid[i] = delta + (i-1) * step;
    }

    return;

}

void compute_simpson_weights(Conf_type *conf){
  int i;
  conf->simps_weight = (double *)calloc(N_tot,sizeof(double));
  conf->simps_weight[1]=step / 3.0;
  conf->simps_weight[N_tot-1]=step / 3.0;
  for (i = 2; i < N_tot-1; i+=2)
  {
    conf->simps_weight[i]= 4*step / 3.0;
  }
  for (i = 3; i < N_tot-1; i+=2)
  {
    conf->simps_weight[i]= 2*step / 3.0;
  }
  
}

void compute_field(Par_type *param, Conf_type *conf){

    int i;
    for (i = 1; i < N_tot; i++)
    {
        param->field[i] = pow(conf->grid[i], BL-1.0) * exp(- 0.5 * BETA * (conf->grid[i]*conf->grid[i]-2*conf->grid[i]));
        param->field_gauss_part[i] = exp(- 0.5 * BETA * (conf->grid[i]*conf->grid[i]-2*conf->grid[i]));
    } 

return;

}


int MakeRandomRegularGraph(Tree_type *tree) {
  int i, j1, j2, index, t, tMax, flag;
  int  numUnchosen;
  
  for (i = 0; i < C * N; i++) {
    tree->unchosen[i] = i % N; 
  }
  numUnchosen = C * N;
  while (numUnchosen) {//while numUnchosen is not 0
    index = FRANDOM * (numUnchosen); 
    j1 = tree->unchosen[index];
    tree->unchosen[index] = tree->unchosen[--numUnchosen];
    tree->unchosen[numUnchosen] = j1;
    t = 0;
    tMax = 10 * numUnchosen;
    do {
      index = FRANDOM * (numUnchosen);
      j2 = tree->unchosen[index];
      t++;
      flag = 0;
      if (j2==j1) flag=1;
      for (i = 0; i < tree->n_neigh[j1]; i++)
	if (tree->ivic[j1][i] == j2) flag = 1;
    } while (flag && t < tMax);
    if (t == tMax) {
      //printf("REDO\n");
      tree->signal=1;
      return 1;
    }
    tree->unchosen[index] = tree->unchosen[--numUnchosen];
    tree->unchosen[numUnchosen] = j2;

    tree->ivic[j1][tree->n_neigh[j1]] = j2;
    tree->ineighpos[j1][tree->n_neigh[j1]]=tree->n_neigh[j2];

    tree->ivic[j2][tree->n_neigh[j2]] = j1;
    tree->ineighpos[j2][tree->n_neigh[j2]]=tree->n_neigh[j1];

    tree->n_neigh[j1]++;
    tree->n_neigh[j2]++;
  }
  for (i = 0; i < N; i++) {
    if (tree->n_neigh[i] != C) fprintf(stderr,"generazione grafo1\n");
  }

  return 0;
}



void compute_integrand(int i, int j, Par_type *param, Conf_type *conf, Tree_type *tree) {

  int k,n_i,n_j;
    for (n_i = 0; n_i < N_tot; n_i++)
    {
    double prod_integr=1.0;
    for (k = 0; k < C; k++)
    {
        int neighbor_index = tree->ivic[i][k];
          
        if (neighbor_index == j)
        {
            continue;
        }
        prod_integr*= conf->eta_hat[3*neighbor_index + tree->ineighpos[i][k]][n_i];
    }

    for (n_j = 0; n_j < N_tot; n_j++)
    {
      param->integrand[n_i][n_j]=param->field[n_i] * param->exp_interaction[n_i][n_j] * prod_integr;
    }

    }
    
    return;
}


void compute_integrand_zi(int i, Par_type *param, Conf_type *conf, Tree_type *tree) {

  int k,n_i;
    for (n_i = 1; n_i < N_tot; n_i++)
    {
    double prod_integr=1.0;
    for (k = 0; k < C; k++)
    {
        int neighbor_index = tree->ivic[i][k];
        prod_integr*= conf->eta_hat[3*neighbor_index + tree->ineighpos[i][k]][n_i];
    }

    param->integrand_zi[n_i]=param->field[n_i] * prod_integr;
    param->integrand_gaussian_part[n_i] = param->field_gauss_part[n_i] * prod_integr;
    }
    
    return;
}




double compute_integral(int n_j,Par_type *param, Conf_type *conf){

  int i;
  double integral=0.0;
      for (i = 1; i < N_tot; i++) //j starts from 1 because I'm integrating from delta, not from 0
      {
        integral+=conf->simps_weight[i]*param->integrand[i][n_j];
      } 

    return integral;
}


double compute_integral_zi(Par_type *param, Conf_type *conf, double *integral_gaussian_part){

  int i;
  double integral=0.0;
  *integral_gaussian_part=0.0;
  for (i = 1; i < N_tot; i++) //j starts from 1 because I'm integrating from delta, not from 0
  {
    integral+=conf->simps_weight[i]*param->integrand_zi[i];
    *integral_gaussian_part += conf->simps_weight[i]*param->integrand_gaussian_part[i];
  } 

  return integral;
}



double computeZij(int i, int j, Par_type *param, Conf_type *conf, Tree_type *tree) //I'm using it like this: double z_ij = computeZ(i, vic[i][j], ...);
{
    int k,n_j;
    double z = 0.0;

    double first_part=0.0;
    double second_part=0.0;
    double prod1=1.0;
    for (k = 0; k < C; k++)
    {
        int neighbor_index = tree->ivic[i][k];
          
        if (neighbor_index == j)
        {
            continue;
        }
        prod1*= conf->eta_hat[3*neighbor_index + tree->ineighpos[i][k]][0];
    }
    
    first_part=pow(delta,BL)*prod1/BL;

    //this is the first part of the update

    compute_integrand(i,j,param,conf,tree);
    for (n_j = 0; n_j < N_tot; n_j++)
    {
      second_part=compute_integral(n_j,param, conf);

      conf->NNetahatij[n_j]=first_part+second_part;
      z += conf->NNetahatij[n_j];
    }

    return z;
}


double computeZi(int i, Par_type *param, Conf_type *conf, Tree_type *tree, double *z_gauss_part) //I'm using it like this: double z_ij = computeZ(i, vic[i][j], ...);
{
    int k;
    double z;

    double first_part=0.0;
    double second_part=0.0, second_part_gauss;
    double prod1=1.0;
    for (k = 0; k < C; k++)
    {
        int neighbor_index = tree->ivic[i][k];
        prod1*= conf->eta_hat[3*neighbor_index + tree->ineighpos[i][k]][0];
    }
    
    first_part=pow(delta,BL)*prod1/BL;

    //this is the first part of the update

    compute_integrand_zi(i,param,conf,tree);
    second_part=compute_integral_zi(param, conf, &second_part_gauss);

    z = first_part+second_part;
    *z_gauss_part = second_part_gauss;

    return z;
}

void BP_iteration(Par_type *param, Conf_type *conf,Tree_type * tree){

  int i,j, n_j,k;
  double z_ij;
  // Compute New Messagges
  for (k = 0; k < N; k++){
    i = (int)(FRANDOM*N);
    for (j = 0; j < C; j++){

	  z_ij = computeZij(i, tree->ivic[i][j], param, conf, tree);
	      
    for (n_j = 0; n_j < N_tot; n_j++){

	    conf->eta_hat[3*i+j][n_j] = DAMP * conf->NNetahatij[n_j]/z_ij + (1 - DAMP) * conf->eta_hat[3*i+j][n_j]; 
	    if(conf->eta_hat[3*i+j][n_j]<1e-320)
	      conf->eta_hat[3*i+j][n_j]=0;
      }
    }
  }
  
  return;
  
}


void get_true_marginals(Par_type *param, Conf_type *conf,Tree_type * tree){

  int i, n_i,k;
  double z_i;
  double z_gauss_part;
  // Compute New Messagges
  for (i = 0; i < N; i++){
    
	  z_i = computeZi(i, param, conf, tree, &z_gauss_part);
	      
    for (n_i = 0; n_i < N_tot; n_i++){

	    conf->eta_i[i][n_i] = param->integrand_zi[n_i]/z_i;
      conf->gauss_i[i][n_i] = param->integrand_gaussian_part[n_i]/z_gauss_part; 
	    if(conf->eta_i[i][n_i]<1e-320)
	      conf->eta_i[i][n_i]=0;
      if(conf->gauss_i[i][n_i]<1e-320)
	      conf->gauss_i[i][n_i]=0;

      }
    
  }
  
  return;
  
}


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


int convergence(Par_type *param, Conf_type *conf, Tree_type *tree, Obs_type *obs, 
                int consecutive_threshold){
  int iter;
  int consecutive_iterations = 0;
  for (iter = 1; iter < NITER + 1; iter++) {
    double diff_means = 0.;
    double diff_variances = 0.;
    BP_iteration(param, conf, tree);
    compute_eta_hat_MOMENTS(obs, conf);

    for (int i = 0; i < N*C; i++)
    {
      if (fabs((obs->mean[i]-obs->previous_mean[i])/obs->previous_mean[i])>diff_means)
      {
        diff_means=fabs((obs->mean[i]-obs->previous_mean[i])/obs->previous_mean[i]);          
      }      
      if (fabs((obs->variance[i]-obs->previous_variance[i])/obs->previous_variance[i])>diff_variances)
      {
        diff_variances=fabs((obs->variance[i]-obs->previous_variance[i])/obs->previous_variance[i]);
      }    
    }
          
    if (diff_means < 0.000001 && diff_variances <0.000001) {
      consecutive_iterations++;
      if (consecutive_iterations > consecutive_threshold) {
        break;
      }
    }
    else {
      consecutive_iterations = 0; 
    }

    for (int i = 0; i < N*C; i++)
    {        
      obs->previous_mean[i]=obs->mean[i];
      obs->previous_variance[i]=obs->variance[i];
    }
  
    if (iter % NPRINT_OBS == 0) {
      fprintf(stderr, "iter %d diff_m=%.6lf diff_v=%.6lf\n", iter, diff_means, diff_variances);
    }
  }
  return iter;

}


void average_true_marginal(Conf_type *conf, Obs_type *obs){
    int n_i,i;
    for (n_i = 0; n_i < N_tot; n_i++)
    {
      obs->av_eta_i[n_i]=0.;
      obs->av_gauss_i[n_i]=0.;
      for (i = 0; i < N; i++)
      {
        obs->av_eta_i[n_i]+=conf->eta_i[i][n_i];
        obs->av_gauss_i[n_i]+=conf->gauss_i[i][n_i];
      }
      obs->av_eta_i[n_i]=obs->av_eta_i[n_i]/N;
      obs->av_gauss_i[n_i]=obs->av_gauss_i[n_i]/N;
      
    }
    return;
}


void print_distr(Conf_type *conf, Obs_type *obs, FILE *file_eta, FILE *file_gauss){

  int i, n_i;
  fprintf(file_eta, "#n_i av(eta) eta_i[0] eta_i[1] ... eta_i[N-1]\n");
  for (n_i = 0; n_i < N_tot; n_i++)
  {
    fprintf(file_eta, "%lf\t%lf", conf->grid[n_i], obs->av_eta_i[n_i]);
    for (i = 0; i < N; i++){
      fprintf(file_eta, "\t%lf", conf->eta_i[i][n_i]);
    }
    fprintf(file_eta, "\n");
  }

  fprintf(file_gauss, "#n_i av(gauss) gauss_i[0] gauss_i[1] ... gauss_i[N-1]\n");
  for (n_i = 0; n_i < N_tot; n_i++)
  {
    fprintf(file_gauss, "%lf\t%lf", conf->grid[n_i], obs->av_gauss_i[n_i]);
    for (i = 0; i < N; i++){
      fprintf(file_gauss, "\t%lf", conf->gauss_i[i][n_i]);
    }
    fprintf(file_gauss, "\n");
  }
  
  return;
  
}


int main(int argc, char *argv[])
{
  Tree_type tree;
  Conf_type conf;
  Par_type param;
  Obs_type obs;
 
  if (argc != 6&&argc!=7) {
    fprintf(stderr, "usage: %s MU T N_tot N damping [seed]\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  
  MU = atof(argv[1]);
  T = atof(argv[2]);
  N_tot = atoi(argv[3]);
  N = atoi(argv[4]);
  DAMP = atof(argv[5]);
  
  if (argc == 7) {
    myrand = (unsigned)atoi(argv[6]);
    if (myrand == 2147483647)
      error("seed must be less than 2147483647");
  } else {
    FILE *devran = fopen("/dev/random","r");
    fread(&myrand, 4, 1, devran);
    fclose(devran);
  }


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


  BETA= 1./T;
  BL=BETA*LAMBDA;
  compute_field(&param,&conf);

  compute_exp_interaction(&param, &conf);

  int consecutive_threshold = 10;
  int conv = 1;
  computeZ_hat_ij_0(&conf);
  int iter = convergence(&param, &conf, &tree, &obs, consecutive_threshold);

  
  get_true_marginals(&param, &conf, &tree);
  compute_true_MOMENTS(&obs, &conf);
  compute_gauss_part_MOMENTS(&obs, &conf);

  printf("%d %d %d %.6lf %.6lf %.6lf %.6lf %.6lf %.6lf %.6lf\n", 
          atoi(argv[9]), iter, conv, obs.mean_AV, obs.variance_AV, 
          obs.true_mean_AV, obs.true_variance_AV,
          obs.gauss_part_mean_AV, obs.gauss_part_variance_AV, obs.gauss_part_skewness_AV);


  average_true_marginal(&conf, &obs);


  char filename_eta[300];
  snprintf(filename_eta, sizeof(filename_eta), "BP_cont_LV_sigma0_true_marginals_N_%d_T_%.3f_mu_%.3f.txt", N, T, MU);
  FILE *fileout_eta = fopen(filename_eta, "w");
  
  char filename_gauss[300];
  snprintf(filename_gauss, sizeof(filename_gauss), "BP_cont_LV_sigma0_true_marginals_gaussian_part_N_%d_T_%.3f_mu_%.3f.txt", N, T, MU);
  FILE *fileout_gauss = fopen(filename_gauss, "w");

  print_distr(&conf, &obs, fileout_eta, fileout_gauss);

  free_all_memory(&param, &conf, &tree, &obs);


  return 0;
}


