typedef struct parameters{ 

  double *field;
  double **exp_interaction;
  double **integrand;

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
  
  return;
  
}

void allocate_conf(Conf_type * conf, Tree_type * tree){
 
  int i;
  
  conf->eta_hat = (double **)calloc(N*C,sizeof(double *));
  for(i=0;i<N*C;i++) conf->eta_hat[i] = (double *)calloc(N_tot,sizeof(double));

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

    const double sqrt_pi = sqrt(M_PI);

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
  for(j = 0; j < N_tot; j++){
    obs->mean[i]+=conf->grid[j]*conf->eta_hat[i][j];
  }

  sum_mean += obs->mean[i]/(N*C);

  for(j = 0; j < N_tot; j++){
    obs->variance[i]+= pow(conf->grid[j]-obs->mean[i],2)*conf->eta_hat[i][j];
  }

  sum_var += obs->variance[i]/(N*C);

}

obs->mean_AV=sum_mean;
obs->variance_AV=sum_var;
    
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
    }
    free(conf->eta_hat);

    
    free(obs->mean);
    free(obs->variance);

for (i = 0; i < N_tot; i++) {
        free(param->exp_interaction[i]);
}
free(param->exp_interaction);


    return;
}

void create_grid(Conf_type *conf){

    conf->grid = (double *)calloc(N_tot,sizeof(double));

    double start = 0.0;
    double end = 1.2;

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
    for (i = 1; i < N_tot; i++) //if BL < 1, how should we consider field[0]?
    {
        param->field[i] = pow(conf->grid[i], BL-1.0) * exp(- 0.5 * BETA * (conf->grid[i]*conf->grid[i]-2*conf->grid[i]));
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
      //printf("RIFARE\n");
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


double compute_integral(int n_j,Par_type *param, Conf_type *conf){

  int i;
  double integral=0.0;
      for (i = 1; i < N_tot; i++) //j da 1 perchè integro da delta, non da 0
      {
        integral+=conf->simps_weight[i]*param->integrand[i][n_j];
      } 

    return integral;
}



double computeZij(int i, int j, Par_type *param, Conf_type *conf, Tree_type *tree) //la richiamo così: double z_ij = computeZ(i, vic[i][j], ...);
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

    //così faccio la prima parte dell'update

    compute_integrand(i,j,param,conf,tree);
    for (n_j = 0; n_j < N_tot; n_j++)
    {
      second_part=compute_integral(n_j,param, conf);

      conf->NNetahatij[n_j]=first_part+second_part;
      z += conf->NNetahatij[n_j];
    }

    return z;
}

void BP_iteration(Par_type *param, Conf_type *conf,Tree_type * tree){

  int i,j, n_i,k;
  double z_ij;
  // Compute New Messagges
  for (k = 0; k < N; k++){
    i=(int)(FRANDOM*N);
    for (j = 0; j < C; j++){

	z_ij = computeZij(i, tree->ivic[i][j], param, conf, tree);
	      
      for (n_i = 0; n_i < N_tot; n_i++){

	conf->eta_hat[3*i+j][n_i] = conf->NNetahatij[n_i]/z_ij; 
	if(conf->eta_hat[3*i+j][n_i]<1e-320)
	  conf->eta_hat[3*i+j][n_i]=0;
      }
    }
  }
  
  return;
  
}


// void compute_true_marginal_MOMENTS(Obs_type * obs){

//   int i,j, n_i;

//   double sum_mean=0;    
//   double sum2_mean=0;
    
//   double sum_var=0;    
//   double sum2_var=0;
  
//   double sum_kurt=0;    
//   double sum2_kurt=0;

//   double sum_mod_kurt=0;    
//   double sum2_mod_kurt=0;  
  
//   for(i = 0; i < N; i++){

//     obs->mean[i]=0; // mean value of n_i for each site i
//     obs->variance[i]=0; 
    
//     for(j = 0; j < N_MAX; j++){
//       obs->mean[i]+=j*obs->p_marg[i][j];
//     }


//     for(j = 0; j < N_MAX; j++){
//       obs->variance[i]+= pow(j-obs->mean[i],2)*obs->p_marg[i][j];
//     }


//     sum_mean += obs->mean[i]/N;
//     sum2_mean += obs->mean[i]*obs->mean[i]/N;

//     sum_var  += obs->variance[i]/N;
//     sum2_var += obs->variance[i]*obs->variance[i]/N;
    
//   }

//   for (n_i = 0; n_i < N_MAX; n_i++)
//   {
//     double sum_MeanP_marg=0;
//     for (i = 0; i < N; i++)
//     {
//       sum_MeanP_marg+= obs->p_marg[i][n_i];
//     }
//     obs->MeanP_marg[n_i]= sum_MeanP_marg/N; 
//   }
  
  
//   obs->mean_AV=sum_mean;
  
//   return;
  
// }