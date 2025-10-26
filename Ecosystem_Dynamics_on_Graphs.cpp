/* PREPROCESSOR INCLUSION */
//#include <cmath>  // Sono già contenuti in altri include ma per fortuna sono anche dotati di include guard
//#include <ctime>  // Forse si possono cancellare
//#include <cstring>
//#include <cstdio> 
//#include <cstdlib>
//#include <chrono>

#include "pietro_toolbox.h"
#include "graph_lv_toolbox.h"
#include "numerical_integration_genLotka-Volterra_toolbox.h"
#include "graph_from_sequence_toolbox.h"
#include "matrices_toolbox.h"

/* PREPROCESSOR DEFINITION */

/* STRUCTURE */

/* TYPE DEFINITION */

/* PROTOTYPES */

/* GLOBAL VARIABLES */

/* MAIN BEGINNING */

int main(int argc, char **argv){
    if(argc != 10) {
        fprintf(stderr,
                "Usage: %s N c mu sigma epsilon T N_ext N_prev_ext N_meas\n",
                argv[0]);
        return MY_FAIL;
    }
    clock_t begin = clock(); //Clock iniziale
    float time_spent;
    char *N_label = argv[1]; // Ricorda argv[0] è il nome dell'eseguibile!
    int N = atoi(N_label);
    char *c_label = argv[2];
    double c = atof(c_label);
    char *mu_label = argv[3];
    double mu = atof(mu_label);
    char *sigma_label = argv[4];
    double sigma = atof(sigma_label);
    char *epsilon_label = argv[5];
    double epsilon = atof(epsilon_label);
    char *T_label = argv[6];
    double T = atof(T_label);
    int N_ext = atoi(argv[7]); 
    int N_previous_ext = atoi(argv[8]);
    int N_meas = atoi(argv[9]);

    char ia_label[] = "Partially_AsymGauss";
    void (*ia_generator)(double*, double*, double, double, double) = &partyally_asymmetric_gaussian;
    //void (*ia_generator)(double*, double*, double, double, double, double, double, double) = &asymmetric_antagonistic_lognormal_extraction;
    //char gr_label[] = "Trivially";
    //double (*gr_generator)(double, double) = &UNG; //&trivial_generator; // Pointer to a function (weird syntax but somehow it has a reason)
    //double gr_value = 1.0; 
    //double growth_rate_interval[2] = {1.0, 1.0};
    //char cc_label[] = "Trivially";
    //double (*cc_generator)(double, double) = &UNG; //&trivial_generator; // Pointer to a function (weird syntax but somehow it has a reason)
    //double carrying_cap_interval[2] = {1.0, 1.0};
    //double cc_value = 1.0;
    void (*save_lv_function)(graph*, double, FILE*) = &save_LotkaVolterra_trajectories_only; // &save_LotkaVolterra_trajectories_only // &save_LotkaVolterra_nothing

    bool sparsity_flag = MY_TRUE;
    double population_factor = 1.0; // Fattore estrazione popolazione iniziale
    double lambda = 1e-06; // Tasso di Migrazione
    double h = 1e-2;
    double t_max = 4000;
    double deltat_save = 2.0;
    double log_factor = 1.1;
    //double extinction_log_threshold = log(pow(10, -6)); // Soglia di estinzione
    //double delta = 0.01; // Soglia per considerare due punti di equilibrio distinti

    int *rds, *rdsp, *fingerprint;
    int *crossindex = NULL;  // It is mandatory to have a NULL or a valid pointer to use realloc, not an uninitialized one
    int kstar, fhs;
    RealMatrix interaction_matrix(N, N);
    //double **average_traj, **std_traj, *traj_time;
    //int max_time_idx = (int) (t_max / deltat_save); // + 1; // The +1 apparently is not needed
    bool graphical, divergence;
    int j, n;  
    int divergence_counter;
    char *name_buffer, *dir_name, *remove_buffer;
    graph* ecosystem;
    time_t my_seed;
    FILE *fp, *fp_RK, *fp_eq, *fp_r, *fp_K, *fp_summary; // Pointer to File
    // FILE *fp_av, *fp_std;

    printf("Stiamo eseguendo: %s\n", argv[0]);
    printf("Orario di inizio: \n");
    system("date");
    name_buffer = my_char_malloc(CHAR_LENGHT); // Memoria liberata esplicitamente nel main
    dir_name = my_char_malloc(CHAR_LENGHT); // Memoria liberata esplicitamente nel main
    remove_buffer = my_char_malloc(CHAR_LENGHT); // Memoria liberata esplicitamente nel main
    //srand48(time(0)); // Inizializzazione generatore numeri casuali, va fatto una volta sola nel main IMPORTANTE.
        
    // Creazione directory mu, sigma, epsilon and T
    snprintf(name_buffer, sizeof(char) * CHAR_LENGHT, "mkdir epsilon_%s_%s_lambda_%g_h_%g", epsilon_label, ia_label, lambda, h); // sizeof(char)=1 so it can be omitted
    system(name_buffer); 
    // Inizializziamo la struttura generale dell'ecosistema
    ecosystem = ecosystem_initialization(N); // Memoria liberata in erase_ecosystem
    //average_traj = my_double_arofars_malloc(N, max_time_idx); // Memoria liberata in my_double_arofars_free
    // = my_double_arofars_malloc(N, max_time_idx); // Memoria liberata in my_double_arofars_free
    //traj_time = my_double_malloc(max_time_idx);
    // Creazione directory
    snprintf(dir_name, CHAR_LENGHT, "epsilon_%s_%s_lambda_%g_h_%g/N_%d_c_%.2f", epsilon_label, ia_label, lambda, h, ecosystem->size, c); // sizeof(char)=1 so it can be omitted
    create_bunch_of_directories(dir_name);
    snprintf(name_buffer, CHAR_LENGHT, "mkdir %s/Equilibrium_Points", dir_name);
    system(name_buffer);
    // Scriviamo su file i parametri del codice
    // Memory Allocation
    rds = my_int_malloc(N);
    rdsp = my_int_calloc(N);
    fingerprint = my_int_malloc(N);
    // I contatori di divergenze e multi_eq vanno resettati ad ogni taglia N dell'ecosistema
    divergence_counter = 0;
    snprintf(name_buffer, sizeof(char) * CHAR_LENGHT, "mkdir Results"); // sizeof(char)=1 so it can be omitted
    system(name_buffer); 
    snprintf(name_buffer, CHAR_LENGHT, "./Results/Lotka-Volterra_epsilon_%s_%s_lambda_%g_h_%g_N_%d_c_%.2f_mu_%s_sigma_%s_T_%s_Equilibrium_Points.txt", epsilon_label, ia_label, lambda, h, ecosystem->size, c, mu_label, sigma_label, T_label);
    fp_summary = my_open_appending_file(name_buffer); // Apertura file scrittura
    fprintf(fp_summary, "# Extraction\tMeasure\tSeed\ttime\tNum_NonEq\tav(#extinctions)\tav(<n_i>)\tstd(<n_i>)\tav(<(n_i - <n_i>)^{2}>)\tstd(<(n_i - <n_i>)^{2}>)\tcomputation time(ms)\n");        
    
    // Estraiamo N_ext volte il grafo e per ognuna di esse evolviamo N_meas volte diverse il sistema modificando ogni volta le condizioni iniziali
    for(n=(N_previous_ext+1); n<=(N_previous_ext+N_ext); n++){
        // Estrazione Matrice di Interazione
        seeds_from_dev_random_and_time(&my_seed); // Generate Random Seed
        //my_seed = 980890644513876229;
        srand48(my_seed); // The random seed is given
        // Generate the sequence
        graphical = MY_FALSE;
        while(not graphical){
            fill_array(rds, (int) c, N);
            //BubbleSort(rds, N); // Sort the degrees
            qsort(rds, N, sizeof(int), compare_desc);
            fhs = rds[0]; // First Hub Size
            if(N > fhs){
                crossindex = my_int_realloc(crossindex, fhs+1);
                kstar = build_crossing_index_table_remove_zeros(crossindex, rds, fhs, N);
                graphical = erdos_gallai_test(rds, N, crossindex, kstar);
            }
        }   // We get a graphical degree sequence
        // GENERATE INTERACTION MATRIX
        //graph_generator(interaction_matrix, rds, rdsp, fingerprint, N, p_c, weights_generator, od_p1, od_p2, epsilon);
        graph_from_degree_sequence(interaction_matrix, N, rds, rdsp, fingerprint, crossindex, kstar, ia_generator, mu, sigma, epsilon);
        // Save to file in order to generate ecosystem
        snprintf(name_buffer, CHAR_LENGHT, "%s/Extractions/Extraction_%d_Interaction_Matrix_Sparse_mu_%s_sigma_%s_T_%s.txt", dir_name, n, mu_label, sigma_label, T_label);
        fp = my_open_writing_file(name_buffer);
        save_matrix_to_file_sparse_ijaij(interaction_matrix, N, fp);
        fclose(fp);
        // Save to file for David
        snprintf(name_buffer, CHAR_LENGHT, "%s/Extractions/Extraction_%d_Interaction_Matrix_David.txt", dir_name, n);
        fp = my_open_writing_file(name_buffer);
        save_matrix_to_file_sparse_ijaijaji(interaction_matrix, N, fp);
        fclose(fp);
        // CREATE ECOSYSTEM
        snprintf(name_buffer, CHAR_LENGHT, "Array_of_ones.txt");
        generate_file_of_ones(name_buffer, N);
        //snprintf(name_buffer, CHAR_LENGHT, "Array_of_ones.txt");
        //snprintf(name_buffer, CHAR_LENGHT, "%s/Extractions/Extraction_%d_Growth_Rates.txt", dir_name, n);
        fp_r = my_open_reading_file(name_buffer);
        //snprintf(name_buffer, CHAR_LENGHT, "Array_of_ones.txt");
        //snprintf(name_buffer, CHAR_LENGHT, "%s/Extractions/Extraction_%d_Carrying_Capacities.txt", dir_name, n);
        fp_K = my_open_reading_file(name_buffer);
        snprintf(name_buffer, CHAR_LENGHT, "%s/Extractions/Extraction_%d_Interaction_Matrix_Sparse_mu_%s_sigma_%s_T_%s.txt", dir_name, n, mu_label, sigma_label, T_label);
        fp = my_open_reading_file(name_buffer);
        load_ecosystem_from_file(ecosystem, sparsity_flag, fp_r, fp_K, fp);
        fclose(fp);
        fclose(fp_r);
        fclose(fp_K);
        //my_double_arofars_reset(average_traj, N, max_time_idx);
        //my_double_arofars_reset(std_traj, N, max_time_idx);
        // MEASURES
        for(j=1; j<=N_meas; j++){
            if(j!=1){ // Use a different seed for measures after the first one
                seeds_from_dev_random_and_time(&my_seed); // Generate Random Seed
                //my_seed = 3939269325854888696;
                srand48(my_seed); // The random seed is given
            }
            // INITIAL CONDITIONS
            // Estraiamo condizioni iniziali
            //snprintf(name_buffer, CHAR_LENGHT, "%s/Initial_Conditions/Log_Initial_Condition_Extraction_%d_Measure_%d_Seed_%ld.txt", dir_name, n, j, my_seed);
            snprintf(name_buffer, CHAR_LENGHT, "%s/Initial_Conditions/Log_Initial_Condition_Extraction_%d_PopFactor_%.2f_mu_%s_sigma_%s_T_%s.txt", dir_name, n, population_factor, mu_label, sigma_label, T_label);
            fp = my_open_writing_file(name_buffer);
            extract_and_save_random_initial_conditions(ecosystem, population_factor, fp);
            //extract_and_save_delta_initial_conditions(ecosystem, population_factor, fp);
            fclose(fp);
            // Apertura file di scrittura
            snprintf(name_buffer, CHAR_LENGHT, "%s/Evolutions/Lotka-Volterra_Extraction_%d_Measure_%d_mu_%s_sigma_%s_T_%s.txt", dir_name, n, j, mu_label, sigma_label, T_label);
            fp_RK = my_open_writing_file(name_buffer); // Apertura file scrittura
            snprintf(name_buffer, CHAR_LENGHT, "%s/Equilibrium_Points/Lotka-Volterra_mu_%s_sigma_%s_T_%s_Extraction_%d_Measure_%d_Equilibrium_Points.txt", dir_name, mu_label, sigma_label, T_label, n, j);
            fp_eq = my_open_writing_file(name_buffer); // Apertura file scrittura
            // Integrazione numerica con Milstein
            //srand48(my_seed); // It is important to provide the seed at this point, so that it is used just for the noise and not for the initial condition
            //Milstein_driver_GLV_demographic_noise_with_single_species_measure(ecosystem, t_max, h, log_factor, j, &divergence, T, lambda, deltat_save, save_lv_function, fp_RK, fp_eq);
            Milstein_driver_GLV_demographic_noise_with_single_species_measure_and_summary(ecosystem, t_max, h, log_factor, j, &divergence, T, lambda, deltat_save, save_lv_function, fp_summary, fp_RK, fp_eq, my_seed, n);
            fflush(fp_RK);
            fclose(fp_RK);
            fflush(fp_eq);
            fclose(fp_eq);
            //Milstein_driver_GLV_demographic_noise_only_with_final_summary(ecosystem, t_max, h, log_factor, j, &divergence, T, lambda, deltat_save, fp_summary, my_seed, n);
            fflush(fp_summary);
            // Cleaning files
            //snprintf(name_buffer, CHAR_LENGHT, "%s/Evolutions/Lotka-Volterra_Extraction_%d_Measure_%d_mu_%s_sigma_%s_T_%s.txt", dir_name, n, j, mu_label, sigma_label, T_label);
            //snprintf(remove_buffer, CHAR_LENGHT, "rm %s", name_buffer);
            system(remove_buffer);
            printf("Extraction %d Measure %d \n", n, j);
            if(divergence == MY_TRUE){
                divergence_counter++;
            }
        } // Fine loop sulle misure
        // Calcolo traiettoria media e stampa su file
        //snprintf(name_buffer, CHAR_LENGHT, "./Results/Average_Trajectories_Lotka_Volterra_epsilon_%s_%s_lambda_%g_h_%g_N_%d_c_%.2f_mu_%s_sigma_%s_T_%s_Ext_%d_Nmeas_%d.txt", epsilon_label, ia_label, lambda, h, ecosystem->size, c, mu_label, sigma_label, T_label, n, N_meas);
        //fp_av = my_open_writing_file(name_buffer); // Apertura file scrittura
        //snprintf(name_buffer, CHAR_LENGHT, "./Results/Standard_Deviation_Trajectories_Lotka_Volterra_epsilon_%s_%s_lambda_%g_h_%g_N_%d_c_%.2f_mu_%s_sigma_%s_T_%s_Ext_%d_Nmeas_%d.txt", epsilon_label, ia_label, lambda, h, ecosystem->size, c, mu_label, sigma_label, T_label, n, N_meas);
        //fp_std = my_open_writing_file(name_buffer); // Apertura file scrittura
        //compute_and_print_average_and_std_trajectories(average_traj, std_traj, traj_time, N_meas, max_time_idx, N, fp_av, fp_std);
        // Pulizia interazioni grafo
        clean_up_graph(ecosystem);
        // Cleaning files
        snprintf(name_buffer, CHAR_LENGHT, "%s/Initial_Conditions/Log_Initial_Condition_Extraction_%d_PopFactor_%.2f_mu_%s_sigma_%s_T_%s.txt", dir_name, n, population_factor, mu_label, sigma_label, T_label);
        snprintf(remove_buffer, CHAR_LENGHT, "rm %s", name_buffer);
        system(remove_buffer);
        snprintf(name_buffer, CHAR_LENGHT, "%s/Extractions/Extraction_%d_Interaction_Matrix_Sparse_mu_%s_sigma_%s_T_%s.txt", dir_name, n, mu_label, sigma_label, T_label);
        snprintf(remove_buffer, CHAR_LENGHT, "rm %s", name_buffer);
        system(remove_buffer);
    } // Fine del ciclo sulle estrazioni
    fclose(fp_summary); // Chiudo file summary
    // Pulizia memoria
    erase_ecosystem(ecosystem); // Pulizia memoria grafo
    //my_double_arofars_free(average_traj, N);
    //my_double_arofars_free(std_traj, N);
    free(rds);
    free(rdsp);
    free(fingerprint);
    free(name_buffer); // Pulizia memoria stringa
    free(dir_name); // Pulizia memoria stringa
    free(remove_buffer);
    // Calcolo tempo impiegato
    clock_t end = clock(); //Clock conclusivo
    time_spent = ((float)(end - begin)) / CLOCKS_PER_SEC; //Calcolo tempo di esecuzione
    printf("\n # Il programma ha impiegato %f secondi.\n # A breve dovrebbe arrivare la mail (se opportunamente richiesto presso i nostri uffici).\n \n", time_spent);
    printf("Orario di conclusione: \n");
    system("date");

    return MY_SUCCESS;
}

/* MAIN - THE END */

/* FUNCTIONS */