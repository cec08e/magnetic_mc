#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf_pow_int.h>
#include <time.h>
#include "h3d_nl.h"


/*
Requirements: GNU Scientific Library (gsl) and CBLAS

gcc -lgsl -lgslcblas h3d_nl.c
ALT: gcc -fPIC -shared -o h3d_nl.so -lgsl -lgslcblas h3d_nl.c

************ N LAYER VERSION -- OFFICIAL, UPDATED *****************
*/


int main(){
  parse_config_file();
  echo_params(stdout);
  build_lattice();
}

void parse_config_file(){
  /* to do: make structure more permissable - whitespace allowance, etc. */
  /* MUST ALLOCATE SPACE TO PARAM LISTS BEFORE ACCEPTING VALUES */

  FILE* fp;
  char config_line[80];
  char param[80];
  char val_list[80], *str_val;
  int val, i;
  float val_d;

  fp = fopen(CONFIG_FILE, "r");
  while(fgets(config_line, 80,fp)){
    if(sscanf(config_line, "%s = %d\n", param, &val) == 2){
      /* Case includes: SIM_NUM, NUM_L, NUM_R, NUM_C, OVER_FLAG, ANNEAL_TIME,
      EQ_TIME, COR_TIME, J_INTRA_P, J_INTER_P, D_INTRA_P, D_INTER_P, K_P*/
      if(!strcmp(param, "SIM_NUM")) SIM_NUM = val;
      else if(!strcmp(param, "NUM_L")){
         NUM_L = val;
         /* Number of layers now known. Allocate space for constant
         interaction strength parameters */
         J_INTRA_CONST = malloc(NUM_L*sizeof(float));
         J_INTER_CONST = malloc(NUM_L*sizeof(float));
         D_INTRA_CONST = malloc(NUM_L*sizeof(float));
         D_INTER_CONST = malloc(NUM_L*sizeof(float));
         K_CONST = malloc(NUM_L*sizeof(float));

       }
      else if(!strcmp(param, "NUM_R")) NUM_R = val;
      else if(!strcmp(param, "NUM_C")) NUM_C = val;
      else if(!strcmp(param, "OVER_FLAG")) OVER_FLAG = val;
      else if(!strcmp(param, "ANNEAL_TIME")) ANNEAL_TIME = val;
      else if(!strcmp(param, "EQ_TIME")) EQ_TIME = val;
      else if(!strcmp(param, "COR_TIME")) COR_TIME = val;
      else if(!strcmp(param, "J_INTRA_P")) J_INTRA_P = val;
      else if(!strcmp(param, "J_INTER_P")) J_INTER_P = val;
      else if(!strcmp(param, "D_INTRA_P")) D_INTRA_P = val;
      else if(!strcmp(param, "D_INTER_P")) J_INTER_P = val;
      else if(!strcmp(param, "K_P")) K_P = val;
    }
    else if(sscanf(config_line, "%s = %f", param, &val_d) == 2){
      /* Case includes: B_CONST, J_INTRA_CONST, J_INTER_CONST, D_INTRA_CONST, D_INTER_CONST, K_CONST*/
      if(!strcmp(param, "INIT_T")) INIT_T = val_d;
      else if(!strcmp(param, "FINAL_T")) FINAL_T = val_d;
      else if(!strcmp(param, "DELTA_T")) DELTA_T = val_d;
      else if(!strcmp(param, "B_CONST")) B_CONST = val_d;
      else if(!strcmp(param, "DELTA_B")) DELTA_B = val_d;
      else if(!strcmp(param, "RADIUS")) RADIUS = val_d;

    }
    else if(sscanf(config_line, "%s = [%s]\n", param, val_list) == 2){
      i=0;
      str_val = strtok(val_list,",");
      if(!strcmp(param, "J_INTRA_CONST")){
        while (str_val!= NULL){
          J_INTRA_CONST[i++] = atof(str_val);
          str_val = strtok (NULL, ",");
        }
      }
      else if(!strcmp(param, "J_INTER_CONST")){
        while (str_val != NULL){
          J_INTER_CONST[i++] = atof(str_val);
          str_val = strtok (NULL, ",");
        }
      }
      else if(!strcmp(param, "D_INTRA_CONST")){
        while (str_val != NULL){
          D_INTRA_CONST[i++] = atof(str_val);
          str_val = strtok (NULL, ",");
        }
      }
      else if(!strcmp(param, "D_INTER_CONST")){
        while (str_val != NULL){
          D_INTER_CONST[i++] = atof(str_val);
          str_val = strtok (NULL, ",");
        }
      }
      else if(!strcmp(param, "K_CONST")){
        while (str_val != NULL){
          K_CONST[i++] = atof(str_val);
          str_val = strtok (NULL, ",");
        }
      }
    }
    else{
      printf("Invalid line encountered in configuration file: %s\n", config_line);
    }

  }
}

void echo_params(FILE *fp){
  int i = 0;
  /* Print the parameter list to fp */
  fprintf(fp, "Simulation number: %d\n", SIM_NUM);
  fprintf(fp, "Number of layers: %d\n", NUM_L);
  fprintf(fp, "Number of rows: %d\n", NUM_R);
  fprintf(fp, "Number of columns: %d\n", NUM_C);
  fprintf(fp, "Initial temperature: %f\n", INIT_T);
  fprintf(fp, "Final temperature: %f\n", FINAL_T);
  fprintf(fp, "Temperature annealing interval dT: %f\n", DELTA_T);
  fprintf(fp, "Initial external field strength B = %f\n", B_CONST);
  fprintf(fp, "Field scanning interval dB = %f\n", DELTA_B);
  fprintf(fp, "Anneal time: %d\n", ANNEAL_TIME);
  fprintf(fp, "Equilibration time: %d\n", EQ_TIME);
  fprintf(fp, "Correlation time: %d\n", COR_TIME);

  fprintf(fp, "J_INTRA_CONST: [");
  for(i = 0; i < NUM_L; i++){
    fprintf(fp, " %f ", J_INTRA_CONST[i]);
  }
  fprintf(fp, "]\n");

  fprintf(fp, "J_INTER_CONST: [");
  for(i = 0; i < NUM_L; i++){
    fprintf(fp, " %f ", J_INTER_CONST[i]);
  }
  fprintf(fp, "]\n");

  fprintf(fp, "D_INTRA_CONST: [");
  for(i = 0; i < NUM_L; i++){
    fprintf(fp, " %f ", D_INTRA_CONST[i]);
  }
  fprintf(fp, "]\n");

  fprintf(fp, "D_INTER_CONST: [");
  for(i = 0; i < NUM_L; i++){
    fprintf(fp, " %f ", D_INTER_CONST[i]);
  }
  fprintf(fp, "]\n");

  fprintf(fp, "K_CONST: [");
  for(i = 0; i < NUM_L; i++){
    fprintf(fp, " %f ", K_CONST[i]);
  }
  fprintf(fp, "]\n");

  fprintf(fp, "J_INTRA_P: %d\n", J_INTRA_P);
  fprintf(fp, "J_INTER_P: %d\n", J_INTER_P);
  fprintf(fp, "D_INTRA_P: %d\n", D_INTRA_P);
  fprintf(fp, "D_INTER_P: %d\n", D_INTER_P);
  fprintf(fp, "K_P: %d\n", K_P);
  fprintf(fp, "RADIUS: %f\n", RADIUS);

}

void build_lattice(){
  /* Allocate memory for the lattice structure */
  /* Allocate mem and set values for the parameter interaction strengths

  To do: define path for custom periodic functions
  defining interaction strengths, initialize D vector
  */
  /* Generate random spins to populate the lattice */

  int i, j, k;
  lattice_t lattice = (spin_t ***)malloc(NUM_L*sizeof(spin_t **));
  J_INTER = (float ***)malloc(NUM_L*sizeof(float **));
  J_INTRA = (float ***)malloc(NUM_L*sizeof(float **));
  D_INTER = (float ***)malloc(NUM_L*sizeof(float **));
  D_INTRA = (float ***)malloc(NUM_L*sizeof(float **));
  K = (float ***)malloc(NUM_L*sizeof(float **));
  B = (float ***)malloc(NUM_L*sizeof(float **));
  for(i = 0; i < NUM_L; i++){
    lattice[i] = (spin_t **)malloc(NUM_R*sizeof(spin_t *));
    J_INTER[i] = (float **)malloc(NUM_R*sizeof(float *));
    J_INTRA[i] = (float **)malloc(NUM_R*sizeof(float *));
    D_INTER[i] = (float **)malloc(NUM_R*sizeof(float *));
    D_INTRA[i] = (float **)malloc(NUM_R*sizeof(float *));
    K[i] = (float **)malloc(NUM_R*sizeof(float *));
    B[i] = (float **)malloc(NUM_R*sizeof(float *));

    for(j = 0; j < NUM_R; j++){
      lattice[i][j] = (spin_t *)malloc(NUM_C*sizeof(spin_t));
      J_INTER[i][j] = (float *)malloc(NUM_C*sizeof(float));
      J_INTRA[i][j] = (float *)malloc(NUM_C*sizeof(float));
      D_INTER[i][j] = (float *)malloc(NUM_C*sizeof(float));
      D_INTRA[i][j] = (float *)malloc(NUM_C*sizeof(float));
      K[i][j] = (float *)malloc(NUM_C*sizeof(float));
      B[i][j] = (float *)malloc(NUM_C*sizeof(float));

    }
  }


  if(!J_INTRA_P){
    for(i = 0; i < NUM_L; i++)
      for(j = 0; j < NUM_R; j++)
        for(k = 0; k < NUM_C; k++)
          J_INTRA[i][j][k] = J_INTRA_CONST[i];
  }

  if(!J_INTER_P){
    for(i = 0; i < NUM_L; i++)
      for(j = 0; j < NUM_R; j++)
        for(k = 0; k < NUM_C; k++)
          J_INTER[i][j][k] = J_INTER_CONST[i];
  }

  if(!D_INTRA_P){
    for(i = 0; i < NUM_L; i++)
      for(j = 0; j < NUM_R; j++)
        for(k = 0; k < NUM_C; k++)
          D_INTRA[i][j][k] = D_INTRA_CONST[i];
  }

  if(!D_INTER_P){
    for(i = 0; i < NUM_L; i++)
      for(j = 0; j < NUM_R; j++)
        for(k = 0; k < NUM_C; k++)
          D_INTER[i][j][k] = D_INTER_CONST[i];
  }

  if(!K_P){
    for(i = 0; i < NUM_L; i++)
      for(j = 0; j < NUM_R; j++)
        for(k = 0; k < NUM_C; k++)
          K[i][j][k] = K_CONST[i];
  }

  rng = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set (rng, time(NULL));
  for(i = 0; i < NUM_L; i++){
    for(j = 0; j < NUM_R; j++){
      for(k = 0; k < NUM_C; k++){
        gen_random_spin(&lattice[i][j][k]);
      }
    }
}

}

void gen_random_spin(spin_t* spin){
    double x1, x2, mag_sq;
    gsl_rng_uniform(rng);

    x1 = gsl_rng_uniform(rng)*2.0 - 1.0;
    x2 = gsl_rng_uniform(rng)*2.0 - 1.0;
    mag_sq = gsl_pow_2(x1) + gsl_pow_2(x2);

    while(mag_sq >= 1){
      x1 = gsl_rng_uniform(rng)*2.0 - 1.0;
      x2 = gsl_rng_uniform(rng)*2.0 - 1.0;
      mag_sq = gsl_pow_2(x1) + gsl_pow_2(x2);
    }

    spin->x = 2.0*x1*sqrt(1.0-mag_sq);
    spin->y = 2.0*x2*sqrt(1.0-mag_sq);
    spin->z = 1.0-2.0*mag_sq;

}
