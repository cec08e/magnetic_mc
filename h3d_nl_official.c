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

/**** To do: manage overrelaxation */


int main(){
  int i = 0;
  double ** results = (double **)malloc(10000*sizeof(double *));
  for(i = 0; i < 10000; i++)
    results[i] = (double *)malloc(8*sizeof(double));

  parse_config_file();
  echo_params(stdout);
  build_lattice();
  M_v_B(results);
  //X_v_T(results);
  //C_v_T(results);

  cleanup();
}

void parse_config_file(){
  /* to do: make structure more permissable - whitespace allowance, etc. */
  /* MUST ALLOCATE SPACE TO PARAM LISTS BEFORE ACCEPTING VALUES */

  /* FIX TYPING ISSUE */

  if(DEBUG)
    printf("Parsing config file...");

  FILE* fp;
  char config_line[80];
  char param[80];
  char val_list[80], *str_val;
  int val, i;
  float val_d;

  fp = fopen(CONFIG_FILE, "r");
  while(fgets(config_line, 80,fp)){

    if(sscanf(config_line, "%s = %f", param, &val_d) == 2){
      /* Case includes: B_CONST, J_INTRA_CONST, J_INTER_CONST, D_INTRA_CONST, D_INTER_CONST, K_CONST*/
      if(!strcmp(param, "INIT_T")) INIT_T = val_d;
      else if(!strcmp(param, "FINAL_T")) FINAL_T = val_d;
      else if(!strcmp(param, "DELTA_T")) DELTA_T = val_d;
      else if(!strcmp(param, "B_CONST")){ printf("val_d = %f\n", val_d); B_CONST = val_d;}
      else if(!strcmp(param, "DELTA_B")) DELTA_B = val_d;
      else if(!strcmp(param, "RADIUS")) RADIUS = val_d;
      else{
        printf("param = %s\n", param);
      }

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
    else if(sscanf(config_line, "%s = %d\n", param, &val) == 2){
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
      else if(!strcmp(param, "DEBUG")) DEBUG = val;
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
  fprintf(fp, "Overrelaxation: %d\n", OVER_FLAG);
  fprintf(fp, "Debugging statements: %d\n", DEBUG);


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

  /* allocate spin vectors for linear algebra */

  if(DEBUG)
    printf("Building lattice...\n");

  int i, j, k;
  lattice = (spin_t ***)malloc(NUM_L*sizeof(spin_t **));
  J_INTER = (float ***)malloc(NUM_L*sizeof(float **));
  J_INTRA = (float ***)malloc(NUM_L*sizeof(float **));
  D_INTER = (float ***)malloc(NUM_L*sizeof(float **));
  D_INTRA = (float ***)malloc(NUM_L*sizeof(float **));
  K = (float ***)malloc(NUM_L*sizeof(float **));
  //B = (float ***)malloc(NUM_L*sizeof(float **));
  for(i = 0; i < NUM_L; i++){
    lattice[i] = (spin_t **)malloc(NUM_R*sizeof(spin_t *));
    J_INTER[i] = (float **)malloc(NUM_R*sizeof(float *));
    J_INTRA[i] = (float **)malloc(NUM_R*sizeof(float *));
    D_INTER[i] = (float **)malloc(NUM_R*sizeof(float *));
    D_INTRA[i] = (float **)malloc(NUM_R*sizeof(float *));
    K[i] = (float **)malloc(NUM_R*sizeof(float *));
    //B[i] = (float **)malloc(NUM_R*sizeof(float *));

    for(j = 0; j < NUM_R; j++){
      lattice[i][j] = (spin_t *)malloc(NUM_C*sizeof(spin_t));
      J_INTER[i][j] = (float *)malloc(NUM_C*sizeof(float));
      J_INTRA[i][j] = (float *)malloc(NUM_C*sizeof(float));
      D_INTER[i][j] = (float *)malloc(NUM_C*sizeof(float));
      D_INTRA[i][j] = (float *)malloc(NUM_C*sizeof(float));
      K[i][j] = (float *)malloc(NUM_C*sizeof(float));
      //B[i][j] = (float *)malloc(NUM_C*sizeof(float));

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

  B = B_CONST;

  rng = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set (rng, time(NULL));
  for(i = 0; i < NUM_L; i++){
    for(j = 0; j < NUM_R; j++){
      for(k = 0; k < NUM_C; k++){
        gen_random_spin(&lattice[i][j][k]);
      }
    }
  }

  spin_vector = gsl_vector_alloc(3);
  point_vector = gsl_vector_alloc(3);
  temp_vector = gsl_vector_calloc(3);
  rot_matrix = gsl_matrix_alloc(3,3);

  delta_vector = gsl_vector_alloc(3);
  neighbor_vector = gsl_vector_alloc(3);
  cross_temp = gsl_vector_alloc(3);
  inter_vector = gsl_vector_alloc(3);

  neighbors[0] = gsl_vector_alloc(3);
  neighbors[1] = gsl_vector_alloc(3);
  neighbors[2] = gsl_vector_alloc(3);
  neighbors[3] = gsl_vector_alloc(3);

  D_vec[0] = gsl_vector_alloc(3);
  D_vec[1] = gsl_vector_alloc(3);
  D_vec[2] = gsl_vector_alloc(3);
  D_vec[3] = gsl_vector_alloc(3);

  gsl_vector_set(D_vec[0], 0, 1);
  gsl_vector_set(D_vec[0], 1, 0);
  gsl_vector_set(D_vec[0], 2, 0);


  gsl_vector_set(D_vec[1], 0, -1);
  gsl_vector_set(D_vec[1], 1, 0);
  gsl_vector_set(D_vec[1], 2, 0);


  gsl_vector_set(D_vec[2], 0, 0);
  gsl_vector_set(D_vec[2], 1, -1);
  gsl_vector_set(D_vec[2], 2, 0);


  gsl_vector_set(D_vec[3], 0, 0);
  gsl_vector_set(D_vec[3], 1, 1);
  gsl_vector_set(D_vec[3], 2, 0);

  /* Calculate initial energy and magnetization */
  magnetization = calc_magnetization(-1);
  energy = calc_energy();

  if(DEBUG){
    printf("Initial magnetization: %f \n", magnetization);
    printf("Initial energy: %f \n", energy);

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

void simulate(int num_sweeps, double T){
  int i, num_or, num_accept = 0;
  for(i = 0; i < num_sweeps; i++){
    num_accept += sweep(T);
    //for(num_or = 0; num_or < OVER_FLAG; num_or++){
    //  overrelax();
    //}
  }
}
/*
void overrelax(){
  spin_t temp;
  int i, j, k;
  for(i = 0; i < NUM_L; i++){
    for(j=0; j < NUM_R; j++)
      for(k = 0; k < NUM_C; k++){
        eff_project(&temp, i, j, k);
        lattice[i][j][k].x = 2*temp.x - lattice[i][j][k].x;
        lattice[i][j][k].y = 2*temp.y - lattice[i][j][k].y;
        lattice[i][j][k].z = 2*temp.z - lattice[i][j][k].z;
      }
  }
}
*/

int sweep(double T){
  int i, num_accept = 0;
  int layer, row, col;
  double delta_E;
  double random_num;
  spin_t temp_spin;

  /* Perform as many MC steps as there are lattice sites */
  for(i=0; i < NUM_L*NUM_R*NUM_C; i++){

    /* Choose a random spin site on a random lattice */
    layer = gsl_rng_get(rng) % NUM_L;
    row = gsl_rng_get(rng) % NUM_R;
    col = gsl_rng_get(rng) % NUM_C;

    /* Generate a perturbation of the spin at the chosen lattice site */
    perturb_spin(&temp_spin, &lattice[layer][row][col]);
    delta_E = calc_delta_E(&temp_spin, &lattice[layer][row][col], layer, row, col);
    random_num = gsl_rng_uniform(rng);
    if (delta_E/T > 600){
      continue;
    }
    else if( !( (delta_E > 0) && (random_num >= gsl_sf_exp(-(1.0/T)*delta_E)) ) ){
      magnetization += temp_spin.z - lattice[layer][row][col].z;    /* Update magnetization */
      energy += delta_E;
      lattice[layer][row][col].x = temp_spin.x;
      lattice[layer][row][col].y = temp_spin.y;
      lattice[layer][row][col].z = temp_spin.z;


      num_accept += 1;
    }
  }

  return num_accept;

}


void perturb_spin(spin_t* temp_spin, spin_t* spin){
  double r, arg;
  double x, y, z, x_new, y_new, u_x, u_y;
  double theta, phi;
  double norm;

  gsl_vector_set(spin_vector, 0, spin->x);
  gsl_vector_set(spin_vector, 1, spin->y);
  gsl_vector_set(spin_vector, 2, spin->z);

  /* Generate random position on the tangent disc */
  r = sqrt(gsl_rng_uniform(rng))*RADIUS;
  arg = gsl_rng_uniform(rng)*2*M_PI;

  /* Express random point as cartesian coordinate in same reference frame as spin */
  x = r*sin(arg);
  y = r*cos(arg);
  z = 0;

  /* Grab spherical coordinates of spin vector */

  theta = acos(spin->z);
  if(spin->x == 0){
    if(spin->y > 0)
      phi = M_PI/2.0;
    else
      phi = 3*M_PI/2.0;
  }
  else
    phi = atan(spin->y/spin->x);

  /* Rotate random point with phi */
  x_new = x*cos(phi) - y*sin(phi);
  y_new = x*sin(phi) + y*cos(phi);

  gsl_vector_set(point_vector, 0, x_new);
  gsl_vector_set(point_vector, 1, y_new);
  gsl_vector_set(point_vector, 2, z);

  /* Now, rotate random point with theta - rotate around y' = -sin(phi), cos(phi) axis */
  u_x = -sin(phi);
  u_y = cos(phi);

  /* Creating rotation matrix */
  gsl_matrix_set(rot_matrix, 0, 0, cos(theta) +  gsl_pow_2(u_x)*(1-cos(theta)));
  gsl_matrix_set(rot_matrix, 0, 1, u_x*u_y*(1-cos(theta)));
  gsl_matrix_set(rot_matrix, 0, 2, u_y*sin(theta));
  gsl_matrix_set(rot_matrix, 1, 0, u_x*u_y*(1-cos(theta)));
  gsl_matrix_set(rot_matrix, 1, 1, cos(theta) + gsl_pow_2(u_y)*(1-cos(theta)));
  gsl_matrix_set(rot_matrix, 1, 2, -u_x*sin(theta));
  gsl_matrix_set(rot_matrix, 2, 0, -u_y*sin(theta));
  gsl_matrix_set(rot_matrix, 2, 1, u_x*sin(theta));
  gsl_matrix_set(rot_matrix, 2, 2,  cos(theta));


  gsl_blas_dgemv(CblasNoTrans,1.0, rot_matrix, point_vector,0.0,temp_vector);
  gsl_vector_add(temp_vector, spin_vector);
  norm = sqrt(gsl_pow_2(gsl_vector_get(temp_vector, 0))
              + gsl_pow_2(gsl_vector_get(temp_vector, 1))
              + gsl_pow_2(gsl_vector_get(temp_vector, 2)));
  gsl_vector_scale(temp_vector, 1/norm);

  temp_spin->x = gsl_vector_get(temp_vector, 0);
  temp_spin->y = gsl_vector_get(temp_vector, 1);
  temp_spin->z = gsl_vector_get(temp_vector, 2);

}

double calc_delta_E(spin_t* temp_spin, spin_t* spin, int layer, int row, int col){
  double delta_dot_neighbor;
  double delta_dot_inter;
  double delta_a;
  double delta_E;
  double delta_D;
  double temp;
  int i;

  /* FIRST TERM */
  /* Calculate change in spin */
  gsl_vector_set(delta_vector, 0, temp_spin->x - spin->x);
  gsl_vector_set(delta_vector, 1, temp_spin->y - spin->y);
  gsl_vector_set(delta_vector, 2, temp_spin->z - spin->z);

  /* Calculate neighbor sum */
  gsl_vector_set(neighbor_vector, 0, lattice[layer][(((row-1)%NUM_R) + NUM_R) % NUM_R][col].x
                                    + lattice[layer][(((row+1)%NUM_R) + NUM_R) % NUM_R][col].x
                                    + lattice[layer][row][(((col-1)%NUM_C) + NUM_C) % NUM_C].x
                                    + lattice[layer][row][(((col+1)%NUM_C) + NUM_C) % NUM_C].x);
  gsl_vector_set(neighbor_vector, 1, lattice[layer][(((row-1)%NUM_R) + NUM_R) % NUM_R][col].y
                                    + lattice[layer][(((row+1)%NUM_R) + NUM_R) % NUM_R][col].y
                                    + lattice[layer][row][(((col-1)%NUM_C) + NUM_C) % NUM_C].y
                                    + lattice[layer][row][(((col+1)%NUM_C) + NUM_C) % NUM_C].y);
  gsl_vector_set(neighbor_vector, 2, lattice[layer][(((row-1)%NUM_R) + NUM_R) % NUM_R][col].z
                                    + lattice[layer][(((row+1)%NUM_R) + NUM_R) % NUM_R][col].z
                                    + lattice[layer][row][(((col-1)%NUM_C) + NUM_C) % NUM_C].z
                                    + lattice[layer][row][(((col+1)%NUM_C) + NUM_C) % NUM_C].z);

  gsl_blas_ddot(delta_vector, neighbor_vector, &delta_dot_neighbor);
  /* END FIRST TERM */

  /* SECOND TERM */
  gsl_vector_set(inter_vector, 0, J_INTER[layer][row][col]*lattice[(layer+1)%NUM_L][row][col].x + J_INTER[(((layer-1)%NUM_L)+NUM_L)%NUM_L][row][col]*lattice[(((layer-1)%NUM_L)+NUM_L)%NUM_L][row][col].x);
  gsl_vector_set(inter_vector, 1, J_INTER[layer][row][col]*lattice[(layer+1)%NUM_L][row][col].y + J_INTER[(((layer-1)%NUM_L)+NUM_L)%NUM_L][row][col]*lattice[(((layer-1)%NUM_L)+NUM_L)%NUM_L][row][col].y);
  gsl_vector_set(inter_vector, 2, J_INTER[layer][row][col]*lattice[(layer+1)%NUM_L][row][col].z + J_INTER[(((layer-1)%NUM_L)+NUM_L)%NUM_L][row][col]*lattice[(((layer-1)%NUM_L)+NUM_L)%NUM_L ][row][col].z);

  gsl_blas_ddot(delta_vector, inter_vector, &delta_dot_inter);
  /* END SECOND TERM */

  /* THIRD TERM */
  /* Calculate anisotropy change */
  delta_a = K[layer][row][col]*(gsl_pow_2(temp_spin->z) - gsl_pow_2(spin->z));
  /* END THIRD TERM */

  /* FOURTH TERM */
  /* TO DO: DELTA D CALCULATION */

  /* NORTH */

  gsl_vector_set(neighbors[0], 0, lattice[layer][(((row-1)%NUM_R) + NUM_R) % NUM_R][col].x);
  gsl_vector_set(neighbors[0], 1, lattice[layer][(((row-1)%NUM_R) + NUM_R) % NUM_R][col].y);
  gsl_vector_set(neighbors[0], 2, lattice[layer][(((row-1)%NUM_R) + NUM_R) % NUM_R][col].z);

  /* SOUTH */
  gsl_vector_set(neighbors[1], 0, lattice[layer][(((row+1)%NUM_R) + NUM_R) % NUM_R][col].x);
  gsl_vector_set(neighbors[1], 1, lattice[layer][(((row+1)%NUM_R) + NUM_R) % NUM_R][col].y);
  gsl_vector_set(neighbors[1], 2, lattice[layer][(((row+1)%NUM_R) + NUM_R) % NUM_R][col].z);

  /* EAST */
  gsl_vector_set(neighbors[2], 0, lattice[layer][row][(((col+1)%NUM_C) + NUM_C) % NUM_C].x);
  gsl_vector_set(neighbors[2], 1, lattice[layer][row][(((col+1)%NUM_C) + NUM_C) % NUM_C].y);
  gsl_vector_set(neighbors[2], 2, lattice[layer][row][(((col+1)%NUM_C) + NUM_C) % NUM_C].z);

  /* WEST */
  gsl_vector_set(neighbors[3], 0, lattice[layer][row][(((col-1)%NUM_C) + NUM_C) % NUM_C].x);
  gsl_vector_set(neighbors[3], 1, lattice[layer][row][(((col-1)%NUM_C) + NUM_C) % NUM_C].y);
  gsl_vector_set(neighbors[3], 2, lattice[layer][row][(((col-1)%NUM_C) + NUM_C) % NUM_C].z);

  delta_D = 0;

  for(i = 0; i < 4; i++){
    cross_product(delta_vector, neighbors[i], cross_temp);
    gsl_blas_ddot(D_vec[i], cross_temp, &temp);
    delta_D += D_INTRA[layer][row][col]*temp;
  }

  /* END FOURTH TERM */

  //delta_E = -J_INTRA*delta_dot_neighbor + J_INTER*delta_dot_inter + delta_a - B_EXT*gsl_vector_get(delta_vector,2);
  delta_E = -J_INTRA[layer][row][col]*delta_dot_neighbor + delta_dot_inter + delta_a + delta_D - B*gsl_vector_get(delta_vector,2);

  return delta_E;

}

void cool_lattice(double T){
  float curr_temp;
  curr_temp = INIT_T;
  while(curr_temp > T){
    if(DEBUG)
      printf("Annealing to T = %f\n", curr_temp);
    fflush(stdout);
    simulate(ANNEAL_TIME, curr_temp);
    curr_temp -= DELTA_T;
  }

}

double calc_energy(){
  int i,j,k;
  spin_t zero_spin;
  zero_spin.x = 0;
  zero_spin.y = 0;
  zero_spin.z = 0;
  for(i=0; i < NUM_L; i++)
      for(j=0; j < NUM_R; j++)
          for(k = 0; k < NUM_C; k++)
              energy += calc_delta_E(&lattice[i][j][k], &zero_spin, i, j, k);

  return energy/2.0;
}


double calc_magnetization(int layer){
      float mag, mag_spin;
      int i,j,k;
      mag = 0.0;
      mag_spin = 0.0;
      if(layer == -1){
        for(i=0; i < NUM_L; i++)
            for(j=0; j < NUM_R; j++)
                for(k = 0; k < NUM_C; k++)
                    mag += lattice[i][j][k].z;
        //mag_spin = mag/(NUM_R*NUM_C*NUM_L);
        mag_spin = mag;
      }
      else{
        for(j=0; j < NUM_R; j++)
            for(k = 0; k < NUM_C; k++)
                mag += lattice[layer][j][k].z;
        mag_spin = mag/(NUM_R*NUM_C);
      }
      return mag_spin;
}

double calc_TC(int layer){
  double solid_angle_sum = 0;
  int i, j, k;
  i = layer;
  for(j =0; j < NUM_R; j++)
    for(k=0; k < NUM_C; k++){
      solid_angle_sum += calc_solid_angle(lattice[i][j][k], lattice[i][j][(((k-1)%NUM_C) + NUM_C) % NUM_C], lattice[i][(j+1)%NUM_R][k]);
      solid_angle_sum += calc_solid_angle(lattice[i][j][k], lattice[i][(j)%NUM_R][(k+1)%NUM_C], lattice[i][(((j-1)%NUM_R) + NUM_R) % NUM_R][(k)%NUM_C]);
    }
  //printf("returning sa\n");
  //fflush(stdout);
  return solid_angle_sum/(4*M_PI);
}

double calc_solid_angle(spin_t n1, spin_t n2, spin_t n3){

  gsl_complex c_temp;
  gsl_vector * n1_vec = gsl_vector_alloc(3);
  gsl_vector * n2_vec = gsl_vector_alloc(3);
  gsl_vector * n3_vec = gsl_vector_alloc(3);
  gsl_vector * n2_cross_n3 = gsl_vector_alloc(3);

  double n1_dot_n2;
  double n2_dot_n3;
  double n3_dot_n1;
  double n1_dot_n2_cross_n3;

  double rho;
  double Omega;


  gsl_vector_set(n1_vec, 0, n1.x);
  gsl_vector_set(n1_vec, 1, n1.y);
  gsl_vector_set(n1_vec, 2, n1.z);

  gsl_vector_set(n2_vec, 0, n2.x);
  gsl_vector_set(n2_vec, 1, n2.y);
  gsl_vector_set(n2_vec, 2, n2.z);

  gsl_vector_set(n3_vec, 0, n3.x);
  gsl_vector_set(n3_vec, 1, n3.y);
  gsl_vector_set(n3_vec, 2, n3.z);

  cross_product(n2_vec, n3_vec, n2_cross_n3);

  gsl_blas_ddot(n1_vec, n2_vec, &n1_dot_n2);
  gsl_blas_ddot(n2_vec, n3_vec, &n2_dot_n3);
  gsl_blas_ddot(n3_vec, n1_vec, &n3_dot_n1);
  gsl_blas_ddot(n1_vec, n2_cross_n3, &n1_dot_n2_cross_n3);

  //printf("n1_dot_n2: %f\n", n1_dot_n2);
  //printf("n2_dot_n3: %f\n", n2_dot_n3);
  //printf("n3_dot_n1: %f\n", n3_dot_n1);

  rho = sqrt(2*(1+n1_dot_n2)*(1+n2_dot_n3)*(1+n3_dot_n1));
  //printf("Rho is %f \n", rho);
  GSL_SET_COMPLEX(&c_temp, (1.0/rho)*(1 + n1_dot_n2 + n2_dot_n3 + n3_dot_n1), (1.0/rho)*n1_dot_n2_cross_n3);
  Omega = 2*GSL_IMAG(gsl_complex_log(c_temp));

  gsl_vector_free(n1_vec);
  gsl_vector_free(n2_vec);
  gsl_vector_free(n3_vec);
  gsl_vector_free(n2_cross_n3);

  return Omega;

}


int C_v_T(double** results){
  /* Calc Heat Capacity from Init temp to final temp, in delta_t intervals
  where C = (<E^2> - <E>^2)/T^{2}

  At each temperature, the energy is measured over 100,000 cor time samples.
  */
  float curr_temp = INIT_T;
  int j = 0, i;
  int num_samples = 10000;
  double e_avg = 0, e_sq_avg = 0;
  while(curr_temp > FINAL_T){
    //if(DEBUG)
    //  printf("Calculating magnetic susceptibility at T = %f\n", curr_temp);
    fflush(stdout);
    simulate(EQ_TIME, curr_temp);
    for(i = 0; i < num_samples; i++){
      simulate(COR_TIME, curr_temp);
      e_avg += energy;
      e_sq_avg += gsl_pow_2(energy);
    }
    //printf("Energy per site at temp T = %f: %f\n", curr_temp, energy/(NUM_L*NUM_R*NUM_C));
    results[j][0] = curr_temp;
    results[j][1] = ((e_sq_avg/num_samples) - (gsl_pow_2(e_avg/num_samples)))/(gsl_pow_2(curr_temp)*NUM_L*NUM_R*NUM_C);

    printf("%f,%f\n", results[j][0], results[j][1]);

    curr_temp -= DELTA_T;
    j++;

    e_avg = 0;
    e_sq_avg = 0;
  }
  return j;
}

int X_v_T(double** results){
  /* Calc magnetic susceptibility from Init temp to final temp, in delta_t intervals
  where C = (<M^2> - <M>^2)/T

  At each temperature, X is measured over 10,000 cor time samples.
  */

  float curr_temp = INIT_T;
  int j = 0, i;
  int num_samples = 10000;
  double m_avg = 0, m_sq_avg = 0;
  while(curr_temp > FINAL_T){
    //if(DEBUG)
    //  printf("Calculating magnetic susceptibility at T = %f\n", curr_temp);
    fflush(stdout);
    simulate(EQ_TIME, curr_temp);
    for(i = 0; i < num_samples; i++){
      simulate(COR_TIME, curr_temp);
      m_avg += fabs(magnetization);
      m_sq_avg += gsl_pow_2(magnetization);
    }

    //printf("|M| per site at temp T = %f: %f\n", curr_temp, fabs(magnetization)/(NUM_L*NUM_R*NUM_C));

    results[j][0] = curr_temp;
    results[j][1] = ((m_sq_avg/num_samples) - (gsl_pow_2(m_avg/num_samples)))/(curr_temp*NUM_L*NUM_R*NUM_C);

    printf("%f,%f\n", results[j][0], results[j][1]);

    curr_temp -= DELTA_T;
    j++;

    m_avg = 0;
    m_sq_avg = 0;
  }
  return j;
}

/* EXPERIMENTS */
int M_v_B(double** results){
    int cor_count = 0;
    int num_samples = 1000;
    int n = 0;
    int sample_counter = 0;
    int i;
    double init_B = B;
    printf("B = %f, init_B = %f, fabs = %f\n", B, init_B, fabs(init_B));

    cool_lattice(FINAL_T);
    while(B < fabs(init_B)){
        //printf("equilibrating\n");
        //fflush(stdout);

        simulate(EQ_TIME, FINAL_T);
        // Measure magnetization
        //printf("recording b\n");
        //fflush(stdout);

        results[sample_counter][0] = B;
        //printf("recording m\n");
        //fflush(stdout);


        results[sample_counter][1] = magnetization;

        //for(i=1; i <= NUM_L; i++)
        //  results[sample_counter][i+1] = calc_magnetization(i-1);
        //printf("recording tc\n");
        //fflush(stdout);


	      results[sample_counter][NUM_L+2] = calc_TC(0);

        //printf("cors\n");
        //fflush(stdout);

        for(cor_count = 1; cor_count < num_samples; cor_count++){
          simulate(COR_TIME, FINAL_T);
          results[sample_counter][1] += magnetization;

          //for(i=0; i <= NUM_L; i++)
          //  results[sample_counter][i+1] += calc_magnetization(i-1);
	        results[sample_counter][NUM_L+2] += calc_TC(0);
        }
        //printf("taking avs\n");
        //fflush(stdout);

        results[sample_counter][1] = results[sample_counter][1]/num_samples;
        //for(i=0; i <= NUM_L; i++)
        //  results[sample_counter][i+1] = results[sample_counter][i+1]/num_samples;
	      results[sample_counter][NUM_L + 2] = results[sample_counter][NUM_L+2]/num_samples;

        if(DEBUG)
          printf("%f,%f,%f\n", B, results[sample_counter][1], results[sample_counter][NUM_L + 2]);

        sample_counter += 1;

        B += DELTA_B;
    }
    while(B > init_B){
        simulate(EQ_TIME, FINAL_T);
        // Measure magnetization
        results[sample_counter][0] = B;
        results[sample_counter][1] = magnetization;
        //for(i=0; i <= NUM_L; i++)
        //  results[sample_counter][i+1] = calc_magnetization(i-1);
	      results[sample_counter][NUM_L+2] = calc_TC(0);

        // Take average over 1000 sweeps
        for(cor_count = 1; cor_count < num_samples; cor_count++){
          simulate(COR_TIME, FINAL_T);
          results[sample_counter][1] += magnetization;
          //for(i=0; i <= NUM_L; i++)
          //  results[sample_counter][i+1] += calc_magnetization(i-1);
	        results[sample_counter][NUM_L+2] += calc_TC(0);
        }

        results[sample_counter][1] = results[sample_counter][1]/num_samples;

        //for(i=0; i <= NUM_L; i++)
        //  results[sample_counter][i+1] = results[sample_counter][i+1]/num_samples;
	      results[sample_counter][NUM_L + 2] = results[sample_counter][NUM_L + 2]/num_samples;

        if(DEBUG)
          printf("%f,%f,%f\n", B, results[sample_counter][1], results[sample_counter][NUM_L + 2]);


        sample_counter += 1;
        B -= DELTA_B;
    }

    return sample_counter;
}

void cross_product(const gsl_vector *u, const gsl_vector *v, gsl_vector *product)
{
        double p1 = gsl_vector_get(u, 1)*gsl_vector_get(v, 2)
                - gsl_vector_get(u, 2)*gsl_vector_get(v, 1);

        double p2 = gsl_vector_get(u, 2)*gsl_vector_get(v, 0)
                - gsl_vector_get(u, 0)*gsl_vector_get(v, 2);

        double p3 = gsl_vector_get(u, 0)*gsl_vector_get(v, 1)
                - gsl_vector_get(u, 1)*gsl_vector_get(v, 0);

        gsl_vector_set(product, 0, p1);
        gsl_vector_set(product, 1, p2);
        gsl_vector_set(product, 2, p3);
}

void cleanup(){
  if(DEBUG)
    printf("Cleaning up...\n");
  gsl_vector_free(spin_vector);
  gsl_vector_free(point_vector);
  gsl_vector_free(temp_vector);
  gsl_matrix_free(rot_matrix);


  gsl_vector_free(delta_vector);
  gsl_vector_free(neighbor_vector);
  gsl_vector_free(inter_vector);
  gsl_vector_free(cross_temp);
  gsl_vector_free(neighbors[0]);
  gsl_vector_free(neighbors[1]);
  gsl_vector_free(neighbors[2]);
  gsl_vector_free(neighbors[3]);
}
