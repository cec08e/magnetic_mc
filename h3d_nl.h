#ifndef H3D_NL_HEADER
#define H3D_NL_HEADER

#define CONFIG_FILE "config.txt"
#define SIM_CONFIG "sim_results/sim_results.txt"   /* Configuration parameter recording file location */


/* Likely to be modified */
int SIM_NUM = 1;       /* Simulation Number */
short NUM_L = 1;           /* Number of layers */
short NUM_R = 20;           /* Number of rows in each lattice layer */
short NUM_C = 20;           /* Number of columns in each lattice layer */
float INIT_T = 3.0;         /* Initial temperature */
float FINAL_T = .0025;      /* Final Temperature */
float DELTA_T = .1;        /* Annealing temperature interval */
float DELTA_B = .01;       /* B sweeping interval */
short OVER_FLAG = 0;       /* Number of overrelaxation sweeps per MC sweep */
int ANNEAL_TIME = 2000;    /* Number of annealing sweeps to perform at each intermediate temperature */
int EQ_TIME = 100000;      /* Number of equilibration sweeps */
int COR_TIME = 10;         /* Number of sweeps between measurements of equilibrated system */


/* Uniform interaction values */
float B_CONST;
float * J_INTRA_CONST;
float * J_INTER_CONST;
float * D_INTRA_CONST;
float * D_INTER_CONST;
float * K_CONST;


/* Periodicity flags - determine whether the interaction parameters should
be overrode by a custom function. If not, a uniform value must be supplied in the
configuration file */
short J_INTRA_P = 0;
short J_INTER_P = 0;
short D_INTRA_P = 0;
short D_INTER_P = 0;
short K_P = 0;

/* Lattice interaction parameters
   All can be accessed via layer,row,col indices to get the parameter value
   for that index. For example, the interlayer coupling between (1,2,3) and (2,2,3)
   is accessible as J_INTER[1][2][3]
*/
float *** J_INTER;     /* Interlayer exchange interaction strength */
float *** J_INTRA;     /* Intralayer exchange interaction strength */
float *** D_INTER;     /* Interlayer DM interaction strength */
float *** D_INTRA;     /* Intralayer DM interaction strength */
float *** K;           /* Anisotropy strength */
float *** B;           /* External field strength */

/* Not likely to be modified */
float RADIUS = .6;         /* Radius of tangent disc in perturbing function */

typedef struct {
  double x;
  double y;
  double z;
} spin_t;             /* Spin type structure - x, y, and z components of spin */

typedef spin_t *** lattice_t;    /* Lattice type */
gsl_rng * rng;

void parse_config_file();
void echo_params();
void build_lattice();
void gen_random_spin(spin_t* spin);

#endif