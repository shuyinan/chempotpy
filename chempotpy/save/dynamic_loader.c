#include <stdio.h>
#include <dlfcn.h>

#define NATOMS 3
#define NSTATES 2
#define DIM 3

// Assume pes takes an int and a double
typedef void (*pes_t)(double*, int*, double*, double*, double*);
typedef void (*dpem_t)(double*, int*, double*, double*);

//typedef void (*pes_t)(double[NATOMS][DIM], int, double[NSTATES], double[NSTATES][NATOMS][DIM], double[NSTATES][NSTATES][NATOMS][DIM]);
//typedef void (*dpem_t)(double[NATOMS][DIM], int, double[NSTATES][NSTATES], double[NSTATES][NSTATES][NATOMS][DIM]);

void call_pes(const char* so_name, double geom[NATOMS][DIM], int igrad, double p[NSTATES], double g[NSTATES][NATOMS][DIM], double d[NSTATES][NSTATES][NATOMS][DIM]) {
  void *handle;
  char *error;
  pes_t pes_func;

  handle = dlopen(so_name, RTLD_LAZY);
  if (!handle) {
    fprintf(stderr, "Cannot open library: %s\n", dlerror());
    return;
  }

  pes_func = (pes_t) dlsym(handle, "pes_");
  if (dlerror() != NULL) {
    fprintf(stderr, "Couldn't find pes function.\n");
    return;
  }

  printf("Input geometry:\n");
  for (int i = 0; i < NATOMS; ++i) {
    for (int j = 0; j < DIM; ++j) {
      printf("%f ", geom[i][j]);  // Note how we access the array; this is due to Fortran's column-major order
    }
    printf("\n");
  }

  printf("REDEFINE Input geometry:\n");
  geom[0][0]=-0.506;
  geom[0][1]=0.0;
  geom[0][2]=-0.094;
  geom[1][0]=-1.789;
  geom[1][1]=0.0;
  geom[1][2]=-0.902;
  geom[2][0]=-1.467;
  geom[2][1]=0.0;
  geom[2][2]=0.0;

  printf("Input geometry:\n");
  for (int i = 0; i < NATOMS; ++i) {
    for (int j = 0; j < DIM; ++j) {
      printf("%f ", geom[i][j]);  // Note how we access the array; this is due to Fortran's column-major order
    }
    printf("\n");
  }

  pes_func(&geom[0][0], &igrad, &p[0], &g[0][0][0], &d[0][0][0][0]);
  //pes_func(geom, igrad, p, g, d);
 


  printf("output surface:\n");
  for (int i = 0; i < NSTATES; ++i) {
    printf("%f ", p[i]);
  }

  dlclose(handle);
}

void call_dpem(const char* so_name, double geom[NATOMS][DIM], int igrad, double u[NSTATES][NSTATES], double ug[NSTATES][NSTATES][NATOMS][DIM]) {
  void *handle;
  dpem_t dpem_func;

  handle = dlopen(so_name, RTLD_LAZY);
  if (!handle) {
    fprintf(stderr, "Cannot open library: %s\n", dlerror());
    return;
  }

  dpem_func = (dpem_t) dlsym(handle, "dpem_");
  if (dlerror() != NULL) {
    fprintf(stderr, "Couldn't find dpem function.\n");
    return;
  }

  dpem_func(&geom[0][0], &igrad, &u[0][0], &ug[0][0][0][0]);
  //dpem_func(geom, igrad, u, ug);

  dlclose(handle);
}

