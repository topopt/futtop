#include <getopt.h>
#include <math.h>
#include <stdalign.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "io.h"
#include "libmultigrid.h"

#ifdef _OPENMP
#include <omp.h>
#endif

double getTime() {
#ifdef _OPENMP
  return omp_get_wtime();
#else
  return (double)clock() / (double)CLOCKS_PER_SEC;
#endif
}

void handle_fut_err(struct futhark_context *fut_ctx, int err) {
  if (err != 0) {
    char *err_string = futhark_context_get_error(fut_ctx);
    printf("error in fut call!! error message:\n%s", err_string);
    free(err_string);
  }
}

void top3dmgcg(const uint_fast32_t nelx, const uint_fast32_t nely,
               const uint_fast32_t nelz, const int iters, const float volfrac,
               const float rmin) {

  int err;
  const uint_fast64_t nelem = nelx * nely * nelz;

  const uint_fast32_t nx = nelx + 1;
  const uint_fast32_t ny = nely + 1;
  const uint_fast32_t nz = nelz + 1;

  const uint_fast64_t ndof = 3 * nx * ny * nz;

  float *F = malloc(sizeof(float) * ndof);
  float *U = malloc(sizeof(float) * ndof);

  for (uint_fast64_t i = 0; i < ndof; i++) {
    F[i] = 0.0;
    U[i] = 0.0;
  }

  float forceMagnitude = -1;

  if (true) { // cantilever
    const int k = 0;
    const double radius = ((double)ny) / 5.0; // snap
    const double radius2 = radius * radius;
    const double center_x = (double)nelx;
    const double center_y = ((double)nely - 1.0) / 2.0;

    int num_elements = 0;
    for (int i = 0; i < nelx; i++) {
      for (int j = 0; j < nely; j++) {
        const double dx = (double)i - center_x;
        const double dy = (double)j - center_y;
        const double dist2 = dx * dx + dy * dy;
        if (dist2 < radius2) {
          num_elements++;
        }
      }
    }

    double nodalForce = forceMagnitude / (4.0 * (double)num_elements);
    int isum = 0;
    for (int i = 0; i < nelx; i++) {
      for (int j = 0; j < nely; j++) {

        const double dx = (double)i - center_x;
        const double dy = (double)j - center_y;
        const double dist2 = dx * dx + dy * dy;

        if (dist2 < radius2) {
          const int nidx1 = (i + 1) * ny * nz + (j + 1) * nz + k;
          const int nidx2 = (i + 1) * ny * nz + j * nz + k;
          const int nidx3 = (i)*ny * nz + (j + 1) * nz + k;
          const int nidx4 = (i)*ny * nz + j * nz + k;
          F[3 * nidx1 + 2] += nodalForce;
          F[3 * nidx2 + 2] += nodalForce;
          F[3 * nidx3 + 2] += nodalForce;
          F[3 * nidx4 + 2] += nodalForce;
        }
      }
    }
  } else if (true) { // MBB Beam and snapthrough with center load

    const int k = nelz;
    const double radius = ((double)ny) / 5.0;
    const double radius2 = radius * radius;
    const double center_x = 0;
    const double center_y = ((double)nely - 1.0) / 2.0;

    float nodalForce = forceMagnitude;
    int isum = 0;
    for (int i = 0; i < nelx; i++) {
      for (int j = 0; j < nely; j++) {

        const double dx = (double)i - center_x;
        const double dy = (double)j - center_y;
        const double dist2 = dx * dx + dy * dy;

        if (dist2 < radius2) {
          const int nidx1 = (i + 1) * ny * nz + (j + 1) * nz + k;
          const int nidx2 = (i + 1) * ny * nz + j * nz + k;
          const int nidx3 = (i)*ny * nz + (j + 1) * nz + k;
          const int nidx4 = (i)*ny * nz + j * nz + k;
          F[3 * nidx1 + 2] += nodalForce;
          F[3 * nidx2 + 2] += nodalForce;
          F[3 * nidx3 + 2] += nodalForce;
          F[3 * nidx4 + 2] += nodalForce;
        }
      }
    }
  }

  float *x = malloc(sizeof(float) * nelem);
  for (uint_fast64_t i = 0; i < nelem; i++) {
    x[i] = volfrac;
  }

  unsigned int loop = 0;
  float change = 1;

#ifdef _OPENMP
  printf(" OpenMP enabled with %d threads\n", omp_get_max_threads());
#endif
  const double totalTime = getTime();

  struct futhark_context_config *fut_cfg = futhark_context_config_new();
  futhark_context_config_set_cache_file(fut_cfg, "./futtop_compile_cache");
  struct futhark_context *fut_ctx = futhark_context_new(fut_cfg);

  // setup forces and displacements on device
  struct futhark_f32_4d *force_gpu =
      futhark_new_f32_4d(fut_ctx, F, nx, ny, nz, 3);
  struct futhark_f32_4d *uold_gpu =
      futhark_new_f32_4d(fut_ctx, U, nx, ny, nz, 3);
  struct futhark_f32_4d *u_gpu;

  // setup densities on device
  struct futhark_f32_3d *x_gpu =
      futhark_new_f32_3d(fut_ctx, x, nelx, nely, nelz);
  struct futhark_f32_3d *xPhys_gpu;
  struct futhark_f32_3d *xnew_gpu;

  err = futhark_entry_densityFilter(fut_ctx, &xPhys_gpu, rmin, x_gpu);
  handle_fut_err(fut_ctx, err);
  err = futhark_context_sync(fut_ctx);
  handle_fut_err(fut_ctx, err);

  struct futhark_opaque_multigridData *mgData = NULL;
  int filenumber = 0;
  int file_freq = 1;

  /* %% START ITERATION */
  while ((change > 1e-2) && (loop < iters)) {

    const double loopTime = getTime();
    double interTime;

    loop++;
    int32_t cgiter;
    float cgres, compliance, vol;

    // construct mg data for design iteration
    interTime = getTime();
    err = futhark_entry_assembleMultigridData(fut_ctx, &mgData, xPhys_gpu);
    handle_fut_err(fut_ctx, err);
    err = futhark_context_sync(fut_ctx);
    handle_fut_err(fut_ctx, err);
    const double assemblyTime = getTime() - interTime;

    // solve state system Ku = f
    interTime = getTime();
    err = futhark_entry_solveSystem(fut_ctx, &u_gpu, &cgres, &cgiter, xPhys_gpu,
                                    mgData, uold_gpu, force_gpu);
    handle_fut_err(fut_ctx, err);
    err = futhark_context_sync(fut_ctx);
    handle_fut_err(fut_ctx, err);
    const double solveTime = getTime() - interTime;

    if (loop % file_freq == 0) {
      err = futhark_values_f32_4d(fut_ctx, u_gpu, U);
      handle_fut_err(fut_ctx, err);
      err = futhark_values_f32_3d(fut_ctx, xPhys_gpu, x);
      handle_fut_err(fut_ctx, err);
      err = futhark_context_sync(fut_ctx);
      handle_fut_err(fut_ctx, err);

      char filename[50];
      sprintf(filename, "designstep_%d.vtr", filenumber);
      filenumber++;
      writeDensityAndDisplacement(nelx, nely, nelz, x, U, filename);
    }

    // clear uneeded data
    err = futhark_free_opaque_multigridData(fut_ctx, mgData);
    handle_fut_err(fut_ctx, err);
    err = futhark_free_f32_4d(fut_ctx, uold_gpu);
    handle_fut_err(fut_ctx, err);
    err = futhark_context_sync(fut_ctx);
    handle_fut_err(fut_ctx, err);
    uold_gpu = u_gpu;

    // compute new x
    interTime = getTime();
    err = futhark_entry_designIteration(fut_ctx, &compliance, &change, &vol,
                                        &xnew_gpu, x_gpu, xPhys_gpu, u_gpu,
                                        volfrac, rmin);
    handle_fut_err(fut_ctx, err);
    err = futhark_context_sync(fut_ctx);
    handle_fut_err(fut_ctx, err);
    const double updateTime = getTime() - interTime;

    // clear uneeded data
    err = futhark_free_f32_3d(fut_ctx, xPhys_gpu);
    handle_fut_err(fut_ctx, err);
    err = futhark_free_f32_3d(fut_ctx, x_gpu);
    handle_fut_err(fut_ctx, err);
    err = futhark_context_sync(fut_ctx);
    handle_fut_err(fut_ctx, err);

    x_gpu = xnew_gpu;
    err = futhark_entry_densityFilter(fut_ctx, &xPhys_gpu, rmin, x_gpu);
    handle_fut_err(fut_ctx, err);
    err = futhark_context_sync(fut_ctx);
    handle_fut_err(fut_ctx, err);

    printf("It.:%4i Obj.:%6.3e Vol.:%6.3f ch.:%4.2f relres: %4.2e iters: %4i\n",
           loop, compliance, vol, change, cgres, cgiter);

    printf("\tassembly: %4.2fs solve: %4.2fs update: %4.2fs total: %6.3fs\n",
           assemblyTime, solveTime, updateTime, getTime() - loopTime);
  }

  printf("End time: %6.3f \n", getTime() - totalTime);

  err = futhark_values_f32_3d(fut_ctx, xPhys_gpu, x);
  handle_fut_err(fut_ctx, err);
  err = futhark_context_sync(fut_ctx);
  handle_fut_err(fut_ctx, err);

  err = futhark_free_f32_4d(fut_ctx, force_gpu);
  handle_fut_err(fut_ctx, err);
  err = futhark_free_f32_3d(fut_ctx, xPhys_gpu);
  handle_fut_err(fut_ctx, err);
  err = futhark_free_f32_4d(fut_ctx, uold_gpu);
  handle_fut_err(fut_ctx, err);
  err = futhark_free_f32_3d(fut_ctx, x_gpu);
  handle_fut_err(fut_ctx, err);
  err = futhark_context_sync(fut_ctx);
  handle_fut_err(fut_ctx, err);

  futhark_context_free(fut_ctx);
  futhark_context_config_free(fut_cfg);

  char filename[50];
  sprintf(filename, "designstep_%d.vtr", filenumber);
  writeDensityAndDisplacement(nelx, nely, nelz, x, F, filename);

  // writeDensity(nelx, nely, nelz, x, "output.vtr");
  free(x);
  free(F);
  free(U);
}

int main(int argc, char *argv[]) {
  int nelx_coarse = 4;
  int nely_coarse = 2;
  int nelz_coarse = 2;
  float volfrac = 0.2;
  float rmin = 1.5;
  int iters = 20;
  int opt;

  while ((opt = getopt(argc, argv, "x:y:z:r:v:i:")) != -1) {
    switch (opt) {
    case 'x':
      nelx_coarse = atoi(optarg);
      break;
    case 'y':
      nely_coarse = atoi(optarg);
      break;
    case 'z':
      nelz_coarse = atoi(optarg);
      break;
    case 'r':
      rmin = atof(optarg);
      break;
    case 'v':
      volfrac = atof(optarg);
      break;
    case 'i':
      iters = atoi(optarg);
      break;
    default:
      fprintf(stderr, "Usage: %s [-xyzrv]\n", argv[0]);
      exit(EXIT_FAILURE);
    }
  }

  int nelx = nelx_coarse * 8;
  int nely = nely_coarse * 8;
  int nelz = nelz_coarse * 8;

  printf("Running topopt with:\n number of coarse els x (-x): %i (fine els = "
         "%i)\n number of coarse els y (-y): %i (fine els = %i)\n number of "
         "coarse els z (-z): %i (fine els = %i)\n total number of elements: "
         "%i\n volume fraction (-v): %f\n filter radius in elements (-r): %f\n "
         "number of design iterations (-i): %i\n\n",
         nelx_coarse, nelx, nely_coarse, nely, nelz_coarse, nelz,
         nelx * nely * nelz, volfrac, rmin, iters);
  top3dmgcg(nelx, nely, nelz, iters, volfrac, rmin);
}
