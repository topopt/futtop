#ifndef IO
#define IO

void writeDensity(const int nelx, const int nely, const int nelz,
                  const float *densityArray, const char *filename);

void writeDensityAndDisplacement(const int nelx, const int nely, const int nelz,
                                 const float *densityArray,
                                 const float *displacementArray,
                                 const char *filename);

void writeDensityAndDisplacementDouble(const int nelx, const int nely,
                                       const int nelz,
                                       const float *densityArray,
                                       const double *displacementArray,
                                       const char *filename);

#endif
