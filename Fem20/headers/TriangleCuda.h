#include "Common.h"

extern void get_triangle_type(TriangleResult* result, ComputeParameters p, int gridSize, int blockSize);
extern double d_u_function(double par_b, double t, double x, double y);
extern double d_v_function(double lbDom, double rbDom,double bbDom, double ubDom, double t, double x, double y );