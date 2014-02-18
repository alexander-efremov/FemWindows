#include <algorithm>    // std::max
#include "cuda.h"
#include "cuda_runtime.h"
#include "../Headers/Common.h"
#include "../Headers/hemi.h"
#include "../Headers/array.h"
#include "math.h"

#define DEBUG

//#define PRINT_VALUES

__device__ double d_u_function_quadrangle(double par_b, double t, double x,
										  double y) {
											  return par_b * y * (1. - y) * (C_pi_device / 2. + atan(-x));
}

__device__ double d_v_function_quadrangle(double lbDom, double rbDom,
										  double bbDom, double ubDom, double t, double x, double y) {
											  return atan(
												  (x - lbDom) * (x - rbDom) * (1. + t) / 10. * (y - ubDom)
												  * (y - bbDom));
}

__device__ void d_quadrAngleType(ComputeParameters& p, Triangle& firstT,
								 Triangle& secondT) {


									 double alpha[2], betta[2], gamma[2], theta[2]; //   -  Vertexes of square. Anticlockwise order from left bottom vertex.
									 double u, v;                           //   -  Items of velocity components.
									 double alNew[2], beNew[2], gaNew[2], thNew[2]; //   -  New positions of vertexes. Vertexes of quadrangle.
									 double vectAlGa[2], vectBeTh[2]; //   -  Vectors: 1) from "alNew" to "gaNew" 2) from "beNew" to "thNew".
									 double a_1LC, b_1LC, c_1LC; //   -  a_1LC * x  +  b_1LC * y  = c_1LC. Line equation through "alNew" and "gaNew".
									 double a_2LC, b_2LC, c_2LC; //   -  a_2LC * x  +  b_2LC * y  = c_2LC. Line equation through "beNew" and "thNew".
									 double AcrP[2];                            //   -  Across point of two lines
									 double vectAlBe[2]; //   -  Vectors for computing vertexes sequence order by vector production.
									 double vectAlTh[2]; //   -  Vectors for computing vertexes sequence order by vector production.
									 double vectBeGa[2]; //   -  Vectors for computing vertexes sequence order by vector production.
									 double vectBeAl[2]; //   -  Vectors for computing vertexes sequence order by vector production.
									 double vectProdOz;                      //   -  Z-axis of vector production.
									 double scalProd;                   //   -  Scalar production of two vectors.

									 //   1. First of all let's compute coordinates of square vertexes.
									 alpha[0] = (p.x[p.i - 1] + p.x[p.i]) / 2.;
									 betta[0] = (p.x[p.i + 1] + p.x[p.i]) / 2.;
									 gamma[0] = (p.x[p.i + 1] + p.x[p.i]) / 2.;
									 theta[0] = (p.x[p.i - 1] + p.x[p.i]) / 2.;

									 alpha[1] = (p.y[p.j] + p.y[p.j - 1]) / 2.;
									 betta[1] = (p.y[p.j] + p.y[p.j - 1]) / 2.;
									 gamma[1] = (p.y[p.j] + p.y[p.j + 1]) / 2.;
									 theta[1] = (p.y[p.j] + p.y[p.j + 1]) / 2.;

									 // TODO remove
#ifdef PRINT_VALUES1
									 //if (p.j ==9 && p.i == 9)
									 printf("cuda i = %d \n %f %f %f %f\ncuda j = %d \n %f %f %f %f\n", p.i,
										 alpha[0], betta[0], gamma[0], theta[0], p.j, alpha[1], betta[1],
										 gamma[1], theta[1]);
#endif

									 //   2. Now let's compute new coordinates on the previous time level of alpha, betta, gamma, theta points.
									 //  alNew.

									 u = d_u_function_quadrangle(p.b, p.currentTimeLevel * p.tau, alpha[0], alpha[1]);
									 v = d_v_function_quadrangle(p.lb, p.rb, p.bb, p.ub, p.currentTimeLevel * p.tau,
										 alpha[0], alpha[1]);
									 alNew[0] = alpha[0] - p.tau * u;
									 alNew[1] = alpha[1] - p.tau * v;

									 //  beNew.
									 u = d_u_function_quadrangle(p.b, p.currentTimeLevel * p.tau, betta[0], betta[1]);
									 v = d_v_function_quadrangle(p.lb, p.rb, p.bb, p.ub, p.currentTimeLevel * p.tau,
										 betta[0], betta[1]);
									 beNew[0] = betta[0] - p.tau * u;
									 beNew[1] = betta[1] - p.tau * v;

									 //  gaNew.

									 u = d_u_function_quadrangle(p.b, p.currentTimeLevel * p.tau, gamma[0], gamma[1]);
									 v = d_v_function_quadrangle(p.lb, p.rb, p.bb, p.ub, p.currentTimeLevel * p.tau,
										 gamma[0], gamma[1]);
									 gaNew[0] = gamma[0] - p.tau * u;
									 gaNew[1] = gamma[1] - p.tau * v;

									 //  thNew.
									 u = d_u_function_quadrangle(p.b, p.currentTimeLevel * p.tau, theta[0], theta[1]);
									 v = d_v_function_quadrangle(p.lb, p.rb, p.bb, p.ub, p.currentTimeLevel * p.tau,
										 theta[0], theta[1]);
									 thNew[0] = theta[0] - p.tau * u;
									 thNew[1] = theta[1] - p.tau * v;

									 //   3.a Let's compute coefficients of first line betweeen "alNew" and "gaNew" points.
									 //   a_1LC * x  +  b_1LC * y  = c_1LC.

									 vectAlGa[0] = gaNew[0] - alNew[0];
									 vectAlGa[1] = gaNew[1] - alNew[1];
									 a_1LC = vectAlGa[1];
									 b_1LC = -vectAlGa[0];
									 c_1LC = vectAlGa[1] * alNew[0] - vectAlGa[0] * alNew[1];

									 //   3.b Let's compute coefficients of second line betweeen "beNew" and "thNew" points.
									 //   a_2LC * x  +  b_2LC * y  = c_2LC.

									 vectBeTh[0] = thNew[0] - beNew[0];
									 vectBeTh[1] = thNew[1] - beNew[1];
									 a_2LC = vectBeTh[1];
									 b_2LC = -vectBeTh[0];
									 c_2LC = vectBeTh[1] * beNew[0] - vectBeTh[0] * beNew[1];

									 //   4. Let's compute coordinates of across point of this two lines.
									 //   Are lines parallel?

									 if (fabs(b_1LC * a_2LC - b_2LC * a_1LC) < 1.e-14) {
										 //   Not checked.



										 //   Pseudo case. Anyway I need to compute some values.

										 //   First triangle.

										 firstT.first[0] = alNew[0];
										 firstT.first[1] = alNew[1];
										 firstT.second[0] = beNew[0];
										 firstT.second[1] = beNew[1];
										 firstT.third[0] = gaNew[0];
										 firstT.third[1] = gaNew[1];

										 //   Vertices of second triagle depends on scalar production.

										 vectAlGa[0] = gaNew[0] - alNew[0];
										 vectAlGa[1] = gaNew[1] - alNew[1];
										 vectBeTh[0] = thNew[0] - beNew[0];
										 vectBeTh[1] = thNew[1] - beNew[1];

										 scalProd = vectAlGa[0] * vectBeTh[0] + vectAlGa[1] * vectBeTh[1];
										 secondT.first[0] = beNew[0];
										 secondT.first[1] = beNew[1];
										 secondT.second[0] = thNew[0];
										 secondT.second[1] = thNew[1];

										 if (scalProd >= 0.) {
											 secondT.third[0] = gaNew[0];
											 secondT.third[1] = gaNew[1];
										 }

										 if (scalProd < 0.) {
											 secondT.third[0] = alNew[0];
											 secondT.third[1] = alNew[1];
										 }

										 return;
									 }

									 AcrP[0] = (b_1LC * c_2LC - b_2LC * c_1LC) / (b_1LC * a_2LC - b_2LC * a_1LC);
									 AcrP[1] = (a_1LC * c_2LC - a_2LC * c_1LC)
										 / (-b_1LC * a_2LC + b_2LC * a_1LC);

									 if (((beNew[1] - AcrP[1]) * (thNew[1] - AcrP[1])) > 0.) {

										 if (((alNew[0] - AcrP[0]) * (gaNew[0] - AcrP[0])) > 0.) {

											 firstT.first[0] = alNew[0];
											 firstT.first[1] = alNew[1];
											 firstT.second[0] = beNew[0];
											 firstT.second[1] = beNew[1];
											 firstT.third[0] = gaNew[0];
											 firstT.third[1] = gaNew[1];

											 //   Second triangle.

											 secondT.first[0] = beNew[0];
											 secondT.first[1] = beNew[1];
											 secondT.second[0] = thNew[0];
											 secondT.second[1] = thNew[1];

											 //   Third vertex computing...

											 vectAlGa[0] = gaNew[0] - alNew[0];
											 vectAlGa[1] = gaNew[1] - alNew[1];

											 vectBeTh[0] = thNew[0] - beNew[0];
											 vectBeTh[1] = thNew[1] - beNew[1];

											 scalProd = vectAlGa[0] * vectBeTh[0] + vectAlGa[1] * vectBeTh[1];

											 if (scalProd >= 0.) {
												 secondT.third[0] = gaNew[0];
												 secondT.third[1] = gaNew[1];
											 }

											 if (scalProd < 0.) {
												 secondT.third[0] = alNew[0];
												 secondT.third[1] = alNew[1];
											 }

											 return;

										 } //   "if(  (  (alNew[0] - AcrP[0])*(gaNew[0] - AcsP[0])  )  >  0.  )".

										 //   Second criterion. Second case.

										 if (((alNew[0] - AcrP[0]) * (gaNew[0] - AcrP[0])) <= 0.) {
											 vectAlBe[0] = beNew[0] - alNew[0];
											 vectAlBe[1] = beNew[1] - alNew[1];
											 vectAlTh[0] = thNew[0] - alNew[0];
											 vectAlTh[1] = thNew[1] - alNew[1];

											 vectProdOz = vectAlBe[0] * vectAlTh[1] - vectAlBe[1] * vectAlTh[0];

											 if (vectProdOz < 0.) {
												 //   The vertex "beNew" is NOT in triangle "alNew - gaNew - thNew".


												 //   Pseudo case. Anyway I need to find some solution. So

												 firstT.first[0] = alNew[0];
												 firstT.first[1] = alNew[1];
												 firstT.second[0] = beNew[0];
												 firstT.second[1] = beNew[1];
												 firstT.third[0] = thNew[0];
												 firstT.third[1] = thNew[1];

												 //   Second triangle.

												 secondT.first[0] = beNew[0];
												 secondT.first[1] = beNew[1];
												 secondT.second[0] = thNew[0];
												 secondT.second[1] = thNew[1];
												 secondT.third[0] = gaNew[0];
												 secondT.third[1] = gaNew[1];

												 return;
											 }

											 if (vectProdOz >= 0.) {
												 //  It's all write. We have a good concave quadrangle.

												 //   Now let's compute all vertices which I need.

												 //   First triangle.

												 firstT.first[0] = alNew[0];
												 firstT.first[1] = alNew[1];
												 firstT.second[0] = beNew[0];
												 firstT.second[1] = beNew[1];
												 firstT.third[0] = thNew[0];
												 firstT.third[1] = thNew[1];

												 //   Second triangle.

												 secondT.first[0] = beNew[0];
												 secondT.first[1] = beNew[1];
												 secondT.second[0] = thNew[0];
												 secondT.second[1] = thNew[1];
												 secondT.third[0] = gaNew[0];
												 secondT.third[1] = gaNew[1];

												 return;
											 }

										 } //   "if(  (  (alNew[0] - AcsP[0])*(gaNew[0] - AcsP[0])  )  <=  0.  )".   //   Last second case of second criterion.

									 } //   end of "if (  (  (beNew[1] - AcrP[1])*(thNew[1] - AcrP[1])  )  >  0.  )"

									 //  Now let's consider SECOND case 5.b "(  (beNew[1] - AcrP[1])*(thNew[1] - AcrP[1])  )  <= 0."

									 if (((beNew[1] - AcrP[1]) * (thNew[1] - AcrP[1])) <= 0.) {
										 if (((alNew[0] - AcrP[0]) * (gaNew[0] - AcrP[0])) > 0.) {
											 //  It means the across point IS NOT between "alNew" and "gaNew" vertices by Ox-axis?

											 //   O.K. the quadrangle IS NOT CONVEX. Is it concave or pseudo? Third criterion.

											 vectBeGa[0] = gaNew[0] - beNew[0];
											 vectBeGa[1] = gaNew[1] - beNew[1];
											 vectBeAl[0] = alNew[0] - beNew[0];
											 vectBeAl[1] = alNew[1] - beNew[1];

											 vectProdOz = vectBeGa[0] * vectBeAl[1] - vectBeGa[1] * vectBeAl[0];

											 if (vectProdOz >= 0.) {

												 //   The quadrangle is concave. First triangle.

												 firstT.first[0] = alNew[0];
												 firstT.first[1] = alNew[1];
												 firstT.second[0] = beNew[0];
												 firstT.second[1] = beNew[1];
												 firstT.third[0] = gaNew[0];
												 firstT.third[1] = gaNew[1];

												 //   Second triangle.

												 secondT.first[0] = alNew[0];
												 secondT.first[1] = alNew[1];
												 secondT.second[0] = thNew[0];
												 secondT.second[1] = thNew[1];
												 secondT.third[0] = gaNew[0];
												 secondT.third[1] = gaNew[1];

												 return;
											 }

											 if (vectProdOz < 0.) {

												 //   This concave quadrangle do has NO write anticlockwise vertices sequence order. It's pseudo.

												 firstT.first[0] = alNew[0];
												 firstT.first[1] = alNew[1];
												 firstT.second[0] = beNew[0];
												 firstT.second[1] = beNew[1];
												 firstT.third[0] = gaNew[0];
												 firstT.third[1] = gaNew[1];

												 //   Second triangle.

												 secondT.first[0] = alNew[0];
												 secondT.first[1] = alNew[1];
												 secondT.second[0] = thNew[0];
												 secondT.second[1] = thNew[1];
												 secondT.third[0] = gaNew[0];
												 secondT.third[1] = gaNew[1];

												 return;
											 }
										 } //   end of "if(  (  (alNew[0] - AcrP[0])*(gaNew[0] - AcsP[0])  )  >  0.  )". First case of second criterion.

										 //   Second criterion. Second case.

										 if (((alNew[0] - AcrP[0]) * (gaNew[0] - AcrP[0])) <= 0.) {
											 //   O.K. the quadrangle is convex. Is it has the same vertices sequence order.

											 vectAlBe[0] = beNew[0] - alNew[0];

											 vectAlBe[1] = beNew[1] - alNew[1];

											 vectAlTh[0] = thNew[0] - alNew[0];

											 vectAlTh[1] = thNew[1] - alNew[1];

											 vectProdOz = vectAlBe[0] * vectAlTh[1] - vectAlBe[1] * vectAlTh[0];

											 if (vectProdOz >= 0.) {

												 //   Convex quadrangle DO HAS WRITE anticlockwise vertices sequence order. It's convex.

												 firstT.first[0] = alNew[0];
												 firstT.first[1] = alNew[1];
												 firstT.second[0] = beNew[0];
												 firstT.second[1] = beNew[1];
												 firstT.third[0] = gaNew[0];
												 firstT.third[1] = gaNew[1];

												 //   Second triangle.

												 secondT.first[0] = alNew[0];
												 secondT.first[1] = alNew[1];
												 secondT.second[0] = thNew[0];
												 secondT.second[1] = thNew[1];
												 secondT.third[0] = gaNew[0];
												 secondT.third[1] = gaNew[1];

												 return;
											 }

											 if (vectProdOz < 0.) {

												 firstT.first[0] = alNew[0];
												 firstT.first[1] = alNew[1];
												 firstT.second[0] = beNew[0];
												 firstT.second[1] = beNew[1];
												 firstT.third[0] = gaNew[0];
												 firstT.third[1] = gaNew[1];
												 secondT.first[0] = alNew[0];
												 secondT.first[1] = alNew[1];
												 secondT.second[0] = thNew[0];
												 secondT.second[1] = thNew[1];
												 secondT.third[0] = gaNew[0];
												 secondT.third[1] = gaNew[1];
												 return;
											 }
										 }
									 }

}

__global__ void get_angle_type_kernel(TriangleResult result, ComputeParameters p) {
	const int offset = hemiGetElementOffset();
	const int stride = hemiGetElementStride();

	for (int opt = offset; opt < result.length; opt += stride) {
		p.i = (opt % result.x_length + 1) ;
		p.j = (opt / result.x_length + 1) + (int) (result.offset / result.x_length);
#ifdef PRINT_VALUES1
		printf("p.i = %d p.j = %d \n", p.i, p.j);
#endif
		d_quadrAngleType(p, result.f[opt], result.s[opt]);
	}
}

TriangleResult get_triangle_type(ComputeParameters p, int gridSize, int blockSize) 
{
	int d_xy_size(0), offset(0), length(0), copy_offset(0), cp_start_index(0), tr_size(0);

	Triangle* first = new Triangle[p.get_inner_matrix_size()];
	Triangle* second = new Triangle[p.get_inner_matrix_size()];
	TriangleResult result = TriangleResult(p);
	double* x = p.x;
	double* y = p.y;
	
	for (int i = 0; i < p.get_part_count(); ++i) 
	{
		offset = i * p.get_chunk_size();
		length = std::min(p.get_inner_chuck_size(), p.size  - offset);
		copy_offset = offset - 2*i;
		tr_size = sizeof(Triangle) * length;
		d_xy_size = sizeof(double) * p.get_chunk_size();
		cp_start_index = (i*p.get_inner_chuck_size()) % p.get_inner_x_size();

		checkCuda(cudaMallocManaged(&result.f, tr_size));
		checkCuda(cudaMallocManaged(&result.s, tr_size));
		checkCuda(cudaMallocManaged(&p.x, d_xy_size));
		checkCuda(cudaMallocManaged(&p.y, d_xy_size));
		memcpy(p.x, &x[cp_start_index], d_xy_size);
		memcpy(p.y, &y[cp_start_index], d_xy_size);

		result.length = length;
		result.setOffset(i);

		get_angle_type_kernel<<<gridSize, blockSize>>>(result, p);
		cudaDeviceSynchronize();


		
		memcpy(&first [copy_offset],  result.f, tr_size);
		memcpy(&second[copy_offset], result.s, tr_size);

		cudaFree(p.x);
		cudaFree(p.y);
		cudaFree(result.f);
		cudaFree(result.s);
	}

	cudaDeviceReset();
	result.f = first;
	result.s = second;
	p.x = x;
	p.y = y;
	return result;
}

