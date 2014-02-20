#include <algorithm>    // std::max
#include "cuda.h"
#include "cuda_runtime.h"
#include "../Headers/Common.h"
#include "../Headers/hemi.h"
#include "../Headers/array.h"
#include "math.h"

#define DEBUG

__constant__ double c_tau;
__constant__ double c_tau_to_current_time_level;
__constant__ double c_lb;
__constant__ double c_rb;
__constant__ double c_ub;
__constant__ double c_bb;
__constant__ double c_b;
__constant__ int c_x_length;
__constant__ int c_length;

__device__ void d_quadrAngleType(Triangle& firstT,
								 Triangle& secondT, double* x, double* y, int i, int j) 
{
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

	//   1. First of all let's compute coordinates of square vertexes.
	alpha[0] = (x[i - 1] + x[i]) / 2.;
	betta[0] = (x[i + 1] + x[i]) / 2.;
	gamma[0] = (x[i + 1] + x[i]) / 2.;
	theta[0] = (x[i - 1] + x[i]) / 2.;

	alpha[1] = (y[j] + y[j - 1]) / 2.;
	betta[1] = (y[j] + y[j - 1]) / 2.;
	gamma[1] = (y[j] + y[j + 1]) / 2.;
	theta[1] = (y[j] + y[j + 1]) / 2.;

	//   2. Now let's compute new coordinates on the previous time level of alpha, betta, gamma, theta points.
	//  alNew.

	u = c_b * alpha[1] * (1. - alpha[1]) * (C_pi_device / 2. + atan(alpha[0]));

	v = atan(
		(alpha[0] - c_lb) * (alpha[0] - c_rb) * (1. + c_tau_to_current_time_level) / 10. * (alpha[1] - c_ub)
		* (alpha[1] - c_bb));
	alNew[0] = alpha[0] - c_tau * u;
	alNew[1] = alpha[1] - c_tau * v;

	//  beNew.

	u = c_b * betta[1] * (1. - betta[1]) * (C_pi_device / 2. + atan(betta[0]));
	v = atan(
		(betta[0] - c_lb) * (betta[0] - c_rb) * (1. + c_tau_to_current_time_level) / 10. * (betta[1] - c_ub)
		* (betta[1] - c_bb));
	beNew[0] = betta[0] - c_tau * u;
	beNew[1] = betta[1] - c_tau * v;

	//  gaNew.


	u = c_b * gamma[1] * (1. - gamma[1]) * (C_pi_device / 2. + atan(gamma[0]));
	v = atan(
		(gamma[0] - c_lb) * (gamma[0] - c_rb) * (1. + c_tau_to_current_time_level) / 10. * (gamma[1] - c_ub)
		* (gamma[1] - c_bb));
	gaNew[0] = gamma[0] - c_tau * u;
	gaNew[1] = gamma[1] - c_tau * v;

	//  thNew.

	u = c_b * theta[1] * (1. - theta[1]) * (C_pi_device / 2. + atan(theta[0]));
	v =  atan(
		(theta[0] - c_lb) * (theta[0] - c_rb) * (1. + c_tau_to_current_time_level) / 10. * (theta[1] - c_ub)
		* (theta[1] - c_bb));
	thNew[0] = theta[0] - c_tau * u;
	thNew[1] = theta[1] - c_tau * v;

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

	if (fabs(b_1LC * a_2LC - b_2LC * a_1LC) < 1.e-14) 
	{
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

		secondT.first[0] = beNew[0];
		secondT.first[1] = beNew[1];
		secondT.second[0] = thNew[0];
		secondT.second[1] = thNew[1];

		if (vectAlGa[0] * vectBeTh[0] + vectAlGa[1] * vectBeTh[1] >= 0.) 
		{
			secondT.third[0] = gaNew[0];
			secondT.third[1] = gaNew[1];
		}
		else if (vectAlGa[0] * vectBeTh[0] + vectAlGa[1] * vectBeTh[1] < 0.) 
		{
			secondT.third[0] = alNew[0];
			secondT.third[1] = alNew[1];
		}
		return;
	}

	AcrP[0] = (b_1LC * c_2LC - b_2LC * c_1LC) / (b_1LC * a_2LC - b_2LC * a_1LC);
	AcrP[1] = (a_1LC * c_2LC - a_2LC * c_1LC)
		/ (-b_1LC * a_2LC + b_2LC * a_1LC);

	if (((beNew[1] - AcrP[1]) * (thNew[1] - AcrP[1])) > 0.) 
	{

		if (((alNew[0] - AcrP[0]) * (gaNew[0] - AcrP[0])) > 0.)
		{

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

			if (vectAlGa[0] * vectBeTh[0] + vectAlGa[1] * vectBeTh[1] >= 0.) {
				secondT.third[0] = gaNew[0];
				secondT.third[1] = gaNew[1];
			} else if (vectAlGa[0] * vectBeTh[0] + vectAlGa[1] * vectBeTh[1] < 0.) {
				secondT.third[0] = alNew[0];
				secondT.third[1] = alNew[1];
			}

			return;

		} //   "if(  (  (alNew[0] - AcrP[0])*(gaNew[0] - AcsP[0])  )  >  0.  )".

		//   Second criterion. Second case.

		if (((alNew[0] - AcrP[0]) * (gaNew[0] - AcrP[0])) <= 0.)
		{
			vectAlBe[0] = beNew[0] - alNew[0];
			vectAlBe[1] = beNew[1] - alNew[1];
			vectAlTh[0] = thNew[0] - alNew[0];
			vectAlTh[1] = thNew[1] - alNew[1];

			if (vectAlBe[0] * vectAlTh[1] - vectAlBe[1] * vectAlTh[0] < 0.) 
			{
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
			else if (vectAlBe[0] * vectAlTh[1] - vectAlBe[1] * vectAlTh[0] >= 0.) 
			{
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

	if (((beNew[1] - AcrP[1]) * (thNew[1] - AcrP[1])) <= 0.) 
	{
		if (((alNew[0] - AcrP[0]) * (gaNew[0] - AcrP[0])) > 0.) 
		{
			//  It means the across point IS NOT between "alNew" and "gaNew" vertices by Ox-axis?

			//   O.K. the quadrangle IS NOT CONVEX. Is it concave or pseudo? Third criterion.

			vectBeGa[0] = gaNew[0] - beNew[0];
			vectBeGa[1] = gaNew[1] - beNew[1];
			vectBeAl[0] = alNew[0] - beNew[0];
			vectBeAl[1] = alNew[1] - beNew[1];

			if (vectBeGa[0] * vectBeAl[1] - vectBeGa[1] * vectBeAl[0] >= 0.)
			{

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
			else if (vectBeGa[0] * vectBeAl[1] - vectBeGa[1] * vectBeAl[0] < 0.) 
			{

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

			if (vectAlBe[0] * vectAlTh[1] - vectAlBe[1] * vectAlTh[0] >= 0.) {

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
			else if (vectAlBe[0] * vectAlTh[1] - vectAlBe[1] * vectAlTh[0] < 0.)
			{
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

__global__ void get_angle_type_kernel(Triangle* firstT, Triangle* secondT, 
									  double *x, 
									  double *y,
									  int offset) 
{
	for (int opt = hemiGetElementOffset(); opt < c_length; opt += hemiGetElementStride()) 
	{
		d_quadrAngleType(firstT[opt], secondT[opt], x, y, (opt % c_x_length + 1), 
			(opt / c_x_length + 1) + (int) (offset / c_x_length));
	}
}

void get_triangle_type(TriangleResult* result, ComputeParameters p, int gridSize, int blockSize) 
{
	int d_xy_size(0), length(0), copy_offset(0), tr_size(0);
	Triangle *first = NULL, *second = NULL;
	double *x = NULL, *y = NULL;
	double tau_to_current_time_level = 0.;
	d_xy_size = sizeof(double) * p.get_chunk_size();
	tr_size = sizeof(Triangle) * p.get_inner_chuck_size();
	length = p.get_inner_chuck_size();

	cudaMallocManaged(&x, d_xy_size);
	cudaMallocManaged(&y, d_xy_size);
	cudaMemcpyToSymbol(c_tau, &p.tau, sizeof(double));
	cudaMemcpyToSymbol(c_lb, &p.lb, sizeof(double));
	cudaMemcpyToSymbol(c_rb, &p.rb, sizeof(double));
	cudaMemcpyToSymbol(c_bb, &p.bb, sizeof(double));
	cudaMemcpyToSymbol(c_ub, &p.ub, sizeof(double));
	cudaMemcpyToSymbol(c_b, &p.b, sizeof(double));
	cudaMemcpyToSymbol(c_x_length, &result->x_length, sizeof(int));
	cudaMemcpyToSymbol(c_length, &length, sizeof(int));

	for (p.reset_time_counter(); p.can_iterate_over_time_level(); p.inc_time_level())
	{
		tau_to_current_time_level = p.currentTimeLevel * p.tau;
		cudaMemcpyToSymbol(c_tau_to_current_time_level, &tau_to_current_time_level, sizeof(double));

		for (int i = 0; i < p.get_part_count(); ++i) 
		{
			copy_offset = i * p.get_inner_chuck_size();

			cudaMallocManaged(&first, tr_size);
			cudaMallocManaged(&second, tr_size);

			memcpy(x, p.x, d_xy_size);
			memcpy(y, p.y, d_xy_size);

			get_angle_type_kernel<<<gridSize, blockSize>>>(first, second, x, y,
				p.get_inner_chuck_size() * i);
			cudaDeviceSynchronize();

			memcpy(&result->f[copy_offset], first, tr_size);
			memcpy(&result->s[copy_offset], second, tr_size);

			cudaFree(first);
			cudaFree(second);
		}
	}
	cudaFree(x);
	cudaFree(y);
	cudaDeviceReset();
}

