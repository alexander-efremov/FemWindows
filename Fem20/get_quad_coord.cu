#include "cuda.h"
#include "cuda_runtime.h"
#include "headers/hemi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "headers/Common.h"

// assert() is only supported // for devices of compute capability 2.0 and higher 
#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 200) 
#undef  assert 
#define assert(arg) 
#endif

#define DEBUG

struct SoATriangle
{
private:
	int size;

public:
	// first triangle
	double* first1;
	double* second1;
	double* third1;

	// second triangle
	double* first2;
	double* second2;
	double* third2;

	SoATriangle(int n)
	{
		size = 2*sizeof(double)*n;
		checkCuda(cudaMallocManaged(&first1, size));
		checkCuda(cudaMallocManaged(&second1, size));
		checkCuda(cudaMallocManaged(&third1,size));

		checkCuda(cudaMallocManaged(&first2, size));
		checkCuda(cudaMallocManaged(&second2, size));
		checkCuda(cudaMallocManaged(&third2, size));
	}

	~SoATriangle()
	{
		cudaFree(first1  );
		cudaFree(second1 );
		cudaFree(third1  );
		cudaFree(first2  );
		cudaFree(second2 );
		cudaFree(third2  );
	}
};

__constant__ double c_tau;
__constant__ double c_tau_to_current_time_level;
__constant__ double c_lb;
__constant__ double c_rb;
__constant__ double c_ub;
__constant__ double c_bb;
__constant__ double c_tau_b;
__constant__ int c_x_length;
__constant__ int c_length;
__constant__ int c_n;


// ядро написано для того, чтобы избежать копирования кусочка x и y на девайс во время вычисления координат
// Посчитаем объем глобальной памяти
// Пусть N  - число элементов внутренней матрицы
//   Тогда для запуска этого ядра необходимо
//   2*N - для хранения X, Y
//   2*N - для хранения координат точки alpha на предыдущем временном слое
//   2*N - для хранения координат точки betta на предыдущем временном слое
//   2*N - для хранения координат точки gamma на предыдущем временном слое
//   2*N - для хранения координат точки theta на предыдущем временном слое
//   итого = 10*N глобальной памяти
// Дополнительно храня 8*N элементов, мы сможем избавиться от копирования x и y на карту 
// и сэкономить 2*N памяти в основном ядре расчетов
// Но самое главное - это сократить нагрузку на регистры в основном ядре, избежав спиллинга регистров в локальную память
// Это позволит запускать на расчет бОльшие сетки. Глобальную память легче маштабировать нежели регистры
__global__ void prepare_data_kernel(double *a_x, double *a_y, 
									double* first1, double* second1, double* third1, 
									double* first2, double* second2, double* third2, 
									double* alpha, double* beta, double* gamma, double* theta)
{
	for (int opt = hemiGetElementOffset(); opt < c_n; opt += hemiGetElementStride()) 
	{
		int i = opt % c_x_length + 1;
		int j = opt / c_x_length + 1;
		double x, y; 

		//   1. First of all let's compute coordinates of square vertexes.
		// Для каждой координаты i, j определяем точки квадрата, описанный вокруг этой точки.
		//   2. Now let's compute new coordinates on the previous time level of alpha, betta, gamma, theta points.
		// определяем координаты данных точек на предыдущем временном слое

		// можно еще соптимизировать константы

		// A
		x = (a_x[i - 1] + a_x[i]) / 2.;
		y = (a_y[j - 1] + a_y[j]) / 2.;
		alpha[2*opt] = x - c_tau_b * y * (1. - y) * (C_pi_device / 2. + atan(-x));;
		alpha[2*opt+1] = y - c_tau * atan((x - c_lb) * (x - c_rb) * c_tau_to_current_time_level * (y - c_ub) * (y - c_bb));

		// B
		x = (a_x[i + 1] + a_x[i]) / 2.;
		y = (a_y[j - 1] + a_y[j]) / 2.;
		beta[2*opt] = x - c_tau_b * y * (1. - y) * (C_pi_device / 2. + atan(-x));
		beta[2*opt+1] = y - c_tau * atan((x - c_lb) * (x - c_rb) * c_tau_to_current_time_level * (y - c_ub) * (y - c_bb));

		// C
		x = (a_x[i + 1] + a_x[i]) / 2.;
		y = (a_y[j + 1] + a_y[j]) / 2.;
		gamma[2*opt] = x - c_tau_b * y * (1. - y) * (C_pi_device / 2. + atan(-x));
		gamma[2*opt+1] = y - c_tau * atan((x - c_lb) * (x - c_rb) * c_tau_to_current_time_level * (y - c_ub) * (y - c_bb));

		// D
		x = (a_x[i - 1] + a_x[i]) / 2.;
		y = (a_y[j + 1] + a_y[j]) / 2.;
		theta[2*opt] = x - c_tau_b * y * (1. - y) * (C_pi_device / 2. + atan(-x));
		theta[2*opt+1] = y - c_tau * atan((x - c_lb) * (x - c_rb) * c_tau_to_current_time_level * (y - c_ub) * (y - c_bb));

		first1[2*opt] = alpha[2*opt];
		first1[2*opt+1] = alpha[2*opt + 1];
		second1[2*opt] = beta[2*opt];
		second1[2*opt+1] = beta[2*opt + 1];
		third1[2*opt] = gamma[2*opt];
		third1[2*opt+1] = gamma[2*opt + 1];

		first2[2*opt] = alpha[2*opt];
		first2[2*opt+1] = alpha[2*opt + 1];
		second2[2*opt] = theta[2*opt];
		second2[2*opt+1] = theta[2*opt + 1];
		third2[2*opt] = gamma[2*opt];
		third2[2*opt+1] = gamma[2*opt + 1]; 
	}
}

__global__ void get_quad_coord(
	double* first1, double* second1, double* third1, 
	double* first2, double* second2, double* third2, 
	const double* alNew, const double* beNew, const double* gaNew, const double* thNew) 
{
	for (int i = hemiGetElementOffset(); i < c_length; i += hemiGetElementStride()) 
	{ 
		//int opt = (i + offset)*2;
		int opt = i*2;

		/*if (i == 192801)
		{
		printf("alnew = %f\n", alNew[opt]);
		printf("alnew = %f\n", alNew[opt+1]);
		}*/

		double vectAlGa[2], vectBeTh[2]; //   -  Vectors: 1) from "alNew" to "gaNew" 2) from "beNew" to "thNew".
		double a_1LC, b_1LC, c_1LC; //   -  a_1LC * x  +  b_1LC * y  = c_1LC. Line equation through "alNew" and "gaNew".
		double a_2LC, b_2LC, c_2LC; //   -  a_2LC * x  +  b_2LC * y  = c_2LC. Line equation through "beNew" and "thNew".
		double AcrP[2];                            //   -  Across point of two lines
		double vectAlBe[2]; //   -  Vectors for computing vertexes sequence order by vector production.
		double vectAlTh[2]; //   -  Vectors for computing vertexes sequence order by vector production.
		double vectBeGa[2]; //   -  Vectors for computing vertexes sequence order by vector production.
		double vectBeAl[2]; //   -  Vectors for computing vertexes sequence order by vector production.

		// Расчитываем коэффициенты диагонали между точкой А и точкой Г

		//   3.a Let's compute coefficients of first line betweeen "alNew" and "gaNew" points.
		//   a_1LC * x  +  b_1LC * y  = c_1LC.
		vectAlGa[0] = gaNew[opt] - alNew[opt];
		vectAlGa[1] = gaNew[opt + 1] - alNew[opt + 1];
		a_1LC = vectAlGa[1];
		b_1LC = -vectAlGa[0];
		c_1LC = vectAlGa[1] * alNew[opt] - vectAlGa[0] * alNew[opt + 1];

		//   3.b Let's compute coefficients of second line betweeen "beNew" and "thNew" points.
		//   a_2LC * x  +  b_2LC * y  = c_2LC.

		vectBeTh[0] = thNew[opt] - beNew[opt];
		vectBeTh[1] = thNew[opt + 1] - beNew[opt + 1];
		a_2LC = vectBeTh[1];
		b_2LC = -vectBeTh[0];
		c_2LC = vectBeTh[1] * beNew[opt] - vectBeTh[0] * beNew[opt + 1];


		AcrP[0] = (b_1LC * c_2LC - b_2LC * c_1LC) / (b_1LC * a_2LC - b_2LC * a_1LC);
		AcrP[1] = (a_1LC * c_2LC - a_2LC * c_1LC)
			/ (-b_1LC * a_2LC + b_2LC * a_1LC);

		if (((beNew[opt + 1] - AcrP[1]) * (thNew[opt + 1] - AcrP[1])) > 0.) 
		{

			if (((alNew[opt] - AcrP[0]) * (gaNew[opt] - AcrP[0])) > 0.)
			{

				first1[opt] = alNew[opt];
				first1[opt+1] = alNew[opt + 1];
				second1[opt] = beNew[opt];
				second1[opt+1] = beNew[opt + 1];
				third1[opt] = gaNew[opt];
				third1[opt+1] = gaNew[opt + 1];

				//   Second triangle.

				first2[opt] = beNew[opt];
				first2[opt+1] = beNew[opt + 1];
				second2[opt] = thNew[opt];
				second2[opt+1] = thNew[opt + 1];

				//   Third vertex computing...

				vectAlGa[0] = gaNew[opt] - alNew[opt];
				vectAlGa[1] = gaNew[opt + 1] - alNew[opt + 1];
				vectBeTh[0] = thNew[opt] - beNew[opt];
				vectBeTh[1] = thNew[opt + 1] - beNew[opt + 1];

				if (vectAlGa[0] * vectBeTh[0] + vectAlGa[1] * vectBeTh[1] >= 0.) {
					third2[opt] = gaNew[opt];
					third2[opt+1] = gaNew[opt + 1];
				} else if (vectAlGa[0] * vectBeTh[0] + vectAlGa[1] * vectBeTh[1] < 0.) {
					third2[opt] = alNew[opt];
					third2[opt+1] = alNew[opt + 1];
				}

				continue;

			} //   "if(  (  (alNew[opt] - AcrP[0])*(gaNew[opt] - AcsP[0])  )  >  0.  )".

			//   Second criterion. Second case.

			if (((alNew[opt] - AcrP[0]) * (gaNew[opt] - AcrP[0])) <= 0.)
			{
				vectAlBe[0] = beNew[opt] - alNew[opt];
				vectAlBe[1] = beNew[opt + 1] - alNew[opt + 1];
				vectAlTh[0] = thNew[opt] - alNew[opt];
				vectAlTh[1] = thNew[opt + 1] - alNew[opt + 1];

				if (vectAlBe[0] * vectAlTh[1] - vectAlBe[1] * vectAlTh[0] < 0.) 
				{
					//   The vertex "beNew" is NOT in triangle "alNew - gaNew - thNew".
					//   Pseudo case. Anyway I need to find some solution. So

					first1[opt] = alNew[opt];
					first1[opt+1] = alNew[opt + 1];
					second1[opt] = beNew[opt];
					second1[opt+1] = beNew[opt + 1];
					third1[opt] = thNew[opt];
					third1[opt+1] = thNew[opt + 1];

					//   Second triangle.

					first2[opt] = beNew[opt];
					first2[opt+1] = beNew[opt + 1];
					second2[opt] = thNew[opt];
					second2[opt+1] = thNew[opt + 1];
					third2[opt] = gaNew[opt];
					third2[opt+1] = gaNew[opt + 1];

					continue;
				}
				else if (vectAlBe[0] * vectAlTh[1] - vectAlBe[1] * vectAlTh[0] >= 0.) 
				{
					//  It's all write. We have a good concave quadrangle.
					//   Now let's compute all vertices which I need.
					//   First triangle.
					first1[opt] = alNew[opt];
					first1[opt+1] = alNew[opt + 1];
					second1[opt] = beNew[opt];
					second1[opt+1] = beNew[opt + 1];
					third1[opt] = thNew[opt];
					third1[opt+1] = thNew[opt + 1];

					//   Second triangle.

					first2[opt] = beNew[opt];
					first2[opt+1] = beNew[opt + 1];
					second2[opt] = thNew[opt];
					second2[opt+1] = thNew[opt + 1];
					third2[opt] = gaNew[opt];
					third2[opt+1] = gaNew[opt + 1];

					continue;
				}

			} //   "if(  (  (alNew[opt] - AcsP[0])*(gaNew[opt] - AcsP[0])  )  <=  0.  )".   //   Last second case of second criterion.

		} //   end of "if (  (  (beNew[opt + 1] - AcrP[1])*(thNew[opt + 1] - AcrP[1])  )  >  0.  )"

		//  Now let's consider SECOND case 5.b "(  (beNew[opt + 1] - AcrP[1])*(thNew[opt + 1] - AcrP[1])  )  <= 0."

		if (((beNew[opt + 1] - AcrP[1]) * (thNew[opt + 1] - AcrP[1])) <= 0.) 
		{
			if (((alNew[opt] - AcrP[0]) * (gaNew[opt] - AcrP[0])) > 0.) 
			{
				//  It means the across point IS NOT between "alNew" and "gaNew" vertices by Ox-axis?

				//   O.K. the quadrangle IS NOT CONVEX. Is it concave or pseudo? Third criterion.

				vectBeGa[0] = gaNew[opt] - beNew[opt];
				vectBeGa[1] = gaNew[opt + 1] - beNew[opt + 1];
				vectBeAl[0] = alNew[opt] - beNew[opt];
				vectBeAl[1] = alNew[opt + 1] - beNew[opt + 1];

				if (vectBeGa[0] * vectBeAl[1] - vectBeGa[1] * vectBeAl[0] >= 0.)
				{

					//   The quadrangle is concave. First triangle.

					first1[opt] = alNew[opt];
					first1[opt+1] = alNew[opt + 1];
					second1[opt] = beNew[opt];
					second1[opt+1] = beNew[opt + 1];
					third1[opt] = gaNew[opt];
					third1[opt+1] = gaNew[opt + 1];

					//   Second triangle.

					first2[opt] = alNew[opt];
					first2[opt+1] = alNew[opt + 1];
					second2[opt] = thNew[opt];
					second2[opt+1] = thNew[opt + 1];
					third2[opt] = gaNew[opt];
					third2[opt+1] = gaNew[opt + 1];

					continue;
				}
				else if (vectBeGa[0] * vectBeAl[1] - vectBeGa[1] * vectBeAl[0] < 0.) 
				{

					//   This concave quadrangle do has NO write anticlockwise vertices sequence order. It's pseudo.

					first1[opt] = alNew[opt];
					first1[opt+1] = alNew[opt + 1];
					second1[opt] = beNew[opt];
					second1[opt+1] = beNew[opt + 1];
					third1[opt] = gaNew[opt];
					third1[opt+1] = gaNew[opt + 1];

					//   Second triangle.

					first2[opt] = alNew[opt];
					first2[opt+1] = alNew[opt + 1];
					second2[opt] = thNew[opt];
					second2[opt+1] = thNew[opt + 1];
					third2[opt] = gaNew[opt];
					third2[opt+1] = gaNew[opt + 1];

					continue;
				}
			} //   end of "if(  (  (alNew[opt] - AcrP[0])*(gaNew[opt] - AcsP[0])  )  >  0.  )". First case of second criterion.

			//   Second criterion. Second case.

			if (((alNew[opt] - AcrP[0]) * (gaNew[opt] - AcrP[0])) <= 0.) {
				//   O.K. the quadrangle is convex. Is it has the same vertices sequence order.

				vectAlBe[0] = beNew[opt] - alNew[opt];

				vectAlBe[1] = beNew[opt + 1] - alNew[opt + 1];

				vectAlTh[0] = thNew[opt] - alNew[opt];

				vectAlTh[1] = thNew[opt + 1] - alNew[opt + 1];

				if (vectAlBe[0] * vectAlTh[1] - vectAlBe[1] * vectAlTh[0] >= 0.) {

					//   Convex quadrangle DO HAS WRITE anticlockwise vertices sequence order. It's convex.

					first1[opt] = alNew[opt];
					first1[opt+1] = alNew[opt + 1];
					second1[opt] = beNew[opt];
					second1[opt+1] = beNew[opt + 1];
					third1[opt] = gaNew[opt];
					third1[opt+1] = gaNew[opt + 1];

					//   Second triangle.

					first2[opt] = alNew[opt];
					first2[opt+1] = alNew[opt + 1];
					second2[opt] = thNew[opt];
					second2[opt+1] = thNew[opt + 1];
					third2[opt] = gaNew[opt];
					third2[opt+1] = gaNew[opt + 1]; 

					continue;
				} 
				else if (vectAlBe[0] * vectAlTh[1] - vectAlBe[1] * vectAlTh[0] < 0.)
				{
					first1[opt] = alNew[opt];
					first1[opt+1] = alNew[opt + 1];
					second1[opt] = beNew[opt];
					second1[opt+1] = beNew[opt + 1];
					third1[opt] = gaNew[opt];
					third1[opt+1] = gaNew[opt + 1];
					first2[opt] = alNew[opt];
					first2[opt+1] = alNew[opt + 1];
					second2[opt] = thNew[opt];
					second2[opt+1] = thNew[opt + 1];
					third2[opt] = gaNew[opt];
					third2[opt+1] = gaNew[opt + 1];
					continue;
				}
			}
		}
	}
}

void convert(TriangleResult* result, SoATriangle* tr, int n)
{
	int k = 0;

	for (int i = 0; i < n; i++)
	{
		result->f[i].first[0] = tr->first1[k];
		result->f[i].first[1] = tr->first1[k+1];
		result->f[i].second[0] = tr->second1[k];
		result->f[i].second[1] = tr->second1[k+1];
		result->f[i].third[0] = tr->third1[k];
		result->f[i].third[1] = tr->third1[k+1];
		result->s[i].first[0] = tr->first2[k];
		result->s[i].first[1] = tr->first2[k+1];
		result->s[i].second[0] = tr->second2[k];
		result->s[i].second[1] = tr->second2[k+1];
		result->s[i].third[0] = tr->third2[k];
		result->s[i].third[1] = tr->third2[k+1];
		k+=2;
	}

}	

float get_quad_coord(TriangleResult* result, ComputeParameters* p, int gridSize, int blockSize) 
{
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	size_t size(0), length(0), n(0);

	double tau_to_current_time_level = (1. + p->currentTimeLevel * p->tau) / 10.;
	double tau_b = p->b * p->tau;
	float elapsedTime;
	n = p->get_inner_matrix_size();
	length = p->get_inner_matrix_size();


	// Start record
	cudaEventRecord(start, 0);

	cudaMemcpyToSymbol(c_tau, &p->tau, sizeof(double));
	cudaMemcpyToSymbol(c_lb,  &p->lb, sizeof(double));
	cudaMemcpyToSymbol(c_rb,  &p->rb, sizeof(double));
	cudaMemcpyToSymbol(c_bb,  &p->bb, sizeof(double));
	cudaMemcpyToSymbol(c_ub,  &p->ub, sizeof(double));
	cudaMemcpyToSymbol(c_tau_b, &tau_b, sizeof(double));
	cudaMemcpyToSymbol(c_x_length, &result->x_length, sizeof(int));
	cudaMemcpyToSymbol(c_length, &length, sizeof(int));
	cudaMemcpyToSymbol(c_tau_to_current_time_level, &tau_to_current_time_level, sizeof(double));
	cudaMemcpyToSymbol(c_n, &n, sizeof(int));

	// сначала расчитаем координаты точек на предыдущем временном слое
	// резать сразу не будем, будем грузить все данные сразу
	double *x = NULL, *y = NULL, *alpha = NULL, *beta = NULL, *gamma = NULL, *theta = NULL;

	//нельзя копировать столько данных, у х длина другая = N
	size = sizeof(double) * p->get_real_x_size();
	checkCuda(cudaMallocManaged(&x, size)    );

	memcpy(x, p->x, size);
	size = sizeof(double) * p->get_real_y_size();
	checkCuda(cudaMallocManaged(&y, size)    );
	memcpy(y, p->y, size);

	size = sizeof(double) * p->get_inner_matrix_size();
	checkCuda(cudaMallocManaged(&alpha, 2*size));
	checkCuda(cudaMallocManaged(&beta,  2*size));
	checkCuda(cudaMallocManaged(&gamma, 2*size));
	checkCuda(cudaMallocManaged(&theta, 2*size));

	SoATriangle *soa = new SoATriangle(p->get_inner_matrix_size());

	prepare_data_kernel<<<gridSize, blockSize>>>(x, y, 
		soa->first1, soa->second1, soa->third1, 
		soa->first2, soa->second2, soa->third2,
		alpha, beta, gamma, theta);

	convert(result, soa, p->get_inner_matrix_size());

	cudaFree(x);
	cudaFree(y);
	cudaFree(alpha);
	cudaFree(beta);
	cudaFree(gamma);
	cudaFree(theta);

	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsedTime, start, stop); // that's our time!
	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	cudaDeviceReset();

	

	delete soa;
	return elapsedTime;
}