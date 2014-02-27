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

__constant__ double c_tau;
__constant__ double c_tau_to_current_time_level;
__constant__ double c_lb;
__constant__ double c_rb;
__constant__ double c_ub;
__constant__ double c_bb;
__constant__ double c_tau_b;
__constant__ double c_pi_half;
__constant__ int c_x_length;
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
__global__ void get_square_coord(double* first1, double* second1, double* third1,
	double* first2, double* second2, double* third2)
{
	for (int opt = hemiGetElementOffset(); opt < c_n; opt += hemiGetElementStride())
	{
		int i = opt % c_x_length + 1;
		int j = opt / c_x_length + 1;
		double x, y, h = 1. / (c_x_length + 1);
		
		// A
		x = (h*(i - 1) + h*i) / 2.;
		y = (h*(j - 1) + h*j) / 2.;
		first1[2 * opt] = first2[2 * opt] = x - c_tau_b * y * (1. - y) * (c_pi_half + atan(-x));
		first1[2 * opt + 1] = first2[2 * opt + 1] = y - c_tau * atan((x - c_lb) * (x - c_rb) * c_tau_to_current_time_level * (y - c_ub) * (y - c_bb));

		// B
		x = (h*(i + 1) + h*i) / 2.;
	    //	y = (h*(j - 1) + h*j) / 2.; // это значение совпадает со значением для предыдущей точки значит его можно не расчитывать
		second1[2 * opt] = x - c_tau_b * y * (1. - y) * (c_pi_half + atan(-x));
		second1[2 * opt + 1] = y - c_tau * atan((x - c_lb) * (x - c_rb) * c_tau_to_current_time_level * (y - c_ub) * (y - c_bb));

		// C
		//x = (a_x[i + 1] + a_x[i]) / 2.; // это значение совпадает со значением для предыдущей точки значит его можно не расчитывать
		y = (h*(j + 1) + h*j) / 2.;
		third1[2 * opt] = third2[2 * opt] = x - c_tau_b * y * (1. - y) * (c_pi_half + atan(-x));
		third1[2 * opt + 1] = third2[2 * opt + 1] = y - c_tau * atan((x - c_lb) * (x - c_rb) * c_tau_to_current_time_level * (y - c_ub) * (y - c_bb));

		// D
		x = (h*(i - 1) + h*i) / 2.;
		//y = (a_y[j + 1] + a_y[j]) / 2.; // это значение совпадает со значением для предыдущей точки значит его можно не расчитывать
		second2[2 * opt] = x - c_tau_b * y * (1. - y) * (c_pi_half + atan(-x));
		second2[2 * opt + 1] = y - c_tau * atan((x - c_lb) * (x - c_rb) * c_tau_to_current_time_level * (y - c_ub) * (y - c_bb));
	}
}

void convert(TriangleResult* result, double* first1, double* second1, double* third1,
	double* first2, double* second2, double* third2, int n)
{
	int k = 0;

	for (int i = 0; i < n; i++)
	{
		result->f[i].first[0] = first1[k];
		result->f[i].first[1] = first1[k + 1];
		result->f[i].second[0] = second1[k];
		result->f[i].second[1] = second1[k + 1];
		result->f[i].third[0] = third1[k];
		result->f[i].third[1] = third1[k + 1];
		result->s[i].first[0] = first2[k];
		result->s[i].first[1] = first2[k + 1];
		result->s[i].second[0] = second2[k];
		result->s[i].second[1] = second2[k + 1];
		result->s[i].third[0] = third2[k];
		result->s[i].third[1] = third2[k + 1];
		k += 2;
	}

}

float get_quad_coord(TriangleResult* result, ComputeParameters* p, int gridSize, int blockSize)
{
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	size_t size(0), n(0);
	double temp(0);
	n = p->get_inner_matrix_size();


	float elapsedTime;
	double *first1 = NULL, *second1 = NULL, *third1 = NULL, *first2 = NULL, *second2 = NULL, *third2 = NULL;

	// Start record
	cudaEventRecord(start, 0);

	cudaMemcpyToSymbol(c_tau, &p->tau, sizeof(double));
	cudaMemcpyToSymbol(c_lb, &p->lb, sizeof(double));
	cudaMemcpyToSymbol(c_rb, &p->rb, sizeof(double));
	cudaMemcpyToSymbol(c_bb, &p->bb, sizeof(double));
	cudaMemcpyToSymbol(c_ub, &p->ub, sizeof(double));
	cudaMemcpyToSymbol(c_n, &n, sizeof(int));
	cudaMemcpyToSymbol(c_x_length, &result->x_length, sizeof(int));

	temp = (1. + p->currentTimeLevel * p->tau) / 10.;
	cudaMemcpyToSymbol(c_tau_to_current_time_level, &temp, sizeof(double));

	temp = p->b * p->tau;
	cudaMemcpyToSymbol(c_tau_b, &temp, sizeof(double));

	temp = C_pi_device / 2.;
	cudaMemcpyToSymbol(c_pi_half, &temp, sizeof(double));
	
	size = 2 * sizeof(double)*n;
	checkCuda(cudaMallocManaged(&first1, size));
	checkCuda(cudaMallocManaged(&second1, size));
	checkCuda(cudaMallocManaged(&third1, size));
	checkCuda(cudaMallocManaged(&first2, size));
	checkCuda(cudaMallocManaged(&second2, size));
	checkCuda(cudaMallocManaged(&third2, size));

	get_square_coord<< <gridSize, blockSize >> >(first1, second1, third1, first2, second2, third2);
	
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsedTime, start, stop);

	convert(result, first1, second1, third1, first2, second2, third2, n);

	cudaFree(first1);
	cudaFree(second1);
	cudaFree(third1);
	cudaFree(first2);
	cudaFree(second2);
	cudaFree(third2);
	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	
	cudaDeviceReset();
	return elapsedTime;
}