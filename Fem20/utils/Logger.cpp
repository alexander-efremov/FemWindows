#include "../headers/Logger.h"

void Logger::WriteComputeResultsToFile()
{
	if (_compute_results == NULL)
		return;

	std::ofstream resultCompute (_filename.c_str(), std::ios_base::app);
	if (!resultCompute.is_open()) return;
	resultCompute << currentDateTime() << std::endl;
	resultCompute << "Count of time steps is " << _compute_results->timeStepCount << std::endl;
	resultCompute << "Size of grid is " << _compute_results->xSize << " x " << _compute_results->ySize << std::endl;
	resultCompute << "CPU time is " << _compute_results->cpuTime << " ms. " << std::endl;
	resultCompute << "OpenMP time is " << _compute_results->openmpTime << " ms. " << std::endl;
	resultCompute << "GPU time is " << _compute_results->gpuTime << " ms. " << std::endl;
	resultCompute << "CPU / GPU = " << _compute_results->cpuTime / _compute_results->gpuTime << " times" << std::endl;
	resultCompute << "CPU / OpenMP = " << _compute_results->cpuTime / _compute_results->openmpTime << " times" << std::endl;
	resultCompute << "OpenMP / GPU = " << _compute_results->openmpTime / _compute_results->gpuTime << " times" << std::endl;
	resultCompute << std::endl;
	resultCompute << std::endl;
	resultCompute.close();
}

void Logger::WriteComputeResultsToConsole()
{
	if (_compute_results == NULL)
		return;
	std::cout << "\nCPU time is " << _compute_results->cpuTime << " ms. " << std::endl;
	std::cout << "OpenMP time is " << _compute_results->openmpTime << " ms. " << std::endl;
	std::cout << "GPU time is " << _compute_results->gpuTime << " ms. " << std::endl;
	std::cout << "CPU / GPU = " << _compute_results->cpuTime / _compute_results->gpuTime << " times" << std::endl;
	std::cout << "CPU / OpenMP = " << _compute_results->cpuTime / _compute_results->openmpTime << " times" << std::endl;
	std::cout << "OpenMP / GPU = " << _compute_results->openmpTime / _compute_results->gpuTime << " times" << std::endl;
}

void Logger::InitComputeResults(ComputeResults *results)
{
	_compute_results = results;
}

void Logger::SetPrecision(int precision)
{
	_precision = precision;
}

void Logger::WriteCpuResultToFile(double *matrix, int n, int m, int timeStepCount)
{
	std::string cpu_filename = "cpu.txt";
	std::stringstream cpu_st;
	cpu_st << n << "_" << m << "_" << timeStepCount << "_" << cpu_filename;
	WriteToFile(cpu_st.str(), matrix, n, m);
}

void Logger::WriteOpenMPResultToFile(double *matrix, int n, int m, int timeStepCount)
{
	std::string openmp_filename = "openmp.txt";
	std::stringstream openmp_st;
	openmp_st << n << "_" << m << "_" << timeStepCount << "_" << openmp_filename;
	WriteToFile(openmp_st.str(), matrix, n, m);
}

void Logger::WriteGpuResultToFile(double *matrix, int n, int m, int timeStepCount)
{
	std::string gpu_filename = "gpu.txt";
	std::stringstream gpu_st;
	gpu_st << n << "_" << m << "_" << timeStepCount << "_" << gpu_filename;
	WriteToFile(gpu_st.str(), matrix, n, m);
}

void Logger::WriteToFile(const std::string &filename, double *matrix, int n, int m)
{
	std::ofstream file (filename.c_str());
	if (!file.is_open()) return;
	for( int k = 0; k < m; k++ )
	{
		for( int j = 0; j < n; j++ )
		{
			int current = n*k + j;
			file << std::fixed << std::setprecision(_precision) << matrix[current] << " ";
		}
		file << std::endl;
	}
	file.close();
}