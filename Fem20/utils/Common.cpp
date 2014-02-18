#include <stdio.h>
#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <iomanip>
#include "../headers/Common.h"

#define LINUX

void print_matrix(const int tempY, const int tempX, const double* array)
{
	for(int i=0;i<tempY;i++)
	{
		for(int j=0;j<tempX;j++)
		{
			printf("%f ", array[tempX*i + j]);
		}
		printf("\n");
	}
}

void print_matrix(const int tempY, const int tempX, const int* array)
{
	for(int i=0;i<tempY;i++)
	{
		for(int j=0;j<tempX;j++)
		{
			printf("%d ", array[tempX*i + j]);
		}
		printf("\n");
	}
}

void print_vector(const int n, const double* vector)
{
	for(int i = 0; i < n; i++)
	{
		printf("%f ", vector[i]);
	}
	printf("\n");
}

// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const std::string currentDateTime()
{
	time_t now = time(0);
	struct tm tstruct;
	char buf[80];
	tstruct = *localtime(&now);
	strftime(buf, sizeof(buf), "%Y-%m-%d_%X", &tstruct);

	return buf;
}

/*void write_graph(double cpu, double cpu_openmp, double gpu, const int size_x, const int size_y, double timeLevel, const std::string &datetime, const std::string &path)
{
	#ifdef LINUX	
	std::ofstream plot_data ("plot_data");
	
	plot_data << "name CPU OpenMP GPU" << std::endl;
	plot_data << "Время ";	plot_data << cpu << " " << cpu_openmp << " " << gpu;
	plot_data.close();

	int upper_bound = (int)std::max(cpu, cpu_openmp);
	upper_bound = (int)std::max((double)upper_bound, gpu);
	upper_bound +=  50;

	std::ostringstream hist_name_str;
	hist_name_str << path << "/" << size_x << "_" << size_y << "_" << timeLevel << "_" << "graph" << "_" << datetime << ".png";		
	std::string hist_name(hist_name_str.str());
	FILE *pipe = popen("gnuplot -persist","w");		
	fprintf(pipe, "set terminal pngcairo  transparent enhanced \n");
	fprintf(pipe, "set output '%s'\n", hist_name.c_str());
	fprintf(pipe, "set title \"Сравнение времени счета CPU, GPU, OpenMP\"\n");
	fprintf(pipe, "set boxwidth 0.9 absolute\n");
	fprintf(pipe, "set style fill   solid 1.00 border lt -1\n");
	fprintf(pipe, "set key inside right top vertical Right noreverse noenhanced autotitles nobox\n");
	fprintf(pipe, "set style histogram clustered gap 5 title  offset character 0, 0, 0\n");
	fprintf(pipe, "set datafile missing '-'\n");
	fprintf(pipe, "set style data histograms\n");
	fprintf(pipe, "set yrange [0:%d]\n", upper_bound);
	fprintf(pipe, "plot './plot_data' using 2:xtic(1) ti col, '' using 3 ti col, '' using 4 ti col\n");		
	pclose(pipe); 
	remove("plot_data");
	#endif
}
*/

void write_csv(double cpu, double cpu_openmp, double gpu, const int size_x, const int size_y, double timeLevel, const std::string &datetime, const std::string &path)
{
	char csv_delimiter = ';';
	std::ostringstream name;
	name << path << "/" << size_x << "_" << size_y << "_" << timeLevel << "_" << "csv" << "_" << datetime << ".dat";
	std::ofstream plot_data (name.str().c_str(), std::ios::out);
	plot_data.imbue(std::locale(std::cout.getloc(), new DecimalSeparator<char>(',')));
	plot_data << size_x << "x" << size_y << csv_delimiter << timeLevel << csv_delimiter << cpu << csv_delimiter << cpu_openmp << csv_delimiter << gpu <<
		csv_delimiter << cpu/cpu_openmp << csv_delimiter << gpu/cpu_openmp << csv_delimiter << cpu/gpu << "\n";
	plot_data.close();
}

void write_info(double cpu, double cpu_openmp, double gpu, const int size_x, const int size_y, double timeLevel)
{
	std::string path = "/home/jane/NVIDIA_GPU_Computing_SDK/C/src/!FemCuda/compute_result";
	std::string datetime = currentDateTime();
	//	write_graph(cpu, cpu_openmp, gpu, size_x, size_y, timeLevel, datetime, path);
	write_csv(cpu, cpu_openmp, gpu, size_x, size_y, timeLevel, datetime, path);
}

void write_openmp_stress_test_info(const std::string &filename, int threadNumber, int n, int m, int tl_number, double time_ms, int ht_on, bool append)
{
	std::ofstream f;
	std::stringstream out;
	out << filename << "_" << threadNumber;
	// если hyperthread включен
	if (ht_on)
		out << "_ht_on";
	else
		out << "_ht_off";

	out << ".log";
	std::string name;
	name = out.str();
	if (append)
	{
		f.open (name.c_str(), std::ios::app);
		f << threadNumber << " ";
		f.width(15);
		f << n << "x" << m;
		f.width(10);
		f << tl_number << " ";
		f.width(10);
		f << time_ms / 1000;
		f << std::endl;
	}
	else
	{
		f.open (name.c_str(), std::ios::out);
		f << "Thread number" << " ";
		f.width(5);
		f << "Size" << " ";
		f.width(10);
		f << "tl number" << " ";
		f.width(10);
		f << "Time, s";
		f << std::endl;
		f << threadNumber << " ";
		f.width(15);
		f << n << "x" << m;;
		f.width(10);
		f << tl_number << " ";
		f.width(10);
		f << time_ms / 1000;
		f << std::endl;
	}

	f.close();
}