#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include "Common.h"

class Logger
{
private:
	std::string _filename;
	ComputeResults* _compute_results;
	int _precision;
	void WriteToFile(const std::string &, double*, int, int);

public:

	Logger()
	{
		_filename = "info.log";
		_precision = 8;
	}

	Logger(std::string filename)
	{
		_filename = filename;
		_precision = 8;
	}

	void SetPrecision(int);
	void InitComputeResults(ComputeResults*);
	void WriteComputeResultsToFile();
	void WriteComputeResultsToConsole();
	void WriteCpuResultToFile(double*, int, int, int);
	void WriteOpenMPResultToFile(double*, int, int, int);
	void WriteGpuResultToFile(double*, int, int, int);
};