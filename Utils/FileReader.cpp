#include "FileReader.h"

double* FileReader::ReadMatrixFromTextFile(const std::string &filename, int n, int m)
{
	double* result = new double[n*m];

	std::ifstream file (filename.c_str());
	if (!file.is_open()) return result;
	for( int k = 0; k < m; k++ )
	{
		for( int j = 0; j < n; j++ )
		{
			int current = n*k + j;
			file >> result[current];
		}
	}
	file.close();
	return result;
}