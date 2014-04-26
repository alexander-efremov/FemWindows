#include "ModelDataProvider.h"
#include "FileReader.h"


double* ModelDataProvider::GetModelData(enum ModelDataType type)
{
	FileReader reader;

	switch (type)
	{
	case Model11:
		return reader.ReadMatrixFromTextFile("11_11_50_cpu_model.txt", 11, 11);
	case Model21:
		return reader.ReadMatrixFromTextFile("21_21_100_cpu_model.txt", 21, 21);
	case Model41:
		return reader.ReadMatrixFromTextFile("41_41_200_cpu_model.txt", 41, 41);
	case Model81:
		return reader.ReadMatrixFromTextFile("81_81_400_cpu_model.txt", 81, 81);
	case Model161:
		return reader.ReadMatrixFromTextFile("161_161_800_cpu_model.txt", 161, 161);
	case Model321:
		return reader.ReadMatrixFromTextFile("321_321_1600_cpu_model.txt", 321, 321);
	case Model641:
		return reader.ReadMatrixFromTextFile("641_641_3200_cpu_model.txt", 641, 641);
	case Model1281:
		return reader.ReadMatrixFromTextFile("1281_1281_6400_cpu_model.txt", 1281, 1281);
	default:
		break;
	}
	return NULL;
}

double* ModelDataProvider::GetModelData1tl(int type)
{
	FileReader reader;

	switch (type)
	{
	case 0:
		return reader.ReadMatrixFromTextFile("11_11_50_cpu_model_1tl.txt", 11, 11);
	}
	return NULL;
}

double* ModelDataProvider::GetModelData(int type)
{
	FileReader reader;

	switch (type)
	{
	case 0:
		return reader.ReadMatrixFromTextFile("11_11_50_cpu_model.txt", 11, 11);
	case 1:
		return reader.ReadMatrixFromTextFile("21_21_100_cpu_model.txt", 21, 21);
	case 2:
		return reader.ReadMatrixFromTextFile("41_41_200_cpu_model.txt", 41, 41);
	case 3:
		return reader.ReadMatrixFromTextFile("81_81_400_cpu_model.txt", 81, 81);
	case 4:
		return reader.ReadMatrixFromTextFile("161_161_800_cpu_model.txt", 161, 161);
	case 5:
		return reader.ReadMatrixFromTextFile("321_321_1600_cpu_model.txt", 321, 321);
	case 6:
		return reader.ReadMatrixFromTextFile("641_641_3200_cpu_model.txt", 641, 641);
	case 7:
		return reader.ReadMatrixFromTextFile("1281_1281_6400_cpu_model.txt", 1281, 1281);
	default:
		break;
	}
	return NULL;
}