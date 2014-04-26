
#include "Timer.h"
#include "gtest/gtest.h"
#include "gtest/gtest-spi.h"
#include "ModelDataProvider.h"
#include <assert.h>     /* assert */
#ifndef COMMON_H_INCLUDED
#include "Common.h"
#endif
#define FULL_TEST false


class TestBase : public testing::Test
{
protected:
	double _accuracy;



	ModelDataProvider _modelDataProvider;

	TestBase()
	{
		_accuracy = 1.0e-8;
		initCompOfGlVar();
		_modelDataProvider = ModelDataProvider();
	}

	virtual ~TestBase()
	{
		// You can do clean-up work that doesn't throw exceptions here.
		memClean();
	}

	// If the constructor and destructor are not enough for setting up
	// and cleaning up each test, you can define the following methods:

	virtual void SetUp()
	{
		// Code here will be called immediately after the constructor (right
		// before each test).
	}

	virtual void TearDown()
	{
		// Code here will be called immediately after each test (right
		// before the destructor).
	}

	int GetSize()
	{
		return C_numOfOXSt * C_numOfOYSt;
	}

	void print_matrix(int n, int m, double* a, int precision = 8)
	{
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < m; ++j)
			{
				int k = i*n + j;
				switch (precision)
				{
				case 1:
					printf("%.1f ", a[k]);
					break;
				case 2:
					printf("%.2f ", a[k]);
					break;
				case 3:
					printf("%.3f ", a[k]);
					break;
				case 4:
					printf("%.4f ", a[k]);
					break;
				case 5:
					printf("%.5f ", a[k]);
					break;
				case 6:
					printf("%.6f ", a[k]);
					break;
				case 7:
					printf("%.7f ", a[k]);
					break;
				case 8:
					printf("%.8f ", a[k]);
					break;
				}
			}
			printf("\n");
		}
	}
};

class cputest : public TestBase
{
protected:

	cputest()
	{
	}

	virtual ~cputest()
	{
	}

	double* GetCpuToLevel(int level)
	{
		return solve_cpu_test(C_par_a, C_par_b, C_lbDom, C_rbDom, C_bbDom,
			C_ubDom, C_tau, C_numOfTSt, masOX, C_numOfOXSt, masOY,
			C_numOfOYSt, level);
	}
};

TEST_F(cputest, CpuTestModel11)
{
	double* data = _modelDataProvider.GetModelData(Model11);
	double* result = GetCpuToLevel(0);

	for (int i = 0; i < GetSize(); i++)
	{
		ASSERT_TRUE(abs(data[i] - result[i]) <= _accuracy) << data[i] << " " << result[i];

	}
}

TEST_F(cputest, CpuTestModel21)
{
	double* data = _modelDataProvider.GetModelData(Model21);
	double* result = GetCpuToLevel(1);

	for (int i = 0; i < GetSize(); i++)
	{
		EXPECT_TRUE(data[i] - result[i] <= _accuracy);
	}
}

TEST_F(cputest, CpuTestModel41)
{
	if (FULL_TEST)
	{
		double* data = _modelDataProvider.GetModelData(Model41);
		double* result = GetCpuToLevel(2);

		for (int i = 0; i < GetSize(); i++)
		{
			EXPECT_TRUE(data[i] - result[i] <= _accuracy);
		}
	}
}

TEST_F(cputest, CpuTestModel81)
{
	if (FULL_TEST)
	{
		double* data = _modelDataProvider.GetModelData(Model81);
		double* result = GetCpuToLevel(3);

		for (int i = 0; i < GetSize(); i++)
		{
			ASSERT_TRUE(data[i] - result[i] <= _accuracy);
		}
	}
}

TEST_F(cputest, CpuTestModel161)
{
	if (FULL_TEST)
	{
		double* data = _modelDataProvider.GetModelData(Model161);
		double* result = GetCpuToLevel(4);

		for (int i = 0; i < GetSize(); i++)
		{
			ASSERT_TRUE(data[i] - result[i] <= _accuracy);
		}
	}
}

TEST_F(cputest, CpuTestModel321)
{
	if (FULL_TEST)
	{
		double* data = _modelDataProvider.GetModelData(Model321);
		double* result = GetCpuToLevel(5);

		for (int i = 0; i < GetSize(); i++)
		{
			ASSERT_TRUE(data[i] - result[i] <= _accuracy);
		}
	}
}

TEST_F(cputest, CpuTestModel641)
{
	if (FULL_TEST)
	{
		double* data = _modelDataProvider.GetModelData(Model641);
		double* result = GetCpuToLevel(6);

		for (int i = 0; i < GetSize(); i++)
		{
			ASSERT_TRUE(data[i] - result[i] <= _accuracy);
		}
	}
}

TEST_F(cputest, CpuTestModel1281)
{
	if (FULL_TEST)
	{
		double* data = _modelDataProvider.GetModelData(Model1281);
		double* result = GetCpuToLevel(7);

		for (int i = 0; i < GetSize(); i++)
		{
			ASSERT_TRUE(data[i] - result[i] <= _accuracy);
		}
	}
}

class gputest : public TestBase
{
protected:

	ModelDataProvider _modelDataProvider;
	gputest()
	{
		_modelDataProvider = ModelDataProvider();
	}

	virtual ~gputest()
	{
	}

	double* GetCpuToLevel(int level)
	{
		return solve_cpu_test(C_par_a, C_par_b, C_lbDom, C_rbDom, C_bbDom,
			C_ubDom, C_tau, C_numOfTSt, masOX, C_numOfOXSt, masOY,
			C_numOfOYSt, level);
	}
};



TEST_F(gputest, get_quad_coord)
{
	const int finishLevel = 1;
	const int startLevel = 0;
	const double error = 1.0e-15;

	for (int level = startLevel; level < finishLevel; ++level)
	{
		std::cout << "level = " << level << std::endl;
		ComputeParameters* p = new ComputeParameters(level, false);
		p->currentTimeLevel = 1;
		TriangleResult* gpu = new TriangleResult(p);
		float t = get_quad_coord(gpu, p);
		printf("gpu time elapsed = %f\n", t);
		double firVfirT[2], secVfirT[2], thiVfirT[2];
		double firVsecT[2], secVsecT[2], thiVsecT[2];

		for (int j = 1; j < p->y_size; j++)
		{
			for (int i = 1; i < p->x_size; i++)
			{
				p->i = i;
				p->j = j;
				quadrAngleType(p->a, p->b, p->lb, p->rb, p->bb, p->ub, p->tau, p->currentTimeLevel, p->i,
					p->x,
					p->get_real_x_size(),
					p->j,
					p->y,
					p->get_real_y_size(), firVfirT, secVfirT, thiVfirT, firVsecT, secVsecT, thiVsecT);
				int c = (p->x_size - 1) * (j - 1) + (i - 1);

				bool b1 = (fabs(gpu->first1[2 * c] - firVfirT[0]) < error) && (fabs(gpu->first1[2 * c + 1] - firVfirT[1]) < error);
				ASSERT_TRUE(b1) << "c = " << c << std::endl;
				bool b2 = (fabs(gpu->second1[2 * c] - secVfirT[0]) < error) && (fabs(gpu->second1[2 * c + 1] - secVfirT[1]) < error);
				ASSERT_TRUE(b2) << "c = " << 2 * c << std::endl;
				bool b3 = (fabs(gpu->third1[2 * c] - thiVfirT[0]) < error) && (fabs(gpu->third1[2 * c + 1] - thiVfirT[1]) < error);
				ASSERT_TRUE(b3) << "c = " << 2 * c << std::endl;

				ASSERT_TRUE(b1&&b2&&b3) << "c = " << c << std::endl;

				bool b4 = (fabs(gpu->first2[2 * c] - firVsecT[0]) < error) && (fabs(gpu->first2[2 * c + 1] - firVsecT[1]) < error);
				ASSERT_TRUE(b4) << "c = " << 2 * c << std::endl;

				bool b5 = (fabs(gpu->second2[2 * c] - secVsecT[0]) < error) && (fabs(gpu->second2[2 * c + 1] - secVsecT[1]) < error);
				ASSERT_TRUE(b5) << "c = " << 2 * c << std::endl;

				bool b6 = (fabs(gpu->third2[2 * c] - thiVsecT[0]) < error) && (fabs(gpu->third2[2 * c + 1] - thiVsecT[1]) < error);
				ASSERT_TRUE(b6) << "c = " << 2 * c << std::endl;

				ASSERT_TRUE(b3&&b5&&b6) << "c = " << 2 * c << std::endl;
			}
		}
		delete p;
		delete gpu;
	}
}

/*TEST_F(gputest, DISABLED_main_test_old)
{
const int finishLevel = 10;
const int startLevel = 0;
const double error = 1.0e-15;

for (int level = startLevel; level < finishLevel; ++level)
{
std::cout << "level = " << level << std::endl;
ComputeParameters* p = new ComputeParameters(level, false);
p->currentTimeLevel = 1;
TriangleResult* gpu = new TriangleResult(p);
float t = get_quad_coord(gpu, p);
printf("gpu time elapsed = %f\n", t);
double firVfirT[2], secVfirT[2], thiVfirT[2];
double firVsecT[2], secVsecT[2], thiVsecT[2];


for (int j = 1; j < p->y_size; j++)
{
for (int i = 1; i < p->x_size; i++)
{
p->i = i;
p->j = j;

quadrAngleType(p->a, p->b, p->lb, p->rb, p->bb, p->ub, p->tau, p->currentTimeLevel, p->i,
p->x,
p->get_real_x_size(),
p->j,
p->y,
p->get_real_y_size(), firVfirT, secVfirT, thiVfirT, firVsecT, secVsecT, thiVsecT ) ;
int c = (p->x_size - 1) * (j - 1) + (i - 1);

bool b1 = (fabs(gpu->first1[2 * c] - firVfirT[0]) < error) && (fabs(gpu->first1[2 * c + 1] - firVfirT[1]) < error);
ASSERT_TRUE(b1) << "c = " << c << std::endl;
bool b2 = (fabs(gpu->second1[2 * c] - secVfirT[0]) < error) && (fabs(gpu->second1[2 * c + 1] - secVfirT[1]) < error);
ASSERT_TRUE(b2) << "c = " << 2 * c << std::endl;
bool b3 = (fabs(gpu->third1[2 * c] - thiVfirT[0]) < error) && (fabs(gpu->third1[2 * c + 1] - thiVfirT[1]) < error);
ASSERT_TRUE(b3) << "c = " << 2 * c << std::endl;

ASSERT_TRUE(b1&&b2&&b3) << "c = " << c << std::endl;

bool b4 = (fabs(gpu->first2[2 * c] - firVsecT[0]) < error) && (fabs(gpu->first2[2 * c + 1] - firVsecT[1]) < error);
ASSERT_TRUE(b4) << "c = " << 2 * c << std::endl;

bool b5 = (fabs(gpu->second2[2 * c] - secVsecT[0]) < error) && (fabs(gpu->second2[2 * c + 1] - secVsecT[1]) < error);
ASSERT_TRUE(b5) << "c = " << 2 * c << std::endl;

bool b6 = (fabs(gpu->third2[2 * c] - thiVsecT[0]) < error) && (fabs(gpu->third2[2 * c + 1] - thiVsecT[1]) < error);
ASSERT_TRUE(b6) << "c = " << 2 * c << std::endl;

ASSERT_TRUE(b3&&b5&&b6) << "c = " << 2 * c << std::endl;
}
}
delete p;
delete gpu;
}
}*/

TEST_F(gputest, get_quad_coord_te)
{
	double time_cpu(-1), time_gpu(0);
	double firVfirT[2], secVfirT[2], thiVfirT[2];
	double firVsecT[2], secVsecT[2], thiVsecT[2];
	int finish_level = 1;
	int start_level = 0;


	for (int level = start_level; level < finish_level; level++)
	{
		printf("level %d \n", level);
		ComputeParameters* p = new ComputeParameters(level, false);
		p->currentTimeLevel = 1;

		std::cout << *p << std::endl;

		printf("Start GPU\n");
		TriangleResult* gpu = new TriangleResult(p);

		time_gpu = get_quad_coord(gpu, p);

		printf("End GPU\n");

		printf("Start CPU\n");

		int y_size;
		int x_size;
		y_size = p->y_size;
		x_size = p->x_size;
		StartTimer();

		for (int j = 1; j < y_size; j++)
		{
			for (int i = 1; i < x_size; i++)
			{
				p->i = i;
				p->j = j;
				quadrAngleType(p->a, p->b, p->lb, p->rb, p->bb, p->ub, p->tau, p->currentTimeLevel, p->i,
					p->x,
					p->get_real_x_size(),
					p->j,
					p->y,
					p->get_real_y_size(), firVfirT, secVfirT, thiVfirT, firVsecT, secVsecT, thiVsecT);
			}
		}

		time_cpu = GetTimer();

		printf("End CPU\n");
		printf("CPU time is = %f\n", time_cpu);
		printf("GPU time is = %f\n", time_gpu);
		printf("CPU/GPU = %f\n", time_cpu / time_gpu);
		delete gpu;
		delete p;
	}
}



/*TEST_F(gputest, DISABLED_main_test)
{
const int finishLevel = 1;
const int startLevel = 0;
const double error = 1.0e-8;

for (int level = startLevel; level < finishLevel; ++level)
{
std::cout << "level = " << level << std::endl;
ComputeParameters* p = new ComputeParameters(level, true);
ASSERT_TRUE(p->result != NULL);
float gpu_time = solve_at_gpu(p);
ASSERT_TRUE(gpu_time != -1);
double* data = _modelDataProvider.GetModelData(level);
print_matrix(p->get_real_x_size(), p->get_real_y_size(), data);
printf("%s\n", "");
print_matrix(p->get_real_x_size(), p->get_real_y_size(), p->result);
printf("%s\n", "Start testing...");
for (int i = 0; i < p->get_real_matrix_size(); i++)
{
ASSERT_TRUE(fabs(data[i] - p->result[i]) <= error) << i << " " << data[i] << " " << p->result[i] << std::endl;
}

delete p;
}
}*/

TEST_F(gputest, main_test_1tl_boundaries)
{
	const int finishLevel = 1;
	const int startLevel = 0;
	const double error = 1.0e-8;

	ComputeParameters* p = new ComputeParameters(0, true);
	ASSERT_TRUE(p->result != NULL);
	float gpu_time = solve_at_gpu(p);
	ASSERT_TRUE(gpu_time != -1);
	double* data = _modelDataProvider.GetModelData1tl(0);
	print_matrix(p->get_real_x_size(), p->get_real_y_size(), data);
	printf("%s\n", "");
	print_matrix(p->get_real_x_size(), p->get_real_y_size(), p->result);
	printf("%s\n", "Start testing...");
	for (int i = 0; i < p->get_real_matrix_size(); i++)
	{
		int n = i % p->get_real_x_size();
		int m = i / p->get_real_y_size();

		// расчет границы
		if (m == 0 || n == 0 || m == p->get_real_y_size() - 1 || n == p->get_real_x_size() - 1)
		{
			ASSERT_TRUE(fabs(data[i] - p->result[i]) <= error) << i << " " << data[i] << " " << p->result[i] << std::endl;
		}
	}

	delete p;
}

TEST_F(gputest, main_test_1tl_inner)
{
	const int finishLevel = 1;
	const int startLevel = 0;
	const double error = 1.0e-8;

	ComputeParameters* p = new ComputeParameters(0, true);
	ASSERT_TRUE(p->result != NULL);
	float gpu_time = solve_at_gpu(p);
	ASSERT_TRUE(gpu_time != -1);
	double* data = _modelDataProvider.GetModelData1tl(0);
	printf("%d\n", p->get_real_x_size());
	print_matrix(p->get_real_x_size(), p->get_real_y_size(), data, 5);
	printf("%s\n", "");
	print_matrix(p->get_real_x_size(), p->get_real_y_size(), p->result, 5);
	printf("%s\n", "Start testing...");
	for (int i = 0; i < p->get_real_matrix_size(); i++)
	{
		int n = i % p->get_real_x_size();
		int m = i / p->get_real_y_size();

		// расчет границы
		if (m == 0 || n == 0 || m == p->get_real_y_size() - 1 || n == p->get_real_x_size() - 1)
		{
			continue;
		}
		ASSERT_TRUE(fabs(data[i] - p->result[i]) <= error) << i << " " << data[i] << " " << p->result[i] << std::endl;
	}

	delete p;
}

TEST_F(gputest, gen_1tl)
{
	const int finishLevel = 1;
	const int startLevel = 0;
	const double error = 1.0e-8;
	double* tl1 = GetCpuToLevel(0);
	//print_matrix(11, 11, tl1);	
}

class CpuVersusGpuFunctionalFemTest : public TestBase
{
protected:

	CpuVersusGpuFunctionalFemTest()
	{
	}

	virtual ~CpuVersusGpuFunctionalFemTest()
	{
	}
};

TEST_F(CpuVersusGpuFunctionalFemTest, AnalyticSolutionEqualsTest)
{
	ComputeParameters* p = new ComputeParameters(0, false);
	double cpu = analytSolut(p->a, p->lb, p->rb, p->bb, p->ub, p->tau, p->x[1],
		p->y[1]);
	double gpu = h_analytSolut(p->tau, p->x[1], p->y[1]);
	ASSERT_EQ(cpu, gpu);
	delete p;
}

TEST_F(CpuVersusGpuFunctionalFemTest, AnalyticSolutionNotEqualsTest)
{
	ComputeParameters* p = new ComputeParameters(0, false);
	double cpu = analytSolut(p->a, p->lb, p->rb, p->bb, p->ub, p->tau, p->x[1],
		p->y[2]);
	double gpu = h_analytSolut(p->tau, p->x[1], p->y[1]);
	ASSERT_NE(cpu, gpu);
	delete p;
}

TEST_F(CpuVersusGpuFunctionalFemTest, d_bottomBound_EqualsTest)
{
	ComputeParameters* p = new ComputeParameters(0, false);
	int currentTimeLevel = 1;
	double t = p->tau * 1;
	int i = 1;

	double cpu = bottonBound(p->a, p->lb, p->rb, p->bb, p->ub, t, p->x[i]);

	p->currentTimeLevel = currentTimeLevel;
	p->i = i;
	double gpu = h_bottomBound(p);

	ASSERT_EQ(cpu, gpu);
	delete p;
}

TEST_F(CpuVersusGpuFunctionalFemTest, d_upperBound_EqualsTest)
{
	ComputeParameters* p = new ComputeParameters(0, false);
	int currentTimeLevel = 1;
	double t = p->tau * 1;
	int i = 1;

	double cpu = upperBound(p->a, p->lb, p->rb, p->bb, p->ub, t, p->x[i]);
	p->currentTimeLevel = currentTimeLevel;
	p->i = i;
	double gpu = h_upperBound(p);
	ASSERT_EQ(cpu, gpu);
	delete p;
}

TEST_F(CpuVersusGpuFunctionalFemTest, d_rightBound_EqualsTest)
{
	ComputeParameters* p = new ComputeParameters(0, false);
	int currentTimeLevel = 1;
	double t = p->tau * 1;
	int j = 1;

	double cpu = rightBound(p->a, p->lb, p->rb, p->bb, p->ub, t, p->y[j]);

	p->currentTimeLevel = currentTimeLevel;
	p->j = j;
	double gpu = h_rightBound(p);
	ASSERT_EQ(cpu, gpu);
	delete p;
}

TEST_F(CpuVersusGpuFunctionalFemTest, d_leftBound_EqualsTest)
{
	ComputeParameters* p = new ComputeParameters(0, false);
	int currentTimeLevel = 1;
	double t = p->tau * 1;
	int j = 1;

	double cpu = leftBound(p->a, p->lb, p->rb, p->bb, p->ub, t, p->y[j]);

	p->currentTimeLevel = currentTimeLevel;
	p->j = j;
	double gpu = h_leftBound(p);

	ASSERT_EQ(cpu, gpu);
	delete p;
}