
#include "Timer.h"
#include "gtest/gtest.h"
#include "gtest/gtest-spi.h"
#include "ModelDataProvider.h"
#include <assert.h>     /* assert */
#ifndef COMMON_H_INCLUDED
#include "Common.h"
#endif
#define FULL_TEST false

class TestBase: public testing::Test
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
};

class CpuTest: public TestBase
{
protected:

	CpuTest()
	{
	}

	virtual ~CpuTest()
	{
	}

	double* GetCpuToLevel(int level)
	{
		return solve_cpu_test(C_par_a, C_par_b, C_lbDom, C_rbDom, C_bbDom,
		                      C_ubDom, C_tau, C_numOfTSt, masOX, C_numOfOXSt, masOY,
		                      C_numOfOYSt, level);
	}
};

TEST_F(CpuTest, CpuTestModel11)
{
	double* data = _modelDataProvider.GetModelData(Model11);
	double* result = GetCpuToLevel(0);

	for (int i = 0; i < GetSize(); i++)
	{
		ASSERT_TRUE(abs(data[i] - result[i]) <= _accuracy) << data[i] << " " << result[i];
		
	}
}

TEST_F(CpuTest, CpuTestModel21)
{
	double* data = _modelDataProvider.GetModelData(Model21);
	double* result = GetCpuToLevel(1);

	for (int i = 0; i < GetSize(); i++)
	{
		EXPECT_TRUE(data[i] - result[i] <= _accuracy);
	}
}

TEST_F(CpuTest, CpuTestModel41)
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

TEST_F(CpuTest, CpuTestModel81)
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

TEST_F(CpuTest, CpuTestModel161)
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

TEST_F(CpuTest, CpuTestModel321)
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

TEST_F(CpuTest, CpuTestModel641)
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

TEST_F(CpuTest, CpuTestModel1281)
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

class GpuTest: public TestBase
{
protected:


	GpuTest()
	{
	}

	virtual ~GpuTest()
	{
	}

	

	double get_cpu_integr_result_internal(ComputeParameters p, double* rhoInPrevTL_asV, Triangle t1,
	                                      Triangle s1)
	{
		double buf_D = integUnderUnunifTr(
			p.a, p.b,
			//
			p.lb, p.rb,
			//
			p.bb, p.ub,
			//
			p.tau, p.currentTimeLevel,
			//
			t1.first, t1.second, t1.third,
			//
			p.x, p.x_size,
			//
			p.y, p.y_size,
			rhoInPrevTL_asV,
			p.i, p.j);


		buf_D += integUnderUnunifTr(
			p.a, p.b, //   -  Analitycal solution parameters.
			//
			p.lb, p.rb, //   -  Left and right boundaries of rectangular domain.
			//
			p.bb, p.ub, //   -  Botton and upper boundaries of rectangular domain.
			//
			p.tau, p.currentTimeLevel, //   -  Index of current time layer.
			//
			s1.first, s1.second, s1.third, //   -  Vertices of first triangle.
			//
			p.x, p.x_size, //   -  Number of OX steps.
			//
			p.y, p.y_size, //   -  Number of OY steps.
			//
			rhoInPrevTL_asV,
			p.i, p.j );

		return buf_D;
	}
};



TEST_F(GpuTest, get_quad_coord)
{
	const int finishLevel = 7;
	const int startLevel = 0;
	const int gridSize = 512;
	const int blockSize = 1024;
	const int needPrint = 0;
	const double error = 1.0e-15;

	for (int level = startLevel; level < finishLevel; ++level)
	{
		std::cout << "level = " << level << std::endl;
		ComputeParameters* p = new ComputeParameters(level, false);
		p->currentTimeLevel = 1;
		TriangleResult* gpu = new TriangleResult(p);
		get_quad_coord(gpu, p, gridSize, blockSize);

		// get cpu data
		Triangle f = Triangle();
		Triangle s = Triangle();

		for (int j = 1; j < p->y_size; j++)
		{
			for (int i = 1; i < p->x_size; i++)
			{
				p->i = i;
				p->j = j;

				int result = h_quadrAngleType(p, f, s);
				int c = (p->x_size - 1) * (j - 1) + (i - 1);

				bool b1 = (fabs(gpu->f[c].first[0] - f.first[0]) < error) && ( fabs(gpu->f[c].first[1] - f.first[1] )< error);
				ASSERT_TRUE( b1 ) << "c = " << c << std::endl << gpu->f[c] << f << std::endl;
				bool b2 = (fabs(gpu->f[c].second[0] -f.second[0]) < error) && ( fabs(gpu->f[c].second[1] - f.second[1] )< error);
				ASSERT_TRUE( b2 ) << "c = " << c << std::endl << gpu->f[c] << f << std::endl;
				bool b3 = (fabs(gpu->f[c].third[0] - f.third[0]) < error) && ( fabs(gpu->f[c].third[1] - f.third[1] )< error);
				ASSERT_TRUE( b3 ) << "c = " << c << std::endl << gpu->f[c] << f << std::endl;


				ASSERT_TRUE( b1&&b2&&b3 ) << "c = " << c << std::endl;

				bool b4 = (fabs(gpu->s[c].first[0] - s.first[0] )< error) && ( fabs(gpu->s[c].first[1] - s.first[1]) < error) ;
				ASSERT_TRUE( b4 ) << "c = " << c << std::endl << gpu->s[c] << std::endl << s << std::endl;

				bool b5 = (fabs(gpu->s[c].second[0] - s.second[0] )< error) && (fabs(gpu->s[c].second[1] - s.second[1] )< error) ;
				ASSERT_TRUE( b5 ) << "c = " << c << std::endl << gpu->s[c] << std::endl << s << std::endl;

				bool b6 = (fabs(gpu->s[c].third[0] - s.third[0]) < error) && ( fabs(gpu->s[c].third[1] - s.third[1]) < error);
				ASSERT_TRUE( b6 ) << "c = " << c << std::endl << gpu->s[c] << std::endl << s << std::endl;

				ASSERT_TRUE( b3&&b5&&b6 ) << "c = " << c << std::endl << gpu->s[c] << std::endl << s << std::endl;
			}
		}
		delete p;
		delete gpu;
	}
	std::cout << "Done!" << std::endl;
	std::cin.get();
}

TEST_F(GpuTest, get_quad_coord_time_elapsed)
{
	int gridSize = 512;
	int blockSize = 1024;
	double time_cpu = -1;
	double time_gpu = -1;
	int finish_level = 8; 
	int start_level = 7;
	int tc = 0;
	
	Triangle f = Triangle();
	Triangle s = Triangle();

	for (int level = start_level; level < finish_level; level++)
	{
		ComputeParameters* p = new ComputeParameters(level, false);
		p->currentTimeLevel = 1;
		
		tc = p->t_count;
		tc = 100;
		//tc = 2;
		//tc = 100;
		

		std::cout << *p << std::endl;

		printf("Start GPU\n");
		TriangleResult* gpu = new TriangleResult(p);
	
		time_gpu =0;
		for (int i = 1; i < tc; i++)
		{
			p->currentTimeLevel = i;
			time_gpu += get_quad_coord(gpu, p, gridSize, blockSize);
		}
		printf("End GPU\n");


		printf("Start CPU\n");

		int y_size;
		int x_size;
		y_size = p->y_size;
		x_size = p->x_size;
		StartTimer();

		for (int i = 1; i < tc; i++)
		{
			p->currentTimeLevel = i;

			for (int j = 1; j < y_size; j++)
			{
				for (int i = 1; i < x_size; i++)
				{
					p->i = i;
					p->j = j;
					int result = h_quadrAngleType(p, f, s);
				}
			}
		}
		time_cpu = GetTimer();
		printf("End CPU\n");
		printf("CPU time is = %f\n", time_cpu);
		printf("GPU time is = %f\n", time_gpu);
		printf("CPU/GPU = %f\n", time_cpu/time_gpu);
		delete gpu;
		delete p;
	}
	//std::cin.get();
}




//основной тест
//TEST_F(GpuTest, DISABLED_d_integUnderUnunifTr)
TEST_F(GpuTest, d_integUnderUnunifTr)
{
	//// параметры размерности задачи
	//int levelStep = 1;
	////int start_level =7;
	//int start_level =0;
	////int start_level =9;
	//int finish_level = 1;

	//// параметры сетки для GPU
	//int grid_size = 512;
	//int block_size = 1024;

	//for (int level = start_level; level < finish_level; level++)
	//{
	//	ComputeParameters p = GetComputeParameters(level);
	//	cout << p << std::endl;

	//	for(p.currentTimeLevel = 1; p.currentTimeLevel < 2; p.currentTimeLevel += levelStep)
	//	{
	//		// init test data
	//		p.print_info();
	//		double time_cpu = -1;
	//		double time_gpu = -1;
	//		double* rho_in_prev_tl = new double [ p.size ];

	//		//   Initial data of rho.
	//		for( int k = 0; k <= p.x_size; k++ )
	//		{
	//			for( int j = 0; j <= p.y_size; j++ )
	//			{
	//				rho_in_prev_tl[ (p.x_size+1)*k + j ] = d_initDataOfSol(p, j, k);
	//			}
	//		}

	//		// ==


	//		// get gpu result of integr
	//		printf("Start GPU triangulation...  ");
	//		TriangleResult* gpu = new TriangleResult(p);
	//		get_triangle_type(gpu, p, grid_size,
	//		                  block_size);
	//		printf("GPU triangulation is done.\n");

	//		int c = -1;
	//		double* cpu_result = new double[p.x_size * p.y_size];

	//		printf("Perform CPU integration...   ");
	//		StartTimer();

	//		for (int j = 1; j < p.y_size; j++)
	//		{
	//			for (int i = 1; i < p.x_size; i++)
	//			{
	//				p.i = i;
	//				p.j = j;
	//				c = (p.x_size - 1) * (j - 1) + (i - 1);

	//				//printf("\n i = %d, j = %d, c is = %d\n", i, j, c);

	//				cpu_result[c] = get_cpu_integr_result_internal(p, rho_in_prev_tl, gpu->f[c], gpu->s[c]);
	//			}
	//		}

	//		time_cpu = GetTimer();
	//		printf("CPU integration is done.\n");

	//		printf("Perform GPU integration...   ");
	//		StartTimer();
	//		double* gpu_integr = GetIntegrGpuResult(p, *gpu, rho_in_prev_tl, grid_size, block_size
	//		                                        , p.get_chunk_size());
	//		time_gpu = GetTimer();
	//		printf("GPU integration is done.\n");


	//		printf("CPU time is = %f ms.\n", time_cpu);
	//		printf("GPU time is = %f ms.\n", time_gpu);
	//		printf("CPU/GPU = %f\n times\n", time_cpu/time_gpu);


	//		// проверям корректность полученных результатов
	//		for (int j = 1; j < p.y_size; j++)
	//		{
	//			for (int i = 1; i < p.x_size; i++)
	//			{
	//				c = (p.x_size - 1) * (j - 1) + (i - 1);

	//				ASSERT_DOUBLE_EQ(cpu_result[c], gpu_integr[c]);
	//			}
	//		}

	//		delete[] cpu_result;
	//		delete[] gpu_integr;
	//		delete[] rho_in_prev_tl;
	//		delete gpu;
	//	}
	//}
}


class CpuVersusGpuFunctionalFemTest: public TestBase
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

TEST_F(CpuVersusGpuFunctionalFemTest, DataInitializationAnalyticSolutionEqualsTest)
{
	ComputeParameters* p = new ComputeParameters(0, false);
	double cpu = initDataOfSol(p->a, p->lb, p->rb, p->bb, p->ub, 1, p->x, 2, p->y);
	double gpu = d_initDataOfSol(p, 1, 2);
	ASSERT_EQ(cpu, gpu);
	delete p;
}

TEST_F(CpuVersusGpuFunctionalFemTest, F_Function_EqualsTest)
{
	ComputeParameters* p = new ComputeParameters(0, false);
	int currentTimeLevel = 1;
	int i = 1;
	int j = 1;
	double cpu = f_function(p->a, p->b, p->lb, p->rb, p->bb, p->ub, p->tau,
	                        currentTimeLevel, i, p->x, p->x_size, j, p->y, p->y_size);
	double gpu = h_f_function(p, currentTimeLevel, i, j);
	ASSERT_EQ(cpu, gpu);
	delete p;
}

TEST_F(CpuVersusGpuFunctionalFemTest, F_Function_NotEqualsTest)
{
	ComputeParameters* p = new ComputeParameters(0, false);
	int currentTimeLevel = 1;
	int i = 1;
	int j = 1;
	double cpu = f_function(p->a, p->b, p->lb, p->rb, p->bb, p->ub, p->tau,
	                        currentTimeLevel, i, p->x, p->x_size, j + 2, p->y, p->y_size);
	double gpu = h_f_function(p, currentTimeLevel, i + 1, j);
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