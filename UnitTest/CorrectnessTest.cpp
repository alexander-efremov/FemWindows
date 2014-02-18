#include "CorrectnessTest.h"
#include "Timer.h"

class FemTest: public testing::Test
{
protected:
	double _accuracy;

	ComputeParameters GetComputeParameters(int level)
	{
		ComputeParameters result;

		result.a = C_par_a;
		result.b = C_par_b;
		result.lb = C_lbDom;
		result.rb = C_rbDom;
		result.bb = C_bbDom;
		result.ub = C_ubDom;
		double value = pow(2., level);
		result.tau = C_tau / value;
		result.t_count = C_numOfTSt * value;

		int n = C_numOfOXSt * value;

		double* arr_x = new double[n + 1];
		double buf_D = (C_rbDom - C_lbDom) / n;
		for (int j = 0; j < n + 1; j++)
		{
			arr_x[j] = C_lbDom + ((double) j) * buf_D;
		}

		result.x = arr_x;
		result.x_size = n;

		int m = C_numOfOYSt * value;
		buf_D = (C_ubDom - C_bbDom) / m;
		double* arr_y = new double[m + 1];
		for (int k = 0; k < m + 1; k++)
		{
			arr_y[k] = C_bbDom + ((double) k) * buf_D;
		}

		result.y = arr_y;
		result.y_size = m;
		result.size = (n + 1) * (m + 1);
		result.result = new double[result.size];
		return result;
	}

	ModelDataProvider _modelDataProvider;

	FemTest()
	{
		_accuracy = 1.0e-8;
		initCompOfGlVar();
		_modelDataProvider = ModelDataProvider();
	}

	virtual ~FemTest()
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

class CpuFemTest: public FemTest
{
protected:

	CpuFemTest()
	{
	}

	virtual ~CpuFemTest()
	{
	}

	double* GetCpuToLevel(int level)
	{
		return solve_cpu_test(C_par_a, C_par_b, C_lbDom, C_rbDom, C_bbDom,
			C_ubDom, C_tau, C_numOfTSt, masOX, C_numOfOXSt, masOY,
			C_numOfOYSt, level);
	}
};

TEST_F(CpuFemTest, CpuTestModel11)
{
	double* data = _modelDataProvider.GetModelData(Model11);
	double* result = GetCpuToLevel(0);

	for (int i = 0; i < GetSize(); i++)
	{
		EXPECT_TRUE(data[i] - result[i] <= _accuracy);
	}
}

TEST_F(CpuFemTest, CpuTestModel21)
{
	double* data = _modelDataProvider.GetModelData(Model21);
	double* result = GetCpuToLevel(1);

	for (int i = 0; i < GetSize(); i++)
	{
		EXPECT_TRUE(data[i] - result[i] <= _accuracy);
	}
}

TEST_F(CpuFemTest, CpuTestModel41)
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

TEST_F(CpuFemTest, CpuTestModel81)
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

TEST_F(CpuFemTest, CpuTestModel161)
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

TEST_F(CpuFemTest, CpuTestModel321)
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

TEST_F(CpuFemTest, CpuTestModel641)
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

TEST_F(CpuFemTest, CpuTestModel1281)
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

class GpuFemTest: public FemTest
{
protected:


	GpuFemTest()
	{
	}

	virtual ~GpuFemTest()
	{
	}

	double* GetGpuToLevel(int level)
	{
		ComputeParameters p = GetComputeParameters(level);
		//solve_cuda_params(p);
		d_solByEqualVolumes(p);
		return p.result;
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

	void TestTriangleTypeInternal(int startLevel, int finishLevel,
		int gridSize, int blockSize, int needPrint)
	{
		for (int level = startLevel; level < finishLevel; ++level)
		{
			cout << "level = " << level << std::endl;
			ComputeParameters p = GetComputeParameters(level);
			p.currentTimeLevel = 1;

			cout << "chunk = " << p.get_chunk_size() << std::endl;
			TriangleResult gpu = get_triangle_type(p, gridSize,
				blockSize);

			// get cpu data
			Triangle f = Triangle();
			Triangle s = Triangle();

			int y_size;
			int x_size;

			y_size = p.y_size;
			x_size = p.x_size;

			for (int j = 1; j < y_size; j++)
			{
				for (int i = 1; i < x_size; i++)
				{
					p.i = i;
					p.j = j;

					int result = h_quadrAngleType(p, f, s);
					int c = (p.x_size - 1) * (j - 1) + (i - 1);

					if (needPrint == 1)
					{
						printf("\n i = %d, j = %d, c is = %d\n", i, j, c);
						cout << "First GPU is " << std::endl << gpu.f[c]
						<< std::endl;
						cout << "Second GPU is " << std::endl << gpu.s[c]
						<< std::endl;
						cout << "First CPU is " << std::endl << f
							<< std::endl;
						cout << "Second CPU is " << std::endl << s
							<< std::endl;
					}

					ASSERT_TRUE(gpu.f[c] == f);
					ASSERT_TRUE(gpu.s[c] == s);
				}
			}
		}
	}
};

TEST_F(GpuFemTest, DISABLED_GpuTestModel11)
{
	//TEST_F(GpuFemTest, GpuTestModel11) {
	double* data = _modelDataProvider.GetModelData(Model11);
	double* result = GetGpuToLevel(0);

	for (int i = 0; i < GetSize(); i++)
	{
		EXPECT_TRUE(data[i] - result[i] <= _accuracy);
	}
}

TEST_F(GpuFemTest, Gpu_h_quadrAngleType)
{

	int finishLevel = 1;
	int startLevel = 0;

	int gridSize = 512;
	int blockSize = 1024;

	int needPrint = 0;
	//needPrint = 1;
	TestTriangleTypeInternal(startLevel, finishLevel, gridSize,
		blockSize, needPrint);
}

TEST_F(GpuFemTest, Gpu_h_quadrAngleType_time_measurment)
{

	int gridSize = 512;
	int blockSize = 1024;
	int level = 0;

	ComputeParameters p = GetComputeParameters(level);
	p.currentTimeLevel = 1;
	int y_size;
	int x_size;
	y_size = p.y_size;
	x_size = p.x_size;

	// get cpu data
	Triangle f = Triangle();
	Triangle s = Triangle();
	double time_cpu = -1;
	double time_gpu = -1;

	printf("Start CPU\n");
	StartTimer();
	cout << p << std::endl;

	for (int t = 1; t <= p.t_count; t++)
	{
		for (int j = 1; j < y_size; j++)
		{
			for (int i = 1; i < x_size; i++)
			{
				p.i = i;
				p.j = j;
				p.currentTimeLevel = t;
				int result = h_quadrAngleType(p, f, s);
				//std::cout  << "TL= " << t << " iteration = " <<  j*y_size + i<< "\n" << flush;
			}
		}
	}

	time_cpu = GetTimer();

	printf("End CPU\n");

	printf("Start GPU\n");

	StartTimer();

	for (int t = 1; t <= p.t_count; t++)
	{
		p.currentTimeLevel = t;
		TriangleResult gpu = get_triangle_type(p, gridSize,
			blockSize);
	}
	time_gpu = GetTimer();
	printf("End GPU\n");

	printf("CPU time is = %f\n", time_cpu);
	printf("GPU time is = %f\n", time_gpu);
	printf("CPU/GPU = %f\n", time_cpu/time_gpu);
}

//основной тест
//TEST_F(GpuFemTest, DISABLED_d_integUnderUnunifTr) {
TEST_F(GpuFemTest, debug_multi_chunks)
{
	int start_level =0;
	int finish_level = 1;

	// параметры сетки для GPU
	int grid_size = 512;
	int block_size = 1024;
	for (int level = start_level; level < finish_level; level++)
	{
		ComputeParameters p = GetComputeParameters(level);
		cout << p << std::endl;
		p.print_info();
		// get gpu result of integr
		printf("Start GPU triangulation...  ");
		TriangleResult gpu = get_triangle_type(p, grid_size, block_size);
		printf("GPU triangulation is done.\n");
	}
}


//основной тест
TEST_F(GpuFemTest, DISABLED_d_integUnderUnunifTr)
	//TEST_F(GpuFemTest, d_integUnderUnunifTr)
{
	// параметры размерности задачи
	int levelStep = 1;
	//int start_level =7;
	int start_level =0;
	//int start_level =9;
	int finish_level = 1;

	// параметры сетки для GPU
	int grid_size = 512;
	int block_size = 1024;

	for (int level = start_level; level < finish_level; level++)
	{
		ComputeParameters p = GetComputeParameters(level);
		cout << p << std::endl;

		for(p.currentTimeLevel = 1; p.currentTimeLevel < 2; p.currentTimeLevel += levelStep)
		{
			// init test data
			p.print_info();
			double time_cpu = -1;
			double time_gpu = -1;
			double* rho_in_prev_tl = new double [ p.size ];

			//   Initial data of rho.
			for( int k = 0; k <= p.x_size; k++ )
			{
				for( int j = 0; j <= p.y_size; j++ )
				{
					rho_in_prev_tl[ (p.x_size+1)*k + j ] = d_initDataOfSol(p, j, k);
				}
			}

			// ==


			// get gpu result of integr
			printf("Start GPU triangulation...  ");

			TriangleResult gpu = get_triangle_type(p, grid_size,
				block_size);
			printf("GPU triangulation is done.\n");

			int c = -1;
			double* cpu_result = new double[p.x_size * p.y_size];

			printf("Perform CPU integration...   ");
			StartTimer();

			for (int j = 1; j < p.y_size; j++)
			{
				for (int i = 1; i < p.x_size; i++)
				{
					p.i = i;
					p.j = j;
					c = (p.x_size - 1) * (j - 1) + (i - 1);

					//printf("\n i = %d, j = %d, c is = %d\n", i, j, c);

					cpu_result[c] = get_cpu_integr_result_internal(p, rho_in_prev_tl, gpu.f[c], gpu.s[c]);
				}
			}

			time_cpu = GetTimer();
			printf("CPU integration is done.\n");

			printf("Perform GPU integration...   ");
			StartTimer();
			double* gpu_integr = GetIntegrGpuResult(p, gpu, rho_in_prev_tl, grid_size, block_size
				, p.get_chunk_size());
			time_gpu = GetTimer();
			printf("GPU integration is done.\n");


			printf("CPU time is = %f ms.\n", time_cpu);
			printf("GPU time is = %f ms.\n", time_gpu);
			printf("CPU/GPU = %f\n times\n", time_cpu/time_gpu);


			// проверям корректность полученных результатов
			for (int j = 1; j < p.y_size; j++)
			{
				for (int i = 1; i < p.x_size; i++)
				{
					c = (p.x_size - 1) * (j - 1) + (i - 1);

					ASSERT_DOUBLE_EQ(cpu_result[c], gpu_integr[c]);
				}
			}

			delete[] cpu_result;
			delete[] gpu_integr;
			delete[] rho_in_prev_tl;
		}
	}
}


class CpuVersusGpuFunctionalFemTest: public FemTest
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
	ComputeParameters p = GetComputeParameters(0);
	double cpu = analytSolut(p.a, p.lb, p.rb, p.bb, p.ub, p.tau, p.x[1],
		p.y[1]);
	double gpu = h_analytSolut(p.tau, p.x[1], p.y[1]);
	ASSERT_EQ(cpu, gpu);
}

TEST_F(CpuVersusGpuFunctionalFemTest, AnalyticSolutionNotEqualsTest)
{
	ComputeParameters p = GetComputeParameters(0);
	double cpu = analytSolut(p.a, p.lb, p.rb, p.bb, p.ub, p.tau, p.x[1],
		p.y[2]);
	double gpu = h_analytSolut(p.tau, p.x[1], p.y[1]);
	ASSERT_NE(cpu, gpu);
}

TEST_F(CpuVersusGpuFunctionalFemTest, DataInitializationAnalyticSolutionEqualsTest)
{
	ComputeParameters p = GetComputeParameters(0);
	double cpu = initDataOfSol(p.a, p.lb, p.rb, p.bb, p.ub, 1, p.x, 2, p.y);
	double gpu = d_initDataOfSol(p, 1, 2);
	ASSERT_EQ(cpu, gpu);
}

TEST_F(CpuVersusGpuFunctionalFemTest, F_Function_EqualsTest)
{
	ComputeParameters p = GetComputeParameters(0);
	int currentTimeLevel = 1;
	int i = 1;
	int j = 1;
	double cpu = f_function(p.a, p.b, p.lb, p.rb, p.bb, p.ub, p.tau,
		currentTimeLevel, i, p.x, p.x_size, j, p.y, p.y_size);
	double gpu = h_f_function(p, currentTimeLevel, i, j);
	ASSERT_EQ(cpu, gpu);
}

TEST_F(CpuVersusGpuFunctionalFemTest, F_Function_NotEqualsTest)
{
	ComputeParameters p = GetComputeParameters(0);
	int currentTimeLevel = 1;
	int i = 1;
	int j = 1;
	double cpu = f_function(p.a, p.b, p.lb, p.rb, p.bb, p.ub, p.tau,
		currentTimeLevel, i, p.x, p.x_size, j + 2, p.y, p.y_size);
	double gpu = h_f_function(p, currentTimeLevel, i + 1, j);
	ASSERT_NE(cpu, gpu);
}

TEST_F(CpuVersusGpuFunctionalFemTest, d_bottomBound_EqualsTest)
{
	ComputeParameters p = GetComputeParameters(0);
	int currentTimeLevel = 1;
	double t = p.tau * 1;
	int i = 1;

	double cpu = bottonBound(p.a, p.lb, p.rb, p.bb, p.ub, t, p.x[i]);

	p.currentTimeLevel = currentTimeLevel;
	p.i = i;
	double gpu = h_bottomBound(p);

	ASSERT_EQ(cpu, gpu);
}

TEST_F(CpuVersusGpuFunctionalFemTest, d_upperBound_EqualsTest)
{
	ComputeParameters p = GetComputeParameters(0);
	int currentTimeLevel = 1;
	double t = p.tau * 1;
	int i = 1;

	double cpu = upperBound(p.a, p.lb, p.rb, p.bb, p.ub, t, p.x[i]);
	p.currentTimeLevel = currentTimeLevel;
	p.i = i;
	double gpu = h_upperBound(p);
	ASSERT_EQ(cpu, gpu);
}

TEST_F(CpuVersusGpuFunctionalFemTest, d_rightBound_EqualsTest)
{
	ComputeParameters p = GetComputeParameters(0);
	int currentTimeLevel = 1;
	double t = p.tau * 1;
	int j = 1;

	double cpu = rightBound(p.a, p.lb, p.rb, p.bb, p.ub, t, p.y[j]);

	p.currentTimeLevel = currentTimeLevel;
	p.j = j;
	double gpu = h_rightBound(p);
	ASSERT_EQ(cpu, gpu);
}

TEST_F(CpuVersusGpuFunctionalFemTest, d_leftBound_EqualsTest)
{
	ComputeParameters p = GetComputeParameters(0);
	int currentTimeLevel = 1;
	double t = p.tau * 1;
	int j = 1;

	double cpu = leftBound(p.a, p.lb, p.rb, p.bb, p.ub, t, p.y[j]);

	p.currentTimeLevel = currentTimeLevel;
	p.j = j;
	double gpu = h_leftBound(p);

	ASSERT_EQ(cpu, gpu);
}