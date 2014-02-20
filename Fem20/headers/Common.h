#include <iostream>
#ifndef COMMON_H_INCLUDED
#define COMMON_H_INCLUDED

#define C_pi_device 3.14159265358979323846264338327


template<typename CharT>
class DecimalSeparator : public std::numpunct<CharT>
{
public:

	DecimalSeparator(CharT Separator)
		: m_Separator(Separator)
	{
	}

protected:

	CharT do_decimal_point()const
	{
		return m_Separator;
	}

private:
	CharT m_Separator;
};

struct ComputeResults
{
	int timeStepCount;
	int xSize;
	int ySize;
	double cpuTime;
	double openmpTime;
	double gpuTime;
};

struct ComputeParameters
{
private:
	int time_i;
	int get_real_x_size()
	{
		return x_size + 1;
	}


	int get_real_y_size()
	{
		return y_size + 1;
	}


public:
	double a;
	double b;
	double lb;
	double rb;
	double bb;
	double ub;
	double tau;

	int size;
	int t_count;
	int x_size;
	int y_size;
	int i;
	int j;
	double *x;
	double *y;
	double *result;
	double currentTimeLevel;

	ComputeParameters() : currentTimeLevel(1), t_count(0)
	{
	}

	int get_inner_chuck_size()
	{
		return get_chunk_size() - 2;
	}

	void reset_time_counter()
	{
		time_i = 1;
	}

	bool can_iterate_over_time_level()
	{
		return time_i <= t_count;
	}

	void inc_time_level()
	{
		time_i++;
	}

	int get_chunk_size()
	{
		int size = get_real_x_size();
		if (size == 11)
		{
			return size;
		}

		if (size == 21)
		{
			return size;
		}
		else if (size == 41)
		{
			return size;
		}
		else if (size == 81)
		{
			return size;
		}
		else if (size == 161)
		{
			return size;
		}
		else if (size == 321)
		{
			return size;
		}
		else if (size == 641)
		{
			return size;
		}
		else if (size == 1281)
		{
			return size;
		}
		else if (size == 2561)
		{
			return size;
		}
		else if (size == 5121)
		{
			return size;
			//return size*569; // 5121 / 9 = 569 - столько ядер запустится для расчета триангуляции
		}
		else if (size == 10242)
		{
			return size;
			//return size*569; // 5121 / 9 = 569 - столько ядер запустится для расчета триангуляции
		}
		return 0;
	}

	// получает размер внутренней матрицы
	int get_inner_matrix_size()
	{
		return (get_real_x_size() - 2) * (get_real_y_size() - 2);
	}

	// получает размер внутренней матрицы
	int get_inner_x_size()
	{
		return get_real_x_size() - 2;
	}


	int get_part_count()
	{
		int chunk_size = get_inner_chuck_size();
		int inner_matrix_size = get_inner_matrix_size();
		return (int) ( inner_matrix_size % chunk_size == 0 ? inner_matrix_size / chunk_size : inner_matrix_size / chunk_size + 1);
	}

	void print_info()
	{
		std::cout << "chunk = " << get_chunk_size() << " current time level " << currentTimeLevel << std::endl;
	}


	friend std::ostream &operator<<( std::ostream &output,
	                                 const ComputeParameters &tr )
	{
		output << "a = " << tr.a << std::endl;
		output << "b = " << tr.b << std::endl;
		output << "lb = " << tr.lb << std::endl;
		output << "rb = " << tr.rb << std::endl;
		output << "bb = " << tr.bb << std::endl;
		output << "ub = " << tr.ub << std::endl;
		output << "size = " << tr.size << std::endl;
		output << "tau = " << tr.tau << std::endl;
		output << "Time levels = " << tr.t_count << std::endl;
		output << "x size = " << tr.x_size << std::endl;
		output << "y size = " << tr.y_size << std::endl;
		return output;
	}
};

struct Triangle
{
	double first[2];
	double second[2];
	double third[2];

	friend std::ostream &operator<<( std::ostream &output,
	                                 const Triangle &tr )
	{
		output << "First Point: ";
		output << "x: " << tr.first[0] << " y: " << tr.first[1] << std::endl;
		output << "Second Point: ";
		output << "x: " << tr.second[0] << " y: " << tr.second[1] << std::endl;
		output << "Third Point: ";
		output << "x: " << tr.third[0] << " y: " << tr.third[1] << std::endl;
		return output;
	}

	friend bool operator==(const Triangle& x, const Triangle& y)
	{
		double const error = 10e-16;
		bool p1 = (x.first[0] - y.first[0] < error) && (x.first[1] - y.first[1] < error);
		bool p2 = (x.second[0] - y.second[0] < error) && (x.second[1] - y.second[1] < error);
		bool p3 = (x.third[0] - y.third[0] < error) && (x.third[1] - y.third[1] < error);
		return p1 && p2 && p3;
	}
};


struct TriangleResult
{
	
private:
	int chunk;

public:
	Triangle* f;
	Triangle* s;

	int length;
	int x_length;
	int y_length;
	int currentTimeLevel;
	int offset;
	TriangleResult(ComputeParameters param)
	{
		currentTimeLevel =param.currentTimeLevel;
		x_length = param.x_size - 1;
		y_length = param.y_size - 1;
		chunk = param.get_inner_chuck_size();
		f = new Triangle[param.get_inner_matrix_size()];
		s = new Triangle[param.get_inner_matrix_size()];
	}

	~TriangleResult()
	{
		delete[] f;
		delete[] s;
	}

	void setOffset(int part_index)
	{
		offset = part_index * chunk;
	}
};

struct VertexData
{
	Triangle* f;
	Triangle* s;
	int* types;
	int length;
	int x_length;
	int y_length;
	int currentTimeLevel;
	int part_number;
	int offset;
	int chunk;
};

void print_vector(const int, const double*);

const std::string currentDateTime();

void print_matrix(const int, const int, const double*);

void print_matrix(const int, const int, const int*);

void write_csv(double, double, double, const int, const int, const double, const std::string&, const std::string&);

void write_info(double, double, double, const int, const int, double);

void write_openmp_stress_test_info(std::string &filename, int threadNumber, int n, int m, int tl_number, double time_ms, int ht_on, bool append);

#endif