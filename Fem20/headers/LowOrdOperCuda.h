#include "Common.h"
extern double h_analytSolut(double t, double x, double y );

extern double h_f_function(ComputeParameters p, const int currentTimeLevel, const int i, const int j);

extern double h_rightBound(ComputeParameters& p);


extern double u_function(double par_b, double t, double x, double y);
extern double v_function(double lbDom, double rbDom, double bbDom, double ubDom, double t, double x, double y );

extern double h_leftBound(ComputeParameters& p);

extern double h_upperBound(ComputeParameters& p);

extern double h_bottomBound(ComputeParameters& p);

extern double d_initDataOfSol(ComputeParameters p, int i, int j);

extern float solve_cuda_params(ComputeParameters p);

extern double d_solByEqualVolumes(ComputeParameters p);

extern int h_quadrAngleType(ComputeParameters& p, Triangle& firstT, Triangle& secondT);

extern int d_quadrAngleType(ComputeParameters& p, Triangle& firstT, Triangle& secondT);

extern double d_integUnderUnunifTr(
	double par_a, //   -  Item of left and right setback (parameter "a" in test).
	double par_b, //   -  Item of second parameter from "u_funcion".
	//
	double lbDom, //   -  Left and right boundaries of rectangular domain.
	double rbDom,
	//
	double bbDom, //   -  Botton and upper boundaries of rectangular domain.
	double ubDom,
	//
	double tau,
	int iCurrTL, //   -  Index of current time layer.
	//
	double * firVer, //   -  First vertex of triangle.
	double * secVer, //   -  Second vertex of triangle.
	double * thiVer, //   -  Third vertex of triangle.
	//
	const double * masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
	int numOfOXSt, //   -  Number of OX steps.
	//
	const double * masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
	int numOfOYSt, //   -  Number of OY steps.
	//
	double * rhoInPrevTL_asV,
	int ii, int jj );

extern double integUnderUnunifTr(
	double par_a, //   -  Item of left and right setback (parameter "a" in test).
	double par_b, //   -  Item of second parameter from "u_funcion".
	//
	double lbDom, //   -  Left and right boundaries of rectangular domain.
	double rbDom,
	//
	double bbDom, //   -  Botton and upper boundaries of rectangular domain.
	double ubDom,
	//
	double tau,
	int iCurrTL, //   -  Index of current time layer.
	//
	double * firVer, //   -  First vertex of triangle.
	double * secVer, //   -  Second vertex of triangle.
	double * thiVer, //   -  Third vertex of triangle.
	//
	const double * masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
	int numOfOXSt, //   -  Number of OX steps.
	//
	const double * masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
	int numOfOYSt, //   -  Number of OY steps.
	//
	double * rhoInPrevTL_asV,
	int ii, int jj ); //!!!!!!!!!!!!!!!!!!!

double* GetIntegrGpuResult(ComputeParameters p, TriangleResult triangles, double *rhoInPrevTL_asV, int gridSize, int blockSize, int chunk);