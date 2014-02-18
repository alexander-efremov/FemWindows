//extern int h_quadrAngleType(ComputeParameters& p, Triangle& firstT,	Triangle& secondT);
extern double solByEqualVolumes(
	double par_a, //   -  Item of left and right setback (parameter "a" in test).
	double par_b, //   -  Item of second parameter from "u_funcion".
	//
	double lbDom, //   -  Left and right boundaries of rectangular domain.
	double rbDom,
	//
	double bbDom, //   -  Bottom and upper boundaries of rectangular domain.
	double ubDom,
	//
	double tau, //   -  Time step.
	int numOfTSt, //   -  A number of time steps.
	//
	double *masOX, //   -  Massive of OX points. Dimension = numOfOXSt +1.
	int numOfOXSt, //   -  Number of OX steps.
	//
	double *masOY, //   -  Massive of OY points. Dimension = numOfOYSt +1.
	int numOfOYSt, //   -  Number of OY steps.
	//
	int numOfSolOrd, //   -  For print only. Solution order which we want to get.
	//
	double *rhoInCurrTL_asV ); //   -  Rho (solution) in Last Time Level which we will compute.

extern double* solve_cpu_test(
	double par_a, //   -  Item of left and right setback (parameter "a" in test).
	double par_b, //   -  Item of second parameter from "u_funcion".
	//
	double lbDom, //   -  Left and right boundaries of rectangular domain.
	double rbDom,
	//
	double bbDom, //   -  Botton and upper boundaries of rectangular domain.
	double ubDom,
	//
	double tau, //   -  Time step.
	int numOfTSt, //   -  A number of time steps.
	//
	double * masOX, //   -  Massive of OX nodes. Dimension = numOfOXSt +1.
	int numOfOXSt, //   -  Number of OX steps.
	//
	double * masOY, //   -  Massive of OY nodes. Dimension = numOfOYSt +1.
	int numOfOYSt, //   -  Number of OY steps.
	int gridStep
);

extern double spaceVolumeInPrevTL(
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
	double iCurrTL, //   -  Index of current time layer. Necessary for velocity.
	//
	int iOfOXN, //   -  Index of current OX node.
	const double * masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
	int numOfOXSt, //   -  Number of OX steps.
	//
	int iOfOYN, //   -  Index of current OY node.
	const double * masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
	int numOfOYSt, //   -  Number of OY steps.	

	double * rhoInPrevTL_asV );

extern double f_function( //   -  It's item of right part of differential equation.
	double par_a, //   -  Item of left and right setback (parameter "a" in test).
	double par_b, //   -  Item of second parameter from "u_funcion".
	//
	double lbDom, //   -  Left and right boundaries of rectangular domain.
	double rbDom,
	//
	double bbDom, //   -  Bottom and upper boundaries of rectangular domain.
	double ubDom,
	//
	double tau,
	int iCurrTL, //   -  Index of current time layer.
	//
	int iOfOXN, //   -  Index of current OX node.
	const double *masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
	int numOfOXSt, //   -  Number of OX steps (segments).
	//
	int iOfOYN, //   -  Index of current OY node.
	const double *masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
	int numOfOYSt ) ; //   -  Number of OY steps (segments).

extern double rightBound(
	double par_a, //   -  Item of left and right setback (parameter "a" in test).
	//
	double lbDom, //   -  Left and right boundaries of rectangular domain.
	double rbDom,
	//
	double bbDom, //   -  Bottom and upper boundaries of rectangular domain.
	double ubDom,
	//
	double t,
	double y );

extern double bottonBound(
	double par_a, //   -  Item of left and right setback (parameter "a" in test).
	//
	double lbDom, //   -  Left and right boundaries of rectangular domain.
	double rbDom,
	//
	double bbDom, //   -  Botton and upper boundaries of rectangular domain.
	double ubDom,
	//
	double t,
	double x );

extern double upperBound(
	double par_a, //   -  Item of left and right setback (parameter "a" in test).
	//
	double lbDom, //   -  Left and right boundaries of rectangular domain.
	double rbDom,
	//
	double bbDom, //   -  Botton and upper boundaries of rectangular domain.
	double ubDom,
	//
	double t,
	double x );

extern double leftBound(
	double par_a, //   -  Item of left and right setback (parameter "a" in test).
	//
	double lbDom, //   -  Left and right boundaries of rectangular domain.
	double rbDom,
	//
	double bbDom, //   -  Botton and upper boundaries of rectangular domain.
	double ubDom,
	//
	double t,
	double y );

extern double integUnderUnunifTr(
	double par_a, //   -  Item of left and right setback (parameter "a" in test).
	double par_b, //   -  Item of second parameter from "u_funcion".
	//
	double lbDom, //   -  Left and right boundaries of rectangular domain.
	double rbDom,
	//
	double bbDom, //   -  Bottom and upper boundaries of rectangular domain.
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
	double * rhoInPrevTL );

extern int quadrAngleType(
	double par_a, //   -  Item of left and right setback (parameter "a" in test).
	double par_b, //   -  Item of second parameter from "u_funcion".
	//
	double lbDom, //   -  Left and right boundaries of rectangular domain.
	double rbDom,
	//
	double bbDom, //   -  bottom and upper boundaries of rectangular domain.
	double ubDom,
	//
	double tau,
	double iCurrTL, //   -  Index of current time layer. Necessary for velocity.
	//
	int iOfOXN, //   -  Index of current OX node.
	const double * masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
	int numOfOXSt, //   -  Number of OX steps.
	//
	int iOfOYN, //   -  Index of current OY node.
	const double * masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
	int numOfOYSt, //   -  Number of OY steps.
	//
	double * firVfirT, //   -  First vertex of first triangle.
	double * secVfirT, //   -  Second vertex of first triangle.
	double * thiVfirT, //   -  Third vertex of first triangle.
	//
	double * firVsecT, //   -  First vertex of second triangle.
	double * secVsecT, //   -  Second vertex of second triangle.
	double * thiVsecT ); //   -  Third vertex of second triangle.

extern double integUnderUnunifTr(
	double par_a, //   -  Item of left and right setback (parameter "a" in test).
	double par_b, //   -  Item of second parameter from "u_funcion".
	//
	double lbDom, //   -  Left and right boundaries of rectangular domain.
	double rbDom,
	//
	double bbDom, //   -  Bottom and upper boundaries of rectangular domain.
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
	double * rhoInPrevTL );

extern
	double itemOfInteg_2SpecType_optimized(
		double Py,
		double Qy,
		//
		double alpha,
		//
		double a,
		double b,
		double betta );

extern double analytSolut(
	double par_a,
	//
	double lbDom, //   -  Left and right boundaries of rectangular domain.
	double rbDom,
	//
	double bbDom, //   -  bottom and upper boundaries of rectangular domain.
	double ubDom,
	//
	double t, double x, double y );


extern void cuda_solve(
	double par_a, //   -  Item of left and right setback (parameter "a" in test).
	double par_b, //   -  Item of second parameter from "u_funcion".
	//
	double lbDom, //   -  Left and right boundaries of rectangular domain.
	double rbDom,
	//
	double bbDom, //   -  bottom and upper boundaries of rectangular domain.
	double ubDom,
	//
	double tau, //   -  Time step.
	int numOfTSt, //   -  A number of time steps.
	//
	double * masOX, //   -  Massive of OX nodes. Dimension = numOfOXSt +1.
	int numOfOXSt, //   -  Number of OX steps.
	//
	double * masOY, //   -  Massive of OY nodes. Dimension = numOfOYSt +1.
	int numOfOYSt, //   -  Number of OY steps.
	//
	bool isTimeStShBeChan, //   -  Is time step should be change?
	bool isGridStShBeChan, //   -  Is grid step should be change?
	//
	int numOfGrStepLayer ); //   -  How many computations with different grid steps we want to make.

extern double initDataOfSol(
	double par_a, //   -  Item of left and right setback (parameter "a" in test).
	//
	double lbDom, //   -  Left and right boundaries of rectangular domain.
	double rbDom,
	//
	double bbDom, //   -  Botton and upper boundaries of rectangular domain.
	double ubDom,
	//
	int iOfOXN, //   -  Index of current OX node.
	const double *masOX, //   -  Massive of abscissa grid steps. Dimension = numOfOxSt +1.
	//
	int iOfOYN, //   -  Index of current OY node.
	const double *masOY ); //   -  Massive of ordinate grid steps. Dimension = numOfOySt +1.

extern double normOfMatrAtL1_asV(
	const double *masOX, //   -  Massive of OX grid nodes. Dimension = dimOX.
	int dimOX,
	//
	const double *masOY, //   -  Massive of OY grid nodes. Dimension = dimOY.
	int dimOY,
	//
	double * mat_asV );

extern double u_function(double par_b, double t, double x, double y);
extern double v_function(double lbDom, double rbDom, double bbDom, double ubDom, double t, double x, double y );