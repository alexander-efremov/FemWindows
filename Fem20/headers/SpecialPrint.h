extern bool printVectorBy3D(
	const char *fileName,
	//
	int n,
	int numOfCurrTimeIter,
	int numOfSolOrd,
	//
	double numOfTimeSteps,
	//
	double *masOx,
	double *VectorOfData,
	//
	int dimOfVect );

extern bool printSurface_asV(
	const char *fileName, //   -  char *fileName,	
	int iCurrTSt, //   -  Index of current time step (layer) IN WHICH we printing solution.
	int numOfSolOrd, //   -  Solution order which we want to print.
	double *masOX, //   -  Massive of abscissa grid points. Dimension = numOfOXSt +1.
	int numOfOXSt, //   -  Number of abscissa grid steps.
	//
	double *masOY, //   -  Massive of ordinate grid points. Dimension = numOfOYSt +1.
	int numOfOYSt, //   -  Number of ordinate grid steps.
	//
	double *matData_asV,
	double *matData_asV2 ) ; //   -  Matrix of data.

extern bool printSurface_asV(
	const char *fileName, //   -  char *fileName,
	int iCurrTSt, //   -  Index of current time step (layer) IN WHICH we printing solution.
	int numOfSolOrd, //   -  Solution order which we want to print.	
	double *masOX, //   -  Massive of abscissa grid points. Dimension = numOfOXSt +1.
	int numOfOXSt, //   -  Number of abscissa grid steps.
	//
	double *masOY, //   -  Massive of ordinate grid points. Dimension = numOfOYSt +1.
	int numOfOYSt, //   -  Number of ordinate grid steps.
	//
	double *matData_asV ); //   -  Matrix of data.