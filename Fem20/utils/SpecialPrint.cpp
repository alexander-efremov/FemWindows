#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../headers/OriginPrint.h" //   -  Data printing to file in special form for software "TecPlot".


using namespace std;

bool printSurface_asV(
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
	double *matData_asV2 ) //   -  Matrix of data.
{
	char strOfName[50];
	char str1[4] = "Nx=", str2[7];
	char str3[3] = "k=", str4[7];
	char str5[5] = ".dat", str6[2];
	strcpy( strOfName, fileName );
	sprintf (str6, "%d", numOfSolOrd);
	strcat( strOfName, str6 );
	strcat( strOfName, str1 );
	sprintf (str2, "%d", numOfOXSt);
	strcat( strOfName, str2 );
	strcat( strOfName, str3 );
	sprintf (str4, "%d", iCurrTSt);
	strcat( strOfName, str4 );
	strcat( strOfName, str5 );
	//	  Data printing.
	printSurfaceByMatrix_asV(
		strOfName,
		masOX, //   -  Massive of abscissa grid points. Dimension = numOfOXSt +1.
		numOfOXSt, //   -  Number of abscissa grid steps.
		masOY, //   -  Massive of ordinate grid points. Dimension = numOfOYSt +1.
		numOfOYSt, //   -  Number of ordinate grid steps.
		matData_asV ); //   -  Matrix of data.
	return true;
}

bool printVectorBy3D(
	const char *fileName,
	//
	int n,
	int numOfCurrTimeIter,
	int numOfSolOrd, //   -  Solution order which we want to get.
	//
	double numOfTimeSteps,
	//
	double *masOx,
	double *VectorOfData,
	//
	int dimOfVect )
{
	char strOfName[50];
	char str1[3] = "N=", str2[7];
	char str3[3] = "k=", str4[7];
	char str5[5] = ".dat", str6[2];
	strcpy( strOfName, fileName );
	sprintf (str6, "%d", numOfSolOrd);
	strcat( strOfName, str6 );
	strcat( strOfName, str1 );
	sprintf (str2, "%d", n);
	strcat( strOfName, str2 );
	strcat( strOfName, str3 );
	sprintf (str4, "%d", numOfCurrTimeIter);
	strcat( strOfName, str4 );
	strcat( strOfName, str5 );
	double *mas_Empty = new double [ 1 ];
	mas_Empty[0] = ( (int)(numOfCurrTimeIter *1000. /numOfTimeSteps) ) /10;
	//	Âûâîä äàííûõ.
	print_tecplot_3d( strOfName,
	                  //
	                  masOx, dimOfVect,
	                  //
	                  mas_Empty, 1,
	                  VectorOfData
	);
	delete[] mas_Empty;
	return true;
}

bool printSurface_asV(
	const char *fileName, //   -  char *fileName,	
	int iCurrTSt, //   -  Index of current time step (layer) IN WHICH we printing solution.
	int numOfSolOrd, //   -  Solution order which we want to print.
	double *masOX, //   -  Massive of abscissa grid points. Dimension = numOfOXSt +1.
	int numOfOXSt, //   -  Number of abscissa grid steps.
	//
	double *masOY, //   -  Massive of ordinate grid points. Dimension = numOfOYSt +1.
	int numOfOYSt, //   -  Number of ordinate grid steps.
	//
	double *matData_asV ) //   -  Matrix of data.
{
	char strOfName[50];
	char str1[4] = "Nx=", str2[7];
	char str3[3] = "k=", str4[7];
	char str5[5] = ".dat", str6[2];
	strcpy( strOfName, fileName );
	sprintf (str6, "%d", numOfSolOrd);
	strcat( strOfName, str6 );
	strcat( strOfName, str1 );
	sprintf (str2, "%d", numOfOXSt);
	strcat( strOfName, str2 );
	strcat( strOfName, str3 );
	sprintf (str4, "%d", iCurrTSt);
	strcat( strOfName, str4 );
	strcat( strOfName, str5 );
	//	  Data printing.
	printSurfaceByMatrix_asV(
		strOfName,
		masOX, //   -  Massive of abscissa grid points. Dimension = numOfOXSt +1.
		numOfOXSt, //   -  Number of abscissa grid steps.
		masOY, //   -  Massive of ordinate grid points. Dimension = numOfOYSt +1.
		numOfOYSt, //   -  Number of ordinate grid steps.
		matData_asV ); //   -  Matrix of data.
	return true;
}