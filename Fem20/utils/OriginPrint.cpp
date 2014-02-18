#include <iostream>
#include <stdio.h>
#include <stdlib.h>

void print_tecplot_3d (
	char *filename,
	double *mas_0x,
	int len_x,
	double *mas_0y,
	int len_y,
	double *mas_0z )
{
	FILE *pfile;

	pfile = fopen( filename, "w");
	if( pfile == NULL )
	{
		std::cout << "\nError can not create-open file " << filename;
		exit(1);
	}


	fprintf(pfile, "TITLE = 'DEM DATA | DEM DATA | DEM DATA | DEM DATA'" );
	fprintf(pfile, "\nVARIABLES = 'X' 'Y' 'E' " );
	fprintf(pfile, "\nZONE T='SubZone'" );
	fprintf(pfile, "\nI=%d J=%d K=%d ZONETYPE=Ordered", len_y, len_x, 1 );
	fprintf(pfile, "\nDATAPACKING=POINT" );
	fprintf(pfile, "\nDT=(SINGLE SINGLE SINGLE )" );

	for (int i=0; i<len_x; i++ )
		for (int j=0; j<len_y; j++ )
			fprintf(pfile, "\n%-30.20g  %-30.20g %-30.20g", mas_0x[i], mas_0y[j], mas_0z[ j*len_x + i ] );

	fclose(pfile);
}

void printSurfaceByMatrix_asV(
	char *filename,
	//
	double *masOX, //   -  Massive of abscissa grid points. Dimension = numOfOXSt +1.
	int numOfOXSt, //   -  Number of abscissa grid steps.
	//
	double *masOY, //   -  Massive of ordinate grid points. Dimension = numOfOYSt +1.
	int numOfOYSt, //   -  Number of ordinate grid steps.
	//
	double *matData_asV ) //   -  Matrix of data.
{
	FILE *pfile;

	pfile = fopen( filename, "w");

	if( pfile == NULL )
	{
		std::cout << "\nError cannot create/open file " << filename;
		exit(1);
	}

	fprintf(pfile, "TITLE = 'DEM DATA | DEM DATA | DEM DATA | DEM DATA'" );
	fprintf(pfile, "\nVARIABLES = 'X' 'Y' 'condensation' " );
	fprintf(pfile, "\nZONE T='SubZone'" );
	fprintf(pfile, "\nI=%d J=%d K=%d ZONETYPE=Ordered", numOfOYSt, numOfOXSt, 1 );
	fprintf(pfile, "\nDATAPACKING=POINT" );
	fprintf(pfile, "\nDT=(SINGLE SINGLE SINGLE )" );

	for ( int i=0; i< numOfOXSt; i++ )
		for ( int j=0; j< numOfOYSt; j++ )
			fprintf(pfile, "\n%-30.20g  %-30.20g %-30.20g", masOX[i], masOY[j], matData_asV[ (numOfOXSt +1)*j + i ] );

	fclose(pfile);
}

void printSurfaceByMatrix_asV_velocity(
	char *filename,
	//
	double *masOX, //   -  Massive of abscissa grid points. Dimension = numOfOXSt +1.
	int numOfOXSt, //   -  Number of abscissa grid steps.
	//
	double *masOY, //   -  Massive of ordinate grid points. Dimension = numOfOYSt +1.
	int numOfOYSt, //   -  Number of ordinate grid steps.
	//
	double *matData_asV, double *matData_asV2) //   -  Matrix of data.
{
	FILE *pfile;

	pfile = fopen( filename, "w");

	if( pfile == NULL )
	{
		std::cout << "\nError cannot create/open file " << filename;
		exit(1);
	}

	fprintf(pfile, "TITLE = 'DEM DATA | DEM DATA | DEM DATA | DEM DATA'" );
	fprintf(pfile, "\nVARIABLES = 'X' 'Y' 'condensation' " );
	fprintf(pfile, "\nZONE T='SubZone'" );
	fprintf(pfile, "\nI=%d J=%d K=%d ZONETYPE=Ordered", numOfOYSt, numOfOXSt, 1 );
	fprintf(pfile, "\nDATAPACKING=POINT" );
	fprintf(pfile, "\nDT=(SINGLE SINGLE SINGLE )" );

	for ( int i=0; i< numOfOXSt; i++ )
		for ( int j=0; j< numOfOYSt; j++ )
			fprintf(pfile, "\n%-30.20g  %-30.20g %-30.20g", masOX[i], masOY[j], matData_asV[ (numOfOXSt +1)*j + i ] );

	fclose(pfile);
}