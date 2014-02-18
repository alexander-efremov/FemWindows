#include <iostream>

using namespace std;

#include <math.h>


const double C_pi = 3.14159265358979323846264338327;


const double C_par_a = 2.; //   -  Solution parameter.

const double C_par_b = 100.; //   -  Item of second parameter from "u_funcion" or "v_funcion".


double C_lbDom = 0., C_rbDom = 1.; //   -  Left and right boundaries of rectangular domain.

const int C_numOfOXSt = 10; //   -  Number of OX steps (segments).

double *masOX; //   -  Massive of OX points. Dimension = C_numOfOXSt +1.


double C_bbDom = 0., C_ubDom = 1.; //   -  Botton and upper boundaries of rectangular domain.

const int C_numOfOYSt = 10; //   -  Number of OY steps (segments).

double *masOY; //   -  Massive of OY points. Dimension = C_numOfOYSt +1.


const double C_StepsRel = 1./5.; //   -  A relation of the time step "C_tau" to the max grid step "maxOfGrSt";

double C_tau; //   -  Time step. Will be computed by relation "C_StepsRel" with small correction.

const double C_timeEnd = 1.; //   -  Finish time.

int C_numOfTSt; //   -  A number of time steps.


void initCompOfGlVar()
{
	double maxOfGrSt = 0.;
	int j;

	//   Massive of OX steps. Dimension = C_numOfOXSt +1.
	masOX = new double[ C_numOfOXSt +1];
	for( j=0; j< C_numOfOXSt +1; j++ )
	{
		//   Regular grid for simplification.
		masOX[j] = C_lbDom + (C_rbDom - C_lbDom) * ( (double)j / C_numOfOXSt );
	}

	//   Let's find out maximum of OX step.
	for( j=0; j< C_numOfOXSt; j++ )
	{
		if( maxOfGrSt < (masOX[j+1] - masOX[j]) )
		{
			maxOfGrSt = (masOX[j+1] - masOX[j]);
		}
	}

	//   Massive of OY steps. Dimension = C_numOfOYSt +1.
	masOY = new double[ C_numOfOYSt +1];
	for( j=0; j< C_numOfOYSt +1; j++ )
	{
		//   Regular grid for simplification.
		masOY[j] = C_bbDom + (C_ubDom - C_bbDom) * ( (double)j / C_numOfOYSt );
	}

	//   Let's find out maximum of OY step.
	for( j=0; j< C_numOfOYSt; j++ )
	{
		if( maxOfGrSt < (masOY[j+1] - masOY[j]) )
		{
			maxOfGrSt = (masOY[j+1] - masOY[j]);
		}
	}

	C_tau = C_StepsRel * maxOfGrSt;

	C_numOfTSt = (int)( C_timeEnd / C_tau );


	if( fabs(C_numOfTSt*C_tau - C_timeEnd) < (C_tau* 1.e-7) )
	{
		//   O.K. This error is error of computer arithmetic.

		return;
	}

	if ( fabs(C_numOfTSt*C_tau - C_timeEnd) >= (C_tau* 1.e-7) )
	{
		//   We need to change number of time steps and time step.

		//cout<<"\nA number of time steps and time step itself have changed."<< flush;

		C_numOfTSt++;

		C_tau = C_timeEnd / C_numOfTSt;

		//cout<<"\nNow time step = "<< C_tau << " and a number of such time steps = " << C_numOfTSt << flush;

		//cout<<"\nSo finish time is "<< C_tau*C_numOfTSt <<flush;

		return;
	}
}

void memClean()
{
	delete[] masOX;
	delete[] masOY;
	return;
}