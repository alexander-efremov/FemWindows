#include <iostream>
#include <math.h>
#include <sstream>
#include <fstream>
#include <iomanip>

#define C_pi_device 3.14159265358979323846264338327

using namespace std;

double h_u_function(double par_b, double t, double x,
                    double y)
{
	return par_b * y * (1. - y) * (C_pi_device / 2. + atan(-x));
}

double h_v_function(double lbDom, double rbDom,
                    double bbDom, double ubDom, double t, double x, double y)
{
	return atan(
		(x - lbDom) * (x - rbDom) * (1. + t) / 10. * (y - ubDom)
		* (y - bbDom));
}


double itemOfInteg_1SpecType(
	double Py,
	double Qy,
	//
	double Gx,
	double Hx,
	//
	double a,
	double b )
{
	double integ;
	integ = (Hx - a)*(Hx - a) - (Gx - a)*(Gx - a);
	integ = integ * ( (Qy - b)*(Qy - b) - (Py - b)*(Py - b) );
	return integ / 4.;
}

double analytSolut(
	double par_a,
	//
	double lbDom, //   -  Left and right boundaries of rectangular domain.
	double rbDom,
	//
	double bbDom, //   -  Botton and upper boundaries of rectangular domain.
	double ubDom,
	//
	double t, double x, double y )
{
	return 1.1 + sin( t * x * y);
}

double initDataOfSol(
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
	const double *masOY ) //   -  Massive of ordinate grid steps. Dimension = numOfOYSt +1.
{
	return analytSolut( par_a, lbDom, rbDom, bbDom, ubDom, 0., masOX[ iOfOXN ], masOY[ iOfOYN ] );
}


double itemOfInteg_2SpecType(
	double Py,
	double Qy,
	//
	double alpha,
	//
	double a,
	double b,
	double betta )
{
	double buf_D, integ;
	//   Computing...
	buf_D = (Qy - alpha) * (a*Qy + b - betta) * (a*Qy + b - betta) * (a*Qy + b - betta);
	buf_D = buf_D - (Py - alpha) * (a*Py + b - betta) * (a*Py + b - betta) * (a*Py + b - betta);
	integ = buf_D / (3. * a);
	buf_D = (a*Qy + b - betta) * (a*Qy + b - betta) * (a*Qy + b - betta) * (a*Qy + b - betta);
	buf_D = buf_D - (a*Py + b - betta) * (a*Py + b - betta) * (a*Py + b - betta) * (a*Py + b - betta);
	return integ - buf_D / (12. *a *a);
}

double integUnderLeftTr_OneCell(
	double par_a, //   -  Solution parameter.
	//
	double lbDom, //   -  Left and right boundaries of rectangular domain.
	double rbDom,
	//
	double bbDom, //   -  Botton and upper boundaries of rectangular domain.
	double ubDom,
	//
	double Py,
	double Qy,
	//
	double a_SL,
	double b_SL,
	double Hx,
	//
	double tau,
	int iCurrTL, //   -  Index of current time layer.
	//
	int * indCurSqOx, //   -  Index of current square by Ox axis.
	int * indCurSqOy, //   -  Index of current square by Oy axis.
	//
	const double * masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
	int numOfOXSt, //   -  Number of OX steps.
	//
	const double * masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
	int numOfOYSt, //   -  Number of OY steps.
	//
	double * rhoInPrevTL_asV )
{
	double hx = masOX[1] - masOX[0];
	double hy = masOY[1] - masOY[0];
	double integ = 0;
	double buf_D, bufInteg_D;
	double rho[2][2];
	double t = tau * (iCurrTL - 1.);
	double x, y;
	if( (indCurSqOx[0] >=0) && (indCurSqOx[1] <=numOfOXSt) )
	{
		if( (indCurSqOy[0] >=0) && (indCurSqOy[1] <=numOfOYSt) )
		{
			rho[0][0] = rhoInPrevTL_asV[ (numOfOXSt +1)*indCurSqOy[0] + indCurSqOx[0] ];
			rho[0][1] = rhoInPrevTL_asV[ (numOfOXSt +1)*indCurSqOy[1] + indCurSqOx[0] ];
			rho[1][0] = rhoInPrevTL_asV[ (numOfOXSt +1)*indCurSqOy[0] + indCurSqOx[1] ];
			rho[1][1] = rhoInPrevTL_asV[ (numOfOXSt +1)*indCurSqOy[1] + indCurSqOx[1] ];
		}
	}
	if( (indCurSqOx[0] < 0) || (indCurSqOx[1] > numOfOXSt) || (indCurSqOy[0] < 0) || (indCurSqOy[1] > numOfOYSt) )
	{
		x = indCurSqOx[0] * hx;
		y = indCurSqOy[0] * hy;
		rho[0][0] = analytSolut( par_a, lbDom, rbDom, bbDom, ubDom, t, x, y );
		x = indCurSqOx[0] * hx;
		y = indCurSqOy[1] * hy;
		rho[0][1] = analytSolut( par_a, lbDom, rbDom, bbDom, ubDom, t, x, y );
		x = indCurSqOx[1] * hx;
		y = indCurSqOy[0] * hy;
		rho[1][0] = analytSolut( par_a, lbDom, rbDom, bbDom, ubDom, t, x, y );
		x = indCurSqOx[1] * hx;
		y = indCurSqOy[1] * hy;
		rho[1][1] = analytSolut( par_a, lbDom, rbDom, bbDom, ubDom, t, x, y );
	}

	//   1.
	buf_D = (Qy - masOY[ indCurSqOy[1] ]) * (Qy - masOY[ indCurSqOy[1] ]) - (Py - masOY[ indCurSqOy[1] ]) * (Py - masOY[ indCurSqOy[1] ]);
	if( (indCurSqOx[1] >= 0) && (indCurSqOy[1] >= 0) )
	{
		buf_D = buf_D * (Hx - masOX[ indCurSqOx[1] ]) * (Hx - masOX[ indCurSqOx[1] ]) /4.;
		bufInteg_D = itemOfInteg_2SpecType( Py, Qy, masOY[ indCurSqOy[1] ], a_SL, b_SL, masOX[ indCurSqOx[1] ] );
	}
	else
	{
		buf_D = buf_D * (Hx - hx * indCurSqOx[1] ) * (Hx - hx * indCurSqOx[1] ) /4.;
		bufInteg_D = itemOfInteg_2SpecType( Py, Qy, hy * indCurSqOy[1], a_SL, b_SL, hx * indCurSqOx[1] );
	}
	buf_D = buf_D - bufInteg_D /2.;
	integ = buf_D * rho[0][0] /hx /hy;
	//   2.
	buf_D = (Qy - masOY[ indCurSqOy[1] ]) * (Qy - masOY[ indCurSqOy[1] ]) - (Py - masOY[ indCurSqOy[1] ]) * (Py - masOY[ indCurSqOy[1] ]);
	if( (indCurSqOx[0] >= 0) && (indCurSqOy[1] >= 0) )
	{
		buf_D = -1. * buf_D * (Hx - masOX[ indCurSqOx[0] ]) * (Hx - masOX[ indCurSqOx[0] ]) /4.;
		bufInteg_D = itemOfInteg_2SpecType( Py, Qy, masOY[ indCurSqOy[1] ], a_SL, b_SL, masOX[ indCurSqOx[0] ] );
	}
	else
	{
		buf_D = -1. * buf_D * (Hx - hx * indCurSqOx[0] ) * (Hx - hx * indCurSqOx[0] ) /4.;
		bufInteg_D = itemOfInteg_2SpecType( Py, Qy, hy * indCurSqOy[1], a_SL, b_SL, hx * indCurSqOx[0] );
	}
	buf_D = buf_D + bufInteg_D /2.;
	integ = integ + buf_D * rho[1][0] /hx /hy;
	//   3.
	buf_D = (Qy - masOY[ indCurSqOy[0] ]) * (Qy - masOY[ indCurSqOy[0] ]) - (Py - masOY[ indCurSqOy[0] ]) * (Py - masOY[ indCurSqOy[0] ]);
	if( (indCurSqOx[1] >= 0) && (indCurSqOy[0] >= 0) )
	{
		buf_D = -1. * buf_D * (Hx - masOX[ indCurSqOx[1] ]) * (Hx - masOX[ indCurSqOx[1] ]) /4.;
		bufInteg_D = itemOfInteg_2SpecType( Py, Qy, masOY[ indCurSqOy[0] ], a_SL, b_SL, masOX[ indCurSqOx[1] ] );
	}
	else
	{
		buf_D = -1. * buf_D * (Hx - hx * indCurSqOx[1] ) * (Hx - hx * indCurSqOx[1] ) /4.;
		bufInteg_D = itemOfInteg_2SpecType( Py, Qy, hy * indCurSqOy[0], a_SL, b_SL, hx * indCurSqOx[1] );
	}
	buf_D = buf_D + bufInteg_D /2.;
	integ = integ + buf_D * rho[0][1] /hx /hy;
	//   4.
	buf_D = (Qy - masOY[ indCurSqOy[0] ]) * (Qy - masOY[ indCurSqOy[0] ]) - (Py - masOY[ indCurSqOy[0] ]) * (Py - masOY[ indCurSqOy[0] ]);
	if( (indCurSqOx[0] >= 0) && (indCurSqOy[0] >= 0) )
	{
		buf_D = buf_D * (Hx - masOX[ indCurSqOx[0] ]) * (Hx - masOX[ indCurSqOx[0] ]) /4.;
		bufInteg_D = itemOfInteg_2SpecType( Py, Qy, masOY[ indCurSqOy[0] ], a_SL, b_SL, masOX[ indCurSqOx[0] ] );
	}
	else
	{
		buf_D = buf_D * (Hx - hx * indCurSqOx[0] ) * (Hx - hx * indCurSqOx[0] ) /4.;
		bufInteg_D = itemOfInteg_2SpecType( Py, Qy, hy * indCurSqOy[0], a_SL, b_SL, hx * indCurSqOx[0] );
	}
	buf_D = buf_D - bufInteg_D /2.;
	return integ + buf_D * rho[1][1] /hx /hy;
}

double integUnderRightTr_OneCell(
	double par_a, //   -  Solution parameter.
	//
	double lbDom, //   -  Left and right boundaries of rectangular domain.
	double rbDom,
	//
	double bbDom, //   -  Botton and upper boundaries of rectangular domain.
	double ubDom,
	//
	double Py,
	double Qy,
	//
	double a_SL,
	double b_SL,
	double Gx,
	//
	double tau,
	int iCurrTL, //   -  Index of current time layer.
	//
	int * indCurSqOx, //   -  Index of current square by Ox axis.
	int * indCurSqOy, //   -  Index of current square by Oy axis.
	//
	const double * masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
	int numOfOXSt, //   -  Number of OX steps.
	//
	const double * masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
	int numOfOYSt, //   -  Number of OY steps.
	//
	double * rhoInPrevTL_asV )
{
	return -1. * integUnderLeftTr_OneCell(
		par_a, //   -  Solution parameter.
		//
		lbDom, rbDom, //   -  Left and right boundaries of rectangular domain.
		//
		bbDom, ubDom, //   -  Botton and upper boundaries of rectangular domain.
		//
		Py, Qy,
		//
		a_SL, b_SL,
		Gx, //   -  double Hx,
		//
		tau, iCurrTL, //   -  Index of current time layer.
		//
		indCurSqOx, //   -  Index of current square by Ox axis.
		indCurSqOy, //   -  Index of current square by Oy axis.
		//
		masOX, numOfOXSt, //   -  Massive of OX steps. Dimension = numOfOXSt +1. Number of OX steps.
		//
		masOY, numOfOYSt, //   -  Massive of OY steps. Dimension = numOfOYSt +1. Number of OY steps.
		//
		rhoInPrevTL_asV );
}

double integUnderRectAng_OneCell(
	double par_a, //   -  Solution parameter.
	//
	double lbDom, //   -  Left and right boundaries of rectangular domain.
	double rbDom,
	//
	double bbDom, //   -  Botton and upper boundaries of rectangular domain.
	double ubDom,
	//
	double Py,
	double Qy,
	//
	double Gx,
	double Hx,
	//
	double tau,
	int iCurrTL, //   -  Index of current time layer.
	//
	int * indCurSqOx, //   -  Index of current square by Ox axis.
	int * indCurSqOy, //   -  Index of current square by Oy axis.
	//
	const double * masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
	int numOfOXSt, //   -  Number of OX steps.
	//
	const double * masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
	int numOfOYSt, //   -  Number of OY steps.
	//
	double * rhoInPrevTL_asV )
{
	//   return ( fabs( (Qy - Py) * (Hx - Gx) ) );
	double hx = masOX[1] - masOX[0];
	double hy = masOY[1] - masOY[0];
	double integ = 0;
	double buf_D;
	double rho[2][2];
	double t = tau * (iCurrTL -1.);
	double x, y;
	if( (indCurSqOx[0] >=0) && (indCurSqOy[0] >=0) )
	{
		rho[0][0] = rhoInPrevTL_asV[ (numOfOXSt +1)*indCurSqOy[0] + indCurSqOx[0] ];
		rho[0][1] = rhoInPrevTL_asV[ (numOfOXSt +1)*indCurSqOy[1] + indCurSqOx[0] ];
		rho[1][0] = rhoInPrevTL_asV[ (numOfOXSt +1)*indCurSqOy[0] + indCurSqOx[1] ];
		rho[1][1] = rhoInPrevTL_asV[ (numOfOXSt +1)*indCurSqOy[1] + indCurSqOx[1] ];
	}
	else
	{
		x = indCurSqOx[0] * hx;
		y = indCurSqOy[0] * hy;
		rho[0][0] = analytSolut( par_a, lbDom, rbDom, bbDom, ubDom, t, x, y );
		x = indCurSqOx[0] * hx;
		y = indCurSqOy[1] * hy;
		rho[0][1] = analytSolut( par_a, lbDom, rbDom, bbDom, ubDom, t, x, y );
		x = indCurSqOx[1] * hx;
		y = indCurSqOy[0] * hy;
		rho[1][0] = analytSolut( par_a, lbDom, rbDom, bbDom, ubDom, t, x, y );
		x = indCurSqOx[1] * hx;
		y = indCurSqOy[1] * hy;
		rho[1][1] = analytSolut( par_a, lbDom, rbDom, bbDom, ubDom, t, x, y );
	}

	if( (indCurSqOx[1] >= 0) && (indCurSqOy[1] >= 0) )
	{
		buf_D = itemOfInteg_1SpecType( Py, Qy, Gx, Hx, masOX[ indCurSqOx[1] ], masOY[ indCurSqOy[1] ] );
	}
	else
	{
		buf_D = itemOfInteg_1SpecType( Py, Qy, Gx, Hx, hx *indCurSqOx[1] , hy * indCurSqOy[1] );
	}
	buf_D = buf_D /hx /hy;
	integ = buf_D * rho[0][0]; //   rhoInPrevTL[ indCurSqOx[0] ][ indCurSqOy[0] ];
	if( (indCurSqOx[0] >= 0) && (indCurSqOy[1] >= 0) )
	{
		buf_D = itemOfInteg_1SpecType( Py, Qy, Gx, Hx, masOX[ indCurSqOx[0] ], masOY[ indCurSqOy[1] ] );
	}
	else
	{
		buf_D = itemOfInteg_1SpecType( Py, Qy, Gx, Hx, hx * indCurSqOx[0] , hy * indCurSqOy[1] );
	}
	buf_D = buf_D /hx /hy;
	integ = integ - buf_D * rho[1][0]; //   rhoInPrevTL[ indCurSqOx[1] ][ indCurSqOy[0] ];
	if( (indCurSqOx[1] >= 0) && (indCurSqOy[0] >= 0) )
	{
		buf_D = itemOfInteg_1SpecType( Py, Qy, Gx, Hx, masOX[ indCurSqOx[1] ], masOY[ indCurSqOy[0] ] );
	}
	else
	{
		buf_D = itemOfInteg_1SpecType( Py, Qy, Gx, Hx, hx * indCurSqOx[1] , hy * indCurSqOy[0] );
	}
	buf_D = buf_D /hx /hy;
	integ = integ - buf_D * rho[0][1]; //   rhoInPrevTL[ indCurSqOx[0] ][ indCurSqOy[1] ];
	if( (indCurSqOx[0] >= 0) && (indCurSqOy[0] >= 0) )
	{
		buf_D = itemOfInteg_1SpecType( Py, Qy, Gx, Hx, masOX[ indCurSqOx[0] ], masOY[ indCurSqOy[0] ] );
	}
	else
	{
		buf_D = itemOfInteg_1SpecType( Py, Qy, Gx, Hx, hx * indCurSqOx[0] , hy * indCurSqOy[0] );
	}
	buf_D = buf_D /hx /hy;
	return integ + buf_D * rho[1][1]; //   rhoInPrevTL[ indCurSqOx[1] ][ indCurSqOy[1] ];
}


double integOfChan_SLRightSd( //   -  The domain is Channel with Slant Line on the right side.
	double par_a, //   -  Solution parameter.
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
	double *bv, int wTrPCI, //   -  Where travel point current (botton vertex) is.
	double *uv, int wTrPNI, //   -  Where travel point next (upper vertex) is.
	//
	int * indCurSqOx, //   -  Index by OX axis where bv and uv are.
	//
	double lb, int * indLB, //   -  Left boundary by Ox. Index by OX axis where lb is.
	//
	int * indCurSqOy, //   -  Index of current square by Oy axis.
	//
	const double * masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
	int numOfOXSt, //   -  Number of OX steps.
	//
	const double * masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
	int numOfOYSt, //   -  Number of OY steps.
	//
	double * rhoInPrevTL_asV )
{
	double mv[2], rv[2]; //   -  Middle and right vertices.
	int wMvI; //   -  Where middle vertex is.
	int indCurSqOxToCh[2]; //   -  Indices of current square by Ox axis to be changed. Under which we want to integrate.
	double h = masOX[1] - masOX[0];
	double a_SL, b_SL; //   -  Coefficients of slant line: x = a_SL *y  +  b_SL.
	double Gx, Hx; //   -  Left boundary for each integration.
	double integ = 0.;
	double buf_D;
	int j;

	//   Let's compute helpful values.

	if( uv[0] <= bv[0] )
	{
		mv[0] = uv[0];
		mv[1] = uv[1];
		wMvI = wTrPNI;
		rv[0] = bv[0];
		rv[1] = bv[1];
	}

	if( uv[0] > bv[0] )
	{
		mv[0] = bv[0];
		mv[1] = bv[1];
		wMvI = wTrPCI;
		rv[0] = uv[0];
		rv[1] = uv[1];
	}


	//   buf_D = fabs(mv[0] - lb) * fabs(uv[1] - bv[1])  +   fabs(uv[1] - bv[1]) * fabs(rv[0] - mv[0]) / 2.;
	//   return  buf_D;


	if( ( fabs(uv[1] - bv[1]) ) <= 1.e-12 )
	{
		//   Computation is impossible. Too smale values. Let's return some approximate value.
		//   buf_D  =  (uv[1] - bv[1])  *  ((uv[0] + bv[0]) /2.  -  lb) * rhoInPrevTL[ indCurSqOx[0] ][ indCurSqOy[0] ];
		return fabs(uv[1] - bv[1]); //   fabs(uv[1] - bv[1]);
	}


	//   First step: from "lb" to "masOX[ indCurSqOx[0] ]" by iteration.
	//   integ  += fabs( mv[0] - lb) * fabs(uv[1] - bv[1]);

	indCurSqOxToCh[0] = indLB[0];
	indCurSqOxToCh[1] = indCurSqOxToCh[0] +1;

	for( j = indLB[0]; j< indCurSqOx[0]; j++ )
	{
		//   If this is first cell we should integrate under rectangle only.
		if( indCurSqOxToCh[0] >= 0 )
		{
			Gx = masOX[ indCurSqOxToCh[0] ];
			Hx = masOX[ indCurSqOxToCh[1] ];
		}


		if( indCurSqOxToCh[0] < 0 )
		{
			Gx = h * indCurSqOxToCh[0];
			Hx = h * indCurSqOxToCh[1];
		}

		if( j == indLB[0] )
		{
			Gx = lb;
		}

		buf_D = integUnderRectAng_OneCell(
			par_a, //   -  Solution parameter.
			//
			lbDom, rbDom, //   -  Left and right boundaries of rectangular domain.
			//
			bbDom, ubDom, //   -  Botton and upper boundaries of rectangular domain.
			//
			bv[1], //   -  double Py,
			uv[1], //   -  double Qy,
			//
			Gx, //   -  double Gx,
			Hx, //   -  double Hx,
			//
			tau, iCurrTL, //   -  Index of current time layer.
			//
			indCurSqOxToCh, //   -  Index of current square by Ox axis.
			indCurSqOy, //   -  Index of current square by Oy axis.
			//
			masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
			numOfOXSt, //   -  Number of OX steps.
			//
			masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
			numOfOYSt, //   -  Number of OY steps.
			//
			rhoInPrevTL_asV );

		integ += buf_D;

		indCurSqOxToCh[0] += 1;
		indCurSqOxToCh[1] = indCurSqOxToCh[0] +1;
	}

	//   Integration. Second step: under [ indCurSqOx[0]; indCurSqOx[1] ] square.

	//   A. Under rectangle.
	if( wMvI == 1 )
	{
		if( indCurSqOx[0] == indLB[0] )
		{
			Gx = lb;
		}

		if( indCurSqOx[0] > indLB[0] )
		{
			if( indCurSqOx[0] >= 0)
			{
				Gx = masOX[ indCurSqOx[0] ];
			}

			if( indCurSqOx[0] < 0)
			{
				Gx = h * indCurSqOx[0];
			}
		}

		buf_D = integUnderRectAng_OneCell(
			par_a, //   -  Solution parameter.
			//
			lbDom, rbDom, //   -  Left and right boundaries of rectangular domain.
			//
			bbDom, ubDom, //   -  Botton and upper boundaries of rectangular domain.
			//
			bv[1], //   -  double Py,
			uv[1], //   -  double Qy,
			//
			Gx, //   -  double Gx,
			mv[0], //   -  double Hx,
			//
			tau, iCurrTL, //   -  Index of current time layer.
			//
			indCurSqOx, //   -  Index of current square by Ox axis.
			indCurSqOy, //   -  Index of current square by Oy axis.
			//
			masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
			numOfOXSt, //   -  Number of OX steps.
			//
			masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
			numOfOYSt, //   -  Number of OY steps.
			//
			rhoInPrevTL_asV );

		integ += buf_D;
	}

	//   B. Under triangle.

	if( ( fabs(uv[1] - bv[1]) ) > 1.e-12 )
	{
		//   integ += fabs(uv[1] - bv[1]) * (rv[0] - mv[0]) /2.;
		//   Coefficients of slant line: x = a_SL *y  +  b_SL.
		a_SL = (uv[0] - bv[0]) / (uv[1] - bv[1]);
		b_SL = bv[0] - a_SL * bv[1];


		//   Integration under one cell triangle.

		if( fabs( a_SL ) > 1.e-12 )
		{
			buf_D = integUnderRightTr_OneCell(
				par_a, //   -  Solution parameter.
				//
				lbDom, rbDom, //   -  Left and right boundaries of rectangular domain.
				//
				bbDom, ubDom, //   -  Botton and upper boundaries of rectangular domain.
				//
				bv[1], //   -  double Py,
				uv[1], //   -  double Qy,
				//
				a_SL,
				b_SL,
				mv[0], //   -  double Gx,
				//
				tau, iCurrTL, //   -  Index of current time layer.
				//
				indCurSqOx, //   -  Index of current square by Ox axis.
				indCurSqOy, //   -  Index of current square by Oy axis.
				//
				masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
				numOfOXSt, //   -  Number of OX steps.
				//
				masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
				numOfOYSt, //   -  Number of OY steps.
				//
				rhoInPrevTL_asV );

			integ += buf_D;
		}
	}
	return integ;
}


double upperBound(
	double par_a, //   -  Item of left and right setback (parameter "a" in test).
	//
	double lbDom, //   -  Left and right boundaries of rectangular domain.
	double rbDom,
	//
	double bbDom, //   -  Botton and upper boundaries of rectangular domain.
	double ubDom,
	//
	double t,
	double x )
{
	return analytSolut( par_a, lbDom, rbDom, bbDom, ubDom, t, x, ubDom );
}

double rightBound(
	double par_a, //   -  Item of left and right setback (parameter "a" in test).
	//
	double lbDom, //   -  Left and right boundaries of rectangular domain.
	double rbDom,
	//
	double bbDom, //   -  Botton and upper boundaries of rectangular domain.
	double ubDom,
	//
	double t,
	double y )
{
	return analytSolut( par_a, lbDom, rbDom, bbDom, ubDom, t, rbDom, y );
}

double leftBound(
	double par_a, //   -  Item of left and right setback (parameter "a" in test).
	//
	double lbDom, //   -  Left and right boundaries of rectangular domain.
	double rbDom,
	//
	double bbDom, //   -  Botton and upper boundaries of rectangular domain.
	double ubDom,
	//
	double t,
	double y )
{
	return analytSolut( par_a,lbDom, rbDom, bbDom, ubDom, t, lbDom, y );
}

double integOfChan_SLLeftSd( //   -  The domain is Channel with Slant Line on the left side.
	double par_a, //   -  Solution parameter.
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
	double *bv, int wTrPCI, //   -  Where travel point current (botton vertex) is.
	double *uv, int wTrPNI, //   -  Where travel point next (upper vertex) is.
	//
	int * indCurSqOx, //   -  Index by OX axis where bv and uv are.
	//
	double rb, int * indRB, //   -  Right boundary by Ox. Index by OX axis where rb is.
	//
	int * indCurSqOy, //   -  Index of current square by Oy axis.
	//
	const double * masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
	int numOfOXSt, //   -  Number of OX steps.
	//
	const double * masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
	int numOfOYSt, //   -  Number of OY steps.
	//
	double * rhoInPrevTL_asV )
{
	double lv[2], mv[2]; //   -  Left and middle vertices.
	int wMvI; //   -  Where middle vertex is.
	int indCurSqOxToCh[2]; //   -  Indices of current square by Ox axis to be changed. Under which we want to integrate.
	double h = masOX[1] - masOX[0];
	double a_SL, b_SL; //   -  Coefficients of slant line: x = a_SL *y  +  b_SL.
	double Gx, Hx; //   -  Left and right boundary for each integration.
	double integ = 0.;
	double buf_D;
	int j;

	//   Let's compute helpful values.

	if( uv[0] <= bv[0] )
	{
		lv[0] = uv[0];
		lv[1] = uv[1];
		mv[0] = bv[0];
		mv[1] = bv[1];
		wMvI = wTrPCI;
	}

	if( uv[0] > bv[0] )
	{
		lv[0] = bv[0];
		lv[1] = bv[1];
		mv[0] = uv[0];
		mv[1] = uv[1];
		wMvI = wTrPNI;
	}

	if( ( fabs(uv[1] - bv[1]) ) <= 1.e-12 )
	{
		//   Computation is impossible. Too smale values. Let's return some approximate value.
		//   buf_D  =  (uv[1] - bv[1])  *  (rb  - (uv[0] + bv[0]) /2.) * rhoInPrevTL[ indCurSqOx[0] ][ indCurSqOy[0] ];
		return fabs(uv[1] - bv[1]); //   fabs(uv[1] - bv[1]);
	}

	//   Integration. First step: under [ indCurSqOx[0]; indCurSqOx[1] ] square.

	//   A. Under triangle.

	if( ( fabs(uv[1] - bv[1]) ) > 1.e-12 )
	{
		//   Coefficients of slant line: x = a_SL *y  +  b_SL.
		a_SL = (uv[0] - bv[0]) / (uv[1] - bv[1]);
		b_SL = bv[0] - a_SL * bv[1];

		//   Integration under one cell triangle.
		if( fabs( a_SL ) > 1.e-12 )
		{
			buf_D = integUnderLeftTr_OneCell(
				par_a, //   -  Solution parameter.
				//
				lbDom, rbDom, //   -  Left and right boundaries of rectangular domain.
				//
				bbDom, ubDom, //   -  Botton and upper boundaries of rectangular domain.
				//
				bv[1], //   -  double Py,
				uv[1], //   -  double Qy,
				//
				a_SL,
				b_SL,
				mv[0], //   -  double Hx,
				//
				tau, iCurrTL, //   -  Index of current time layer.
				//
				indCurSqOx, //   -  Index of current square by Ox axis.
				indCurSqOy, //   -  Index of current square by Oy axis.
				//
				masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
				numOfOXSt, //   -  Number of OX steps.
				//
				masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
				numOfOYSt, //   -  Number of OY steps.
				//
				rhoInPrevTL_asV );


			integ += buf_D;
		}
	}


	//   B. Under rectangle. Need to be cheking.

	if( wMvI == 1 )
	{
		if( indCurSqOx[0] == indRB[0] )
		{
			Hx = rb;
		}

		if( indCurSqOx[0] < indRB[0] )
		{
			if( indCurSqOx[1] >= 0)
			{
				Hx = masOX[ indCurSqOx[1] ];
			}

			if( indCurSqOx[1] < 0)
			{
				Hx = h * indCurSqOx[1];
			}
		}

		buf_D = integUnderRectAng_OneCell(
			par_a, //   -  Solution parameter.
			//
			lbDom, rbDom, //   -  Left and right boundaries of rectangular domain.
			//
			bbDom, ubDom, //   -  Botton and upper boundaries of rectangular domain.
			//
			bv[1], //   -  double Py,
			uv[1], //   -  double Qy,
			//
			mv[0], //   -  double Gx,
			Hx, //   -  double Hx,
			//
			tau, iCurrTL, //   -  Index of current time layer.
			//
			indCurSqOx, //   -  Index of current square by Ox axis.
			indCurSqOy, //   -  Index of current square by Oy axis.
			//
			masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
			numOfOXSt, //   -  Number of OX steps.
			//
			masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
			numOfOYSt, //   -  Number of OY steps.
			//
			rhoInPrevTL_asV );

		integ += buf_D;
	}

	//   Second step: from "masOX[ indCurSqOx[1] ]" to "rb" by iteration.


	indCurSqOxToCh[0] = indCurSqOx[0] +1;
	indCurSqOxToCh[1] = indCurSqOxToCh[0] +1;

	for( j = indCurSqOx[0] +1; j< indRB[0] +1; j++ )
	{
		//   If this is first cell we should integrate under triangle only.

		if( indCurSqOxToCh[1] > 0 )
		{
			Gx = masOX[ indCurSqOxToCh[0] ];
			Hx = masOX[ indCurSqOxToCh[1] ];
		}


		if( indCurSqOxToCh[1] <= 0 )
		{
			Gx = h * indCurSqOxToCh[0];
			Hx = h * indCurSqOxToCh[1];
		}


		if( j == indRB[0] )
		{
			Hx = rb;
		}


		buf_D = integUnderRectAng_OneCell(
			par_a, //   -  Solution parameter.
			//
			lbDom, rbDom, //   -  Left and right boundaries of rectangular domain.
			//
			bbDom, ubDom, //   -  Botton and upper boundaries of rectangular domain.
			//
			bv[1], //   -  double Py,
			uv[1], //   -  double Qy,
			//
			Gx, //   -  double Gx,
			Hx, //   -  double Hx,
			//
			tau, iCurrTL, //   -  Index of current time layer.
			//
			indCurSqOxToCh, //   -  Index of current square by Ox axis.
			indCurSqOy, //   -  Index of current square by Oy axis.
			//
			masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
			numOfOXSt, //   -  Number of OX steps.
			//
			masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
			numOfOYSt, //   -  Number of OY steps.
			//
			rhoInPrevTL_asV );

		integ += buf_D;


		indCurSqOxToCh[0] += 1;
		indCurSqOxToCh[1] = indCurSqOxToCh[0] +1;
	}

	return integ;
}


double integUnderRigAngTr_BottLeft(
	double par_a, //   -  Solution parameter.
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
	double *bv,
	double *uv,
	//
	const double * masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
	int numOfOXSt, //   -  Number of OX steps.
	//
	const double * masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
	int numOfOYSt, //   -  Number of OY steps.
	//
	double * rhoInPrevTL_asV )
{
	double trPC[2]; //   -  Travel point current;
	int wTrPCI = 0; //   -  Where travel point current is?
	double trPN[2]; //   -  Travel point next;
	int wTrPNI = 0; //   -  Where travel point next is?
	double ang; //   -  Angle of slant line. Should be greater zero.
	int indCurSqOx[2], indCurSqOy[2]; //   -  Index of current square by Ox and Oy axes.
	int indRB[2]; //   -  Index of right boundary.
	double distOx, distOy; //   -  Distance to near Ox and Oy straight lines.
	bool isTrDone = false; //   -  Is travel done.
	double hx = masOX[1] - masOX[0];
	double hy = masOY[1] - masOY[0];
	double integOfBottTr = 0.; //   -  Value which we are computing.
	double buf_D;
	//   Initial data.
	trPC[0] = bv[0];
	trPC[1] = bv[1];
	if( ( fabs(bv[0] - uv[0]) ) < 1.e-12 )
	{
		//   This triangle has very small width. I guess further computation isn't correct.
		return fabs(bv[0] - uv[0]);
	}
	ang = (uv[1] - bv[1]) / (bv[0] - uv[0]);
	if( fabs(ang) < 1.e-12 )
	{
		//   This triangle has very small height. I guess further computation isn't correct.
		return fabs(ang);
	}
	indCurSqOx[0] = (int)( (trPC[0] - 1.e-14) /hx); //   -  If trPC[0] is in grid edge I want it will be between in the left side of indCurSqOx[1].
	if( (trPC[0] - 1.e-14) <= 0 )
	{
		indCurSqOx[0] -= 1; //   -  The case when "trPC[0]" ia negative.
	}
	indCurSqOx[1] = indCurSqOx[0] +1; //   -  It's important only in rare case then trPC is in grid edge.
	indRB[0] = indCurSqOx[0];
	indRB[1] = indRB[0] +1;
	indCurSqOy[0] = (int)( (trPC[1] + 1.e-14) /hy); //   -  If trPC[1] is in grid edge I want it will be between indCurSqOx[0] and indCurSqOx[1].
	if( (trPC[1] + 1.e-14) <= 0 )
	{
		indCurSqOy[0] -= 1; //   -  The case when "trPC[0]" ia negative.
	}
	indCurSqOy[1] = indCurSqOy[0] +1; //   -  It's important only in rare case then trPC is in grid edge.
	if( indCurSqOx[0] >= 0)
	{
		distOx = trPC[0] - masOX[ indCurSqOx[0] ];
	}
	if( indCurSqOx[0] < 0 )
	{
		distOx = fabs( trPC[0] - hx * indCurSqOx[0] );
	}
	if( indCurSqOy[1] >= 0 )
	{
		distOy = masOY[ indCurSqOy[1] ] - trPC[1];
	}
	if( indCurSqOy[1] < 0 )
	{
		distOy = fabs( hy * indCurSqOy[1] - trPC[1] );
	}
	do
	{
		//   a. First case.
		if( (distOy /distOx) <= ang )
		{
			//   Across with straight line parallel Ox axis.
			wTrPNI = 1;
			if( indCurSqOy[1] >= 0)
			{
				trPN[1] = masOY[ indCurSqOy[1] ];
			}
			if( indCurSqOy[1] < 0)
			{
				trPN[1] = hy * indCurSqOy[1];
			}
			trPN[0] = bv[0] - (trPN[1] - bv[1]) /ang;
		}
		//   b. Second case.
		if( (distOy /distOx) > ang )
		{
			//   Across with straight line parallel Oy axis.
			wTrPNI = 2;
			if( indCurSqOx[0] >= 0 )
			{
				trPN[0] = masOX[ indCurSqOx[0] ];
			}
			if( indCurSqOx[0] < 0 )
			{
				trPN[0] = hx * indCurSqOx[0];
			}
			trPN[1] = bv[1] - ang * (trPN[0] - bv[0]);
		}
		//   c. Cheking.
		if( trPN[0] < (uv[0] + 1.e-14) )
		{
			trPN[0] = uv[0];
			trPN[1] = uv[1];
			isTrDone = true;
			wTrPNI = 0;
		}
		//   d. Integration.
		buf_D = integOfChan_SLLeftSd( //   -  The domain is Channel with Slant Line on the left side.
			par_a, //   -  Solution parameter.
			//
			lbDom, rbDom,
			//
			bbDom, ubDom,
			//
			tau, iCurrTL, //   -  Index of current time layer.
			//
			trPC, wTrPCI, //   -  double *bv,
			trPN, wTrPNI, //   -  double *uv,
			//
			indCurSqOx, //   -  Indices where trPC and trPN are.
			//
			bv[0], indRB, //   -  double rb  =  Right boundary by Ox.
			//
			indCurSqOy, //   -  Index of current square by Oy axis.
			//
			masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
			numOfOXSt, //   -  Number of OX steps.
			//
			masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
			numOfOYSt, //   -  Number of OY steps.
			//
			rhoInPrevTL_asV );
		integOfBottTr = integOfBottTr + buf_D;
		//   e. Updating.
		if( isTrDone == false )
		{
			//   We will compute more. We need to redefine some values.
			wTrPCI = wTrPNI;
			trPC[0] = trPN[0];
			trPC[1] = trPN[1];
			if( wTrPNI == 1)
			{
				indCurSqOy[0] += 1;
				indCurSqOy[1] += 1;
			}
			if( wTrPNI == 2)
			{
				indCurSqOx[0] -= 1;
				indCurSqOx[1] -= 1;
			}
			if( indCurSqOx[0] >= 0)
			{
				distOx = trPC[0] - masOX[ indCurSqOx[0] ];
			}
			if( indCurSqOx[0] < 0)
			{
				distOx = fabs( trPC[0] - hx * indCurSqOx[0] );
			}
			if( indCurSqOy[1] >= 0 )
			{
				distOy = masOY[ indCurSqOy[1] ] - trPC[1];
			}
			if( indCurSqOy[1] < 0 )
			{
				distOy = fabs( hy * indCurSqOy[1] - trPC[1] );
			}
		}
	} while( !isTrDone );
	return integOfBottTr;
}

double integUnderRigAngTr_BottRight(
	double par_a, //   -  Solution parameter.
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
	double *bv,
	double *uv,
	//
	const double * masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
	int numOfOXSt, //   -  Number of OX steps.
	//
	const double * masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
	int numOfOYSt, //   -  Number of OY steps.
	//
	double * rhoInPrevTL_asV )
{
	double trPC[2]; //   -  Travel point current;
	int wTrPCI = 0; //   -  Where travel point current is?
	double trPN[2]; //   -  Travel point next;
	int wTrPNI = 0; //   -  Where travel point next is?
	double ang; //   -  Angle of slant line. Should be greater zero.
	int indCurSqOx[2], indCurSqOy[2]; //   -  Index of current square by Ox and Oy axes.
	int indLB[2]; //   -  Index of left boundary.
	double distOx, distOy; //   -  Distance to near Ox and Oy straight lines.
	bool isTrDone = false; //   -  Is travel done.
	double hx = masOX[1] - masOX[0];
	double hy = masOY[1] - masOY[0];
	double integOfBottTr = 0.; //   -  Value which we are computing.
	double buf_D;


	trPC[0] = bv[0];
	trPC[1] = bv[1];
	if( ( fabs(bv[0] - uv[0]) ) < 1.e-12 ) return fabs(bv[0] - uv[0]);

	ang = (uv[1] - bv[1]) / (uv[0] - bv[0]);
	if( fabs(ang) < 1.e-12 ) return fabs(ang);

	indCurSqOx[0] = (int)( (trPC[0] + 1.e-14) /hx); //   -  If trPC[0] is in grid edge I want it will be between in the right side.

	if( (trPC[0] + 1.e-14) <= 0 ) indCurSqOx[0] -= 1; //   -  The case when "trPC[0]" ia negative.

	indCurSqOx[1] = indCurSqOx[0] +1; //   -  It's important only in rare case then trPC is in grid edge.
	indLB[0] = indCurSqOx[0];
	indLB[1] = indLB[0] +1;
	indCurSqOy[0] = (int)( (trPC[1] + 1.e-14) /hy); //   -  If trPC[1] is in grid edge I want it will be in the upper side.
	if( (trPC[1] + 1.e-14) <= 0 )
	{
		indCurSqOy[0] -= 1; //   -  The case when "trPC[0]" ia negative.
	}
	indCurSqOy[1] = indCurSqOy[0] +1; //   -  It's important only in rare case then trPC is in grid edge.

	if( indCurSqOx[1] >=0 )
	{
		distOx = fabs( masOX[ indCurSqOx[1] ] - trPC[0] );
	}
	if( indCurSqOx[1] < 0 )
	{
		distOx = fabs( hx * indCurSqOx[1] - trPC[0] );
	}
	if( indCurSqOy[1] >=0 )
	{
		distOy = fabs( masOY[ indCurSqOy[1] ] - trPC[1] );
	}
	if( indCurSqOy[1] < 0 )
	{
		distOy = fabs( hy * indCurSqOy[1] - trPC[1] );
	}
	do
	{
		//   a. First case.
		if( (distOy /distOx) <= ang )
		{
			//   Across with straight line parallel Ox axis.
			wTrPNI = 1;
			if( indCurSqOy[1] >=0 )
			{
				trPN[1] = masOY[ indCurSqOy[1] ];
			}
			if( indCurSqOy[1] < 0 )
			{
				trPN[1] = hy * indCurSqOy[1];
			}
			trPN[0] = bv[0] + (trPN[1] - bv[1]) /ang;
		}
		//   b. Second case.
		if( (distOy /distOx) > ang )
		{
			//   Across with straight line parallel Oy axis.
			wTrPNI = 2;
			if( indCurSqOx[1] >= 0 )
			{
				trPN[0] = masOX[ indCurSqOx[1] ];
			}
			if( indCurSqOx[1] < 0 )
			{
				trPN[0] = hx * indCurSqOx[1];
			}
			trPN[1] = bv[1] + ang * (trPN[0] - bv[0]);
		}
		//   c. Cheking.
		if( trPN[0] > (uv[0] - 1.e-14) )
		{ //   -  Without "fabs"!!!
			trPN[0] = uv[0];
			trPN[1] = uv[1];
			isTrDone = true;
			wTrPNI = 0;
		}
		//   d. Integration.
		buf_D = integOfChan_SLRightSd( //   -  The domain is Channel with Slant Line on the Right side.
			par_a, //   -  Solution parameter.
			//
			lbDom, rbDom,
			//
			bbDom, ubDom,
			//
			tau, iCurrTL, //   -  Index of current time layer.
			//
			trPC, wTrPCI, //   -  double *bv,
			trPN, wTrPNI, //   -  double *uv,
			//
			indCurSqOx, //   -  Indices where trPC and trPN are.
			//
			bv[0], indLB, //   -  double lb  =  Left boundary by Ox.
			//
			indCurSqOy, //   -  Index of current square by Oy axis.
			//
			masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
			numOfOXSt, //   -  Number of OX steps.
			//
			masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
			numOfOYSt, //   -  Number of OY steps.
			//
			rhoInPrevTL_asV );
		integOfBottTr = integOfBottTr + buf_D;
		//   e. Updating.
		if( isTrDone == false )
		{
			//   We will compute more. We need to redefine some values.
			wTrPCI = wTrPNI;
			trPC[0] = trPN[0];
			trPC[1] = trPN[1];
			if( wTrPNI == 1)
			{
				indCurSqOy[0] += 1;
				indCurSqOy[1] += 1;
			}
			if( wTrPNI == 2)
			{
				indCurSqOx[0] += 1;
				indCurSqOx[1] += 1;
			}
			if( indCurSqOx[1] >=0 )
			{
				distOx = fabs( masOX[ indCurSqOx[1] ] - trPC[0] );
			}
			if( indCurSqOx[1] < 0 )
			{
				distOx = fabs( hx * indCurSqOx[1] - trPC[0] );
			}
			if( indCurSqOy[1] >=0 )
			{
				distOy = fabs( masOY[ indCurSqOy[1] ] - trPC[1] );
			}
			if( indCurSqOy[1] < 0 )
			{
				distOy = fabs( hy * indCurSqOy[1] - trPC[1] );
			}
		}
	} while( !isTrDone );
	return integOfBottTr;
}

double integUnderBottTr(
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
	double * LvBt, //   -  Left, Right and Botton vertices of Botton triangle.
	double * RvBt, //   -  Left, Right and Botton vertices of Botton triangle.
	double * BvBt, //   -  Left, Right and Botton vertices of Botton triangle.
	//
	const double * masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
	int numOfOXSt, //   -  Number of OX steps.
	//
	const double * masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
	int numOfOYSt, //   -  Number of OY steps.
	//
	double * rhoInPrevTL_asV,
	int ii, int jj ) // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
{
	double integOfBottTr;
	double buf_D;
	//   Three ways are possible.
	//   1.
	if( BvBt[0] <= LvBt[0] )
	{
		buf_D = integUnderRigAngTr_BottRight(
			par_a, lbDom, rbDom, bbDom, ubDom, tau, iCurrTL,
			//
			BvBt, RvBt, masOX, numOfOXSt, masOY, numOfOYSt, rhoInPrevTL_asV );
		integOfBottTr = buf_D;
		buf_D = integUnderRigAngTr_BottRight(
			par_a, lbDom, rbDom, bbDom, ubDom, tau, iCurrTL,
			//
			BvBt, LvBt, masOX, numOfOXSt, masOY, numOfOYSt, rhoInPrevTL_asV );
		integOfBottTr = integOfBottTr - buf_D;

		//		printf("Bv<Lv: i= %d, j= %d      res= %le",ii,jj,integOfBottTr);  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		return integOfBottTr;
	}
	//   2.
	if( (BvBt[0] > LvBt[0]) && (BvBt[0] < RvBt[0]) )
	{
		buf_D = integUnderRigAngTr_BottLeft(
			par_a, lbDom, rbDom, bbDom, ubDom, tau, iCurrTL,
			//
			BvBt, LvBt, masOX, numOfOXSt, masOY, numOfOYSt, rhoInPrevTL_asV );
		integOfBottTr = buf_D;

		buf_D = integUnderRigAngTr_BottRight(
			par_a, lbDom, rbDom, bbDom, ubDom, tau, iCurrTL,
			//
			BvBt, RvBt, masOX, numOfOXSt, masOY, numOfOYSt, rhoInPrevTL_asV );
		integOfBottTr = integOfBottTr + buf_D;

		//		printf("Bv>Lv & Bv<Rv: i= %d, j= %d      res= %le",ii,jj,integOfBottTr);   // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		return integOfBottTr;
	}
	//   3.
	if( BvBt[0] >= RvBt[0] )
	{
		buf_D = integUnderRigAngTr_BottLeft(
			par_a, lbDom, rbDom, bbDom, ubDom, tau, iCurrTL,
			//
			BvBt, LvBt, masOX, numOfOXSt, masOY, numOfOYSt, rhoInPrevTL_asV );
		integOfBottTr = buf_D;
		buf_D = integUnderRigAngTr_BottLeft(
			par_a, lbDom, rbDom, bbDom, ubDom, tau, iCurrTL,
			//
			BvBt, RvBt, masOX, numOfOXSt, masOY, numOfOYSt, rhoInPrevTL_asV );
		integOfBottTr = integOfBottTr - buf_D;

		//		printf("Bv>Rv: i= %d, j= %d      res= %le",ii,jj,integOfBottTr);     // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		return integOfBottTr;
	}
	return integOfBottTr;
}

double integUnderRigAngTr_UppLeft(
	double par_a, //   -  Solution parameter.
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
	double *bv,
	double *uv,
	//
	const double * masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
	int numOfOXSt, //   -  Number of OX steps.
	//
	const double * masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
	int numOfOYSt, //   -  Number of OY steps.
	//
	double * rhoInPrevTL_asV )
{
	//   return ( fabs( (uv[1] - bv[1]) * (bv[0] - uv[0]) /2.) );
	double trPC[2]; //   -  Travel point current;
	int wTrPCI = 0; //   -  Where travel point current is?
	double trPN[2]; //   -  Travel point next;
	int wTrPNI = 0; //   -  Where travel point next is?
	double ang; //   -  Angle of slant line. Should be greater zero.
	int indCurSqOx[2], indCurSqOy[2]; //   -  Index of current square by Ox and Oy axes.
	int indRB[2]; //   -  Index of right boundary.
	double distOx, distOy; //   -  Distance to near Ox and Oy straight lines.
	bool isTrDone = false; //   -  Is travel done.
	double hx = masOX[1] - masOX[0];
	double hy = masOY[1] - masOY[0];
	double integOfUppTr = 0.; //   -  Value which we are computing.
	double buf_D;
	//   Initial data.
	trPC[0] = bv[0];
	trPC[1] = bv[1];
	if( ( fabs(bv[0] - uv[0]) ) < 1.e-12 ) return fabs(bv[0] - uv[0]);

	ang = (uv[1] - bv[1]) / (uv[0] - bv[0]);
	if( fabs(ang) < 1.e-12 ) return fabs(ang);

	//   The follow equations are quite important.
	indCurSqOx[0] = (int)( (trPC[0] + 1.e-14) /hx); //   -  If trPC[0] is in grid edge I want it will be in the right side.
	if( (trPC[0] + 1.e-14) <= 0 )
	{
		indCurSqOx[0] -= 1; //   -  The case when "trPC[0]" ia negative.
	}
	indCurSqOx[1] = indCurSqOx[0] +1; //   -  It's important only in rare case then trPC is in grid edge.
	indCurSqOy[0] = (int)( (trPC[1] + 1.e-14) /hy); //   -  If trPC[1] is in grid edge I want it will be in the upper square.
	if( (trPC[1] + 1.e-14) <= 0 )
	{
		indCurSqOy[0] -= 1; //   -  The case when "trPC[0]" ia negative.
	}
	indCurSqOy[1] = indCurSqOy[0] +1;
	indRB[0] = (int)( (uv[0] - 1.e-14) /hy); //   -  If uv[0] is in grid edge I want it will be in the left side.
	if( (uv[0] - 1.e-14) <= 0 )
	{
		indRB[0] -= 1; //   -  The case when "trPC[0]" ia negative.
	}
	indRB[1] = indRB[0] +1;
	if( indCurSqOx[1] >= 0)
	{
		distOx = masOX[ indCurSqOx[1] ] - trPC[0];
	}
	if( indCurSqOx[1] < 0)
	{
		distOx = fabs( hx * indCurSqOx[1] - trPC[0] );
	}
	if( indCurSqOy[1] >= 0 )
	{
		distOy = masOY[ indCurSqOy[1] ] - trPC[1];
	}
	if( indCurSqOy[1] < 0 )
	{
		distOy = fabs( hy * indCurSqOy[1] - trPC[1] );
	}
	do
	{
		//   a. First case.
		if( (distOy /distOx) <= ang )
		{
			//   Across with straight line parallel Ox axis.
			wTrPNI = 1;
			if( indCurSqOy[1] >= 0 )
			{
				trPN[1] = masOY[ indCurSqOy[1] ];
			}
			if( indCurSqOy[1] < 0 )
			{
				trPN[1] = hy * indCurSqOy[1];
			}
			trPN[0] = bv[0] + (trPN[1] - bv[1]) /ang;
		}
		//   b. Second case.
		if( (distOy /distOx) > ang )
		{
			//   Across with straight line parallel Oy axis.
			wTrPNI = 2;
			if( indCurSqOx[1] >= 0 )
			{
				trPN[0] = masOX[ indCurSqOx[1] ];
			}
			if( indCurSqOx[1] < 0 )
			{
				trPN[0] = hx * indCurSqOx[1];
			}
			trPN[1] = bv[1] + ang * (trPN[0] - bv[0]);
		}
		//   c. Cheking.
		if( trPN[0] > (uv[0] - 1.e-14) )
		{
			trPN[0] = uv[0];
			trPN[1] = uv[1];
			isTrDone = true;
			wTrPNI = 0;
		}
		//   d. Integration.
		buf_D = integOfChan_SLLeftSd( //   -  The domain is Channel with Slant Line on the left side.
			par_a, //   -  Solution parameter.
			//
			lbDom, rbDom,
			//
			bbDom, ubDom,
			//
			tau, iCurrTL, //   -  Index of current time layer.
			//
			trPC, wTrPCI, //   -  double *bv,
			trPN, wTrPNI, //   -  double *uv,
			//
			indCurSqOx, //   -  Indices where trPC and trPN are.
			//
			uv[0], indRB, //   -  double rb  =  Right boundary by Ox.
			//
			indCurSqOy, //   -  Index of current square by Oy axis.
			//
			masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
			numOfOXSt, //   -  Number of OX steps.
			//
			masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
			numOfOYSt, //   -  Number of OY steps.
			//
			rhoInPrevTL_asV );
		integOfUppTr = integOfUppTr + buf_D;
		//   e. Updating.
		if( isTrDone == false )
		{
			//   We will compute more. We need to redefine some values.
			wTrPCI = wTrPNI;
			trPC[0] = trPN[0];
			trPC[1] = trPN[1];
			if( wTrPNI == 1)
			{
				indCurSqOy[0] += 1;
				indCurSqOy[1] += 1;
			}
			if( wTrPNI == 2)
			{
				indCurSqOx[0] += 1;
				indCurSqOx[1] += 1;
			}
			if( indCurSqOx[1] >= 0)
			{
				distOx = fabs( masOX[ indCurSqOx[1] ] - trPC[0] );
			}
			if( indCurSqOx[1] < 0)
			{
				distOx = fabs( hx * indCurSqOx[1] - trPC[0] );
			}
			if( indCurSqOy[1] >= 0 )
			{
				distOy = fabs( masOY[ indCurSqOy[1] ] - trPC[1] );
			}
			if( indCurSqOy[1] < 0 )
			{
				distOy = fabs( hy * indCurSqOy[1] - trPC[1] );
			}
		}
	} while( !isTrDone );
	return integOfUppTr;
}

double integUnderRigAngTr_UppRight(
	double par_a, //   -  Solution parameter.
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
	double *bv,
	double *uv,
	//
	const double * masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
	int numOfOXSt, //   -  Number of OX steps.
	//
	const double * masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
	int numOfOYSt, //   -  Number of OY steps.
	//
	double * rhoInPrevTL_asV )
{
	//   return ( fabs( (uv[1] - bv[1]) * (bv[0] - uv[0]) /2.) );
	double trPC[2]; //   -  Travel point current;
	int wTrPCI = 0; //   -  Where travel point current is?
	double trPN[2]; //   -  Travel point next;
	int wTrPNI = 0; //   -  Where travel point next is?
	double ang; //   -  Angle of slant line. Should be greater zero.
	int indCurSqOx[2], indCurSqOy[2]; //   -  Index of current square by Ox and Oy axes.
	int indLB[2]; //   -  Index of left boundary.
	double distOx, distOy; //   -  Distance to near Ox and Oy straight lines.
	bool isTrDone = false; //   -  Is travel done.
	double hx = masOX[1] - masOX[0];
	double hy = masOY[1] - masOY[0];
	double integOfUppTr = 0.; //   -  Value which we are computing.
	double buf_D;
	//   Initial data.
	trPC[0] = bv[0];
	trPC[1] = bv[1];
	if( ( fabs(bv[0] - uv[0]) ) < 1.e-12 )
	{
		//   This triangle has very small width. I guess further computation isn't correct.
		return fabs(bv[0] - uv[0]);
	}
	ang = (uv[1] - bv[1]) / (bv[0] - uv[0]);
	if( fabs(ang) < 1.e-12 )
	{
		//   This triangle has very small height. I guess further computation isn't correct.
		return fabs(ang);
	}
	indCurSqOx[0] = (int)( (trPC[0] - 1.e-14) /hx); //   -  If trPC[0] is in grid edge I want it will be between in the left side.
	if( (trPC[0] - 1.e-14) <= 0 )
	{
		indCurSqOx[0] -= 1; //   -  The case when "trPC[0]" ia negative.
	}
	indCurSqOx[1] = indCurSqOx[0] +1; //   -  It's important only in rare case then trPC is in grid edge.
	indLB[0] = (int)( (uv[0] + 1.e-14) /hx);
	if( (uv[0] + 1.e-14) <=0 )
	{
		indLB[0] -= 1; //   -  The case when "trPC[0]" ia negative.
	}
	indLB[1] = indLB[0] +1;
	indCurSqOy[0] = (int)( (trPC[1] + 1.e-14) /hy); //   -  If trPC[1] is in grid edge I want it will be in the upper side.
	if( (trPC[1] + 1.e-14) <= 0 )
	{
		indCurSqOy[0] -= 1; //   -  The case when "trPC[0]" ia negative.
	}
	indCurSqOy[1] = indCurSqOy[0] +1; //   -  It's important only in rare case then trPC is in grid edge.
	if( indCurSqOx[0] >= 0 )
	{
		distOx = fabs( trPC[0] - masOX[ indCurSqOx[0] ] );
	}
	if( indCurSqOx[0] < 0 )
	{
		distOx = fabs( trPC[0] - hx * indCurSqOx[0] );
	}
	if( indCurSqOy[1] >= 0 )
	{
		distOy = fabs( masOY[ indCurSqOy[1] ] - trPC[1] );
	}
	if( indCurSqOy[1] < 0 )
	{
		distOy = fabs( hy * indCurSqOy[1] - trPC[1] );
	}
	do
	{
		//   a. First case.
		if( (distOy /distOx) <= ang )
		{
			//   Across with straight line parallel Ox axis.
			wTrPNI = 1;
			if( indCurSqOy[1] >= 0 )
			{
				trPN[1] = masOY[ indCurSqOy[1] ];
			}
			if( indCurSqOy[1] < 0 )
			{
				trPN[1] = hy * indCurSqOy[1];
			}
			trPN[0] = bv[0] - (trPN[1] - bv[1]) /ang;
		}
		//   b. Second case.
		if( (distOy /distOx) > ang )
		{
			//   Across with straight line parallel Oy axis.
			wTrPNI = 2;
			if( indCurSqOx[0] >= 0 )
			{
				trPN[0] = masOX[ indCurSqOx[0] ];
			}
			if( indCurSqOx[0] < 0 )
			{
				trPN[0] = hx * indCurSqOx[0];
			}
			trPN[1] = bv[1] - ang * (trPN[0] - bv[0]);
		}
		//   c. Cheking.
		if( trPN[0] < (uv[0] + 1.e-14) )
		{
			trPN[0] = uv[0];
			trPN[1] = uv[1];
			isTrDone = true;
			wTrPNI = 0;
		}
		//   d. Integration.
		buf_D = integOfChan_SLRightSd( //   -  The domain is Channel with Slant Line on the Right side.
			par_a, //   -  Solution parameter.
			//
			lbDom, rbDom,
			//
			bbDom, ubDom,
			//
			tau, iCurrTL, //   -  Index of current time layer.
			//
			trPC, wTrPCI, //   -  double *bv,
			trPN, wTrPNI, //   -  double *uv,
			//
			indCurSqOx, //   -  Indices where trPC and trPN are.
			//
			uv[0], indLB, //   -  double lb  =  Left boundary by Ox.
			//
			indCurSqOy, //   -  Index of current square by Oy axis.
			//
			masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
			numOfOXSt, //   -  Number of OX steps.
			//
			masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
			numOfOYSt, //   -  Number of OY steps.
			//
			rhoInPrevTL_asV );
		integOfUppTr = integOfUppTr + buf_D;
		//   e. Updating.
		if( isTrDone == false )
		{
			//   We will compute more. We need to redefine some values.
			wTrPCI = wTrPNI;
			trPC[0] = trPN[0];
			trPC[1] = trPN[1];
			if( wTrPNI == 1)
			{
				indCurSqOy[0] += 1;
				indCurSqOy[1] += 1;
			}
			if( wTrPNI == 2)
			{
				indCurSqOx[0] -= 1;
				indCurSqOx[1] -= 1;
			}
			if( indCurSqOx[0] >= 0 )
			{
				distOx = fabs( trPC[0] - masOX[ indCurSqOx[0] ] );
			}
			if( indCurSqOx[0] < 0 )
			{
				distOx = fabs( trPC[0] - hx * indCurSqOx[0] );
			}
			if( indCurSqOy[1] >= 0 )
			{
				distOy = fabs( masOY[ indCurSqOy[1] ] - trPC[1] );
			}
			if( indCurSqOy[1] < 0 )
			{
				distOy = fabs( hy * indCurSqOy[1] - trPC[1] );
			}
		}
	} while(!isTrDone);
	return integOfUppTr;
}

double integUnderUpperTr(
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
	double * LvUt, //   -  Left, Right and Upper vertices of Upper triangle.
	double * RvUt, //   -  Left, Right and Upper vertices of Upper triangle.
	double * UvUt, //   -  Left, Right and Upper vertices of Upper triangle.
	//
	const double * masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
	int numOfOXSt, //   -  Number of OX steps.
	//
	const double * masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
	int numOfOYSt, //   -  Number of OY step
	//
	double * rhoInPrevTL_asV )
{
	double integOfUppTr;
	double buf_D;
	//   Three ways are possible.
	//   1.
	if( UvUt[0] <= LvUt[0] )
	{
		buf_D = integUnderRigAngTr_UppRight(
			par_a, lbDom, rbDom, bbDom, ubDom, tau, iCurrTL,
			//
			RvUt, UvUt, masOX, numOfOXSt, masOY, numOfOYSt, rhoInPrevTL_asV );
		integOfUppTr = buf_D;
		buf_D = integUnderRigAngTr_UppRight(
			par_a, lbDom, rbDom, bbDom, ubDom, tau, iCurrTL,
			//
			LvUt, UvUt, masOX, numOfOXSt, masOY, numOfOYSt, rhoInPrevTL_asV );
		integOfUppTr = integOfUppTr - buf_D;
		return integOfUppTr;
	}
	//   2.
	if( (UvUt[0] > LvUt[0]) && (UvUt[0] < RvUt[0]) )
	{
		buf_D = integUnderRigAngTr_UppLeft(
			par_a, lbDom, rbDom, bbDom, ubDom, tau, iCurrTL,
			//
			LvUt, UvUt, masOX, numOfOXSt, masOY, numOfOYSt, rhoInPrevTL_asV );
		integOfUppTr = buf_D;

		buf_D = integUnderRigAngTr_UppRight(
			par_a, lbDom, rbDom, bbDom, ubDom, tau, iCurrTL,
			//
			RvUt, UvUt, masOX, numOfOXSt, masOY, numOfOYSt, rhoInPrevTL_asV );
		integOfUppTr = integOfUppTr + buf_D;
		return integOfUppTr;
	}
	//   3.
	if( UvUt[0] >= RvUt[0] )
	{
		buf_D = integUnderRigAngTr_UppLeft(
			par_a, lbDom, rbDom, bbDom, ubDom, tau, iCurrTL,
			//
			LvUt, UvUt, masOX, numOfOXSt, masOY, numOfOYSt, rhoInPrevTL_asV );
		integOfUppTr = buf_D;
		buf_D = integUnderRigAngTr_UppLeft(
			par_a, lbDom, rbDom, bbDom, ubDom, tau, iCurrTL,
			//
			RvUt, UvUt, masOX, numOfOXSt, masOY, numOfOYSt, rhoInPrevTL_asV );
		integOfUppTr = integOfUppTr - buf_D;
		return integOfUppTr;
	}
	return integOfUppTr;
}

double integUnderUnunifTr(
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
	int ii, int jj ) //!!!!!!!!!!!!!!!!!!!
{
	double bv[2], mv[2], uv[2]; //   -  Botton, middle and upper vertices of triangle.
	bool isFirVUsed = false;
	bool isSecVUsed = false;
	bool isThiVUsed = false;
	bool is1VUsed, is2VUsed, is3VUsed;
	double a_LC, b_LC, c_LC; //   -  Coefficients of line betweeen "bv" and "uv" vertices.
	double ap[2]; //   -  Across point of line through "bv" to "uv" and "y == mv[1]"
	double LvBt[2], RvBt[2], BvBt[2]; //   -  Left, Right and Botton vertices of Botton triangle.
	double integOfBottTr; //   -  Item of integral under Botton triangle.
	double LvUt[2], RvUt[2], UvUt[2]; //   -  Left, Right and Upper vertices of Upper triangle.
	double integOfUppTr; //   -  Item of integral under Upper triangle.
	double integ = 0.; //   -  Item which I'm computing.
	//   1. I need to understand which vertex is botton, middle and upper.
	bv[1] = firVer[1];
	bv[0] = firVer[0];
	isFirVUsed = true;
	if( bv[1] > secVer[1] )
	{
		bv[1] = secVer[1];
		bv[0] = secVer[0];
		isFirVUsed = false;
		isSecVUsed = true;
	}
	if( bv[1] > thiVer[1] )
	{
		bv[1] = thiVer[1];
		bv[0] = thiVer[0];
		isFirVUsed = false;
		isSecVUsed = false;
		isThiVUsed = true;
	}
	uv[1] = masOY[0]; //   -  The minimum possible value.
	is1VUsed = false;
	is2VUsed = false;
	is3VUsed = false;
	if( (uv[1] < firVer[1]) && (isFirVUsed == false) )
	{
		uv[1] = firVer[1];
		uv[0] = firVer[0];
		is1VUsed = true;
	}
	if( (uv[1] < secVer[1]) && (isSecVUsed == false) )
	{
		uv[1] = secVer[1];
		uv[0] = secVer[0];
		is2VUsed = true;
		is1VUsed = false;
	}
	if( (uv[1] < thiVer[1]) && (isThiVUsed == false) )
	{
		uv[1] = thiVer[1];
		uv[0] = thiVer[0];
		is3VUsed = true;
		is2VUsed = false;
		is1VUsed = false;
	}
	//   Dangerous.
	if( (isFirVUsed == false) && (is1VUsed == false) )
	{
		mv[1] = firVer[1];
		mv[0] = firVer[0];
	}
	if( (isSecVUsed == false) && (is2VUsed == false) )
	{
		mv[1] = secVer[1];
		mv[0] = secVer[0];
	}
	if( (isThiVUsed == false) && (is3VUsed == false) )
	{
		mv[1] = thiVer[1];
		mv[0] = thiVer[0];
	}
	//   2. I want to compute across point.
	//   2.a Let's compute line coefficients betweeen "bv" and "uv" vertices.
	//   a_LC * x  +  b_LC * y  = c_LC.
	a_LC = uv[1] - bv[1];
	b_LC = bv[0] - uv[0];
	c_LC = (bv[0] - uv[0])*bv[1] + (uv[1] - bv[1])*bv[0];
	//   2.b Across point.
	ap[1] = mv[1];
	if( fabs(a_LC) < 1.e-12 )
	{
		//   This triangle has very small height. I guess further computation isn't correct.
		return 1.e-12;
	}
	ap[0] = (c_LC - b_LC * ap[1]) /a_LC;

	//	printf("i= %d, j= %d : ap[0]= %le      mv[0]= %le \n",ii,jj, ap[0], mv[0]); // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	//   3. There the middle vertex relativly straight line is? Two ways are possible.
	if( mv[0] < ap[0] )
	{
		//   Left, Right and Botton vertices of Botton triangle.
		LvBt[0] = mv[0];
		LvBt[1] = mv[1];
		RvBt[0] = ap[0];
		RvBt[1] = ap[1];
		BvBt[0] = bv[0];
		BvBt[1] = bv[1];
		integOfBottTr = integUnderBottTr(
			par_a, par_b,
			//
			lbDom, rbDom, //   -  Left and right boundaries of rectangular domain.
			//
			bbDom, ubDom, //   -  Botton and upper boundaries of rectangular domain.
			//
			tau, iCurrTL, //   -  Index of current time layer.
			//
			LvBt, RvBt, BvBt, //   -  Left, Right and Botton vertices of Botton triangle.
			//
			masOX, numOfOXSt, //   -  Number of OX steps.
			//
			masOY, numOfOYSt, //   -  Number of OY steps.
			//
			rhoInPrevTL_asV,
			ii, jj ); // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		integ = integOfBottTr;

		//		printf("m<a:   i= %d, j= %d : integ= %le \n",ii,jj, integ); // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		//   Left, Right and Upper vertices of Upper triangle.
		LvUt[0] = mv[0];
		LvUt[1] = mv[1];
		RvUt[0] = ap[0];
		RvUt[1] = ap[1];
		UvUt[0] = uv[0];
		UvUt[1] = uv[1];
		integOfUppTr = integUnderUpperTr(
			par_a, par_b,
			//
			lbDom, rbDom, //   -  Left and right boundaries of rectangular domain.
			//
			bbDom, ubDom, //   -  Botton and upper boundaries of rectangular domain.
			//
			tau, iCurrTL, //   -  Index of current time layer.
			//
			LvUt, RvUt, UvUt, //   -  Left, Right and Botton vertices of Upper triangle.
			//
			masOX, numOfOXSt, //   -  Number of OX steps.
			//
			masOY, numOfOYSt, //   -  Number of OY steps.
			//
			rhoInPrevTL_asV );
		integ = integ + integOfUppTr;
		return integ;
	}
	if( mv[0] >= ap[0] )
	{
		//   Left, Right and Botton vertices of Botton triangle.
		LvBt[0] = ap[0];
		LvBt[1] = ap[1];
		RvBt[0] = mv[0];
		RvBt[1] = mv[1];
		BvBt[0] = bv[0];
		BvBt[1] = bv[1];
		integOfBottTr = integUnderBottTr(
			par_a, par_b,
			//
			lbDom, rbDom, //   -  Left and right boundaries of rectangular domain.
			//
			bbDom, ubDom, //   -  Botton and upper boundaries of rectangular domain.
			//
			tau, iCurrTL, //   -  Index of current time layer.
			//
			LvBt, RvBt, BvBt, //   -  Left, Right and Botton vertices of Botton triangle.
			//
			masOX, numOfOXSt, //   -  Number of OX steps.
			//
			masOY, numOfOYSt, //   -  Number of OY steps.
			//
			rhoInPrevTL_asV,
			ii, jj ); // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		integ = integOfBottTr;

		//		printf("m>a:   i= %d, j= %d : integ= %le \n",ii,jj, integ);  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		//   Left, Right and Upper vertices of Upper triangle.
		LvUt[0] = ap[0];
		LvUt[1] = ap[1];
		RvUt[0] = mv[0];
		RvUt[1] = mv[1];
		UvUt[0] = uv[0];
		UvUt[1] = uv[1];
		integOfUppTr = integUnderUpperTr(
			par_a, par_b,
			//
			lbDom, rbDom, //   -  Left and right boundaries of rectangular domain.
			//
			bbDom, ubDom, //   -  Botton and upper boundaries of rectangular domain.
			//
			tau, iCurrTL, //   -  Index of current time layer.
			//
			LvUt, RvUt, UvUt, //   -  Left, Right and Botton vertices of Upper triangle.
			//
			masOX, numOfOXSt, //   -  Number of OX steps.
			//
			masOY, numOfOYSt, //   -  Number of OY steps.
			//
			rhoInPrevTL_asV );
		return integ + integOfUppTr;
	}
	return integ;
}

double bottonBound(
	double par_a, //   -  Item of left and right setback (parameter "a" in test).
	//
	double lbDom, //   -  Left and right boundaries of rectangular domain.
	double rbDom,
	//
	double bbDom, //   -  Botton and upper boundaries of rectangular domain.
	double ubDom,
	//
	double t,
	double x )
{
	return analytSolut( par_a, lbDom, rbDom, bbDom, ubDom, t, x, bbDom );
}

double u_function(double par_b, double t, double x, double y)
{
	return par_b * y * (1.-y) * ( C_pi_device /2. + atan( -x ) );
}

double v_function(
	double lbDom, //   -  Left and right boundaries of rectangular domain.
	double rbDom,
	//
	double bbDom, //   -  Botton and upper boundaries of rectangular domain.
	double ubDom,
	//
	double t, double x, double y )
{
	return atan((x - lbDom) * (x - rbDom) * (1.+t) /10. * (y - ubDom) * (y - bbDom));
}

double f_function( //   -  It's item of right part of differential equation.
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
	int iOfOXN, //   -  Index of current OX node.
	const double *masOX, //   -  Massive of OX steps. Dimension = numOfOXSt
	int numOfOXSt, //   -  Number of OX steps (segments).
	//
	int iOfOYN, //   -  Index of current OY node.
	const double *masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
	int numOfOYSt ) //   -  Number of OY steps (segments).
{
	double t = tau * iCurrTL;
	double x = masOX[ iOfOXN ];
	double y = masOY[ iOfOYN ];
	double arg_v = (x - lbDom) * (x - rbDom) * (1.+t) /10. * (y - ubDom) * (y - bbDom);
	double rho, dRhoDT, dRhoDX, dRhoDY;
	double u, duDX;
	double v, dvDY;
	rho = analytSolut(par_a, lbDom, rbDom, bbDom, ubDom, t, x, y );
	dRhoDT = x * y * cos( t*x*y );
	dRhoDX = t * y * cos( t*x*y );
	dRhoDY = t * x * cos( t*x*y );
	u = u_function(par_b, t, x, y );
	duDX = -par_b * y * (1.-y) / ( 1. + x * x );
	v = v_function(lbDom, rbDom, bbDom, ubDom, t, x, y );
	dvDY = (x - lbDom) * (x - rbDom) * (1.+t) /10. * (y - bbDom + y - ubDom);
	dvDY = dvDY / ( 1. + arg_v * arg_v );
	return dRhoDT + rho * duDX + u * dRhoDX + rho * dvDY + v * dRhoDY;
}

int quadrAngleType(
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
	//
	double * firVfirT, //   -  First vertex of first triangle.
	double * secVfirT, //   -  Second vertex of first triangle.
	double * thiVfirT, //   -  Third vertex of first triangle.
	//
	double * firVsecT, //   -  First vertex of second triangle.
	double * secVsecT, //   -  Second vertex of second triangle.
	double * thiVsecT ) //   -  Third vertex of second triangle.

{
	int qAngType = 0;

	double alpha[2], betta[2], gamma[2], theta[2]; //   -  Vertexes of square. Anticlockwise order from left botton vertex.
	double u, v; //   -  Items of velocity components.
	double alNew[2], beNew[2], gaNew[2], thNew[2]; //   -  New positions of vertexes. Vertexes of quadrangle.
	double vectAlGa[2], vectBeTh[2]; //   -  Vectors: 1) from "alNew" to "gaNew" 2) from "beNew" to "thNew".
	double a_1LC, b_1LC, c_1LC; //   -  a_1LC * x  +  b_1LC * y  = c_1LC. Line equation through "alNew" and "gaNew".
	double a_2LC, b_2LC, c_2LC; //   -  a_2LC * x  +  b_2LC * y  = c_2LC. Line equation through "beNew" and "thNew".
	double AcrP[2]; //   -  Across point of two lines
	double vectAlBe[2]; //   -  Vectors for computing vertexes sequence order by vector production.
	double vectAlTh[2]; //   -  Vectors for computing vertexes sequence order by vector production.
	double vectBeGa[2]; //   -  Vectors for computing vertexes sequence order by vector production.
	double vectBeAl[2]; //   -  Vectors for computing vertexes sequence order by vector production.
	double vectProdOz; //   -  Z-axis of vector production.
	double scalProd; //   -  Scalar production of two vectors.

	//   1. First of all let's compute coordinates of square vertexes.

	//  OX:

	if( iOfOXN == 0 )
	{
		alpha[0] = masOX[ iOfOXN ];
		betta[0] = ( masOX[iOfOXN] + masOX[iOfOXN +1] ) /2.;
		gamma[0] = ( masOX[iOfOXN] + masOX[iOfOXN +1] ) /2.;
		theta[0] = masOX[ iOfOXN ];
	}

	if( iOfOXN == numOfOXSt )
	{
		alpha[0] = ( masOX[iOfOXN -1] + masOX[iOfOXN] ) /2.;
		betta[0] = masOX[ iOfOXN ];
		gamma[0] = masOX[ iOfOXN ];
		theta[0] = ( masOX[iOfOXN -1] + masOX[iOfOXN] ) /2.;
	}

	if( (iOfOXN > 0) && (iOfOXN < numOfOXSt) )
	{
		alpha[0] = ( masOX[iOfOXN -1] + masOX[iOfOXN] ) /2.;
		betta[0] = ( masOX[iOfOXN +1] + masOX[iOfOXN] ) /2.;
		gamma[0] = ( masOX[iOfOXN +1] + masOX[iOfOXN] ) /2.;
		theta[0] = ( masOX[iOfOXN -1] + masOX[iOfOXN] ) /2.;
	}

	//  OY:
	if( iOfOYN == 0 )
	{
		alpha[1] = masOY[ iOfOYN ];
		betta[1] = masOY[ iOfOYN ];
		gamma[1] = ( masOY[iOfOYN] + masOY[ iOfOYN +1] ) /2.;
		theta[1] = ( masOY[iOfOYN] + masOY[ iOfOYN +1] ) /2.;
	}

	if( iOfOYN == numOfOYSt )
	{
		alpha[1] = ( masOY[iOfOYN] + masOY[ iOfOYN -1] ) /2.;
		betta[1] = ( masOY[iOfOYN] + masOY[ iOfOYN -1] ) /2.;
		gamma[1] = masOY[ iOfOYN ];
		theta[1] = masOY[ iOfOYN ];
	}

	if( (iOfOYN > 0) && (iOfOYN < numOfOYSt) )
	{
		alpha[1] = ( masOY[iOfOYN] + masOY[ iOfOYN -1] ) /2.;
		betta[1] = ( masOY[iOfOYN] + masOY[ iOfOYN -1] ) /2.;
		gamma[1] = ( masOY[iOfOYN] + masOY[ iOfOYN +1] ) /2.;
		theta[1] = ( masOY[iOfOYN] + masOY[ iOfOYN +1] ) /2.;
	}

	//   2. Now let's compute new coordinates on the previous time level of alpha, betta, gamma, theta points.
	//  alNew.

	u = u_function( par_b, iCurrTL * tau, alpha[0], alpha[1] );
	v = v_function( lbDom, rbDom, bbDom, ubDom, iCurrTL * tau, alpha[0], alpha[1] );
	alNew[0] = alpha[0] - tau * u;
	alNew[1] = alpha[1] - tau * v;

	//  beNew.
	u = u_function( par_b, iCurrTL * tau, betta[0], betta[1] );
	v = v_function( lbDom, rbDom, bbDom, ubDom, iCurrTL * tau, betta[0], betta[1] );
	beNew[0] = betta[0] - tau * u;
	beNew[1] = betta[1] - tau * v;
	
	//  gaNew.

	u = u_function( par_b, iCurrTL * tau, gamma[0], gamma[1] );
	v = v_function( lbDom, rbDom, bbDom, ubDom, iCurrTL * tau, gamma[0], gamma[1] );
	gaNew[0] = gamma[0] - tau * u;
	gaNew[1] = gamma[1] - tau * v;

	//  thNew.
	u = u_function( par_b, iCurrTL * tau, theta[0], theta[1] );
	v = v_function( lbDom, rbDom, bbDom, ubDom, iCurrTL * tau, theta[0], theta[1] );
	thNew[0] = theta[0] - tau * u;
	thNew[1] = theta[1] - tau * v;

	//   3.a Let's compute coefficients of first line betweeen "alNew" and "gaNew" points.
	//   a_1LC * x  +  b_1LC * y  = c_1LC.

	vectAlGa[0] = gaNew[0] - alNew[0];
	vectAlGa[1] = gaNew[1] - alNew[1];
	a_1LC = vectAlGa[1];
	b_1LC = -vectAlGa[0];
	c_1LC = vectAlGa[1] * alNew[0] - vectAlGa[0] * alNew[1];

	//   3.b Let's compute coefficients of second line betweeen "beNew" and "thNew" points.
	//   a_2LC * x  +  b_2LC * y  = c_2LC.

	vectBeTh[0] = thNew[0] - beNew[0];
	vectBeTh[1] = thNew[1] - beNew[1];
	a_2LC = vectBeTh[1];
	b_2LC = -vectBeTh[0];
	c_2LC = vectBeTh[1] * beNew[0] - vectBeTh[0] * beNew[1];

	//   4. Let's compute coordinates of across point of this two lines.
	//   Are lines parallel?

	if( fabs(b_1LC*a_2LC - b_2LC*a_1LC) < 1.e-14 )
	{
		//   Not checked.

		qAngType = 0.;

		//   Pseudo case. Anyway I need to compute some values.

		//   First triangle.

		firVfirT[0] = alNew[0];
		firVfirT[1] = alNew[1];
		secVfirT[0] = beNew[0];
		secVfirT[1] = beNew[1];
		thiVfirT[0] = gaNew[0];
		thiVfirT[1] = gaNew[1];


		//   Vertices of second triagle depends on scalar production.

		vectAlGa[0] = gaNew[0] - alNew[0];
		vectAlGa[1] = gaNew[1] - alNew[1];
		vectBeTh[0] = thNew[0] - beNew[0];
		vectBeTh[1] = thNew[1] - beNew[1];


		scalProd = vectAlGa[0] * vectBeTh[0] + vectAlGa[1] * vectBeTh[1];
		firVsecT[0] = beNew[0];
		firVsecT[1] = beNew[1];
		secVsecT[0] = thNew[0];
		secVsecT[1] = thNew[1];


		if( scalProd >= 0. )
		{
			thiVsecT[0] = gaNew[0];
			thiVsecT[1] = gaNew[1];
		}

		if( scalProd < 0. )
		{
			thiVsecT[0] = alNew[0];
			thiVsecT[1] = alNew[1];
		}

		return qAngType;
	}


	AcrP[0] = ( b_1LC*c_2LC - b_2LC*c_1LC ) / ( b_1LC*a_2LC - b_2LC*a_1LC );
	AcrP[1] = ( a_1LC*c_2LC - a_2LC*c_1LC ) / (-b_1LC*a_2LC + b_2LC*a_1LC );

	if ( ( (beNew[1] - AcrP[1])*(thNew[1] - AcrP[1]) ) > 0. )
	{
		if( ( (alNew[0] - AcrP[0])*(gaNew[0] - AcrP[0]) ) > 0. )
		{
			qAngType = 0;
			firVfirT[0] = alNew[0];
			firVfirT[1] = alNew[1];
			secVfirT[0] = beNew[0];
			secVfirT[1] = beNew[1];
			thiVfirT[0] = gaNew[0];
			thiVfirT[1] = gaNew[1];

			//   Second triangle.

			firVsecT[0] = beNew[0];
			firVsecT[1] = beNew[1];
			secVsecT[0] = thNew[0];
			secVsecT[1] = thNew[1];


			//   Third vertex computing...

			vectAlGa[0] = gaNew[0] - alNew[0];
			vectAlGa[1] = gaNew[1] - alNew[1];

			vectBeTh[0] = thNew[0] - beNew[0];
			vectBeTh[1] = thNew[1] - beNew[1];

			scalProd = vectAlGa[0] * vectBeTh[0] + vectAlGa[1] * vectBeTh[1];

			if( scalProd >= 0. )
			{
				thiVsecT[0] = gaNew[0];
				thiVsecT[1] = gaNew[1];
			}

			if( scalProd < 0. )
			{
				thiVsecT[0] = alNew[0];
				thiVsecT[1] = alNew[1];
			}

			return qAngType;
		} //   "if(  (  (alNew[0] - AcrP[0])*(gaNew[0] - AcsP[0])  )  >  0.  )".


		//   Second criterion. Second case.

		if( ( (alNew[0] - AcrP[0])*(gaNew[0] - AcrP[0]) ) <= 0. )
		{
			vectAlBe[0] = beNew[0] - alNew[0];
			vectAlBe[1] = beNew[1] - alNew[1];
			vectAlTh[0] = thNew[0] - alNew[0];
			vectAlTh[1] = thNew[1] - alNew[1];

			vectProdOz = vectAlBe[0] * vectAlTh[1] - vectAlBe[1] * vectAlTh[0];

			if( vectProdOz < 0. )
			{
				//   The vertex "beNew" is NOT in triangle "alNew - gaNew - thNew".

				qAngType = 0;

				//   Pseudo case. Anyway I need to find some solution. So

				firVfirT[0] = alNew[0];
				firVfirT[1] = alNew[1];
				secVfirT[0] = beNew[0];
				secVfirT[1] = beNew[1];
				thiVfirT[0] = thNew[0];
				thiVfirT[1] = thNew[1];

				//   Second triangle.

				firVsecT[0] = beNew[0];
				firVsecT[1] = beNew[1];
				secVsecT[0] = thNew[0];
				secVsecT[1] = thNew[1];
				thiVsecT[0] = gaNew[0];
				thiVsecT[1] = gaNew[1];

				return qAngType;
			}

			if( vectProdOz >= 0. )
			{
				//  It's all write. We have a good concave quadrangle.

				//   Now let's compute all vertices which I need.

				qAngType = 2;

				//   First triangle.

				firVfirT[0] = alNew[0];
				firVfirT[1] = alNew[1];
				secVfirT[0] = beNew[0];
				secVfirT[1] = beNew[1];
				thiVfirT[0] = thNew[0];
				thiVfirT[1] = thNew[1];

				//   Second triangle.

				firVsecT[0] = beNew[0];
				firVsecT[1] = beNew[1];
				secVsecT[0] = thNew[0];
				secVsecT[1] = thNew[1];
				thiVsecT[0] = gaNew[0];
				thiVsecT[1] = gaNew[1];

				return qAngType;
			}
		} //   "if(  (  (alNew[0] - AcsP[0])*(gaNew[0] - AcsP[0])  )  <=  0.  )".   //   Last second case of second criterion.
	} //   end of "if (  (  (beNew[1] - AcrP[1])*(thNew[1] - AcrP[1])  )  >  0.  )"

	//  Now let's consider SECOND case 5.b "(  (beNew[1] - AcrP[1])*(thNew[1] - AcrP[1])  )  <= 0."

	if ( ( (beNew[1] - AcrP[1])*(thNew[1] - AcrP[1]) ) <= 0. )
	{
		if( ( (alNew[0] - AcrP[0])*(gaNew[0] - AcrP[0]) ) > 0. )
		{
			//  It means the across point IS NOT between "alNew" and "gaNew" vertices by Ox-axis?

			//   O.K. the quadrangle IS NOT CONVEX. Is it concave or pseudo? Third criterion.

			vectBeGa[0] = gaNew[0] - beNew[0];
			vectBeGa[1] = gaNew[1] - beNew[1];
			vectBeAl[0] = alNew[0] - beNew[0];
			vectBeAl[1] = alNew[1] - beNew[1];

			vectProdOz = vectBeGa[0] * vectBeAl[1] - vectBeGa[1] * vectBeAl[0];

			if( vectProdOz >= 0. )
			{
				qAngType = 2;

				//   The quadrangle is concave. First triangle.

				firVfirT[0] = alNew[0];
				firVfirT[1] = alNew[1];
				secVfirT[0] = beNew[0];
				secVfirT[1] = beNew[1];
				thiVfirT[0] = gaNew[0];
				thiVfirT[1] = gaNew[1];

				//   Second triangle.

				firVsecT[0] = alNew[0];
				firVsecT[1] = alNew[1];

				secVsecT[0] = thNew[0];
				secVsecT[1] = thNew[1];

				thiVsecT[0] = gaNew[0];
				thiVsecT[1] = gaNew[1];

				return qAngType;
			}

			if( vectProdOz < 0. )
			{
				qAngType = 0;

				//   This concave quadrangle do has NO write anticlockwise vertices sequence order. It's pseudo.

				firVfirT[0] = alNew[0];
				firVfirT[1] = alNew[1];

				secVfirT[0] = beNew[0];
				secVfirT[1] = beNew[1];

				thiVfirT[0] = gaNew[0];
				thiVfirT[1] = gaNew[1];

				//   Second triangle.

				firVsecT[0] = alNew[0];
				firVsecT[1] = alNew[1];

				secVsecT[0] = thNew[0];
				secVsecT[1] = thNew[1];

				thiVsecT[0] = gaNew[0];
				thiVsecT[1] = gaNew[1];

				return qAngType;
			}
		} //   end of "if(  (  (alNew[0] - AcrP[0])*(gaNew[0] - AcsP[0])  )  >  0.  )". First case of second criterion.


		//   Second criterion. Second case.

		if( ( (alNew[0] - AcrP[0])*(gaNew[0] - AcrP[0]) ) <= 0. )
		{
			//   O.K. the quadrangle is convex. Is it has the same vertices sequence order.

			vectAlBe[0] = beNew[0] - alNew[0];

			vectAlBe[1] = beNew[1] - alNew[1];

			vectAlTh[0] = thNew[0] - alNew[0];

			vectAlTh[1] = thNew[1] - alNew[1];

			vectProdOz = vectAlBe[0] * vectAlTh[1] - vectAlBe[1] * vectAlTh[0];

			if( vectProdOz >= 0. )
			{
				qAngType = 1;

				//   Convex quadrangle DO HAS WRITE anticlockwise vertices sequence order. It's convex.

				firVfirT[0] = alNew[0];
				firVfirT[1] = alNew[1];

				secVfirT[0] = beNew[0];
				secVfirT[1] = beNew[1];

				thiVfirT[0] = gaNew[0];
				thiVfirT[1] = gaNew[1];

				//   Second triangle.

				firVsecT[0] = alNew[0];
				firVsecT[1] = alNew[1];

				secVsecT[0] = thNew[0];
				secVsecT[1] = thNew[1];

				thiVsecT[0] = gaNew[0];
				thiVsecT[1] = gaNew[1];

				return qAngType;
			}

			if( vectProdOz < 0. )
			{
				qAngType = 0;
				firVfirT[0] = alNew[0];
				firVfirT[1] = alNew[1];
				secVfirT[0] = beNew[0];
				secVfirT[1] = beNew[1];
				thiVfirT[0] = gaNew[0];
				thiVfirT[1] = gaNew[1];
				firVsecT[0] = alNew[0];
				firVsecT[1] = alNew[1];
				secVsecT[0] = thNew[0];
				secVsecT[1] = thNew[1];
				thiVsecT[0] = gaNew[0];
				thiVsecT[1] = gaNew[1];
				return qAngType;
			}
		} //   end of "if(  (  (alNew[0] - AcrP[0])*(gaNew[0] - AcsP[0])  )  <=  0.  )". //   Second case of second criterion.
	}
	return qAngType;
}

double normOfMatrAtL1_asV(
	const double *masOX, //   -  Massive of OX grid nodes. Dimension = dimOX.
	int dimOX,
	//
	const double *masOY, //   -  Massive of OY grid nodes. Dimension = dimOY.
	int dimOY,
	//
	double * mat_asV )
{
	double norm = 0.;
	double hx = masOX[1] - masOX[0];
	double hy = masOY[1] - masOY[0];

	for( int k=1; k< dimOY -1; k++ )
	{
		for( int j=1; j< dimOX -1; j++ )
		{
			norm += fabs( mat_asV[ dimOX*k + j ] );
		}
	}
	return hx * hy * norm;
}

double spaceVolumeInPrevTL(
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
	double * rhoInPrevTL_asV )
{
	double firVfirT[2], secVfirT[2], thiVfirT[2]; //   -  First, second and third vertices of first triangle.
	double firVsecT[2], secVsecT[2], thiVsecT[2]; //   -  First, second and third vertices of second triangle.
	int qAngType; //   -  Type of quadrangle: 0 - pseudo; 1 - convex; 2 - concave;
	double buf_D;
	//   Let's understand what type of quadrangle we have.
	qAngType = quadrAngleType(
		par_a, par_b, //   -  Analitycal solution parameters.
		//
		lbDom, rbDom, //   -  Left and right boundaries of rectangular domain.
		//
		bbDom, ubDom, //   -  Botton and upper boundaries of rectangular domain.
		//
		tau, iCurrTL, //   -  Time data. Necessary for velocity.
		//
		iOfOXN, masOX, numOfOXSt, //   -  OX data.
		//
		iOfOYN, masOY, numOfOYSt, //   -  OY data.
		//
		firVfirT, secVfirT, thiVfirT, //   -  Vertices of first triangle.
		//
		firVsecT, secVsecT, thiVsecT ); //   -  Vertices of second triangle.

	//	printf("i= %d, j=%d, 1: %le    2: %le    3: %le \n", iOfOXN, iOfOYN, firVfirT[1], secVfirT[1], thiVfirT[1]);

	if( qAngType != 1 ) return -1.;

	buf_D = integUnderUnunifTr(
		par_a, par_b, //   -  Analitycal solution parameters.
		//
		lbDom, rbDom, //   -  Left and right boundaries of rectangular domain.
		//
		bbDom, ubDom, //   -  Botton and upper boundaries of rectangular domain.
		//
		tau, iCurrTL, //   -  Index of current time layer.
		//
		firVfirT, secVfirT, thiVfirT, //   -  Vertices of first triangle.
		//
		masOX, numOfOXSt, //   -  Number of OX steps.
		//
		masOY, numOfOYSt, //   -  Number of OY steps.
		rhoInPrevTL_asV,
		iOfOXN, iOfOYN );

	//	printf("i= %d, j=%d, sq: %le \n", iOfOXN, iOfOYN, buf_D);

	//	printf("********************************************************************* \n");

	return buf_D + integUnderUnunifTr(
		par_a, par_b, //   -  Analitycal solution parameters.
		//
		lbDom, rbDom, //   -  Left and right boundaries of rectangular domain.
		//
		bbDom, ubDom, //   -  Botton and upper boundaries of rectangular domain.
		//
		tau, iCurrTL, //   -  Index of current time layer.
		//
		firVsecT, secVsecT, thiVsecT, //   -  Vertices of first triangle.
		//
		masOX, numOfOXSt, //   -  Number of OX steps.
		//
		masOY, numOfOYSt, //   -  Number of OY steps.
		//
		rhoInPrevTL_asV,
		iOfOXN, iOfOYN );
}

double solByEqualVolumes(
	double par_a, //   -  Item of first initial or boundary data parameter.
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
	double *rhoInCurrTL_asV ) //   -  Rho (solution) in Current Time Layer which we will compute.
{
	double *rhoInPrevTL_asV; //   -  Rho (solution) in Previous Time Layer which we have computed.
	double spVolInPrevTL; //   -  Space volume of rho in previous time layer.
	double RPInCurrTL; //   -  Right part in current time layer.
	int iCurrTL; //   -  Index of current time Layer IN WHICH we computing solution.
	int iOfOXN; //   -  Index of current OX node.
	int iOfOYN; //   -  Index of current OY node.
	int iOfThr; //   -  Through OXY plane index.
	double buf_D; //   -  Buffer. Only.
	int j,k;


	//   New memory.
	rhoInPrevTL_asV = new double [ (numOfOXSt +1) * (numOfOYSt +1) ];
	//   Initial data of rho.
	for( k = 0; k < numOfOYSt + 1; k++ )
	{
		for( j = 0; j < numOfOXSt + 1; j++ )
		{
			rhoInPrevTL_asV[ (numOfOXSt +1)*k + j ] = initDataOfSol( par_a, lbDom, rbDom, bbDom, ubDom, j, masOX, k, masOY );
		}
	}

	bool canPrint = true;
	if (canPrint)
	{
		cout<<"\r \t \t \t \t \t \t \t \t \r";
		cout << "SchOrd = " << numOfSolOrd << ", Nx = " << numOfOXSt;
		cout << ", count of time levels " << numOfTSt;
	}

	for( iCurrTL = 1; iCurrTL < numOfTSt + 1; iCurrTL++ )
	{
		if (canPrint)
		{
			cout<<"\r \t \t \t \t \t \t \t \t \r";
			cout << "SchOrd = " << numOfSolOrd << ", Nx = " << numOfOXSt;
			cout << ", indexOfCurrTL = " << iCurrTL << " ( " << numOfTSt << " ) " << flush;
		}


		//   If we know solution on the boundary we can use it.
		for( iOfOXN = 0; iOfOXN < numOfOXSt + 1; iOfOXN++ )
		{
			//   Bottom boundary.
			rhoInCurrTL_asV[ iOfOXN ] = bottonBound( par_a, lbDom, rbDom, bbDom, ubDom, tau*iCurrTL, masOX[ iOfOXN ] );

			//   Upper boundary.
			rhoInCurrTL_asV[ (numOfOXSt + 1)*numOfOYSt + iOfOXN ] = upperBound( par_a, lbDom, rbDom, bbDom, ubDom, tau*iCurrTL, masOX[ iOfOXN ] );
		}

		for( iOfOYN = 0; iOfOYN < numOfOYSt + 1; iOfOYN++ )
		{
			rhoInCurrTL_asV[ (numOfOXSt + 1)*iOfOYN ] = leftBound( par_a, lbDom, rbDom, bbDom, ubDom, tau*iCurrTL, masOY[ iOfOYN ] );
			rhoInCurrTL_asV[ (numOfOXSt + 1)*iOfOYN + numOfOXSt ] = rightBound( par_a, lbDom, rbDom, bbDom, ubDom, tau*iCurrTL, masOY[ iOfOYN ] );
		}


		//   Enumeration from first unknown element to last one.
		for( iOfOYN = 1; iOfOYN < numOfOYSt; iOfOYN++ )
		{
			for( iOfOXN = 1; iOfOXN < numOfOXSt; iOfOXN++ )
			{
				spVolInPrevTL = spaceVolumeInPrevTL(
					par_a, par_b, //   -  Items of parameters.
					//
					lbDom, rbDom, //   -  Left and right boundaries of rectangular domain.
					//
					bbDom, ubDom, //   -  Bottom and upper boundaries of rectangular domain.
					//
					tau, iCurrTL, //   -  Time data. Necessary for velocity.
					//
					iOfOXN, masOX, numOfOXSt, //   -  OX data.
					//
					iOfOYN, masOY, numOfOYSt, //   -  OY data.
					rhoInPrevTL_asV );

				buf_D = (masOX[iOfOXN +1] - masOX[iOfOXN -1]) /2.;
				spVolInPrevTL = spVolInPrevTL / buf_D;

				buf_D = (masOY[iOfOYN +1] - masOY[iOfOYN -1]) /2.;
				spVolInPrevTL = spVolInPrevTL / buf_D;

				RPInCurrTL = f_function(par_a, par_b, lbDom, rbDom, bbDom, ubDom, tau, iCurrTL, iOfOXN,
				                        masOX, //   -  Massive of OX steps. Dimension = numOfOXSt +1.
				                        numOfOXSt, //   -  Number of OX steps.
				                        //
				                        iOfOYN, //   -  Index of current OY node.
				                        masOY, //   -  Massive of OY steps. Dimension = numOfOYSt +1.
				                        numOfOYSt ); //   -  Number of OY steps.

				rhoInCurrTL_asV[ (numOfOXSt + 1)*iOfOYN + iOfOXN ] = spVolInPrevTL;
				rhoInCurrTL_asV[ (numOfOXSt + 1)*iOfOYN + iOfOXN ] += tau * RPInCurrTL;
			}
		}

		for( iOfThr = 0; iOfThr < (numOfOXSt +1) * (numOfOYSt +1); iOfThr++ )
			rhoInPrevTL_asV[ iOfThr ] = rhoInCurrTL_asV[ iOfThr ];
	}


	delete[] rhoInPrevTL_asV;

	return 0;
}

double* solve_cpu_test(
	double par_a,
	double par_b,
	//
	double lbDom,
	double rbDom,
	//
	double bbDom,
	double ubDom,
	//
	double tau,
	int numOfTSt,
	//
	double * masOX,
	int numOfOXSt,
	//
	double * masOY,
	int numOfOYSt,
	int gridStep)
{
	double varTau = tau;
	int varNumOfTSt = numOfTSt;
	double * varMasOX = NULL;
	int varNumOfOXSt = numOfOXSt;
	double * varMasOY = NULL;
	int varNumOfOYSt = numOfOYSt;
	double *rhoInCurrTL_asV = NULL;
	double anSol;
	
	double buf_D;
	int j,k;
	bool bul;

	//   New time step.
	varTau = tau / pow(2., gridStep);
	varNumOfTSt = numOfTSt * pow(2., gridStep);
	//   New absciss grid steps.
	varNumOfOXSt = numOfOXSt * pow(2., gridStep);
	//   New absciss grid.
	varMasOX = new double [ varNumOfOXSt +1 ];
	buf_D = (rbDom - lbDom) / varNumOfOXSt;
	for( j=0; j< varNumOfOXSt+1; j++ )
	{
		varMasOX[j] = lbDom + ( (double)j ) * buf_D;
	}
	//   New ordinate grid steps.
	varNumOfOYSt = numOfOYSt * pow(2., gridStep);
	//   New absciss grid.
	varMasOY = new double [ varNumOfOYSt +1 ];
	buf_D = (ubDom - bbDom) / varNumOfOYSt;
	for( k=0; k< varNumOfOYSt+1; k++ )
	{
		varMasOY[k] = bbDom + ( (double)k ) * buf_D;
	}
	rhoInCurrTL_asV = new double [ (varNumOfOXSt +1) * (varNumOfOYSt +1) ];

	//   Computation of solution.
	solByEqualVolumes(
		par_a, par_b,
		//
		lbDom, rbDom,
		//
		bbDom, ubDom,
		//
		varTau, //   -  Time step.
		varNumOfTSt, //   -  A number of time steps.
		//
		varMasOX, //   -  Massive of abscissa grid points. Dimension = varNumOfOxGrSt +1.
		varNumOfOXSt, //   -  Variable number of abscissa grid steps.
		//
		varMasOY, //   -  Massive of ordinate grid points. Dimension = varNumOfOyGrSt +1.
		varNumOfOYSt, //   -  Variable number of ordinate grid steps.
		//
		1, //   -  For print only. Solution order which we want to get.
		//
		rhoInCurrTL_asV ); //   -  Rho (solution) in Current (Last) Time Level.

	
	delete[] varMasOX;
	delete[] varMasOY;
	cout << std::endl;
	return rhoInCurrTL_asV;
}