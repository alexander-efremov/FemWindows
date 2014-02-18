#include <algorithm>
#include <iostream>
#include "Const.h" //   -  Helpful used consts.
#include "OriginPrint.h" //   -  Printing to file for.
#include "SpecialPrint.h" //   -  Printing to file for.
#include "LowOrdOper.h" //   -  Operator on low order accuracy.
#include "LowOrdOperCuda.h" //   -  Operator on low order accuracy CUDA

using namespace std;

int menu()
{
	int maxValue = 2;
	int choice = -1;

	cout << "Input your choice: " << endl;

	while (true)
	{
		cout << "1. Fem CPU" << endl;
		cout << "2. Fem GPU" << endl;
		cout << "0. Exit" << endl;
		cout << ">> ";
		cin >> choice;
		if (choice >= 0 && choice <= maxValue )
		{
			cout << endl;
			return choice;
		}
		else
		{
			cout << "Unknown value" << endl;
			cout << endl;
		}
	}
	return 0;
}


int main(int argc, char *argv[])
{
	cout << "Welcome to Git" << endl;
	int numOfGrStepLayer = 1;
	while (true)
	{
		int ch = menu();
		switch (ch)
		{
			case 1:
				cout << "Compute on CPU" << endl;
				initCompOfGlVar();
				solve_cpu_test(C_par_a, C_par_b, C_lbDom, C_rbDom, C_bbDom, C_ubDom, C_tau, C_numOfTSt, masOX, C_numOfOXSt, masOY, C_numOfOYSt, numOfGrStepLayer);
				memClean();
				break;
			case 2:
				cout << "Compute on GPU" << endl;

				break;
			default:
				return 0;
		}
	}

	return 0;
}