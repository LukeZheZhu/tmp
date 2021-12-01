#include <iostream>
#include <windows.h>
#include "original.h"
#include "pre_op.h"
using namespace std;

void original(double **X, double **Y) {
	double **Z;
	int i, j, k;
	DWORD t_s, t_run;

	t_run = 0;
	Z = (double**)malloc(m * sizeof(double *));
	for (i = 0; i<m; i++) {
		Z[i] = (double *)malloc(s * sizeof(double));
		memset(Z[i], 0, s * sizeof(double));
	}

	t_s = ::GetTickCount();
	for (i = 0; i < m; i++) {
		for (j = 0; j < s; j++) {
			for (k = 0; k < n; k++) {
				Z[i][j] += X[i][k] * Y[k][j];
			}
		}
	}
	t_run = ::GetTickCount() - t_s;
	cout << "Original:" << t_run << "ms" << endl;
	//mat_show(Z, m, s);
	/*for (i = 0; i < m; i++) { delete[] Z[i];}
	delete[] Z;*/
}