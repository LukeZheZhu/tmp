#include <iostream>
#include <cstring>
#include "original.h"
#include "pre_op.h"
using namespace std;

void original(double **X, double **Y) {
	double **Z;
	int i, j, k;
    struct timespec ts_start, ts_end;
    double ts;

    Z = (double**)malloc(m * sizeof(double *));
	for (i = 0; i<m; i++) {
		Z[i] = (double *)malloc(s * sizeof(double));
		memset(Z[i], 0, s * sizeof(double));
	}

    clock_gettime(CLOCK_MONOTONIC, &ts_start);
	for (i = 0; i < m; i++) {
		for (j = 0; j < s; j++) {
			for (k = 0; k < n; k++) {
				Z[i][j] += X[i][k] * Y[k][j];
			}
		}
	}
    clock_gettime(CLOCK_MONOTONIC, &ts_end);
    ts = (ts_end.tv_sec - ts_start.tv_sec) * 1000.0 +
         (ts_end.tv_nsec - ts_start.tv_nsec) / 1000000.0;
	cout << "Original time: " << ts << " msecond" << endl;

	//mat_show(Z, m, s);
	/*for (i = 0; i < m; i++) { delete[] Z[i];}
	delete[] Z;*/
}
