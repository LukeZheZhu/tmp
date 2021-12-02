#include <iostream>
#include <cstring>
#include "Lei.h"

#include "pre_op.h"
using namespace std;

double **P_left_mul( int *P,  double *P_R,  double **mat, int rows, int cols) {
	int i, j;
	double **result;

	result = (double**)malloc(rows * sizeof(double *));
	for (i = 0; i<rows; i++) {
		result[i] = (double *)malloc(cols * sizeof(double));
		memset(result[i], 0, cols * sizeof(double));
	}

	for (i = 0; i < rows; i++) {
		for (j = 0; j < cols; j++) {
			result[i][j] = mat[P[i]][j] * P_R[i];
		}
	}
	return result;
}

double **P_right_mul( int *P,  double *P_R,  double **mat, int rows, int cols) {
	int i, j;
	double **result;

	result = (double**)malloc(rows * sizeof(double *));
	for (i = 0; i<rows; i++) {
		result[i] = (double *)malloc(cols * sizeof(double));
		memset(result[i], 0, cols * sizeof(double));
	}

	for (j = 0; j < cols; j++) {
		for (i = 0; i < rows; i++) {
			result[i][P[j]] = mat[i][j] * P_R[j];
		}
	}
	return result;
}

double **P_inv_left_mul( int *P,  double *P_R,  double **mat, int rows, int cols) {
	int i, j;
	double **result;

	result = (double**)malloc(rows * sizeof(double *));
	for (i = 0; i<rows; i++) {
		result[i] = (double *)malloc(cols * sizeof(double));
		memset(result[i], 0, cols * sizeof(double));
	}

	for (i = 0; i < rows; i++) {
		for (j = 0; j < cols; j++) {
			result[P[i]][j] = mat[i][j] / P_R[i];
		}
	}
	return result;
}

double **P_inv_right_mul( int *P,  double *P_R,  double **mat, int rows, int cols) {
	int i, j;
	double **result;

	result = (double**)malloc(rows * sizeof(double *));
	for (i = 0; i<rows; i++) {
		result[i] = (double *)malloc(cols * sizeof(double));
		memset(result[i], 0, cols * sizeof(double));
	}

	for (j = 0; j < cols; j++) {
		for (i = 0; i < rows; i++) {
			result[i][j] = mat[i][P[j]] / P_R[j];
		}
	}
	return result;
}


void Lei(double **X, double **Y) {
	int i, j, k;
    struct timespec ts_start, ts_end;
    double ts;

/*######
key_Gen
######*/
	 int *P1 = new int[m];
	 double *P1_R = new double[m];
	 int *P2 = new int[n];
	 double *P2_R = new double[n];
	 int *P3 = new int[s];
	 double *P3_R = new double[s];

	for (i = 0; i < m; i++) { P1[i] = i; }
	for (i = 0; i < n; i++) { P2[i] = i; }
	for (i = 0; i < s; i++) { P3[i] = i; }

	// gen P
    clock_gettime(CLOCK_MONOTONIC, &ts_start);
	int swap_tmp;
	for (i = 0; i < m; i++) {
		j = rand() % m;
		swap_tmp = P1[i];
		P1[i] = P1[j];
		P1[j] = swap_tmp;

		P1_R[i] = rand() % DOMAIN_MAX + 1;
	}

	for (i = 0; i < n; i++) {
		j = rand() % n;
		swap_tmp = P2[i];
		P2[i] = P2[j];
		P2[j] = swap_tmp;

		P2_R[i] = rand() % DOMAIN_MAX + 1;
	}

	for (i = 0; i < s; i++) {
		j = rand() % s;
		swap_tmp = P3[i];
		P3[i] = P3[j];
		P3[j] = swap_tmp;

		P3_R[i] = rand() % DOMAIN_MAX + 1;
	}

    clock_gettime(CLOCK_MONOTONIC, &ts_end);
    ts = (ts_end.tv_sec - ts_start.tv_sec) * 1000.0 +
         (ts_end.tv_nsec - ts_start.tv_nsec) / 1000000.0;
	cout << "key_Gen time: " << ts << " msecond" << endl;

/*######
data_Enc
特殊矩阵相乘时：
i遍历P的size，再进一步找规律
X's dimension: m*n
######*/
    clock_gettime(CLOCK_MONOTONIC, &ts_start);

    double **tmp;

	// X_enc = P1 * X * P2_inv
    double **X_enc;
	// P1 * X
	tmp = P_left_mul(P1, P1_R, X, m, n);
	// (P1 * X) * P2_inv
	X_enc = P_inv_right_mul(P2, P2_R, tmp, m, n);

	// Y_enc = P2 * Y * P3_inv
    double **Y_enc;
	// P2 * Y
	tmp = P_left_mul(P2, P2_R, Y, n, s);
	// (P2 * Y) * P3_inv
	Y_enc = P_inv_right_mul(P3, P3_R, tmp, n, s);

    clock_gettime(CLOCK_MONOTONIC, &ts_end);
    ts = (ts_end.tv_sec - ts_start.tv_sec) * 1000.0 +
         (ts_end.tv_nsec - ts_start.tv_nsec) / 1000000.0;
	cout << "Data_Enc time: " << ts << "msecond" << endl;

/*######
cloud_MMC
######*/
    clock_gettime(CLOCK_MONOTONIC, &ts_start);

    double **Z_enc;
	Z_enc = (double**)malloc(m * sizeof(double *));
	for (i = 0; i<m; i++) {
		Z_enc[i] = (double *)malloc(s * sizeof(double));
		memset(Z_enc[i], 0, s * sizeof(double));
	}

	for (i = 0; i < m; i++) {
		for (j = 0; j < s; j++) {
			for (k = 0; k < n; k++) {
				Z_enc[i][j] += X_enc[i][k] * Y_enc[k][j];
			}
		}
	}

    clock_gettime(CLOCK_MONOTONIC, &ts_end);
    ts = (ts_end.tv_sec - ts_start.tv_sec) * 1000.0 +
         (ts_end.tv_nsec - ts_start.tv_nsec) / 1000000.0;
	cout << "Cloud_MMC time: " << ts << "msecond" << endl;

/*######
result_Dec
######*/
    clock_gettime(CLOCK_MONOTONIC, &ts_start);

	// Z = P1_inv * Z_enc * P3
    double **Z;
	// P1_inv * Z_enc
	tmp = P_inv_left_mul(P1, P1_R, Z_enc, m, s);
	// (P1_inv * Z_enc) * P3
	Z = P_right_mul(P3, P3_R, tmp, m, s);

    clock_gettime(CLOCK_MONOTONIC, &ts_end);
    ts = (ts_end.tv_sec - ts_start.tv_sec) * 1000.0 +
         (ts_end.tv_nsec - ts_start.tv_nsec) / 1000000.0;
	cout << "Result_Dec time: " << ts << "msecond" << endl;

/*######
result_verification
######*/
    clock_gettime(CLOCK_MONOTONIC, &ts_start);

	bool flag = true;
	for (k = 0; k < l; k++) {
		 double r[s];
		for (i = 0; i < s; i++) { r[i] = rand() % 2; }
		 double verify_tmp[n] = { 0 };
		 double verify_1[m] = { 0 };
		for (i = 0; i < n; i++) {
			for (j = 0; j < s; j++) {
				verify_tmp[i] += Y[i][j] * r[j];
			}
		}
		for (i = 0; i < m; i++) {
			for (j = 0; j < n; j++) {
				verify_1[i] += X[i][j] * verify_tmp[j];
			}
		}

		 double verify_2[m] = { 0 };
		for (i = 0; i < m; i++) {
			for (j = 0; j < s; j++) {
				verify_2[i] += Z[i][j] * r[j];
			}
		}

		for (i = 0; i < m; i++) {
			if (verify_1[i] != verify_2[i]) {
				//flag = false;
				//break;
			}
		}
		//if (!flag)
			//break;
	}
	//if (!flag)
		//cout << "False result";
	//else
		//cout << "Verification Succeed";

    clock_gettime(CLOCK_MONOTONIC, &ts_end);
    ts = (ts_end.tv_sec - ts_start.tv_sec) * 1000.0 +
         (ts_end.tv_nsec - ts_start.tv_nsec) / 1000000.0;
	cout << "Result_Verify time: " << ts << "msecond" << endl;

	//mat_show(Z, m, s);
}
