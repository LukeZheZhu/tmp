#include <iostream>
#include <cstring>
#include "Our.h"
#include "pre_op.h"
using namespace std;

//M: rows * rows
double **M_left_mul( int *P,  double *P_R,  double *H,  double **mat, int rows, int cols) {
	int i, j;
	double *col_sum;
	double **result;

	result = (double**)malloc(rows * sizeof(double *));
	for (i = 0; i<rows; i++) {
		result[i] = (double *)malloc(cols * sizeof(double));
		memset(result[i], 0, cols * sizeof(double));
	}

	col_sum = (double *)malloc(cols * sizeof(double));
	memset(col_sum, 0, cols * sizeof(double));
	for (j = 0; j < cols; j++) {
		for (i = 0; i < rows; i++) {
			col_sum[j] += mat[i][j];
		}
	}

	for (i = 0; i < rows; i++) {
		for (j = 0; j < cols; j++) {
			// P * mat
			result[i][j] += mat[P[i]][j] * P_R[i];
			// (P * mat) + H * mat
			result[i][j] += H[i] * col_sum[j];
		}
	}
	return result;
}

//M: cols * cols
double **M_right_mul( int *P,  double *P_R,  double *H,  double **mat, int rows, int cols) {
	int i, j;
	double *row_sum;
	double **result;

	result = (double**)malloc(rows * sizeof(double *));
	for (i = 0; i<rows; i++) {
		result[i] = (double *)malloc(cols * sizeof(double));
		memset(result[i], 0, cols * sizeof(double));
	}

	row_sum = (double *)malloc(rows * sizeof(double));
	memset(row_sum, 0, rows * sizeof(double));
	for (i = 0; i < rows; i++) {
		for (j = 0; j < cols; j++) {
			row_sum[i] += H[j] * mat[i][j];
		}
	}

	for (j = 0; j < cols; j++) {
		for (i = 0; i < rows; i++) {
			// mat * P
			result[i][P[j]] += mat[i][j] * P_R[j];
			// (mat * P) + mat * H
			result[i][j] += row_sum[i];
		}
	}
	return result;
}

//M: rows * rows
double **M_inv_left_mul( int *P,  double *P_R,  double *H,  double **mat, int rows, int cols) {
	int i, j;
	double **result, **Pinv_X, **tmp;

	// P_inv * X
	Pinv_X = (double**)malloc(rows * sizeof(double *));
	for (i = 0; i<rows; i++) {
		Pinv_X[i] = (double *)malloc(cols * sizeof(double));
		memset(Pinv_X[i], 0, cols * sizeof(double));
	}

	for (i = 0; i < rows; i++) {
		for (j = 0; j < cols; j++) {
			Pinv_X[P[i]][j] = mat[i][j] / P_R[i];
		}
	}

	// H * (P_inv * X)
	double *col_sum;
	col_sum = (double *)malloc(cols * sizeof(double));
	memset(col_sum, 0, cols * sizeof(double));

	for (j = 0; j < cols; j++) {
		for (i = 0; i < rows; i++) {
			col_sum[j] += Pinv_X[i][j];
		}
	}

	//tmp = P_inv * (H * (P_inv * X))
	tmp = (double**)malloc(rows * sizeof(double *));
	for (i = 0; i<rows; i++) {
		tmp[i] = (double *)malloc(cols * sizeof(double));
		memset(tmp[i], 0, cols * sizeof(double));
	}

	for (i = 0; i < rows; i++) {
		for (j = 0; j < cols; j++) {
			//result -> tmp
			//mat -> (H * (P_inv * X))[i][j] = H[i]*col_sum[j]
			tmp[P[i]][j] = H[i] * col_sum[j] / P_R[i];
		}
	}

	//tr(H * P_inv)
	double tr = 0;
	for (i = 0; i < rows; i++) {
		// mat[i][j] = H[i]
		tr += H[i] / P_R[i];
	}
	double coe = 1 / (1 + tr);

	//result
	result = (double**)malloc(rows * sizeof(double *));
	for (i = 0; i<rows; i++) {
		result[i] = (double *)malloc(cols * sizeof(double));
		memset(result[i], 0, cols * sizeof(double));
	}

	for (i = 0; i < rows; i++) {
		for (j = 0; j < cols; j++) {
			result[i][j] = Pinv_X[i][j] - coe * tmp[i][j];
		}
	}

	return result;
}

//M: cols * cols
double **M_inv_right_mul( int *P,  double *P_R,  double *H,  double **mat, int rows, int cols) {
	int i, j;
	double **result, **X_Pinv, **tmp;

	// X * P_inv
	X_Pinv = (double**)malloc(rows * sizeof(double *));
	for (i = 0; i<rows; i++) {
		X_Pinv[i] = (double *)malloc(cols * sizeof(double));
		memset(X_Pinv[i], 0, cols * sizeof(double));
	}

	for (j = 0; j < cols; j++) {
		for (i = 0; i < rows; i++) {
			X_Pinv[i][j] = mat[i][P[j]] / P_R[j];
		}
	}

	// (X * P_inv) * H
	double *row_sum;
	row_sum = (double *)malloc(rows * sizeof(double));
	memset(row_sum, 0, rows * sizeof(double));
	for (i = 0; i < rows; i++) {
		for (j = 0; j < cols; j++) {
			row_sum[i] += H[j] * X_Pinv[i][j];
		}
	}

	//tmp = ((X * P_inv) * H) * P_inv
	tmp = (double**)malloc(rows * sizeof(double *));
	for (i = 0; i<rows; i++) {
		tmp[i] = (double *)malloc(cols * sizeof(double));
		memset(tmp[i], 0, cols * sizeof(double));
	}

	for (j = 0; j < cols; j++) {
		for (i = 0; i < rows; i++) {
			//result -> tmp
			//mat -> ((X * P_inv) * H)[i][j] = row_sum[i]
			tmp[i][j] = row_sum[i] / P_R[j];
		}
	}

	//tr(H * P_inv)
	double tr = 0;
	for (i = 0; i < cols; i++) {
		// mat[i][j] = H[i]
		tr += H[i] / P_R[i];
	}
	double coe = 1 / (1 + tr);

	// result
	result = (double**)malloc(rows * sizeof(double *));
	for (i = 0; i<rows; i++) {
		result[i] = (double *)malloc(cols * sizeof(double));
		memset(result[i], 0, cols * sizeof(double));
	}

	for (i = 0; i < rows; i++) {
		for (j = 0; j < cols; j++) {
			result[i][j] = X_Pinv[i][j] - coe * tmp[i][j];
		}
	}

	return result;
}

void Our(double **X, double **Y) {
	int i, j, k;
    struct timespec ts_start, ts_end;
    double ts;

/*######
key_Gen
######*/
	 int *P1 = new int[m];
	 double *P1_R = new double[m];
	 double *H1 = new double[m];
	 int *P2 = new int[n];
	 double *P2_R = new double[n];
	 double *H2 = new double[n];
	 int *P3 = new int[s];
	 double *P3_R = new double[s];
	 double *H3 = new double[s];

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
		H1[i] = rand() % DOMAIN_MAX + 1;
	}
	for (i = 0; i < n; i++) {
		j = rand() % n;
		swap_tmp = P2[i];
		P2[i] = P2[j];
		P2[j] = swap_tmp;

		P2_R[i] = rand() % DOMAIN_MAX + 1;
		H2[i] = rand() % DOMAIN_MAX + 1;
	}

	for (i = 0; i < s; i++) {
		j = rand() % s;
		swap_tmp = P3[i];
		P3[i] = P3[j];
		P3[j] = swap_tmp;

		P3_R[i] = rand() % DOMAIN_MAX + 1;
		H3[i] = rand() % DOMAIN_MAX + 1;
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

	// X_enc = M1 * X * M2_inv
    double **X_enc;
	// M1 * X
	tmp = M_left_mul(P1, P1_R, H1, X, m, n);
	// (M1 * X) * M2_inv
	X_enc = M_inv_right_mul(P2, P2_R, H2, tmp, m, n);

	// Y_enc = M2 * Y * M3_inv
	 double **Y_enc;
	// M2 * Y
	tmp = M_left_mul(P2, P2_R, H2, Y, n, s);
	// (M2 * Y) * M3_inv
	Y_enc = M_inv_right_mul(P3, P3_R, H3, tmp, n, s);

    clock_gettime(CLOCK_MONOTONIC, &ts_end);
    ts = (ts_end.tv_sec - ts_start.tv_sec) * 1000.0 +
         (ts_end.tv_nsec - ts_start.tv_nsec) / 1000000.0;
	cout << "Data_Enc time: " << ts << " msecond" << endl;

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
    ts = (ts_end.tv_sec - ts_start.tv_sec) * 1000 + (ts_end.tv_nsec - ts_start.tv_nsec)/1000000;
	cout << "cloud_mmc time: " << ts << " msecond" << endl;

/*######
result_Dec
######*/
    clock_gettime(CLOCK_MONOTONIC, &ts_start);

	// Z = M1_inv * Z_enc * M3
    double **Z;
	// M1_inv * Z_enc
	tmp = M_inv_left_mul(P1, P1_R, H1, Z_enc, m, s);
	// (M1_inv * Z_enc) * M3
	Z = M_right_mul(P3, P3_R, H3, tmp, m, s);

    clock_gettime(CLOCK_MONOTONIC, &ts_end);
    ts = (ts_end.tv_sec - ts_start.tv_sec) * 1000.0 +
         (ts_end.tv_nsec - ts_start.tv_nsec) / 1000000.0;
	cout << "Result_Gen time: " << ts << " msecond" << endl;

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
	cout << "Result_Verify time: " << ts << " msecond" << endl;

	//mat_show(Z, m, s);
}
