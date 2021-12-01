#include <iostream>
#include "pre_op.h"
using namespace std;


void mat_show(double **mat, int rows, int cols) {
	int i, j;
	for (i = 0; i<rows; i++) {
		for (j = 0; j<cols; j++) {
			cout << mat[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl;
}