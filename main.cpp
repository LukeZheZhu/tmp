#include <iostream>
#include <time.h>
#include "Lei.h"
#include "Our.h"
#include "original.h"
#include "pre_op.h"
using namespace std;

int main()
{
    double **X, **Y;
	int i, j, k;

	X = (double **)malloc(m * sizeof(double *));
	for (i = 0; i<m; i++) {X[i] = (double *)malloc(n * sizeof(double));}
	Y = (double**)malloc(n * sizeof(double *));
	for (i = 0; i<n; i++) { Y[i] = (double *)malloc(s * sizeof(double)); }

	// Ëæ»úÖÖ×Ó
	srand(200);

	for (i = 0; i<m; i++) { for (j = 0; j<n; j++) { X[i][j] = rand() % DOMAIN_MAX; } }
	for (i = 0; i<n; i++) { for (j = 0; j<s; j++) { Y[i][j] = rand() % DOMAIN_MAX; } }


    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    clock_gettime(CLOCK_MONOTONIC, &ts);


	srand(time(NULL));

	cout << endl << endl << endl;
	cout << "############Original############" << endl;
	original(X, Y);

	cout << endl << endl << endl;
	cout << "############Original############" << endl;
	original(X, Y);

	cout << endl << endl << endl;
	cout << "###########Lei scheme###########" << endl;
	Lei(X, Y);

	cout << endl << endl << endl;
	cout << "###########Our scheme###########" << endl;
	Our(X, Y);


//	getchar();

	return 0;

}
