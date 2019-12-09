/*
The project requires a max of three methods which optimize
on the naive solution. Two methods will be commented out at
all times to allow for the testing to work. Otherwise it
would not run since they will all be called dgemm, and have
the same arguments.

3 OF THESE METHODS WILL BE IMPLEMENTED:
BLOCKING/TILING - 40PTS
PRE-FETCHING - 4OPTS
SSE INSTRUCTIONS - 40PTS
UNROLLING - 40PTS
PADDING MATRICES - 4OPT
*/

void dgemm( int m, int n, float *A, float *C ) {
	//UNROLLING
	int i, j, k, l;
	for (j = 0; j < m; j++) {
		for (k = 0; k < n - 3; k += 4) {
			for (i = 0; i < m; i++) {
				int index = i + j * m;
				C[index] += A[i + k * m] * A[j + k * m];
				C[index] += A[i + (k + 1) * m] * A[j + (k + 1) * m];
				C[index] += A[i + (k + 2) * m] * A[j + (k + 2) * m];
				C[index] += A[i + (k + 3) * m] * A[j + (k + 3) * m];
			}
		}
		for (l = k; l < n; l++) {
			float a = A[j + l * m];
			for (i = 0; i < m; i ++) {
				C[i + j * m] += A[i + l * m] * a;
			}
		}
	}
}

/*
void dgemm( int m, int n, float *A, float *C ) {
// FILL-IN  
}

void dgemm( int m, int n, float *A, float *C ) {
// FILL-IN 
}
*/
