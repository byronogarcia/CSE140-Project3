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

/*
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
*/
/*
//DOESNT WORK TOO WELL YET, CACHE BLOCKING BELOW
void dgemm( int m, int n, float *A, float *C ) {
int i, j, k, l;
    int leftOvers = n - (n % 4);
    for (l = 0; l < leftOvers; l += 4) {
        for (j = 0; j < m; j++) {
            for (i = 0; i < m; i++) {
                float r = C[i+j*m];
                for (k = l; k < l + 4; k++) {
                    r += A[j+k*m] *A[i+k*m];
                }
                C[i+j*m] = r;
            }
        }
    }
    for (k = leftOvers; k < n; k++) {
        for (j = 0; j < m; j++) {
            for (i = 0; i < m; i++) {
                C[i+j*m] += A[j+k*m] * A[i+k*m];
            }
        }
    }
}

Consider A,B,C to be N-by-N matrices of b-by-b subblocks where
b=n / N is called the block size 
pseudocode for blocked (Tiled)
	   for i = 1 to N
 	      for j = 1 to N
       	{read block C(i,j) into fast memory}
       	for k = 1 to N
           	       {read block A(i,k) into fast memory}
           	       {read block B(k,j) into fast memory}
          	        C(i,j) = C(i,j) + A(i,k) * B(k,j) {do a matrix multiply on blocks}
      	 {write block C(i,j) back to slow memory}



*/

void dgemm(int m, int n, float *A, float *C ) { //PRE-FETCH
    //We will re-order the operations so that pre-fetch will be easier to predict
    //I couldn't get faster results using __builtin_prefetch()
    for(int k = 0; k < n; k++){
        for(int i = 0; i < m; i++){
            float* statA = &A[i + k*m]; //fetch the stationary address
            float* moveA = &A[k*m]; //fetch the moving address, we reordered so that we can increment by 1!
            for(int j = 0; j < m; j++){
                C[i+j*m] += (*statA) * (*moveA++); //prefetch the next moveA
            }
        }
    }
}
