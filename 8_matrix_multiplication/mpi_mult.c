/*This code implements parallel vector multiplication
 *                   y = A*x,    (A is m by n)
 * with MPI. Each process gets a block of y and x 
 * as well as a block of A. Each process gets n/comm_sz
 * rows of A and we use block distribution. In order for all
 * process to get a copy of the entire array x, we will
 * use the MPI_Allgather() function.
 */

#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <stdlib.h>

float error(float y_paral[],
	    float y_seq[],
	    int m){
	float error = 0.0f;
	float sum = 0.0f;
	for (int i = 0; i < m; i++){
		error = y_paral[i]-y_seq[i];
		sum += (error * error);
	}

	return sum;

}/*error*/

void matrix_mult(float* local_A,
		float*  local_x,
		float*  local_y,
		int     local_m,
		int     local_n,
		int     n){

	//allocate memory for the vector x
	float* x = NULL;
	x = (float*)malloc(n*sizeof(float));

	//Initialize my local block of A
	for (int i = 0; i < local_m; i++){
		for(int j = 0; j < n; j++){
			local_A[i*n + j] = rand() / (float)RAND_MAX;
		}
	}
	
	//Each process will initialize a local block of x
	for (int j = 0; j < local_n; j++){
		local_x[j] = rand() / (float)RAND_MAX;
	}

	//Processes use MPI_Allgather() to make sure every process has x
	MPI_Allgather(local_x, local_n, MPI_FLOAT, x, local_n, MPI_FLOAT, MPI_COMM_WORLD); 

	//perform a local matrix multiplication
	for (int i = 0; i < local_m; i++){
		local_y[i] = 0.0f;
		for (int j = 0; j < n; j++){
			local_y[i] += local_A[i*n + j]*x[j];
		}
	}
		
	free(x);

}/*matrix_mult*/

int main(){

	int m = 1024;
	int n = 2048;
	
	int my_rank;
	int comm_sz;
	
	float* x = NULL;
	float* y = NULL; 
	float* A = NULL;
	float* local_y = NULL;
	float* local_x = NULL;
	float* local_A = NULL;

	MPI_Init(NULL, NULL);

	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	//Important: we assume m and n are multiples of comm_sz.
	int local_m = m/comm_sz; 
	int local_n = n/comm_sz; 
	
	local_y = (float*)malloc(local_m * sizeof(float));
	local_x = (float*)malloc(local_n * sizeof(float));
	local_A = (float*)malloc( (local_m * n) * sizeof(float));
	
	//call matrix_mult()
	matrix_mult(local_A, local_x, local_y, local_m, local_n, n);

	if(my_rank == 0){
		
		float* y_seq = NULL;

		A = malloc(n*m*sizeof(float));
		y = malloc(m*sizeof(float));
		y_seq = (float*) malloc(m*sizeof(float));	
		x = malloc(n*sizeof(float));
			
		//Processes call MPI_Gather() to gather their blocks in process 0. We need A for sequential check.
		MPI_Gather(local_A, (m*n)/comm_sz, MPI_FLOAT, A, (m*n)/comm_sz, MPI_FLOAT, 0, MPI_COMM_WORLD);
		
		//Processes call MPI_Gather() to gather the local_y into y in process 0.
		MPI_Gather(local_y, local_m, MPI_FLOAT, y, local_m, MPI_FLOAT, 0,  MPI_COMM_WORLD);

		//Processes call MPI_Gather() to gather the local_x into x in process 0.
		MPI_Gather(local_x, local_n, MPI_FLOAT, x, local_n, MPI_FLOAT, 0, MPI_COMM_WORLD);

		//Perform sequential multiplication.
		for (int i = 0; i < m; i++){
			y_seq[i] = 0.0f;
			for (int j = 0; j < n; j++){
				y_seq[i] += A[i*n + j]*x[j];
			}	

		}

		float err = error(y, y_seq, m);
		printf("The error is %f: \n", err);

		free(y_seq);		

	} else {

		MPI_Gather(local_A, (m*n)/comm_sz, MPI_FLOAT, A, (m*n)/comm_sz, MPI_FLOAT, 0, MPI_COMM_WORLD);

		MPI_Gather(local_y, local_m, MPI_FLOAT, y, local_m, MPI_FLOAT, 0, MPI_COMM_WORLD);

		MPI_Gather(local_x, local_n, MPI_FLOAT, x, local_n, MPI_FLOAT, 0, MPI_COMM_WORLD);

	}

	free(x);
	free(y);
	free(A);
	free(local_y);
        free(local_x);
        free(local_A);

	MPI_Finalize();
	return 0;

}
