#include <stdio.h>
#include <stdlib.h>
#include <mpi.h> 
#include <math.h>


void vector_add(
		float vector_a[],
		float vector_b[], 
		float vector_c[],
		int n
		){

	for (int i = 0; i < n; i++){
		vector_c[i] = vector_a[i] + vector_b[i];
	}	

}/*vector_add*/


float distance(float a[], float b[], int n){
	float result = 0.0f;
	for (int i = 0; i < n; i++){
		float diff = a[i] - b[i];
		result += (diff*diff);
	}

	return result;
		
 }/*distance*/


int main(){

	int my_rank;
	int comm_sz;
	int n = 1024;

	MPI_Init(NULL,NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	int local_n = n/comm_sz;


	float* vector_a = NULL;
	float* vector_b = NULL;
	float* vector_c = NULL;

	float local_a[local_n];
	float local_b[local_n];
	float local_c[local_n];

	if (my_rank == 0){
		float vector_a[n];
		float vector_b[n];
		float vector_c[n];
		float vector_d[n]; //for comparing 

		//Initialize 
		for(int i = 0; i < n; i++){
			vector_a[i] = rand()/(float)RAND_MAX;
			vector_b[i] = rand()/(float)RAND_MAX;
			
		} 
		
		vector_add(vector_a, vector_b, vector_d, n);

		MPI_Scatter(vector_a, local_n, MPI_FLOAT, local_a, local_n, MPI_FLOAT, 0, MPI_COMM_WORLD);
		MPI_Scatter(vector_b, local_n, MPI_FLOAT, local_b, local_n, MPI_FLOAT, 0, MPI_COMM_WORLD);

	        vector_add(local_a, local_b, local_c, local_n);	

		//gather the arrays
		MPI_Gather(local_c, local_n, MPI_FLOAT, vector_c, local_n, MPI_FLOAT, 0, MPI_COMM_WORLD);
		
		//compute the error
		float error = distance(vector_c, vector_d, n);

		printf("The error is %f \n", error);

	} else {

		//call to MPI_Scatter 
		MPI_Scatter(vector_a, local_n, MPI_FLOAT, local_a, local_n, MPI_FLOAT, 0, MPI_COMM_WORLD);
		MPI_Scatter(vector_b, local_n, MPI_FLOAT, local_b, local_n, MPI_FLOAT, 0, MPI_COMM_WORLD);

		vector_add(local_a, local_b, local_c, local_n);

		//call to MPI_Gather;
		MPI_Gather(local_c, local_n, MPI_FLOAT, vector_c, local_n, MPI_FLOAT, 0, MPI_COMM_WORLD);
	}


	MPI_Finalize();

	return 0;

}
