#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

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
	//float* vector_c = NULL;

	float local_a[local_n];
	float local_b[local_n];
	float local_c[local_n];

	if (my_rank == 0){
		float vector_a[n];
		float vector_b[n];
		float vector_c[n];//we'll need later for gather

		for(int i = 0; i < n; i++){
			vector_a[i] = rand()/(float)RAND_MAX;
			vector_b[i] = rand()/(float)RAND_MAX;
		} 
		printf("The sum in process 2 should be \n");
		for(int i = 2*local_n; i < 3*local_n; i++){
			vector_c[i] = vector_a[i] + vector_b[i];
			printf("%f ", vector_c[i]);
		}	

		printf("\n");

		MPI_Scatter(vector_a, local_n, MPI_FLOAT, local_a, local_n, MPI_FLOAT, 0, MPI_COMM_WORLD);
		MPI_Scatter(vector_b, local_n, MPI_FLOAT, local_b, local_n, MPI_FLOAT, 0, MPI_COMM_WORLD);
		
	}else {
	
		MPI_Scatter(vector_a, local_n, MPI_FLOAT, local_a, local_n, MPI_FLOAT, 0, MPI_COMM_WORLD);
		MPI_Scatter(vector_b, local_n, MPI_FLOAT, local_b, local_n, MPI_FLOAT, 0, MPI_COMM_WORLD);
	}

	vector_add(local_a, local_b, local_c, local_n);

	if(my_rank == 2){
		printf("The array of sums assigned to process %d is \n ", my_rank);
		for (int i = 0; i < local_n; i++){
			printf("%f ", local_c[i]);
		}
		printf("\n");

	}	

	MPI_Finalize();

	return 0;

}
