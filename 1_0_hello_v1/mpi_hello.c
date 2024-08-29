/*Hello world code with MPI. We have process 0 receive greetings from other processes.
 Process 0 will then print it.*/

#include <stdio.h>
#include <string.h>
#include <mpi.h>

const int MAX_STRING = 100; 

int main(){

	char greeting[MAX_STRING];

	int my_rank; 
	int comm_size;


	MPI_Init(NULL, NULL);

	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	if (my_rank != 0){
		//store a string message to the buffer
		sprintf(greeting,"Hello from process %d of %d!", my_rank, comm_size);

		//send it to process 0
		MPI_Send(greeting, strlen(greeting) + 1, MPI_CHAR, 0, 1, MPI_COMM_WORLD);
	} else {
		printf("Hello from process %d of %d!\n", my_rank, comm_size);

		//receive messages
		for (int i = 1; i < comm_size; i++){
			MPI_Recv(greeting, MAX_STRING, MPI_CHAR, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			printf("%s \n", greeting);
		}
	}

	MPI_Finalize();

	return 0;
}
