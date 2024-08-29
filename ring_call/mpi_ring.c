/*In this code, I would like to demonstrate 
 a ring communication among process.*/

#include <stdio.h>
#include <mpi.h>

int main(){
	int comm_size;
	int my_rank;

	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	
	int message;

	MPI_Send(&my_rank, 1, MPI_INT, (my_rank + 1)%comm_size, 0, MPI_COMM_WORLD);
	MPI_Recv(&message, 1, MPI_INT, (my_rank - 1 + comm_size)%comm_size, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	//process 1 now have process 1's rank. Hell next send it process 2. 
	 
	int message_recv;
	
	for (int i = 1; i < comm_size; i++){
		MPI_Send(&message, 1, MPI_INT, (my_rank + 1)%comm_size, 1, MPI_COMM_WORLD);
		MPI_Recv(&message_recv, 1, MPI_INT, (my_rank - 1 + comm_size)%comm_size, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		message = message_recv;
	}
	
	printf("I am process %d of %d. The message I got at the end is: %d \n", my_rank, comm_size, message);

	MPI_Finalize();


	return 0;


}
