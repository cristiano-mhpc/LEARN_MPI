/*In this version, we added the function get_input to get 
  process 0 to read a, b, and n from stdin. */

#include <stdio.h>
#include <mpi.h>

float f(float x){
	return x*x + 4.0;

}/*f*/


void get_input(int my_rank, 
		int comm_size,
		float* a_p,
		float* b_p,
		int* n_p){

	if (my_rank == 0){
		printf("Enter a, b, and n. Choose an n that is multiple of comm_size:\n");
		scanf( "%f %f %d" , a_p, b_p, n_p);

		for(int i = 1; i < comm_size; i++){
			MPI_Send(a_p, 1, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
			MPI_Send(b_p, 1, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
			MPI_Send(n_p, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
		}

		} else {
			MPI_Recv(a_p, 1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(b_p, 1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(n_p, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);	
	
	}

}/*get_input*/

float Trap(float left_endpoint, 
	   float right_endpoint, 
	   int trap_count, 
	   float base){
		
	float approx = (f(left_endpoint) + f(right_endpoint))/2.0;

	for (int i = 1; i < trap_count; i++){
		float x = left_endpoint + i*base;
	      	approx += f(x);	
	}

	return base*approx;

}/*Trap*/


int main(){
	float a, b;
       	int n;

	int comm_size, my_rank;

	MPI_Init(NULL,NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	
	get_input(my_rank, comm_size, &a, &b, &n);

	float h = (b-a)/n;

	int local_n = n/comm_size; // we assume n is a multiple of comm_size 

	float local_a = a + local_n * my_rank*h;
	float local_b = local_a + local_n*h;


	float my_approx = Trap(local_a, local_b, local_n, h);

	if (my_rank != 0){

		MPI_Send(&my_approx, 1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);

	} else { /*my_rank = 0 */
		
		float total;

		for (int i = 1; i < comm_size; i++){
			MPI_Recv(&total, 1, MPI_FLOAT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			my_approx += total;
		}
		printf("With n = %d number of trapezoids and %d processes, we approximatre the integral as %f\n", n, comm_size, my_approx);
	}



	MPI_Finalize();

	return 0;



}
















