/*In this version, we used the MPI function MPI_Allreduce.
 To perform the global sum. We also use MPI_Bcast() to
 optimize distribution of input data.  */

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
	}

	MPI_Bcast(a_p, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(b_p, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(n_p, 1, MPI_INT, 0, MPI_COMM_WORLD);

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
	
	float total = 0.0f;

	float my_approx = Trap(local_a, local_b, local_n, h);
	
	MPI_Allreduce(&my_approx, &total, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);

	if (my_rank == 1){
		printf("I am process %d. With n = %d number of trapezoids and %d processes, we approximate the integral as %f\n", my_rank, n, comm_size, total);
	}


	MPI_Finalize();

	return 0;



}
















