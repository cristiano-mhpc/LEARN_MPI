/*We use MPI to parallelize the implementation of the trapezoidal rule.
 Each process will get a semi interval over which to calculate its
 contribution. They will then send their contribution to processs 0.
 Process 0 will add up all the contributions. */

#include <stdio.h>
#include <mpi.h>

float f(float x){
	return x*x + 4.0;

}/*f*/

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

	float a = -5.0, b = 5.0;
       	int n = 1024;
	float h = (b-a)/n;

	int comm_size, my_rank;

	MPI_Init(NULL,NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	
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
		printf("With n = %d number of trapezoids, we approximatre the integral as %f\n", n, my_approx);
	}



	MPI_Finalize();

	return 0;



}
