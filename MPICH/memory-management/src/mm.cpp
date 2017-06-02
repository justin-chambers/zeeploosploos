/**
 @file: mm.cpp
 @author: Justin Chambers
 @date: 20170314
 @version: 1.0
 @brief: snippet - demonstrate array passing in MPI
*/

#include <iostream>
#include <cstring>
#include <cstdlib>
#include <mpi.h>

#define CONTROLLER		0

/**
 @brief:        main
 @details:      pass an array to a node, have node write stuff in part of that
                array, then return it, and store it.
 @param:        argc - number of command line arguments
                argv - command line strings
 @return:       0 upon completion
 @exception:    none
*/

int main(int argc, char ** argv) {
    	// variables
    	int m=10, n=10, p, row, offset;
        int * buffer = NULL;
	int * final = NULL;
    	int world_size, world_rank;

    	//MPI init...
    	MPI_Init(NULL, NULL);
    	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	MPI_Status status;
	p = world_size - 1;

        // get partition size...
        if(n % p == 0) offset = (int)(n/p);
        else offset = (int)(n/p)+1;

        // init buffers...
	int chunk_size = offset*m;
        buffer = new int [chunk_size];
	final = new int [chunk_size*p];

        // init array row members (m of them...)
        for(int i = 0; i < (chunk_size); ++i) buffer[i] = 0;
	for(int i = 0; i < (chunk_size*p); ++i) final[i] = 0;

    	if(world_rank == CONTROLLER) {
		// CONTROLLER stuff:

		// send rows to nodes...
		int p_i;
		for(p_i=1, row=0; p_i<=p;++p_i,row=row+offset) {
			MPI_Send(&(final[row]),chunk_size,MPI_INT,p_i,p_i,MPI_COMM_WORLD);
            	}
		// receive calculated pixels from nodes...
		for(int p_i=1, c_row = 0; p_i<=p;++p_i, c_row = c_row+offset) {
			MPI_Recv(&(final[c_row*m]),chunk_size,MPI_INT,p_i,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
		}

    		// display final results...
		std::cout << world_rank << ": final display...\n";
		for(int i = 0; i < (m*n); ++i) {
			std::cout << final[i] << " ";
			if(((i+1) % m) == 0) std::cout << "\n";
		}
		delete [] buffer;
		delete [] final;
		buffer = NULL;
		final = NULL;
    	} else {
            // NODE stuff:

            // get row starting location...
    		MPI_Recv(&(buffer[0]),chunk_size,MPI_INT,CONTROLLER,world_rank,MPI_COMM_WORLD, &status);
            // do stuff...
		row = (world_rank-1)*offset;
//std::cout << world_rank << ": row = " << row << "\noffset = " << offset << std::endl;
            for(int i=0; i<offset;++i) {
                for(int j=0; j<m; ++j) {
//std::cout << world_rank << ": buffer[" << i*m+j << "] = " << (row+i)*m+j << " ";
			buffer[i*m+j] = (row+i)*m+j;
                }
            }

            // send stuff back to controller
    		MPI_Send(&(buffer[0]),chunk_size,MPI_INT,CONTROLLER,world_rank,MPI_COMM_WORLD);
    	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
    	return 0;
}
