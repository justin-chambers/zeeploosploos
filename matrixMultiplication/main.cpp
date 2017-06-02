/**
 @file: mm_para.cpp
 @version: 1.0
 @brief: Main program for matrix multiplication (parallel)
*/

#include "Utilities.h"
//#include <mpi.h>

#define ROW_MAX		10000
#define UPPER		100
#define CONTOLLER	0
#define DISPLAY_ARRAY	0
#define WRITE_A		1
#define WRITE_B		1
#define WRITE_C		1

/**
 @brief: function prototypes...
*/

void multiply_matrices(int * A,
                       const int mA,
                       const int nA,
                       int * B,
                       const int mB,
                       const int nB,
                       int * C );

bool is_perfect_square(int x, int& t);

void matrix_to_tile(	int * M,
                        int size_of_m,
                        int * T,
                        int size_of_t,
                        int src_id,
                        int sqrt_p);

void tile_to_matrix(	int * M,
                        int size_of_m,
                        int * T,
                        int size_of_t,
                        int src_id,
                        int sqrt_p);

void init_hilbert_matrix(float ** M, int n);

void shift_rows_left(int n,
                     int * M,
                     int * snd,
                     int * rcv,
                     int src_id,
                     int w_size);

void shift_cols_up(int n,
                  int * M,
                  int * snd,
                  int * rcv,
                  int src_id,
                  int shft_amt,
                  int w_size);

void flatten(int n, float M[n][n], float buf[]);

/**
 @brief:        main
 @details:      reads in a matrix size (N x N only), randomly fills two matrices, multiplies them, stores the
                result in a third matrix.
 @param:        argc -	the number of commandline arguments
                argv -	pointer to the arguments
 			argv[0] is the name of the binary
			argv[1] is the dimension of the n x n matrices A, B to be multiplied
			argv[2] is the number of tiles (i.e. processors) to break the matrix into
			argv[3] is the name of the output file to be written (debug only)
 @return:       none
 @exception:    none
*/

int main(int argc, char ** argv) {
	// variables
	int * A;
	int * B;
	int * C;
	int * buffer;
	int world_rank = 0, world_size;
	int sqrt_p = 1, tile_size = 1;
	int i = 0, j = 0;
	double start = 0, delta = 0;

	Utilities util;

	// handle command line
	util.handle_argc(argc,5);

	// init MPI
	//MPI_Init(&argc, &argv);
	//MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	//MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	// get the number of processors...
	int p = util.get_positive_integer_from_cmd_line(argv[1]);

	// get the problem size...
	int n = util.get_array_size_from_file(argv[2]);

	// fast exit conditions...
	// number of processors is not a perfect square...
	if(!is_perfect_square(p, sqrt_p)) {
		std::cout << world_rank << ": number of processors is not a perfect square. Exiting.\n";
		//MPI_Finalize();
		return 1;
	} else if((n%sqrt_p)!=0) { // bad problem size...
		std::cout << world_rank << ": " << n << " x " << n
			<< " matrix cannot be partitioned evenly into "
			<< sqrt_p << " x " << sqrt_p << " tiles. Exiting.\n";
		//MPI_Finalize();
		return 1;
	}

	// update the tile size...
	tile_size = n / sqrt_p;

	// allocate memory...
	A = new int[n*n];
	B = new int[n*n];
	C = new int[n*n];
	int * tileA = new int[tile_size*tile_size];


    if(world_rank == CONTOLLER) {
        // read in matrices...
        util.read_file_to_one_dee_integer_array(argv[2], A, n*n);
        util.read_file_to_one_dee_integer_array(argv[3], B, n*n);

        // multiply them...
        multiply_matrices(A, n, n, B, n, n, C);

        // write to file...
        std::string output_fn = "c_out.txt";
        util.write_flattened_two_dee_integer_array_to_file(output_fn.c_str(), C, n, n);

        /*

        // disperse matrices A, B to other nodes...
        for(int i=0;i<p;++i) {
            // tile matrix...
            matrix_to_tile(A, n, tileA, tile_size, i, sqrt_p);

            // display tile...
            std::cout << "Tile for processor " << i << std::endl;
            util.display_flattened_two_dee_array(tileA,tile_size,tile_size);

            // place tiles in to C...
            tile_to_matrix(C, n, tileA, tile_size, i, sqrt_p);
        }

        std::cout << "Restored Matrix...\n";
        util.display_flattened_two_dee_array(C,n,n);

         */

        // start timer...
        //start = MPI_Wtime();
    }

	// cannon's algorithm...

	// align the tiles for multiplication...
	//int shift_amt = (int)(world_rank/sqrt_p);
	//for(i=0;i<shift_amt;++i) shift_rows_left(A, buffer, tile_size, world_rank);
	//for(i=0;i<shift_amt;++i) shift_cols_up(B, buffer, tile_size, world_rank);

	// multiply and shift the remaining tiles...
	/*
	for(i=1;i<sqrt_p-1;++i) {
		multiply_matrices(A, tile_size, tile_size,
			B, tile_size, tile_size,
			C);
		shift_rows_left(A, buffer, tile_size, world_rank);
		shift_cols_up(B, buffer, tile_size, world_rank);
	}*/

	// stop timer...
	//delta = MPI_Wtime() - start;

	// log results...
	//std::cout << world_size << "," << p << "," << n << "," << delta << std::endl;

#if DISPLAY_ARRAY
	for(i=0;i<n;++i) {
		for(j=0;j<n;++j) {
			std::cout << C[i*n+j] << " ";
 		}
		std::cout << std::endl;
	}
#endif

#if WRITE_A
	std::string output_a = "a_from_memory.txt";
	util.write_flattened_two_dee_integer_array_to_file(output_a.c_str(), A, n, n);
#endif

#if WRITE_B
	std::string output_b = "b_from_memory.txt";
	util.write_flattened_two_dee_integer_array_to_file(output_b.c_str(), B, n, n);
#endif

#if WRITE_C
	std::string output_c = argv[4];
	util.write_flattened_two_dee_integer_array_to_file(output_c.c_str(), C, n, n);
#endif

	// finalize
	//MPI_Finalize();

	return 0;
}

/**
 @brief:        multiply_matrices
 @details:      multiplies two matrices and stores the result in a third matrix.
 @param:        A   -   pointer to the first matrix
                mA  -   count of rows in A
                nA  -   count of columns in A
                B   -   pointer to the seconds matrix
                mB  -   count of rows in B
                nB  -   count of columns in B
                C   -   pointer to the third matrix
 @return:       none
 @exception:    none
*/

void multiply_matrices(int * A,
                       const int mA,
                       const int nA,
                       int * B,
                       const int mB,
                       const int nB,
                       int * C ) {
	if(mA == nB) {
		int i = 0, j = 0, k = 0;
		for(i=0; i<mA;++i) { // rows of A
			for(j=0; j<nB;++j) { // columns of B
			//compute the cell value
                C[i*mA+j] = 0;
				for(k=0; k<nB;++k) {
					C[i*mA+j] += A[i*mA+k] * B[k*nB+j];
				}
			}
		}
	}
}

/**
 @brief:        is_perfect_square
 @details:      determines whether a number is a perfect square
 @param:        x   -   integer to be tested
 @return:       true, if x is a perfect square
                false, otherwise
 @exception:    none
*/

bool is_perfect_square(int x, int& t) {
	t = 1;
	while(t*t < x) {
		++t;
	}
	if(t*t != x) return false;
	else return true;
}

/**
 @brief:        matrix_to_tile
 @details:      converts a portion of a large n x n matrix into a
                smaller t x t tile
 @param:        M   -   matrix storage pointer
                size_of_m - dimensions of M
                T   -   tile storage pointer
                size_of_t - dimensions of T
                src_id - the id of the receiving processor
 @return:       none
 @exception:    none
*/

void matrix_to_tile(	int * M,
                        int size_of_m,
                        int * T,
                        int size_of_t,
                        int src_id,
                        int sqrt_p) {
    // if matrix can be sliced evenly into this tile...
    if(size_of_m % size_of_t == 0) {
        int tile_size = size_of_t;
        int row_loc = src_id * tile_size + (src_id/sqrt_p) * size_of_m * (tile_size-1);
        for(int i = 0; i < tile_size; ++i) {
            for(int j = 0; j < tile_size; ++j) {
                T[i*tile_size+j] = M[row_loc+j];
            }
            row_loc = row_loc + size_of_m;
        }
    }
}


/**
 @brief:        tile_to_matrix
 @details:      places a t x t tile into a larger n x n matrix
 @param:        M   -   matrix storage pointer
                size_of_m - dimensions of M
                T   -   tile storage pointer
                size_of_t - dimensions of T
                src_id - the id of the receiving processor
 @return:       none
 @exception:    none
*/

void tile_to_matrix(	int * M,
                        int size_of_m,
                        int * T,
                        int size_of_t,
                        int src_id,
                        int sqrt_p) {
    // if matrix can be sliced evenly into this tile...
    if(size_of_m % size_of_t == 0) {
        int tile_size = size_of_t;
        int row_loc = src_id * tile_size + (src_id/sqrt_p) * size_of_m * (tile_size-1);
        for(int i = 0; i < tile_size; ++i) {
            for(int j = 0; j < tile_size; ++j) {
                M[row_loc+j] = T[i*tile_size+j];
            }
            row_loc = row_loc + size_of_m;
        }
    }
}

/**
 @brief:        init_hilbert_matrix
 @details:      generates a n x n Hilbert matrix and stores it
 @param:        M   -   matrix storage pointer
                n   -   dimensions of M
 @return:       none
 @exception:    none
*/
/*
void init_hilbert_matrix(float ** M, int n) {
	int i = 1, j = 1;
	for(i;i<=n;++i) {
		for(j;j<=n;++j) {
			M[i][j] = 1/(i+j-1);
		}
	}
}

/**
 @brief:        shift_rows_left
 @details:      flattens the matrix rows into a send buffer, sends the rows to the
                "left" processor, receives flattened rows from "right" processor,
                unflattens the receive buffer, and stores the result
 @param:        M   -   matrix storage pointer
                snd -   pointer to send buffer
                rcv -   pointer to receive buffer
                src_id - the id of the sending processor
                n   -   dimensions of M
                w_size -    the total number processors
 @return:       none
 @exception:    none
*/
/*
void shift_rows_left(int n,
                     float M[n][n],
                     float snd[],
                     float rcv[],
                     int src_id,
                     int w_size) {

	// determine which processors are "left" and "right"...
	int snd_to = src_id - 1;
	int rcv_from = src_id + 1;
	if(snd_to < 0) snd_to = w_size - 1;
	if(rcv_from > w_size) rcv_from = 0;

	// flatten rows...
	flatten(n, M, snd);

	// send and receive rows...
	MPI_Sendrecv(snd, n*n, MPI_FLOAT, snd_to, snd_to,
				rcv, n*n, MPI_FLOAT, rcv_from, rcv_from,
				MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	// unflatten rows...
	unflatten(n, M, rcv);
}

/**
 @brief:        shift_cols_up
 @details:      flattens the matrix columns into a send buffer, sends the columns
                to the "above" processor, receives flattened columns from "below"
                processor, unflattens the receive buffer, and stores the result
 @param:        M   -   matrix storage pointer
                snd -   pointer to send buffer
                rcv -   pointer to receive buffer
                src_id - the id of the sending processor
                n   -   dimensions of M
 @return:       none
 @exception:    none
*/
/*
void shift_cols_up(int n,
                     float M[n][n],
                     float snd[],
                     float rcv[],
                     int src_id,
                     int shft_amt,
                     int w_size) {
	int dst_id = src_id - 1;
	if(src_id < 0) src_id = w_size - 1;

	// determine which processors are "above" and "below"...
	int snd_to = src_id - shift_amt;
	int rcv_from = src_id + shift_amt;
	if(snd_to < 0) snd_to = w_size - shift_amt;
	if(rcv_from > w_size) rcv_from = 0;

	// flatten columns...
	flatten(n, M, snd);

	// send and receive columns...
	MPI_Sendrecv(snd, n*n, MPI_FLOAT, snd_to, snd_to,
		rcv, n*n, MPI_FLOAT, rcv_from, rcv_from,
		MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	// unflatten columns...
	unflatten(n, M, rcv);
}

void flatten(int n, float M[n][n], float buf[]) {
	int i = 0;
	for(i;i<n;++i) {
		memcpy(buf[i*n], M[i][0], sizeof(float)*n);
	}
}

void unflatten(int n, float M[n][n], float buf[]) {
	int i = 0;
	for(i;i<n;++i) {
		memcpy(M[i][0], buf[i*n], sizeof(float)*n);
	}
}
*/