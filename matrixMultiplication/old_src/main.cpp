/**
 @file: mm_seq.cpp
 @author: Justin Chambers
 @date: 20170410
 @version: 1.0
 @brief: Main program for matrix multiplication (sequential)
*/

#include "Utilities.h"

#define ROW_MAX     1000000
#define COL_MAX     1000000
#define UPPER       100
#define CONTOLLER   0
#define DISPLAY_ARRAY   0
#define WRITE_OUTPUT    1

void multiply_matrices(int ** A,
                       const int mA,
                       const int nA,
                       int ** B,
                       const int mB,
                       const int nB,
                       int ** C );

/**
 @brief:        main
 @details:      reads in a matrix size (N x N only), randomly fills two matrices, multiplies them, stores the
                result in a third matrix.
 @param:        argc -	the number of commandline arguments
                argv -	pointer to the arguments
 			argv[0] is the name of the binary
			argv[1] is the size of the matrix
			argv[2] is the number of tiles (later, when we parallelize) to break the matrix into
			argv[3] is the name of the output file to be written (debug only)
 @return:       none
 @exception:    none
*/

int main(int argc, char ** argv) {
    int ** A = NULL;
    int ** B = NULL;
    int ** C = NULL;
    int world_rank = 0;
    int i = 0, j = 0;

    Utilities util;

    // handle command line
    util.handle_argc(argc,4);

    // get the problem size
    int n = util.get_positive_integer_from_cmd_line(argv[1]);
    int p = util.get_positive_integer_from_cmd_line(argv[2]);

    // allocate memory
    A = new int * [ROW_MAX];
    B = new int * [ROW_MAX];
    C = new int * [ROW_MAX];

    for(i=0;i<n;++i) {
        A[i] = new int[COL_MAX];
        B[i] = new int[COL_MAX];
        C[i] = new int[COL_MAX];
    }

    // initialize data...
    std::srand(world_rank);
    for(i=0;i<n;++i) {
        for(j=0;j<n;++j) {
            A[i][j] = rand() % UPPER;
            B[i][j] = rand() % UPPER;
        }
    }

    // multiply matrices
    multiply_matrices(A, n, n, B, n, n, C);

#if DISPLAY_ARRAY
    for(i=0;i<n;++i) {
        for(j=0;j<n;++j) {
            std::cout << C[i][j] << " ";
            //std::cout << A[i][j] << " ";
        }
        std::cout << std::endl;
    }
#endif

#if WRITE_OUTPUT
    std::string output_fn = argv[3];
    util.write_two_dee_integer_array_to_file(argv[3], C, n, n);
#endif

    // finalize

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

void multiply_matrices(int * A[],
                       const int mA,
                       const int nA,
                       int * B[],
                       const int mB,
                       const int nB,
                       int * C[] ) {
    if(mA == nB) {
        int i = 0, j = 0, k = 0;
        for(i=0; i<mA;++i) { // rows of A
            for(j=0; j<nB;++j) { // columns of B
                //compute the index
                for(k=0; k<mB;++k) { // rows of B
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }
    }
}