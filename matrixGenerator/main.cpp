#include "Utilities.h"
#include <vector>

#define MIN_MATRIX_DIMENSION    120
#define MAX_MATRIX_DIMENSION    10000
#define MIN_SAMPLE_DIMENSION    60
#define MAX_SAMPLE_DIMENSION    120
#define GENERATE_SAMPLE_SET     1
#define GENERATE_FULL_SET       0

/**
 * @brief main
 * @description reads file names and dimensions, generates matrices, writes them to file
 * @param argc
 * @param argv - pointer to arg string
 *             - argv[1] is the integer representing the n x n dimensions of argv[2], argv[3]
 *             - argv[2] is the name of the first matrix to be generated
 *             - argv[3] is the name of the second matrix to be generated
 * @return 0
 */

void generateMatrix(int * M, int dimM, int starting_val);

int main(int argc, char ** argv) {
    Utilities util;
    util.handle_argc(argc, 4);
    int dim = 0;
    std::vector<int> dimensions;

#if GENERATE_FULL_SET
    for(int i = MIN_MATRIX_DIMENSION; i <= MAX_MATRIX_DIMENSION; i = i+60) {
        dimensions.push_back(i);
    }
#endif

#if GENERATE_SAMPLE_SET
    for(int i = 3; i <= 4; ++i) {
        dimensions.push_back(i);
    }

    /*
    for(int i = MIN_SAMPLE_DIMENSION; i <= MAX_SAMPLE_DIMENSION; i = i+60) {
        dimensions.push_back(i);
    }
     */
#endif

    for(int i = 0; i < (int)dimensions.size(); ++i) {
        dim = dimensions[i];
        int * A = new int[dim*dim];
        int * B = new int[dim*dim];
        generateMatrix(A, dim, 1);
        generateMatrix(B, dim, dim*dim+1);
        std::string matA = "a";
        std::string matB = "b";
        matA = matA + std::to_string(dim) + ".txt";
        matB = matB + std::to_string(dim) + ".txt";
        util.write_flattened_two_dee_integer_array_to_file(matA.c_str(), A, dim, dim);
        util.write_flattened_two_dee_integer_array_to_file(matB.c_str(), B, dim, dim);
    }
    return 0;
}

void generateMatrix(int * M, int dimM, int starting_val) {
    for(int i = 0; i < dimM*dimM; ++i) {
        M[i] = starting_val;
        ++starting_val;
    }
}