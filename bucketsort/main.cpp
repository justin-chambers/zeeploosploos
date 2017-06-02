#include "Utilities.h"

#define UPPER   1000
#define LOWER   0

int compare(const void*a, const void*b);
void partition(const int x, int * bucket, const int size);

/**
 @brief:        main
 @details:      reads in a file with unsorted numbers, bucketsorts the numbers, writes the sorted numbers to a file.
                bucketsorting is timed so that it can be compared to parallelized version.
 @param:        argc - the number of commandline arguments
                argv - pointer to the arguments
                       argv[0] is the name of the binary
                       argv[1] is the name of the file to be read in
                       argv[2] is the number of buckets to be used (in the parallelized version, this is also the number
                               of processors)
 @return:       none
 @exception:    none
*/

int main(int argc, char ** argv) {
    // handle the command line...
    Utilities util;
    util.handle_argc(argc, 3);

    // variables...
    std::string output_fn;
    int n = util.get_array_size_from_file((const char *)argv[1]);
    int bucket_count = util.get_positive_integer_from_cmd_line(argv[2]);
    int bucket_size = 0;
    int determinant = 0;
    int * buffer = NULL;
    int ** buckets = NULL;

    // allocate memory to buffer and buckets...
    buffer = new int[n];
    if(bucket_count > n) {
        std::cout << "Warning: number of buckets (" << bucket_count << ") exceeds the count of numbers in the file ("
                                                               << n << "). Fixing...\n";
        bucket_count = n;
    }
    if(n % bucket_count == 0) bucket_size = (n/bucket_count);
    else bucket_size = (n/bucket_count)+1;
    buckets = new int * [bucket_count];
    for(int i=0; i<bucket_count;++i) buckets[i] = new int[bucket_size];

    // initialize buffer...
    for(int i=0; i<n; ++i) buffer[i] = LOWER-1;

    // initialize buckets...
    for(int i=0; i<bucket_count; ++i) {
        for(int j=0; j<bucket_size; ++j) {
            buckets[i][j] = LOWER-1;
        }
    }

    // read in the numbers...
    util.read_file_to_one_dee_integer_array((const char *)argv[1], buffer, n);

    // start timer...

    // partition into buckets...
    determinant = UPPER/bucket_count;
    int x = 0;
    for(int i=0; i<n; ++i) {
        x = buffer[i];
        partition(x,buckets[x/determinant], bucket_size);
    }

    // sort buckets...
    for(int i=0; i<bucket_count; ++i) {
        std::qsort(buckets[i], (size_t)bucket_size, sizeof(int), compare);
    }

    // stop timer...

    // merge buckets with buffer...
    for(int i=0; i<bucket_count; ++i) {
        std::memcpy(buffer+(bucket_size*i),buckets[i],(bucket_size*sizeof(int)));
    }

    std::cout << "Merged\n";
    for(int i=0;i<n;++i) std::cout << buffer[i] <<std::endl;

    // write buffer to file...
    output_fn = argv[1];
    output_fn = output_fn + "_sorted.txt";
    util.write_one_dee_integer_array_to_file(output_fn.c_str(), buffer, n);

    // display bucket_count, n and computation time...
    std::cout << bucket_count << "," << n << ", (MPI_Wtime)" << std::endl;
    return 0;
}

/**
 @brief:        compare
 @details:      compares two integers and determines which one is smaller (for sorting)
 @param:        a - pointer to the first integer
                b - pointer to the second integer
 @return:       integer signal used to determine the order in which a, b are to sorted
 @exception:    none
*/

int compare(const void*a, const void*b)
{
    return (*(int*)a-*(int*)b);
}

/**
 @brief:        partition
 @details:      places the number into the next open spot in the partition
 @param:        x - the integer to be placed
                bucket - pointer to partition into which x is to be placed
 @return:       none
 @exception:    none
*/

void partition(const int x, int * bucket, const int size) {
    int empty_flag = LOWER-1;
    int i=0;
    while(bucket[i] > empty_flag) {
        ++i;
    }
    if(i<size) bucket[i] = x;
    else std::cout << "Warning: " << x << " was not partitioned because the given partition is full.\n";
}