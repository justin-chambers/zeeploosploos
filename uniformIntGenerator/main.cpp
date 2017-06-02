#include "Utilities.h"

#define MAX_SIZE	1000
#define MIN_SIZE    0

void place_rand(const int x, int* dst, const int dst_size);
void place_lin(const int x, int* dst, const int dst_size);

int main(int argc, char ** argv) {
    // handle command line...
    Utilities util;
    util.handle_argc(argc, 3);

    // variables...
    std::string out_fn = "output";
    int * buffer = NULL;
    int ** subintervals = NULL;
    int i = 0;
    int number_of_points = util.get_positive_integer_from_cmd_line(argv[1]);
    int uni_distr = util.get_positive_integer_from_cmd_line(argv[2]);
    long random_num = 0;

    // allocate memory for buffer, subintervals...
    buffer = new int[number_of_points];
    for(int i=0; i<number_of_points; ++i) buffer[i] = MIN_SIZE-1;
    subintervals = new int * [uni_distr];
    int interval_size = 0;
    if((number_of_points % uni_distr) == 0) interval_size = number_of_points/uni_distr;
    else interval_size = number_of_points/uni_distr + 1;
    for(int i=0; i<uni_distr; ++i) {
        subintervals[i] = new int[interval_size];
    }

    // generate uniformly distributed pseudo-random numbers...
    int determinant = MAX_SIZE/uni_distr;
    std::srand(std::time(NULL));
    for(i=0;i<uni_distr;i++){
        for(int j=0; j<interval_size; ++j) {
            random_num = rand();
            subintervals[i][j] = (int)(random_num % determinant) + (determinant*i);
        }
    }

    // randomly place numbers in final buffer...
    for(int i=0; i<uni_distr; ++i){
        for(int j=0; j<interval_size; ++j) {
            place_rand(subintervals[i][j],buffer,number_of_points);
        }
    }

    // write buffer to file...
    out_fn = out_fn + "_" + std::to_string(number_of_points);
    out_fn = out_fn + "_" + std::to_string(uni_distr) + ".txt";
    util.write_one_dee_integer_array_to_file(out_fn.c_str(),buffer,number_of_points);
    return 0;
}

void place_rand(const int x, int* dst, const int dst_size) {
    long rand_loc = rand();
    int loc = (int)(rand_loc % dst_size);
    int count = 1;
    while((dst[loc]>MIN_SIZE-1) && (count < dst_size)) {
        rand_loc = random();
        loc = (int)(rand_loc % dst_size);
        ++count;
    }
    if(count==dst_size) {
        std::cout << "Warning: " << x << " could not be placed randomly. Using linear placement...\n";
        place_lin(x,dst,dst_size);
    } else {
        dst[loc] = x;
    }
}

void place_lin(const int x, int* dst, const int dst_size) {
    int loc = 0;
    while(dst[loc] > MIN_SIZE-1) {
        ++loc;
    }
    if(loc==dst_size) std::cout << "Warning: " << x << " could not be placed linearly as buffer is full.\n";
    else dst[loc] = x;
}