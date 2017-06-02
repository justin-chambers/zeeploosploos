#include "Utilities.h"

#define MAX_SIZE	1000

int main(int argc, char ** argv) {
    Utilities util;
    util.handle_argc(argc, 3);
    int * buffer = NULL;
    int i = 0;
    int number_of_points = util.get_positive_integer_from_cmd_line(argv[1]);
    long random_num = 0;
    buffer = new int[number_of_points+1];
    buffer[0] = number_of_points;
    //std::cout << number_of_points << std::endl;
    for(i=0;i<number_of_points;i++){
        random_num = random();
        buffer[i+1] = (int)(random_num % MAX_SIZE);
        //std::cout << (random_num % MAX_SIZE) << std::endl;
    }
    util.write_one_dee_integer_array_to_file(argv[2],buffer,number_of_points+1);
    return 0;
}