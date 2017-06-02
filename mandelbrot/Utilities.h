#ifndef MANDELBROT_UTILITIES_H
#define MANDELBROT_UTILITIES_H

#import <stdexcept>
#include <iostream>
#include <cstring>
#include <string>

class Utilities {
public:
    void handle_argc(int argc, int valid_argc);
    int get_positive_integer_from_cmd_line(char * arg);
    void display_one_dee_array(int * arr, int n);
    void display_two_dee_array(int ** arr, int m, int n);
    void display_one_dee_array(unsigned char * arr, int n);
    void initialize_one_dee_array(unsigned char *arr,
                                  int m,
                                  unsigned char init_val);
    void display_two_dee_array(unsigned char ** arr, int m, int n);
    void initialize_two_dee_array(unsigned char **arr,
                                             int m,
                                             int n,
                                             unsigned char init_val);
};

#endif //MANDELBROT_UTILITIES_H
