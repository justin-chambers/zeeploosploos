/**
 @file: Utilities.h
 @author: Justin Chambers
 @date: 20170310
 @version: 1.0
 @brief: Specification file for Utilites class
*/

#ifndef MANDELBROT_UTILITIES_H
#define MANDELBROT_UTILITIES_H

#include <iostream>
#include <cstring>
#include <string>
#include <cstdlib>
#include <stdexcept>

class Utilities {
public:
    void handle_argc(int argc, int valid_argc);
    int get_positive_integer_from_cmd_line(char * arg);
    void display_one_dee_array(int * arr, int n);
    void display_two_dee_array(int ** arr, int m, int n);
    void display_one_dee_array(unsigned char * arr, int n);
    void display_two_dee_array(unsigned char ** arr, int m, int n);
};

#endif //MANDELBROT_UTILITIES_H
