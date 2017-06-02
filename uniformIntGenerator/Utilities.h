/**
 @file: Utilities.h
 @author: Justin Chambers
 @date: 20170310
 @version: 1.0
 @brief: Specification file for Utilities class
*/

#ifndef CPP_UTILITIES_H
#define CPP_UTILITIES_H

#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <stdexcept>
#include <exception>

class Utilities {
public:
    void handle_argc(int argc, int valid_argc);
    int get_positive_integer_from_cmd_line(char * arg);
    void display_one_dee_array(int * arr, int n);
    void display_two_dee_array(int ** arr, int m, int n);
    void display_one_dee_array(unsigned char * arr, int n);
    void display_two_dee_array(unsigned char ** arr, int m, int n);
    int get_array_size_from_file(const char * filename);
    void read_file_to_one_dee_integer_array(const char * filename, int * arr, int n);
    void write_one_dee_integer_array_to_file(const char * filename, int * arr, int n);
};

#endif //CPP_UTILITIES_H

