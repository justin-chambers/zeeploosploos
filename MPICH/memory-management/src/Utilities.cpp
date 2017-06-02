/**
 @file: Utilities.cpp
 @author: Justin Chambers
 @date: 20170310
 @version: 1.0
 @brief: Implementation file for Utilites class
*/

#include "Utilities.h"

/**
 @brief:        handle_argc
 @details:      Determines whether or not the program has the proper command line args
 @param:        argc - the number of arguments provided by the user
                valid_argc - the number of arguments necessary
 @return:       none
 @exception:    standard error thrown if arc and valid_argc are not equal
*/

void Utilities::handle_argc(int argc, int valid_argc) {
    if(argc != valid_argc) std::cerr << "Program takes " << valid_argc << " ARGS. Exiting.\n";
}

/**
 @brief:        get_positive_integer_from_cmd_line
 @details:      converts a given command line to integer and tests whether it is
                positive or not. the number is returned if it is positive, else
                the function throws an error.
 @param:        arg - pointer to c-string
 @return:       i   - the converted integer
 @exception:    logic error is throw if a) conversion fails or b) number is
                not positive
*/

int Utilities::get_positive_integer_from_cmd_line(char * arg) {
    int i = 0;
    try {
        i = atoi(arg);
    } catch (std::logic_error& e) {
        std::cerr << e.what() << std::endl;
    }
    if(i<=0) std::cerr << "Argument must be a positive integer. " << arg << " was given. Exiting.";
    return i;
}

/**
 @brief:        display_one_dee_array
 @details:      displays the contents of a one-dimensional array of integers
 @param:        arr - pointer to the integer array
                n   - the size of the array
 @return:       none
 @exception:    none
*/

void Utilities::display_one_dee_array(int *arr, int n) {
    for(int i=0;i<n;++i) std::cout << arr[i] << " ";
    std::cout << "\n";
}

/**
 @brief:        display_two_dee_array
 @details:      displays the contents of a two-dimensional array of integers
 @param:        arr - pointer to the integer array
                n   - the size of the array
 @return:       none
 @exception:    none
*/

void Utilities::display_two_dee_array(int **arr, int m, int n) {
    for(int i=0;i<m;++i) display_one_dee_array(arr[i],n);
}

/**
@brief:        display_one_dee_array
@details:      displays the contents of a one-dimensional array of unsigned chars
@param:        arr - pointer to the integer array
               n   - the size of the array
@return:       none
@exception:    none
*/

void Utilities::display_one_dee_array(unsigned char *arr, int n) {
    for(int i=0;i<n;++i) std::cout << arr[i] << " ";
    std::cout << "\n";
}

/**
@brief:        display_two_dee_array
@details:      displays the contents of a two-dimensional array of unsigend chars
@param:        arr - pointer to the integer array
               n   - the size of the array
@return:       none
@exception:    none
*/

void Utilities::display_two_dee_array(unsigned char **arr, int m, int n) {
    for(int i=0;i<m;++i) display_one_dee_array(arr[i],n);
}
