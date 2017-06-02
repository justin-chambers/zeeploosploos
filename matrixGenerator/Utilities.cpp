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
    int j = 0;
    for(j=0;j<n;++j) std::cout << arr[j] << " ";
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
    int i = 0;
    for(i=0;i<m;++i) display_one_dee_array(arr[i],n);
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
    int i = 0;
    for(i=0;i<n;++i) std::cout << arr[i] << " ";
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
    int i = 0;
    for(i=0;i<m;++i) display_one_dee_array(arr[i],n);
}

/**
@brief:        get_array_size_from_file
@details:      opens a file, gets the size of the array from the file, closes the file
@param:        filename - pointer to the file
@return:       n - size of the array
@exception:    unable to open filename, n is not positive
*/

int Utilities::get_array_size_from_file(const char *filename) {
    std::fstream fin;
    try {
        fin.open(filename, std::ios::in);
    }
    catch (std::exception& e) {
        std::cout << e.what() << std::endl;
    }
    int n;
    fin >> n;
    fin.close();
    if(n<=0) std::cerr << "Array size must be positive. " << n << " was given.\n";
    return n;
}

/**
@brief:        read_file_to_one_dee_integer_array
@details:      reads the contents of a file (first entry is the count of subsequent entries) into an integer array
@param:        filename - pointer to the name of the file to read
               arr   - pointer to the integer array
               n     - size of the integer array
@return:       none
@exception:    unable to open file stream object
*/

void Utilities::read_file_to_one_dee_integer_array(const char *filename, int *arr, int n) {
    std::fstream fin;
    try {
        fin.open(filename, std::ios::in);
    }
    catch (std::exception& e) {
        std::cout << e.what() << std::endl;
    }
    int skip;
    fin >> skip;
    for(int i = 0; i < n; ++i) fin >> arr[i];
    fin.close();
}

/**
@brief:        write_one_dee_integer_array_to_file
@details:      write the contents of an integer array (first entry is the count of subsequent entries) into a file
@param:        filename - pointer to the name of the file to written
               arr   - pointer to the integer array
               n     - size of the integer array
@return:       none
@exception:    unable to open file stream object
*/

void Utilities::write_one_dee_integer_array_to_file(const char *filename, int *arr, int n) {
    std::fstream fout;
    try {
        fout.open(filename, std::ios::out);
    }
    catch (std::exception& e) {
        std::cout << e.what() << std::endl;
    }
    fout << n << std::endl;
    for(int i=0; i < n; ++i) fout << arr[i] << "\n";
    fout.close();
}

/**
@brief:        write_two_dee_integer_array_to_file
@details:      write the contents of an integer array (first two entries are the dimensions) into a file
@param:        filename - pointer to the name of the file to written
               arr   - pointer to the integer array
               n     - size of the integer array
@return:       none
@exception:    unable to open file stream object
*/

void Utilities::write_two_dee_integer_array_to_file(const char * filename, int ** arr, int m, int n) {
    std::fstream fout;
    try {
        fout.open(filename, std::ios::out);
    }
    catch (std::exception& e) {
        std::cout << e.what() << std::endl;
    }
    fout << m << " " << n << std::endl;
    int i = 0, j = 0;
    for(i = 0; i < m; ++i) {
        for (j = 0; j < n; ++j) {
            fout << arr[i][j] << " ";
        }
        fout << "\n";
    }
    fout.close();
}

/**
@brief:        write_two_dee_integer_array_to_file
@details:      write the contents of an integer array (first two entries are the dimensions) into a file
@param:        filename - pointer to the name of the file to written
               arr   - pointer to the integer array
               n     - size of the integer array
@return:       none
@exception:    unable to open file stream object
*/

void Utilities::write_flattened_two_dee_integer_array_to_file(const char *filename, int *arr, int m, int n) {
    std::fstream fout;
    try {
        fout.open(filename, std::ios::out);
    }
    catch (std::exception& e) {
        std::cout << e.what() << std::endl;
    }
    if(m == n) {
        fout << m << std::endl;
    } else {
        fout << m << " " << n << std::endl;
    }
    int i = 0, j = 0;
    for(i = 0; i < m; ++i) {
        for (j = 0; j < n; ++j) {
            fout << arr[i*m+j] << " ";
        }
        fout << "\n";
    }
    fout.close();
}
