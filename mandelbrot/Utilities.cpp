#include "Utilities.h"

void Utilities::handle_argc(int argc, int valid_argc) {
    if(argc != valid_argc) std::cerr << "Program takes " << valid_argc << " ARGS. Exiting.\n";
}

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

void Utilities::display_one_dee_array(int *arr, int n) {
    for(int i=0;i<n;++i) std::cout << arr[i] << " ";
    std::cout << "\n";
}

void Utilities::display_two_dee_array(int **arr, int m, int n) {
    for(int i=0;i<m;++i) display_one_dee_array(arr[i],n);
}

void Utilities::display_one_dee_array(unsigned char *arr, int n) {
    for(int i=0;i<n;++i) std::cout << arr[i] << " ";
    std::cout << "\n";
}

void Utilities::display_two_dee_array(unsigned char **arr, int m, int n) {
    for(int i=0;i<m;++i) display_one_dee_array(arr[i],n);
}

void Utilities::initialize_two_dee_array(unsigned char **arr,
                                         int m,
                                         int n,
                                         unsigned char init_val) {
    arr = new unsigned char * [n];
    for(int i = 0; i < n; ++i) arr[i] = new unsigned char[m];
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < m; ++j) {
            arr[i][j] = init_val;
        }
    }
}

void Utilities::initialize_one_dee_array(unsigned char *arr, int m, unsigned char init_val) {
    arr = new unsigned char [m];
    for(int i = 0; i < m; ++i) arr[i] = init_val;
}
