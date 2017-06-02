#include <iostream>

/**
 * @brief main
 * @description This program takes an exponent, raises 2 by that power, then sums up the digits of that number.
 * @param argc - the count of command line arguments
 * @param argv - the command strings
 *               [0] program name
 *               [1] exponent
 * @return 0 upon success
 */

int main(int argc, char ** argv) {
    int power = 0, sum = 0;
    long int x = 1;
    power = atoi(argv[1]);
    for(int i = 0; i < power; ++i) {
        x *= 2;
    }
    long int quotient = x;
    while(quotient > 0) {
        sum += quotient % 10;
        quotient /= 10;
    }
    std::cout << "x = " << x << std::endl;
    std::cout << "sum of digits = " << sum << std::endl;
    return 0;
}