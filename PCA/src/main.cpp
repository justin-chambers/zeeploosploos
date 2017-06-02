#include "PCAApplication.h"

int main( int argc, char *argv[] )
{
    bool validTestValue = true;
    char * argvPtr = argv[2];

    // Command line parsing...
    if( argv[2][0] == '-' && argv[2][1] == 't' && argv[2][2] == '\0' )
    {
        // Valid training argument detected...
        PCAApplication pcaTrainingMode(argv, argc, 1);
    } else
    {
        while( *argvPtr != '\0' && validTestValue )
        {
            if( !isdigit(*argvPtr) )
            {
                validTestValue = false;
            }
            ++argvPtr;
        }

        if( validTestValue )
        {
            std::cout << "PCA EXPERIMENT MODE:\n"
            << "======================================================================================\n\n"
            << " Number of experiments: " << argv[2] << "\n\n"
            << "======================================================================================\n\n";

            std::string experimentBound;
            int maxBound = atoi(argv[2]);
            for( int i = maxBound; i >= 0; --i )
            {
                experimentBound = std::to_string(i);
                argv[2] = (char *)experimentBound.c_str();
                std::cout << "Experiment #" << experimentBound << "...\n";
                PCAApplication pcaTest(argv, argc, 0);
            }
        } else
        {
            PCAApplication usage;
        }
    }

    return 0;
}