//Precompiler Directives

//Libraries
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <dirent.h>
#include <armadillo>

//Classes
#include "GrayscaleImage.h"
#include "EigenImage.h"
#include "EigenspaceWeights.h"

//Free Function Prototypes
void DisplayUsage();
std::vector<GrayscaleImage> InitializeDatabase( std::string & directoryName );
std::vector<EigenImage> InitializeEigenDatabase( std::string & directoryName );
std::vector<GrayscaleImage> InitializePhiDatabase( std::vector<GrayscaleImage>& aRawDB,
                                                  GrayscaleImage& anAvgImage );
std::vector<EigenImage> InitializeEigenDatabase( arma::vec& anEigValVec,
                                                 arma::mat& anEigVecMatrix,
                                                 GrayscaleImage& anAvgImg );
std::vector<EigenspaceWeights> InitializeEigenspaceWeightsDB ( const std::string & directoryName );
void ReadPGMToObject( const char * fname, GrayscaleImage & obj );
void WriteObjectToPGM( const char * fname, const GrayscaleImage & src );
void WriteDBtoPGMs( const std::vector<GrayscaleImage>& src, const int bound );
void WriteDBtoPGMs( const std::vector<EigenImage>& src, const int bound );
void WriteWeightsForImg( const char *, arma::vec& weightVec );
void ComputeAvgImage( const std::vector<GrayscaleImage> & src,
                      GrayscaleImage & dst );

//Main Function Implementation
int main( int argc, char *argv[] )
{
    // Variables
    std::vector<GrayscaleImage> imageDB;
    std::vector<EigenImage> eigenFDB;
    std::vector<GrayscaleImage> phiDB;
    std::vector<EigenspaceWeights> omegaDB;
    double * avgImgVec = NULL;
    double retainPercent = 0, eigenValSum = 0;
    int eigenVectorsToUse = 0;

    // Command-line branching...
    if( argc > 3 )
    {
        DisplayUsage();
        return 0;
    } else
    {
        std::string tmp = argv[1];

        // Step 1 & 2 ... Capture face data and store as vectors (and scale data)...
        std::cout << "Preprocessing: Attempting to create image database (this may take some time)...";
        imageDB = InitializeDatabase(tmp);

        std::cout << "done.\nDatabase formed from directory " << argv[1] << ".\n";

        if(argc == 3)
        {
            if(argv[2][0] == '-' && argv[2][1] == 't' && argv[2][2] == '\0')
            {
                std::cout << "Training mode: Generating the eigenspace (this may take some time)...";

                // Step 3... Compute & write average face...
                GrayscaleImage avgImage;
                avgImage.SetName("eigenspace/__average_face.pgm");
                avgImage.SetWidth(imageDB[0].GetWidth());
                avgImage.SetHeight(imageDB[0].GetHeight());
                ComputeAvgImage(imageDB, avgImage);
                WriteObjectToPGM(avgImage.GetName().c_str(),avgImage);

                // Step 4... Subtract every face from the average face...
                phiDB = InitializePhiDatabase(imageDB, avgImage);


                // Spot check: can we rebuild imageDB_0 from avgImg + phiDB_0?
                //arma::colvec compositeVec(avgImage.GetScaledValVec());
                //arma::colvec phiVec(phiDB[0].GetScaledValVec());
                //arma::colvec originalImgVec(imageDB[0].GetScaledValVec());
                //compositeVec = compositeVec + phiVec;
                // Yes...

                // Step 5... Compute covariance matrix C...

                // Populate A
                arma::mat A((avgImage.GetHeight()*avgImage.GetWidth()),1, arma::fill::zeros);
                for(int i = 0; i < (int)phiDB.size(); ++i)
                {
                    arma::colvec tmp(phiDB[i].GetScaledValVec());
                    A.insert_cols(i, tmp);
                }
                A.shed_col(phiDB.size());

                // C = (1/M) * A^T * A
                arma::mat C = (1/(double)phiDB.size()) * (A.t() * A);

                // Solve and sort C (note: ascending order here, need descending)...
                arma::vec eigVal;
                arma::mat eigVec;
                eig_sym(eigVal, eigVec, C);

                // Generate eigenFaces (sorts to descending order here)...
                arma::mat U = A*eigVec;
                eigenFDB = InitializeEigenDatabase(eigVal, U, avgImage);

                std::cout << "done.\nWhat percent of image information would you like to preserve? (ex. 90.2): ";
                std::cin >> retainPercent;
                while(retainPercent > 100 || retainPercent < 0)
                {
                    std::cout << "\nPlease select a value between 0 and 100: ";
                    std::cin >> retainPercent;
                }

                // Calculate # of eigenFaces to use...
                if( retainPercent < 100 )
                {
                    eigenValSum = arma::sum(eigVal);
                    eigenValSum =  ( retainPercent/100 ) * eigenValSum;
                    while( eigenValSum > 0 && eigenVectorsToUse < (int)eigVec.n_cols )
                    {
                        eigenValSum -= eigVal.at((int)eigVal.n_rows - eigenVectorsToUse - 1);
                        ++eigenVectorsToUse;
                    }

                    std::cout << eigenVectorsToUse << " of " << eigVec.n_cols << " eigenfaces selected.\n"
                              << "Finalizing eigenspace (this may take some time)...";

                } else
                {
                    eigenVectorsToUse = (int)eigVec.n_cols;
                }

                // Update U (i.e. sort U, normalize U, etc.)...
                U.shed_cols(0, U.n_cols-1);
                for ( int i = 0; i < eigenVectorsToUse; ++i )
                {
                    double eigenvalueCheck = eigenFDB[i].GetEigVal();
                    arma::colvec tmpEV(eigenFDB[i].GetScaledValVec());
                    U.insert_cols(i,tmpEV);
                }
                U = arma::normalise(U);

                // Spot check: does || u_I || = 1?
                //arma::colvec normalizedUVec(U.col(0));
                //double normVal = arma::norm(normalizedUVec);
                // Yes...

                // Calculate & store weights for DB...for each face, Omega_i = U^T * Phi_i
                for( int i = 0; i < (int)phiDB.size(); ++i )
                {
                    arma::colvec phi(phiDB[i].GetScaledValVec());
                    arma::colvec omega_i = U.t()*phi;
                    omega_i = arma::normalise(omega_i);
                    WriteWeightsForImg(imageDB[i].GetName().c_str(), omega_i);

                    // Spot check: can we rebuild phiDB_i from sum(omega_k * u_k) + avgImage?
                    /*
                    arma::colvec mean(avgImage.GetScaledValVec());
                    arma::colvec phiHat_i = U * Omega_i;
                    phiHat_i = phiHat_i + mean;
                    GrayscaleImage rebuilt(avgImage);
                    rebuilt.SetName("rebuiltImage.pgm");
                    std::vector<double> tmpScaledVals;
                    for( int i = 0; i < (int)phiHat_i.n_rows; ++i )
                    {
                        tmpScaledVals.push_back(phiHat_i.at(i));
                    }
                    rebuilt.SetScaledValues(tmpScaledVals);
                    rebuilt.ConvertScaledtoRaw();
                    WriteObjectToPGM(rebuilt.GetName().c_str(), rebuilt);
                    */
                    // pretty good...

                    // Verify norm(I - I_hat) is very small...
                    if( retainPercent == 100 )
                    {
                        // build I_hat =
                        arma::vec I_hat = omega_i*U;

                        // norm(I -I_hat) should be very small...
                        double dffs = arma::norm(phi - I_hat);
                        std::cout << dffs << std::endl;
                    }
                }

                // Write eigenfaces to files...
                WriteDBtoPGMs(eigenFDB, eigenVectorsToUse);

                std::cout << "done.\n";

            } else {
                std::cout << "Error: Invalid training flag detected. Exiting.\n";
                DisplayUsage();
                return 0;
            }
        } else
        {
            //Testing mode...
            std::cout << "Testing mode: Loading the eigenspace (this may take some time)...";


            // Build average image...
            GrayscaleImage avgImage;
            ReadPGMToObject("eigenspace/__average_face.pgm", avgImage);
            avgImage.ScaleRawValues();

            // Build phiDB...
            phiDB = InitializePhiDatabase(imageDB, avgImage);

            // Build eigenFDB...
            std::string tmpEigenDirName = "eigenspace";
            eigenFDB = InitializeEigenDatabase(tmpEigenDirName);
            
            // Build U...
            arma::mat U((avgImage.GetHeight()*avgImage.GetWidth()),1, arma::fill::zeros);
            for(int i = 0; i < (int)eigenFDB.size(); ++i)
            {
                arma::colvec tmp(eigenFDB[i].GetScaledValVec());
                U.insert_cols(i, tmp);
            }
            U.shed_col(eigenFDB.size());
            U = arma::normalise(U);

            // Build omegaDB (with corresponding training DB image name)...
            omegaDB = InitializeEigenspaceWeightsDB( "eigenspace" );

            double tolerance = 0;

            std::cout << "done.\nAttempting to match images. Set tolerance: ";
            std::cin >> tolerance;

            // Step 2: Phi = sum(w_j * u_j) = U (2880 x K) * Omega_i (K * 1)

            // Project each image onto eigenspace...
            for( int i = 0; i < (int)phiDB.size(); ++i )
            {
                double difs_i = INFINITY, minDifs = INFINITY;
                arma::colvec phiVec(phiDB[i].GetScaledValVec());
                std::string matchedImageName;

                // Compare image across training data...
                for( int j = 0; j < (int)omegaDB.size(); ++j )
                {
                    difs_i = INFINITY;
                    arma::colvec omega_j(omegaDB[j].GetWeights()); // is the training image wt vector
                    //omega_j = arma::normalise(omega_j); // put weights between 0 and 1...
                    arma::colvec omega_i = U.t() * phiVec; // is the projection of test image in eigenspace
                    omega_i = arma::normalise(omega_i); // put weights between 0 and 1...
                    difs_i = arma::norm(omega_i - omega_j);
                    if(difs_i < minDifs)
                    {
                        minDifs = difs_i;
                        matchedImageName = omegaDB[j].GetImageName();
                    }
                }
                if(minDifs <= tolerance)
                {
                    std::cout << phiDB[i].GetName() << " has been recognized. Closest matching training image: "
                              << matchedImageName << ".\n";
                } else
                {
                    std::cout << phiDB[i].GetName() << " has not been recognized.\n";
                }
            }

            std::cout << "Done.\n";
        }
    }
    return 0;
}


//Free Function Implementation
void DisplayUsage()
{
    std::cout << "Usage: ./PCA <directoryName> -t [training flag, optional]\n";
    std::cout << "Example: ./PCA ../faceDataDirectory -t\n";
}

std::vector<GrayscaleImage> InitializeDatabase(std::string & directoryName)
{
    std::vector<GrayscaleImage> tmpDB;
    GrayscaleImage rawData;
    struct dirent * directoryEntityPtr = NULL;
    DIR *d = NULL;
    std::string filePath, fileName;

    // Validate the directory
    d = opendir( directoryName.c_str() );
    if( d == NULL )
    {
        perror("Error: Couldn't open directory. Exiting.");
        exit(2);
    }

    // Filename loop
    while( directoryEntityPtr = readdir(d) ) {

        filePath = directoryName + "/" + directoryEntityPtr->d_name;
        fileName = directoryEntityPtr->d_name;

        // process valid file names only
        if(fileName.size() > 4)
        {
            std::string fileType;
            fileType = fileName.substr(fileName.size()-3);
            if(fileType == "pgm")
            {
                ReadPGMToObject(filePath.c_str(), rawData);
                rawData.SetName(fileName);
                rawData.ScaleRawValues();
                tmpDB.push_back(rawData);
            }
        }
    }

    closedir(d);
    return tmpDB;
}

std::vector<EigenImage> InitializeEigenDatabase( std::string & directoryName )
{
    std::vector<EigenImage> tmpDB;
    GrayscaleImage rawData;
    struct dirent * directoryEntityPtr = NULL;
    DIR *d = NULL;
    std::string filePath, fileName;

    // Validate the directory
    d = opendir( directoryName.c_str() );
    if( d == NULL )
    {
        perror("Error: Couldn't open directory. Exiting.");
        exit(2);
    }

    // Filename loop
    while( directoryEntityPtr = readdir(d) ) {

        filePath = directoryName + "/" + directoryEntityPtr->d_name;
        fileName = directoryEntityPtr->d_name;

        // process valid file names only
        if(fileName.size() > 4)
        {
            std::string fileType, fileID;
            fileType = fileName.substr(fileName.size()-3);
            fileID = fileName.substr(0,5);
            if( (fileType == "pgm") && (fileID != "__ave") )
            {
                int tmpRank = std::atoi(fileID.c_str());
                ReadPGMToObject(filePath.c_str(), rawData);
                rawData.ScaleRawValues();
                EigenImage tmpEigen(rawData);
                tmpEigen.SetName(fileName);
                tmpDB.push_back(tmpEigen);
            }
        }
    }

    closedir(d);
    return tmpDB;
}

std::vector<GrayscaleImage> InitializePhiDatabase( std::vector<GrayscaleImage> &aRawDB,
                                                  GrayscaleImage &anAvgImage )
{
    std::vector<GrayscaleImage> tmpDB;
    int * gammaRaw = NULL;
    int * psiRaw = NULL;
    int arrLength = anAvgImage.GetWidth()*anAvgImage.GetHeight();

    psiRaw = anAvgImage.GetRawValues();
    for( int i = 0; i < (int)aRawDB.size(); ++i)
    {
        GrayscaleImage tmpImg(aRawDB[i]);

        gammaRaw = tmpImg.GetRawValues();
        for( int j = 0; j < arrLength; ++j )
        {
            gammaRaw[j] -= psiRaw[j];
        }
        tmpImg.SetRawValues((const int *&)gammaRaw);
        tmpDB.push_back(tmpImg);
    }
    return tmpDB;
}

std::vector<EigenImage> InitializeEigenDatabase( arma::vec& anEigValVec,
                                                 arma::mat& anEigVecMatrix,
                                                 GrayscaleImage& anAvgImg )
{
    int j = 0;
    std::vector<EigenImage> tmpEDB;
    std::string tmpName;

    for( int i = (int)anEigVecMatrix.n_cols-1; i >= 0; --i )
    {
        if( j+1 < 10 )
            tmpName = "0000";
        else if( j+1 < 100 && j+1 >= 10 )
            tmpName = "000";
        else if( j+1 < 1000 && j+1 >= 100 )
            tmpName = "00";
        else if( j+1 < 10000 && j+1 >= 1000 )
            tmpName = "0";
        tmpName = tmpName + std::to_string(j+1) + "_ef.pgm";
        EigenImage eImg(anAvgImg);
        eImg.SetName(tmpName);
        eImg.SetEigVal(anEigValVec.at(i));
        eImg.SetEigVec(anEigVecMatrix.col(i));
        eImg.NormalizeValuesByEigVec();
        tmpEDB.push_back(eImg);
        ++j;
    }
    return tmpEDB;
}

std::vector<EigenspaceWeights> InitializeEigenspaceWeightsDB ( const std::string & directoryName )
{
    std::vector<EigenspaceWeights> tmpDB;
    EigenspaceWeights rawData;
    struct dirent * directoryEntityPtr = NULL;
    DIR *d = NULL;
    std::string filePath, fileName;

    // Validate the directory
    d = opendir( directoryName.c_str() );
    if( d == NULL )
    {
        perror("Error: Couldn't open directory. Exiting.");
        exit(2);
    }

    // Filename loop
    while( directoryEntityPtr = readdir(d) ) {

        filePath = directoryName + "/" + directoryEntityPtr->d_name;
        fileName = directoryEntityPtr->d_name;

        // process valid file names only
        if(fileName.size() > 4)
        {
            std::string fileType;
            fileType = fileName.substr(fileName.size()-3);
            if( fileType == "csv" )
            {
                std::fstream fs;
                fs.open( filePath.c_str(), std::ios::in );

                if(fs)
                {
                    // parse number of eigenspace weights...
                    std::size_t parseEnd = fileName.find("_e");
                    std::size_t parseStart = fileName.find("___");
                    std::string esSubString = fileName.substr(parseStart+3,parseEnd-3);
                    int eigenspaceSize = std::atoi(esSubString.c_str());

                    EigenspaceWeights tmpESW;
                    std::string tmpImgName, tmpRawWtStr;
                    std::vector<double> tmpWts;
                    char c;
                    double tmpWtVal = 0;
                    while( fs.get(c) ) {
                        if( fs.fail() )
                        {
                            std::cout << "failure.\n";
                        }
                        tmpImgName = "";
                        tmpWts.clear();
                        while ( c != ',' ) {
                            tmpImgName += c;
                            fs.get(c);
                        }
                        for (int i = 0; i < eigenspaceSize; ++i) {
                            tmpRawWtStr = "";
                            fs.get(c);
                            while (c != ',') {
                                tmpRawWtStr += c;
                                fs.get(c);
                            }
                            tmpWtVal = std::atof(tmpRawWtStr.c_str());
                            tmpWts.push_back(tmpWtVal);
                        }
                        tmpESW.SetImageName(tmpImgName);
                        tmpESW.SetWeights(tmpWts);
                        tmpDB.push_back(tmpESW);
                        fs.get(c);
                    }
                } else {
                    std::cout << "Error: Could not open " << fileName << ". Exiting...\n";
                    exit(0);
                }
            }
        }
    }
    closedir(d);
    return tmpDB;
}

void ReadPGMToObject( const char * fname, GrayscaleImage & obj )
{
    std::ifstream infile;
    std::string tmpName = fname;
    char header[100];
    char *headerPtr = header;
    long int width, height;
    unsigned char *imageBuffer = NULL;
    int * imageVec = NULL;

    infile.open( fname, std::ios::in | std::ios::binary );

    //Check file opened properly
    if(!infile)
    {
        std::cout << "Cannot open file " << fname << "." << std::endl;
        return;
    }

    //Check that type is P5 = PGM
    infile.getline(header,100,'\n');
    if ( (header[0]!=80) ||    /* 'P' */
         (header[1]!=53) )     /* '5' */
    {
        std::cout << "Image " << fname << " is not PGM" << std::endl;
        return;
    }

    //Reposition pointer to capture dimensions (fa_H pgms have no comments)
    ++headerPtr;
    ++headerPtr;

    //Get width and height
    width = strtol(headerPtr, &headerPtr,0);
    height = atoi(headerPtr);

    //Read in rest of data as char array
    imageBuffer = (unsigned char *) new unsigned char [width*height];

    infile.read( reinterpret_cast<char *>(imageBuffer),
                 (width*height)*sizeof(unsigned char) );

    if(infile.fail())
    {
        std::cout << "Image " << fname << " has wrong size" << std::endl;
        return;
    }

    //close file
    infile.close();

    //convert raw data
    imageVec = new int[width*height];
    for(int i = 0; i < (int)width*height; ++i)
    {
        imageVec[i] = (int)imageBuffer[i];
    }

    //store to object
    obj.SetWidth(width);
    obj.SetHeight(height);
    obj.SetRawValues((const int *&)imageVec);
    obj.SetName(tmpName);

    //clean up
    delete [] imageBuffer;
    imageBuffer = NULL;

    return;
}

void WriteObjectToPGM( const char * fname, const GrayscaleImage & src )
{
    int length = src.GetWidth()*src.GetHeight();
    unsigned char * charImage;
    std::ofstream outfile;

    charImage = (unsigned char *) new unsigned char [length];

    int i = 0;
    for(i=0; i<length; ++i)
    {
        charImage[i] = (unsigned char)src.GetRawValueAtIndex(i);
    }

    outfile.open(fname, std::ios::out | std::ios::binary);

    if (!outfile) {
        std::cout << "Can't open file: " << fname << std::endl;
        exit(1);
    }

    //write header info
    outfile << "P5 " << src.GetWidth() << " " << src.GetHeight()
            << " 255" << std::endl;

    //write byte data
    outfile.write( reinterpret_cast<char *>(charImage),
                   (src.GetWidth()*src.GetHeight())*sizeof(unsigned char) );

    if (outfile.fail()) {
        std::cout << "Can't write image " << fname << std::endl;
        exit(0);
    }

    outfile.close();

    delete [] charImage;
    charImage = NULL;
}

void WriteDBtoPGMs( const std::vector<GrayscaleImage>& src, const int bound )
{
    if( bound <= (int)src.size() )
    {
        for( int i = 0; i < bound; ++i )
        {
            WriteObjectToPGM( src[i].GetName().c_str(), src[i] );
        }
    } else
    {
        std::cout << "Error: illegal bound passed as argument. Write failed.\n";
    }
}

void WriteDBtoPGMs( const std::vector<EigenImage>& src, const int bound )
{
    std::string tmpName;
    if( bound <= (int)src.size() )
    {
        for( int i = 0; i < bound; ++i )
        {
            std::string tmpName = "eigenspace/";
            tmpName += src[i].GetName();
            WriteObjectToPGM( tmpName.c_str(), src[i] );
        }
    } else
    {
        std::cout << "Error: illegal bound passed as argument. Write failed.\n";
    }
}

void WriteWeightsForImg( const char * imgName, arma::vec& weightVec )
{
    std::string tmpFileName = "eigenspace/";
    tmpFileName = tmpFileName + "___" +  std::to_string(weightVec.n_rows) + "_eigvecspace_wts.csv";

    std::fstream fs;
    double wt = 0;
    fs.open(tmpFileName, std::ios::out | std::ios::app );
    if(fs)
    {
        fs << imgName << ",";
        for( int i = 0; i < (int)weightVec.n_rows; ++i )
        {
            wt = weightVec.at(i);
            fs << wt << ",";
        }
        fs << "\n";
    }
    fs.close();
}

void ComputeAvgImage( const std::vector<GrayscaleImage> & src,
                      GrayscaleImage & dst )
{
    int length = src[0].GetWidth()*src[0].GetHeight();
    int i = 0, tmpRawVal = 0;
    double * tmpAvgVec = new double[length];

    for ( i = 0; i < length; ++i )
    {
        tmpAvgVec[i] = 0;
    }

    // Sum raw numbers..
    for( i = 0; i < (int)src.size(); ++i )
    {
        for( int j = 0; j < length; ++j )
        {
            tmpRawVal = src[i].GetRawValueAtIndex(j);
            tmpAvgVec[j] += tmpRawVal;
        }
    }

    // Divide by total...
    for( i = 0; i < length; ++i )
    {
        tmpAvgVec[i] /= (double)src.size();
    }

    // Store scaled values to object...
    dst.SetRawValues((const double *&) tmpAvgVec);

}