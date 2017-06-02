#ifndef PCAAPPLICATION_H
#define PCAAPPLICATION_H

//Libraries
#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <cstdlib>
#include <dirent.h>
#include <armadillo>

//Classes
#include "GrayscaleImage.h"
#include "EigenImage.h"
#include "EigenspaceWeights.h"

class PCAApplication {
public:
    PCAApplication() { DisplayUsage(); };
    PCAApplication(char** argv, const int argc, bool isTrainingInstance);
    ~PCAApplication() {};

private:
    // Variables
    std::vector<GrayscaleImage> imageDB;
    std::vector<EigenImage> eigenFDB;
    std::vector<GrayscaleImage> phiDB;
    std::vector<EigenspaceWeights> omegaDB;
    std::string matchFileName;
    std::vector<std::string> imagesMatched;

    // Operations
    void TrainingSequence();
    void TestingSequence(const int gallerySize);

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
    void WriteTestResult( const std::string& anImgName );
    void WritePerformanceToFile( const std::string& aFileName,
                                                 const int aGallerySize,
                                                 const int numMatches,
                                                 const int totalNumTestImages );
    void ComputeAvgImage( const std::vector<GrayscaleImage> & src,
                          GrayscaleImage & dst );
    int GetImageID( const std::string& src );
    void UpdateExistingPotentialMatches( const double aDifs,
                                         const int anID,
                                         std::map<int, double>& srcMap );
    bool IsMatch( const int anID, std::map<double, int> srcMap );
};

#endif