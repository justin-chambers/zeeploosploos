#include "EigenImage.h"

GrayscaleImage EigenImage::ExtractImage() const
{
    GrayscaleImage tmp;
    tmp = *this;
    return tmp;
}

void EigenImage::NormalizeValuesByEigVec()
{
    std::vector<double> tmpValVec;
    for( int i = 0; i < this->eigenVector.n_rows; ++i )
    {
        tmpValVec.push_back(eigenVector.at(i));
    }
    this->SetScaledValues(tmpValVec); // also updates raw values and range...
    this->ScaleRawValues(); // scales values between 0 and 1
}