#ifndef EIGENIMAGE_H
#define EIGENIMAGE_H

#include <armadillo>
#include "GrayscaleImage.h"

class EigenImage: public GrayscaleImage
{
public:
    EigenImage(): GrayscaleImage(), eigenValue(0) {};
    EigenImage(GrayscaleImage& obj): GrayscaleImage(obj), eigenValue(0) {};
    ~EigenImage() {};

    double GetEigVal() const { return eigenValue; };
    arma::vec GetEigVec() const { return eigenVector; };
    GrayscaleImage ExtractImage() const;

    void SetEigVal(const double aVal) { eigenValue = aVal; };
    void SetEigVec(const arma::vec& aVec) { eigenVector = aVec; };

    void NormalizeValuesByEigVec();

private:
    double eigenValue;
    arma::vec eigenVector;
};


#endif // EIGENIMAGE_H
