#ifndef EIGENSPACEWEIGHTS_H
#define EIGENSPACEWEIGHTS_H

#include <string>
#include <vector>

class EigenspaceWeights {
public:
    EigenspaceWeights() {};
    EigenspaceWeights( const std::string& aName,
                       const std::vector<double>& aWeightVec ):
            imageName(aName), weights(aWeightVec) {};
    EigenspaceWeights( const EigenspaceWeights& obj )
    {
        *this = obj;
    }
    EigenspaceWeights& operator=( const EigenspaceWeights& rhs )
    {
        if( this == &rhs )
            return *this;
        imageName = rhs.GetImageName();
        weights = rhs.GetWeights();
        return *this;
    }

    std::string GetImageName() const { return imageName; };
    std::vector<double> GetWeights() const { return weights; };

    void SetImageName( const std::string& aName ) { imageName = aName; };
    void SetWeights( const std::vector<double>& aWeightVec ) { weights = aWeightVec; };

private:
    std::string imageName;
    std::vector<double> weights;
};


#endif //EIGENSPACEWEIGHTS_H
