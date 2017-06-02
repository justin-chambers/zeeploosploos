#ifndef GRAYSCALEIMAGE_H
#define GRAYSCALEIMAGE_H

#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>

class GrayscaleImage {
public:
    //Constructors & Deconstructor
    GrayscaleImage();
    GrayscaleImage(const int aWidth,
                   const int aHeight,
                   const int *& aValueArray,
                   const std::string aName);
    GrayscaleImage(const GrayscaleImage& obj);
    ~GrayscaleImage();

    //Overloaded Operators
    GrayscaleImage& operator=(const GrayscaleImage& rightSide);

    // Class Operations
    static GrayscaleImage AddFaces(const GrayscaleImage& f1, const GrayscaleImage f2);
    static GrayscaleImage SubtractFaces(const GrayscaleImage& f1, const GrayscaleImage f2);

    //Accessors
    int GetWidth() const { return width; };
    int GetHeight() const { return height; };
    int GetMaxRawValue() const { return maxRaw; };
    int GetMinRawValue() const { return minRaw; };
    double GetMaxScaleValue() const { return maxScale; };
    double GetMinScaleValue() const { return minScale; };

    int * GetRawValues() const;
    std::vector<double> GetRawValVec() const;
    int GetRawValueAtIndex(const int anIndex) const { return rawValues[anIndex]; };

    double * GetScaledValues() const;
    std::vector<double> GetScaledValVec() const;
    double GetScaledValueAtIndex(const int anIndex) const { return scaledValues[anIndex]; };

    std::string GetName() const { return fileName; };

    //Mutators
    void SetWidth(const int aWidth) { width = aWidth; };
    void SetHeight(const int aHeight) { height = aHeight; };
    void SetRawMaxValue(int aVal) { maxRaw = aVal; };
    void SetRawMinValue(int aVal) { minRaw = aVal; };
    void SetScaledMaxValue(double aVal) { maxScale = aVal; };
    void SetScaledMinValue(double aVal) { minScale = aVal; };

    void SetRawValues(const int *& aValueArray);
    void SetRawValues(const double *& aValueArray);

    void SetScaledValues(const double *& aValueArray);
    void SetScaledValues(const std::vector<double>& aValVec);
    void SetScaledValueAtIndex(const int anIndex, const double aVal) { scaledValues[anIndex] = aVal; };

    void SetName(const std::string aName) { fileName = aName; };

    // PCA-Specific Operations
    void ScaleRawValues();
    void ConvertScaledtoRaw();

private:
    //Class shared members

    //Instance members
    int width;
    int height;
    int maxRaw;
    int minRaw;
    double maxScale;
    double minScale;
    int * rawValues;
    double * scaledValues;
    std::string fileName;

    //Utility function prototypes
    int GetArrayLength( const int *& src );
    void CopyValues( const int *& src, int *& dst );
    void CopyValues( const double *& src, double *& dst );
    void ComputeRawRange();
    void ComputeScaledRange();

};

#endif //GRAYSCALEIMAGE_H
