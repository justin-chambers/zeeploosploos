#include "GrayscaleImage.h"

GrayscaleImage::GrayscaleImage():width(0) , height(0), maxRaw(0), maxScale(0),
                                 minRaw(0), minScale(0), rawValues(NULL),
                                 scaledValues(NULL){ }

GrayscaleImage::GrayscaleImage( const int aWidth,
                                const int aHeight,
                                const int *& aValueArray,
                                const std::string aName )
        :width(aWidth),height(aHeight),rawValues(NULL), scaledValues(NULL)
{
    if( (width*height) == GetArrayLength(aValueArray) )
    {
        rawValues = new int[width*height];
        CopyValues(aValueArray, rawValues);
    }
    ScaleRawValues();
    fileName = aName;
}

GrayscaleImage::GrayscaleImage(const GrayscaleImage &obj)
{
    *this = obj;
}

GrayscaleImage::~GrayscaleImage()
{
    if( rawValues != NULL )
        delete [] rawValues;
    rawValues = NULL;

    if( scaledValues != NULL )
    {
        delete [] scaledValues;
    }
    scaledValues = NULL;
}

GrayscaleImage& GrayscaleImage::operator=(const GrayscaleImage &rightSide)
{
    if(this == &rightSide)
        return *this;
    this->width = rightSide.GetWidth();
    this->height = rightSide.GetHeight();
    this->minScale = rightSide.GetMinScaleValue();
    this->maxScale = rightSide.GetMaxScaleValue();
    this->minRaw = rightSide.GetMinRawValue();
    this->maxRaw = rightSide.GetMaxRawValue();
    this->rawValues = rightSide.GetRawValues();
    this->scaledValues = rightSide.GetScaledValues();
    this->fileName = rightSide.GetName();
    return *this;
}

int GrayscaleImage::GetArrayLength( const int *& src )
{
    int length = 0;
    while(*src < 256)
    {
        ++length;
        ++src;
    }
    return length;
}

int * GrayscaleImage::GetRawValues() const
{
    int * copyValues = new int[width*height];
    for(int i = 0; i < width*height; ++i)
    {
        copyValues[i] = rawValues[i];
    }
    return copyValues;
}

std::vector<double> GrayscaleImage::GetRawValVec() const
{
    std::vector<double> tmpVec;
    for(int i = 0; i < height*width; ++i)
    {
        tmpVec.push_back((double)rawValues[i]);
    }
    return tmpVec;
}

double * GrayscaleImage::GetScaledValues() const {
    if(scaledValues != NULL)
    {
        double * copyValues = new double[width*height];
        for(int i = 0; i < width*height; ++i)
        {
            copyValues[i] = scaledValues[i];
        }
        return copyValues;
    } else {
        return NULL;
    }
}

std::vector<double> GrayscaleImage::GetScaledValVec() const
{
    std::vector<double> tmpVec;
    for(int i = 0; i < height*width; ++i)
    {
        tmpVec.push_back(scaledValues[i]);
    }
    return tmpVec;
}

void GrayscaleImage::SetRawValues(const int *& aValueArray)
{
    CopyValues(aValueArray, rawValues);
    ScaleRawValues();
}

void GrayscaleImage::SetRawValues(const double *& aValueArray)
{
    if(rawValues != NULL)
    {
        delete [] rawValues;
    }

    rawValues = new int[height*width];

    for( int i = 0; i < height*width; ++i )
    {
        rawValues[i] = (int)(round(aValueArray[i]));
    }

    ComputeRawRange();
    ScaleRawValues();
}

void GrayscaleImage::SetScaledValues(const double *& aValueArray)
{
    CopyValues(aValueArray, scaledValues);

}

void GrayscaleImage::SetScaledValues(const std::vector<double>& aValVec)
{
    if( height*width == (int)aValVec.size() )
    {
        for( int i = 0; i < height*width; ++i )
        {
            scaledValues[i] = aValVec[i];
        }
        ComputeScaledRange();
        ConvertScaledtoRaw();
    }
}

void GrayscaleImage::CopyValues( const int *& src,
                                 int *& dst )
{
    int length = width*height;

    if( dst != NULL ) {
        delete[] dst;
    }
    dst = new int[length];

    for(int i = 0; i < length; ++i)
    {
        dst[i] = src[i];
    }

    ComputeRawRange();
}

void GrayscaleImage::CopyValues( const double *& src,
                                 double *& dst )
{
    int length = width*height;

    if(dst != NULL)
        delete [] dst;
    dst = new double[length];

    for(int i = 0; i < length; ++i)
    {
        dst[i] = src[i];
    }

    ComputeScaledRange();
    ConvertScaledtoRaw();
}

void GrayscaleImage::ScaleRawValues()
{
    int divisor = maxRaw - minRaw;
    if(scaledValues != NULL)
    {
        delete [] scaledValues;
    }
    scaledValues = new double[height*width];
    for( int i = 0; i < height*width; ++i )
    {
        scaledValues[i] = 0;
        scaledValues[i] = (double)(rawValues[i] - minRaw);
        scaledValues[i] /= divisor;
    }
    ComputeScaledRange();
}

void GrayscaleImage::ConvertScaledtoRaw()
{
    if( rawValues != NULL )
    {
        delete [] rawValues;
    }

    rawValues = new int[height*width];

    for( int i = 0; i < height*width; ++i )
    {
        rawValues[i] = 0;
        rawValues[i] = (int)( (scaledValues[i] - minScale) * 255 / ( maxScale - minScale ) );
    }
    ComputeRawRange();
}

GrayscaleImage GrayscaleImage::AddFaces(const GrayscaleImage &f1, const GrayscaleImage f2)
{
    GrayscaleImage addF;
    std::string tmpName = "__comp_";

    tmpName = tmpName + f1.GetName() + f2.GetName();
    addF.SetName(tmpName);

    if( ( f1.GetWidth() == f2.GetWidth() ) && ( f1.GetHeight() == f2.GetHeight() ) )
    {
        int length = f1.GetHeight()*f1.GetWidth();
        double * tmpSVal = new double[length];

        addF.SetWidth(f1.GetWidth());
        addF.SetHeight(f1.GetHeight());

        for( int i = 0; i < length; ++i )
        {
            tmpSVal[i] = 0;
            tmpSVal[i] = f1.GetScaledValueAtIndex(i) + f2.GetScaledValueAtIndex(i);
        }

        addF.SetScaledValues( (const double *&) tmpSVal );
        addF.ConvertScaledtoRaw();

    }

    return addF;
}

GrayscaleImage GrayscaleImage::SubtractFaces(const GrayscaleImage &f1, const GrayscaleImage f2)
{
    GrayscaleImage addF;
    std::string tmpName = "__comp_";

    tmpName = tmpName + f1.GetName() + f2.GetName();
    addF.SetName(tmpName);

    if( ( f1.GetWidth() == f2.GetWidth() ) && ( f1.GetHeight() == f2.GetHeight() ) )
    {
        int length = f1.GetHeight()*f1.GetWidth();
        double * tmpSVal = new double[length];

        addF.SetWidth(f1.GetWidth());
        addF.SetHeight(f1.GetHeight());

        for( int i = 0; i < length; ++i )
        {
            tmpSVal[i] = 0;
            tmpSVal[i] = f1.GetScaledValueAtIndex(i) - f2.GetScaledValueAtIndex(i);
        }

        addF.SetScaledValues( (const double *&) tmpSVal );
        addF.ConvertScaledtoRaw();

    }

    return addF;
}

void GrayscaleImage::ComputeRawRange()
{
    maxRaw = minRaw = rawValues[0];
    for(int i = 1; i < height*width; ++i)
    {
        if(rawValues[i] > maxRaw)
        {
            maxRaw = rawValues[i];
        } else if(rawValues[i] < minRaw)
        {
            minRaw = rawValues[i];
        }
    }
}

void GrayscaleImage::ComputeScaledRange()
{
    maxScale = minScale = scaledValues[0];
    for(int i = 1; i < height*width; ++i)
    {
        if(scaledValues[i] > maxScale)
        {
            maxScale = scaledValues[i];
        } else if(scaledValues[i] < minScale)
        {
            minScale = scaledValues[i];
        }
    }
}