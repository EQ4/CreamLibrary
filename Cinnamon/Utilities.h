#ifndef UTILITIES_H
#define UTILITIES_H

#include <cmath>
#define PI 3.14159265359

inline void CopyAndPad (float* dest, const float* src, const int size)
{
    ::memcpy(dest     , src, size*sizeof(float));
    ::memset(dest+size,   0, size*sizeof(float));
}



inline void Multpiply(float* dest, const float* re1, const float* im1, const float* re2, const float* im2, const int size)
{
    dest[0]    = re1[0]*re2[0];
    dest[size] = re1[size]*re2[size];
    for (int i = 1 ; i < size ; i ++)
    {
        dest[i]      = re1[i]*re2[i] - im1[i]*im2[i];
        dest[size+i] = re1[i]*im2[i] + im1[i]*re2[i];
    }
}

#endif//UTILITIES_H
