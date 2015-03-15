#include "OptimFFT.h"
#include "Utilities.h"
#include <cmath>

#define PI 3.14159265359

OptimFFT::OptimFFT(const int size):
    _bufSize                                    ( size       ),
    _computeSize                                (_bufSize<<1 ),
    _complexSize(audiofft::AudioFFT::ComplexSize(_computeSize))
{
    _buffer1 = new float[_computeSize];
    _buffer2 = new float[_computeSize];

    _fft     = new audiofft::AudioFFT();
    _result1 = new ComplexBuffer(_complexSize);
    _result2 = new ComplexBuffer(_complexSize);
    _fft->init(_computeSize);
}

OptimFFT::~OptimFFT()
{
    if (_buffer1)
    {
        delete[] _buffer1;
    }
    if (_buffer2)
    {
        delete[] _buffer2;
    }

    if (_fft)
    {
        _fft->init(0);
        delete _fft;
    }
    if (_result1)
    {
        delete _result1;
    }
    if (_result2)
    {
        delete _result2;
    }
}

void OptimFFT::FFT(float* destRe, float* destIm, const float* segment)
{
    CopyAndPad(_buffer1,  segment      , _bufSize                  );
    _fft->fft (_buffer1, _result1->re(), _result1->im()            );
    ::memcpy  (destRe  , _result1->re(), _complexSize*sizeof(float));
    ::memcpy  (destIm  , _result1->im(), _complexSize*sizeof(float));
}

inline void Radix2(float* reald, float* imagd, const float* real1, const float* real2, const float* imag1, const float* imag2, const int size)
{
    double tr[2], ti[2];
    for (int i = 0 ; i < size ; i += 2)
    {
        double wr = cos((float)i*2.0*PI/(double)size);
        double wi = sin((float)i*2.0*PI/(double)size);
        
        tr[0] = (real1[i] + real2[i]);
        ti[0] = (imag1[i] + imag2[i]);
        tr[1] = (real1[i+1] - real2[i+1]) * wr - ((imag1[i+1] - imag2[i+1]) * wi);
        ti[1] = (imag1[i+1] - imag2[i+1]) * wr + ((real1[i+1] - real2[i+1]) * wi);
        
        reald[i*2] = tr[0];
        reald[i*2+2] = tr[1];
        imagd[i*2] = ti[0];
        imagd[i*2+2] = ti[1];
    }

    //dest[size<<1]   = src1[size] - src2[size];
}

inline void EvenPoints(float* dest, const float* src1, const float* src2, const int size)
{
    for (int i = 0 ; i < size ; i += 2)
    {
        dest[i*2  ] = src1[i   ] + src2[i  ];
        dest[i*2+2] = src1[i+1 ] - src2[i+1];
    }
    dest[size<<1]   = src1[size] + src2[size];
}

inline void OddPoints2(float* dest, const float* src1, const float* src2, const int size)
{
    for (int i = 0 ; i < size ; i += 2)
    {
        double wr = cos((float)i*2.0*PI/(double)size);
        double wi = sin((float)i*2.0*PI/(double)size);
        dest[i*2  ] = src1[i   ] + src2[i  ] * wr;
        dest[i*2+2] = src1[i+1 ] - src2[i+1] * wi;
    }
    dest[size<<1]   = src1[size] + src2[size];
}

inline void OddPoints (float* rea, float* ima, const float* real1, const float* real2, const float* imag1, const float* imag2, const int size)
{    
    for (int i = 0 ; i < size ; i += 2)
    {
        double wr = cos((float)i*2.0*PI/(double)size);
        double wi = sin((float)i*2.0*PI/(double)size);
        double t1 = wr * real2[i] - wi * imag2[i];
        double t2 = wi * real2[i] + wr * imag2[i];
        rea[i*2] = real1[i] - t1;//(real1[i] - real2[i]) * wr - ((imag1[i] - imag2[i]) * wi);
        ima[i*2] = imag1[i] - t2;
        rea[i*2+2] = real1[i+1] + t1;//(imag1[i+1] - imag2[i+1]) * wr + ((real1[i+1] - real2[i+1]) * wi);
        ima[i*2+2] = imag1[i+1] + t2;
        
    }
    //dest[size<<1]   = src1[size] - src2[size];
}

void OptimFFT::EvenFFT(float* destRe, float* destIm, const float* segment)
{
    CopyAndPad(_buffer1,  segment         , _bufSize                );
    CopyAndPad(_buffer2,  segment+_bufSize, _bufSize                );
    
    _fft->fft (_buffer1, _result1->re()   , _result1->im()          );
    _fft->fft (_buffer2, _result2->re()   , _result2->im()          );
    
    EvenPoints( destRe , _result1->re()   , _result2->re(), _bufSize);
    EvenPoints( destIm , _result1->im()   , _result2->im(), _bufSize);
}

void OptimFFT::OddFFT(float* destRe, float* destIm, const float* segment)
{
    CopyAndPad(_buffer1,  segment         , _bufSize                );
    CopyAndPad(_buffer2,  segment+_bufSize, _bufSize                );

    _fft->fft (_buffer1, _result1->re()   , _result1->im()          );
    _fft->fft (_buffer2, _result2->re()   , _result2->im()          );
    /*
    OddPoints2( destRe , _result1->re()   , _result2->re(), _bufSize);
    OddPoints2( destIm , _result1->im()   , _result2->im(), _bufSize);
    */
    OddPoints(destRe, destIm, _result1->re(), _result2->re(), _result1->im(), _result2->im(), _bufSize);
    //OddPoints(destIm, _result1->im(), _result2->im(),  _result1->re(), _result2->re(), _bufSize);
}

