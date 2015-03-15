#ifndef OPTIMFFT_H
#define OPTIMFFT_H

#include "../Caramel/AudioFFT.h"
#include "ComplexBuffer.h"

class OptimFFT
{
public:
    OptimFFT(const int size);
    ~OptimFFT();

    void     FFT(float* destRe, float* destIm, const float* segment);
    void  OddFFT(float* destRe, float* destIm, const float* segment);
    void EvenFFT(float* destRe, float* destIm, const float* segment);

private:
    const int           _bufSize;
    const int           _computeSize;
    const int           _complexSize;

    float*              _buffer1;
    float*              _buffer2;

    audiofft::AudioFFT* _fft;
    ComplexBuffer*      _result1;
    ComplexBuffer*      _result2;
};

#endif//OPTIMFFT_H
