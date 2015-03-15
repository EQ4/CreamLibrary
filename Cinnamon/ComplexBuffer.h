#ifndef COMPLEX_BUFFER_H
#define COMPLEX_BUFFER_H

#include <cstring>

class ComplexBuffer
{
private:
    float* _re  ;
    float* _im  ;
    int    _size;

public:
    ComplexBuffer(int initSize)
    {
        if (initSize > 0)
        {
            _size = initSize;
            _re   = new float[_size];
            _im   = new float[_size];
        }
        setZero();
    }

    ~ComplexBuffer()
    {
        clear();
    }

    void clear()
    {
        setZero();
        if (_re)
        {
            delete[] _re;
        }
        if (_im)
        {
            delete[] _im;
        }
        _size   = 0;
    }

    void setZero()
    {
        ::memset(_re, 0, _size*sizeof(float));
        ::memset(_im, 0, _size*sizeof(float));
    }

    float* re()
    {
        return _re;
    }

    float* im()
    {
        return _im;
    }
};

#endif//COMPLEX_BUFFER_H
