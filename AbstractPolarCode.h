#ifndef ABSTRACTPOLARCODE_H
#define ABSTRACTPOLARCODE_H

#include "AbstractChannel.h"

class AbstractPolarCode
{
public:
    AbstractPolarCode(int m, int k, vector<int> frozen_bits = vector<int>(), bool debug = false);

    vector<bool> encode(const vector<bool> &message);
    virtual vector<bool> decode(const vector<vector<double> > &encoder_estimates) = 0;


protected:
    int m_n;
    int m_m;
    int m_k;
    bool m_debug;
    vector<int> m_frozen_bits;
};

#endif // ABSTRACTPOLARCODE_H
