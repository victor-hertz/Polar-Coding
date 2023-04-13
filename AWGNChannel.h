#ifndef AWGNCHANNEL_H
#define AWGNCHANNEL_H

#include "AbstractChannel.h"

class AWGNChannel : public AbstractChannel
{
public:
    AWGNChannel(double snr = 0);

    vector<double> transmit(const vector<bool> &message) const;
    vector<vector<double>> estimate(const vector<double> &transmission) const;

    string get_name() const;

private:
    double get_variance() const;
};

#endif // AWGNCHANNEL_H
