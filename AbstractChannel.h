#ifndef ABSTRACTCHANNEL_H
#define ABSTRACTCHANNEL_H

#include "Utils.h"

class AbstractChannel
{
public:
    AbstractChannel(double noise);

    virtual vector<double> transmit(const vector<bool> &message) const = 0;
    virtual vector<vector<double>> estimate(const vector<double> &transmission) const = 0;

    virtual string get_name() const = 0;

    void set_noise(double noise);
    double get_noise() const;

protected:
    double m_noise;
};

#endif // ABSTRACTCHANNEL_H
