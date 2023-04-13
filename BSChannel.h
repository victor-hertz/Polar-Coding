#ifndef BSCHANNEL_H
#define BSCHANNEL_H

#include "AbstractChannel.h"

class BSChannel : public AbstractChannel
{
public:
    BSChannel(double p_fail = 0);

    vector<double> transmit(const vector<bool> &message) const;
    vector<vector<double>> estimate(const vector<double> &transmission) const;

    string get_name() const;

private:
    bool transmit_bool(bool input) const;
};

#endif // BSCHANNEL_H
