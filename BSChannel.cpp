#include "BSChannel.h"

BSChannel::BSChannel(double p_fail): AbstractChannel(p_fail)
{
}

vector<double> BSChannel::transmit(const vector<bool> &message) const
{
    vector<double> t;

    for(int i = 0; i < message.size(); i++)
        t.push_back(double(transmit_bool(message[i])));

    return t;
}

vector<vector<double>> BSChannel::estimate(const vector<double> &transmission) const
{
    vector<vector<double>> estimate;
    for(int i = 0; i < transmission.size(); i++)
    {
        if(transmission[i] == 0)
            estimate.push_back(vector<double>({1-m_noise, m_noise, 0}));
        else if(transmission[i] == 1)
            estimate.push_back(vector<double>({m_noise, 1-m_noise, 0}));
        else
        {
            cout<<"Undefined symbol in transmission: "<<transmission[i]<<endl;
            exit(EXIT_FAILURE);
        }
    }

    return estimate;
}

bool BSChannel::transmit_bool(bool input) const
{
    float r = static_cast <float> (arc4random_uniform(RAND_MAX)) / static_cast <float> (RAND_MAX);

    if(r < m_noise)
        return !input;

    return input;
}

string BSChannel::get_name() const
{
    return "BSC";
}
