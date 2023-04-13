#include "AWGNChannel.h"

AWGNChannel::AWGNChannel(double snr): AbstractChannel(snr)
{

}

vector<double> AWGNChannel::transmit(const vector<bool> &message) const
{
    default_random_engine generator(random_device{}());
    normal_distribution<double> distribution{0,sqrt(get_variance())};

    vector<double> t;

    for(int i = 0; i < message.size(); i++)
    {
        double signal = distribution(generator);
        t.push_back(1.0f-2.0f*double(message[i])+signal);
    }

    return t;
}

vector<vector<double>> AWGNChannel::estimate(const vector<double> &transmission) const
{
    double variance = get_variance();
    double p = 1/sqrt(2*M_PI*variance);

    vector<vector<double>> estimate;

    for(int i = 0; i < transmission.size(); i++)
        estimate.push_back({exp(-pow(abs(transmission[i]-1),2)/(2*variance)), exp(-pow(abs(transmission[i]+1),2)/(2*variance)), log(p)});

    return estimate;
}

double AWGNChannel::get_variance() const
{
    return 1/exp(m_noise*log(10)/10.f);
}

string AWGNChannel::get_name() const
{
    return "AWGN";
}

