#include "AbstractChannel.h"

AbstractChannel::AbstractChannel(double noise)
{
    m_noise = noise;
}

void AbstractChannel::set_noise(double noise)
{
    m_noise = noise;
}

double AbstractChannel::get_noise() const
{
    return m_noise;
}
