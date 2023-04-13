#include "AbstractPolarCode.h"

AbstractPolarCode::AbstractPolarCode(int m, int k, vector<int> frozen_bits, bool debug)
{
    m_n = pow(2, m);
    m_m = m;
    m_k = k;
    m_frozen_bits = frozen_bits;
    m_debug = debug;
}


vector<bool> AbstractPolarCode::encode(const vector<bool> &message)
{
    if(message.size() != m_k)
    {
        cout<<"Message's size was expected to be "<<m_k<<" but received "<<message.size()<<endl;
        exit(EXIT_FAILURE);
    }

    if(m_frozen_bits.size() < m_n - m_k)
    {
        cout<<"Frozen bits not defined"<<endl;
        exit(EXIT_FAILURE);
    }

    vector<bool> message_e;
    int n_half = m_n;
    int j = 0;

    for(int i = 0; i < m_n; i++)
    {
        bool is_frozen = find(m_frozen_bits.begin(), m_frozen_bits.end(), i) != m_frozen_bits.end();

        if(!is_frozen)
        {
            message_e.push_back(message[j]);
            j++;
        }

        else
            message_e.push_back(0);
    }

    message_e = reverse_bits(message_e, m_m);

    for(int j = 0; j < m_m; j++)
    {
        n_half /= 2;

        for(int k = 0; k < m_n; k+= n_half*2)
        {
            for(int i = 0; i < n_half; i++)
            {
                int u = i+k;
                message_e[u] = message_e[u] ^ message_e[u+n_half];
            }
        }
    }

    return message_e;
}
