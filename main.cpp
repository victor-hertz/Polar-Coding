#include "BSChannel.h"
#include "AWGNChannel.h"
#include "PolarCode.h"

int main()
{
    // Parameters
    int m = 7;
    int n = pow(2,m);
    double rate = 0.5;
    int k = round(n*rate);

    int trials_frozen = 1e4;
    int trials_decoding = 1e4;

    // Display
    cout<<"Code parameters - n: "<<n<<"  k: "<<k<<"  Rate: "<<float(k)/float(n)<<endl;
    cout<<endl;

    // Execution (example)
    AWGNChannel channel;
    PolarCode code(m, k, {}, 0);

    vector<double> noise = {0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5};
    vector<vector<int>> frozen_bits;
    vector<vector<int>> active_bits;
    double progress = 0;

    for(int i = 0; i < noise.size(); i++)
    {
        channel.set_noise(noise[i]);
        frozen_bits.push_back(code.estimate_frozen_bits(channel, trials_frozen, false, false));
        active_bits.push_back(code.get_active_bits());
        cout<<"Noise: "<<noise[i]<<endl;

        double curr_progress = float(i+1)/float(noise.size());
        if(curr_progress >= progress+0.1)
        {
            progress += 0.1;
            cout<<progress*100<<"%"<<endl;
        }
    }

    cout<<endl<<"Frozen bits:"<<endl;
    print_list(frozen_bits);

    cout<<endl<<"Active bits:"<<endl;
    print_list(active_bits);

    return 0;
}


