#include "Utils.h"

bool random_bool()
{
   static const int shift = static_cast<int>(std::log2(RAND_MAX));
   return (arc4random() >> shift) & 1;
}

vector<bool> random_bool_list(int size)
{
    vector<bool> l;

    for(int i = 0; i < size; i++)
        l.push_back(random_bool());

    return l;
}

vector<bool> zero_bool_list(int size)
{
    vector<bool> l;
    bool a = false;

    for(int i = 0; i < size; i++)
    {
        l.push_back(a);
        a = !a;
    }

    return l;
}

vector<int> random_integer_list(int max, int size)
{
    vector<int> vec(max, 0);
    iota(begin(vec), end(vec), 0);
    mt19937 rng(std::random_device{}());
    shuffle(begin(vec), end(vec), rng);

    vec.resize(size);

    return vec;
}

int reverse_bits(int number, int m)
{
     int ans = 0;

     for(int i = m-1; i >= 0; i--)
     {
        ans |= (number & 1) <<i;
        number>>=1;
     }

     return ans;
  }

vector<bool> reverse_bits(vector<bool> &l, int m)
{
    vector<bool> r;

    for(int i = 0; i < l.size(); i++)
        r.push_back(l[reverse_bits(i, m)]);

    return r;
}

int nb_differences(const vector<bool> &a, const vector<bool> &b)
{
    if(a.size() != b.size())
    {
        cout<<"Elements not of same size: "<<a.size()<<" and "<<b.size()<<endl;
        exit(EXIT_FAILURE);
    }

    int d = 0;
    for(int i = 0; i < a.size(); i++)
    {
        if(a[i] != b[i])
            d++;
    }

    return d;
}

vector<int> pos_difference(const vector<bool> &a, const vector<bool> &b)
{
    vector<int> pos_dif;
    if(a.size() != b.size())
    {
        cout<<"Elements not of same size: "<<a.size()<<" and "<<b.size()<<endl;
        exit(EXIT_FAILURE);
    }

    for(int i = 0; i < a.size(); i++)
    {
        if(a[i] != b[i])
            pos_dif.push_back(i);
    }

    return pos_dif;
}

vector<vector<double>> print_stats(vector<vector<double>> lik)
{
    vector<vector<double>> lik_std(lik[0].size(), vector<double>(3, 0));

    for(int i = 0; i < lik[0].size(); i++)
    {
        vector<double> l;

        for(int j = 0; j < lik.size(); j++)
        {
            l.push_back(lik[j][i]);
            lik_std[i][0] += lik[j][i]/lik.size();
        }

        sort(l.begin(), l.end());

        vector<double>::iterator quant1 =l.begin();
        advance(quant1, ((float)l.size()) *0.01);
        lik_std[i][1] = lik_std[i][0]-*quant1;
        vector<double>::iterator quant2 =l.begin();
        advance(quant2, ((float)l.size()) *0.99);
        lik_std[i][2] = *quant2-lik_std[i][0];
    }

    print_list(lik_std);

    return lik_std;
}

void print_comp(vector<vector<double>> lik_std, vector<vector<double>> lik_par, vector<int> pos_error)
{
    double avg_dist = 0;
    double avg_sep = 0;

    for(int i = 0; i < lik_std.size();i++)
    {
        for(int j = 0; j < lik_std[0].size(); j++)
        {
            if(abs(lik_std[i][j] - lik_par[i][j]) > 0.5)
            {
                avg_dist += double(j-pos_error[i]);
                break;
            }
        }

        for(int j = 0; j < lik_std[0].size(); j++)
            avg_sep += (lik_par[i][j]/lik_std[i][j]-1)/double(lik_std[0].size());
    }

    cout<<"Avg distance (steps): "<<avg_dist/lik_std.size()<<endl;
    cout<<"Avg separation (%): "<<avg_sep/lik_std.size()*100<<endl;

}

void print_list(vector<bool> l)
{
    for(int i = 0; i < l.size(); i++)
        cout<<l[i]<<" ";
    cout<<endl;
}

void print_list(vector<double> l)
{
    cout<<"[";

    for(int i = 0; i < l.size(); i++)
    {
        cout<<l[i];

        if(i != l.size()-1)
            cout<<", ";
    }

    cout<<"]";

    cout<<endl;
}

void print_list(vector<int> l)
{
    cout<<"[";

    for(int i = 0; i < l.size(); i++)
    {
        cout<<l[i];

        if(i != l.size()-1)
            cout<<", ";
    }

    cout<<"]";

    cout<<endl;
}

void print_list(vector<unsigned long> l)
{
    cout<<"[";

    for(int i = 0; i < l.size(); i++)
    {
        cout<<l[i];

        if(i != l.size()-1)
            cout<<", ";
    }

    cout<<"]";

    cout<<endl;
}


void print_list(vector<vector<double>> l)
{
    if(l.size() > 1)
        cout<<"[";

    for(int i = 0; i < l.size(); i++)
    {
        cout<<"[";

        for(int j = 0; j < l[i].size(); j++)
        {
            cout<<l[i][j];
            if(j != l[i].size()-1)
                cout<<", ";
        }

        cout<<"]";

        if(i != l.size()-1)
            cout<<", ";
    }

    if(l.size() > 1)
        cout<<"] ";
    cout<<endl;
}

void print_list(vector<vector<int>> l)
{
    if(l.size() > 1)
        cout<<"[";

    for(int i = 0; i < l.size(); i++)
    {
        cout<<"[";

        for(int j = 0; j < l[i].size(); j++)
        {
            cout<<l[i][j];
            if(j != l[i].size()-1)
                cout<<", ";
        }

        cout<<"]";

        if(i != l.size()-1)
            cout<<", ";
    }

    if(l.size() > 1)
        cout<<"] ";
    cout<<endl;
}
