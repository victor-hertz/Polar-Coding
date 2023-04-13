#include "PolarCode.h"

PolarCode::PolarCode(int m, int k, vector<int> frozen_bits, bool debug): AbstractPolarCode(m, k, frozen_bits, debug)
{
    m_decoder_estimates_h = vector<vector<GraphFlow>>(m_n, vector<GraphFlow>(m_m+1, GraphFlow()));
    m_likelihoods=vector<double>(m_n);
}


vector<bool> PolarCode::decode(const vector<vector<double>> &encoder_estimates)
{
    return decode(encoder_estimates, true);
}

vector<bool> PolarCode::decode(const vector<vector<double> > &message_est, bool compute_l)
{
    m_iterations = 0;

    if(m_frozen_bits.size() < m_n - m_k)
    {
        cout<<"Frozen bits not defined"<<endl;
        exit(EXIT_FAILURE);
    }

    vector<bool> message_dec = vector<bool>(m_n, 0);

    initialize_data(message_est);

    for(int i = 0; i < m_n; i++)
    {
        int i_r = reverse_bits(i, m_m);

        if(find(m_frozen_bits.begin(), m_frozen_bits.end(), i) != m_frozen_bits.end())
        {
            if(compute_l)
            {
                vector<double> marginal = compute_block(0, i_r, DIRECTION::LEFT);
                m_likelihoods[i] = marginal[2]+log(marginal[0]);
            }

            m_decoder_estimates_h[i_r][0].f = ZERO_MESSAGE;
        }

        else
        {
            auto start = chrono::high_resolution_clock::now();

            vector<double> marginal = compute_block(0, i_r, DIRECTION::LEFT);

            auto stop = chrono::high_resolution_clock::now();

            m_time += chrono::duration_cast<chrono::microseconds>(stop-start).count();

            if(marginal[0] > marginal[1])
            {
                if(compute_l)
                    m_likelihoods[i] = marginal[2]+log(marginal[0]);

                m_decoder_estimates_h[i_r][0].f = ZERO_MESSAGE;
            }
            else
            {
                if(compute_l)
                    m_likelihoods[i] = marginal[2]+log(marginal[1]);

                message_dec[i] = 1;
                m_decoder_estimates_h[i_r][0].f = ONE_MESSAGE;
            }
        }

        update_block(0, i_r, DIRECTION::RIGHT);
    }

    vector<bool> message_dec_m;

    for(int i = 0; i < m_n; i++)
    {
        if(find(m_frozen_bits.begin(), m_frozen_bits.end(), i) == m_frozen_bits.end())
            message_dec_m.push_back(message_dec[i]);
    }

    return message_dec_m;
}

vector<bool> PolarCode::parity_decode(const vector<vector<double>> &message_est, const vector<bool> &message, vector<bool> random_bits,  vector<int> bit_flips, bool parity_check)
{
    if(random_bits.size() == 0)
        random_bits = random_bool_list(m_n);

    int froz_it = 0;
    bool random_frozen = false;
    m_iterations = 0;

    if(m_frozen_bits.size() < m_n - m_k)
    {
        cout<<"Frozen bits not defined"<<endl;
        exit(EXIT_FAILURE);
    }

    vector<bool> message_dec = vector<bool>(m_n, 0);

    initialize_data(message_est);

    for(int i = 0; i < m_n; i++)
    {
        int i_r = reverse_bits(i, m_m);

        if(find(m_frozen_bits.begin(), m_frozen_bits.end(), i) != m_frozen_bits.end())
        {
            froz_it++;
            vector<double> marginal = compute_block(0, i_r, DIRECTION::LEFT);

            if(!random_frozen || (random_frozen & random_bits[i]) || !parity_check)
            {
                m_likelihoods[i] = marginal[2]+log(marginal[0]);
                m_decoder_estimates_h[i_r][0].f = ZERO_MESSAGE;
            }
            else
            {
                m_likelihoods[i] = marginal[2]+log(marginal[1]);
                m_decoder_estimates_h[i_r][0].f = ONE_MESSAGE;
            }
        }

        else
        {
            auto start = chrono::high_resolution_clock::now();

            vector<double> marginal = compute_block(0, i_r, DIRECTION::LEFT);

            auto stop = chrono::high_resolution_clock::now();

            m_time += chrono::duration_cast<chrono::microseconds>(stop-start).count();

            bool bit_flip = find(bit_flips.begin(), bit_flips.end(), i) != bit_flips.end();

            if(((marginal[0] > marginal[1]) & !bit_flip ) | (marginal[0] < marginal[1]) & bit_flip)
            {
                m_likelihoods[i] = marginal[2]+log(marginal[0]);
                m_decoder_estimates_h[i_r][0].f = ZERO_MESSAGE;
            }
            else
            {
                m_likelihoods[i] = marginal[2]+log(marginal[1]);
                message_dec[i] = 1;
                m_decoder_estimates_h[i_r][0].f = ONE_MESSAGE;
            }

            if(message_dec[i] != message[i-froz_it])
            {
                random_frozen = true;
            }
        }

        update_block(0, i_r, DIRECTION::RIGHT);
    }

    vector<bool> message_dec_m;

    for(int i = 0; i < m_n; i++)
    {
        if(find(m_frozen_bits.begin(), m_frozen_bits.end(), i) == m_frozen_bits.end())
            message_dec_m.push_back(message_dec[i]);
    }

    return message_dec_m;
}

vector<bool> PolarCode::parity_decode_alg(const vector<vector<double>> &message_est, const vector<bool> &message, const vector<vector<double>> &success_lik, vector<bool> random_bits, vector<int> fixed_candidates, bool alg_avg, bool parity_check, bool debug)
{
    m_candidate_likelihoods.clear();
    m_candidates.clear();
    m_iterations_parity=0;

    if(random_bits.size()==0)
        random_bits=random_bool_list(m_n);
    vector<bool> message_dec;

    message_dec = parity_decode(message_est, message, random_bits, {}, parity_check);
    m_candidate_likelihoods.push_back(m_likelihoods);
    double curr_lik;

    if(alg_avg)
    {
        curr_lik = 0;
        for(int k=0; k < m_n; k++)
            curr_lik += m_likelihoods[k]/m_n;
    }

    else
        curr_lik = m_likelihoods[m_n-1];

    if(debug)
    {
        cout<<"Frozen bits:"<<endl;
        print_list(m_frozen_bits);
        cout<<"Likelihoods:"<<endl;
        print_list(m_likelihoods);
        cout<<"Current likelihood: "<<curr_lik<<endl;
    }

    bool improv = true;
    vector<int> candidates_f;

    while((improv && fixed_candidates.empty()) || (m_iterations_parity == 0 && !fixed_candidates.empty()))
    {
        m_iterations_parity++;
        improv=false;

        vector<int> candidates;
        int candidate_best;
        vector<double> likelihood_best;

        if(fixed_candidates.empty())
        {
            for(int i = 0; i < m_n; i++ )
            {
                if((find(m_frozen_bits.begin(), m_frozen_bits.end(), i) == m_frozen_bits.end()) && (find(candidates_f.begin(), candidates_f.end(), i) == candidates_f.end()))
                    candidates.push_back(i);
            }
        }

        else
            candidates=fixed_candidates;

        if(debug)
        {
            cout<<"Candidates: "<<endl;
            print_list(candidates);
        }

        for(int j=0; j < candidates.size(); j++)
        {
            vector<int> candidates_union = candidates_f;
            candidates_union.push_back(candidates[j]);
            vector<bool> message_dec_t = parity_decode(message_est, message, random_bits, candidates_union, parity_check);

            double avg;

            if(alg_avg)
            {
                avg=0;

                for(int k=0; k < m_n; k++)
                    avg+= m_likelihoods[k]/m_n;
            }

            else
                avg=m_likelihoods[m_n-1];

            if(debug)
                cout<<"j: "<<candidates[j]<<" - Avg Likelihood: "<<avg <<endl;
            if(avg > curr_lik)
            {
                if(debug)
                    cout<<"Found better"<<endl;

                improv=true;
                curr_lik = avg;
                likelihood_best = m_likelihoods;
                message_dec=message_dec_t;
                candidate_best = candidates[j];
            }
        }

        if(improv)
        {
            m_candidate_likelihoods.push_back(likelihood_best);
            m_candidates.push_back(candidate_best);
            candidates_f.push_back(candidate_best);
        }
    }

    if(debug)
        cout<<"Iterations: "<<m_iterations_parity<<endl;

    return message_dec;
}

vector<int> PolarCode::estimate_frozen_bits(const AbstractChannel &channel, int trials, bool enable_log, bool show_success)
{
    m_frozen_bits.clear();
    m_active_bits.clear();
    m_bits_success = vector<double>(m_n, 0);
    int k_temp = m_k;
    m_k = m_n;

    for(int j = 0; j < trials; j++)
    {
        vector<bool> message = random_bool_list(m_k);
        vector<bool> message_encod = encode(message);
        vector<double> message_trans = channel.transmit(message_encod);
        vector<vector<double>> encoder_estimates = channel.estimate(message_trans);

        initialize_data(encoder_estimates);

        for(int i = 0; i < m_n; i++)
        {
            int i_r = reverse_bits(i, m_m);

            vector<double> marginal = compute_block(0, i_r, DIRECTION::LEFT);

            if(message[i] == 0)
            {
                if(enable_log)
                    m_bits_success[i] += log(marginal[0]);
                else
                    m_bits_success[i] += marginal[0];

                m_decoder_estimates_h[i_r][0].f = ZERO_MESSAGE;
            }

            else
            {
                if(enable_log)
                    m_bits_success[i] += log(marginal[1]);
                else
                    m_bits_success[i] += marginal[1];

                m_decoder_estimates_h[i_r][0].f = ONE_MESSAGE;
            }

            update_block(0, i_r, DIRECTION::RIGHT);
        }
    }

    m_k = k_temp;

    if(show_success)
        print_list(m_bits_success);

    auto sorted_indices = tag_sort(m_bits_success);

    for(int i = 0; i < m_n-m_k; i++)
        m_frozen_bits.push_back(sorted_indices[i]);

    for(int i = m_n-m_k; i < m_n; i++)
        m_active_bits.push_back(sorted_indices[i]);

    return m_frozen_bits;
}

vector<double> PolarCode::evaluate_ideal_likelihood(const vector<vector<double>> &encoder_estimates, vector<bool> message, vector<int> error_bits)
{
    int f = 0;

    for(int k = 0; k < error_bits.size(); k++)
        message[error_bits[k]] = !message[error_bits[k]];

    initialize_data(encoder_estimates);

    for(int i = 0; i < m_n; i++)
    {
        int i_r = reverse_bits(i, m_m);

        vector<double> marginal = compute_block(0, i_r, DIRECTION::LEFT);
        bool frozen = find(m_frozen_bits.begin(), m_frozen_bits.end(), i) != m_frozen_bits.end();

        if(frozen)
            f++;

        if(frozen | (message[max(i-f,0)] == 0))
        {
            m_likelihoods[i] = marginal[2]+log(marginal[0]);
            m_decoder_estimates_h[i_r][0].f = ZERO_MESSAGE;
        }

        else
        {
            m_likelihoods[i] = marginal[2]+log(marginal[1]);
            m_decoder_estimates_h[i_r][0].f = ONE_MESSAGE;
        }

        update_block(0, i_r, DIRECTION::RIGHT);
    }

    return m_likelihoods;
}

void PolarCode::evaluate_ideal_likelihoods(const AbstractChannel &channel, int decoding_trials, int frozen_trials, vector<int> error_bits, vector<int> frozen_bits)
{
    cout<<"Evaluating ideal likelihoods"<<endl<<endl;
    cout<<"Code: n: "<<m_n<<"  k: "<<m_k<<"  Rate: "<<float(m_k)/float(m_n)<<endl;
    cout<<"Channel: "<<channel.get_name()<<endl;
    cout<<"Decoding trials: "<<decoding_trials<<endl;
    if(!frozen_bits.empty())
    {
        cout<<"User defined frozen bits"<<endl;
        frozen_trials = 1;
        m_frozen_bits=frozen_bits;
    }
    else
        cout<<"Frozen trials: "<<frozen_trials<<endl;
    cout<<"Noise level: "<<channel.get_noise()<<endl;

    if(error_bits.size() >0)
    {
        cout<<"Error bits:"<<endl;
        print_list(error_bits);
    }

    if(frozen_bits.empty())
        estimate_frozen_bits(channel, frozen_trials, false, false);

    vector<vector<double>> lik(decoding_trials, vector<double>(m_n, 0));

    double progress = 0;

    for(int j = 0; j < decoding_trials; j++)
    {
        vector<bool> message = random_bool_list(m_k);
        vector<bool> message_encod = encode(message);
        vector<double> message_trans = channel.transmit(message_encod);
        vector<vector<double>> encoder_estimates = channel.estimate(message_trans);

        evaluate_ideal_likelihood(encoder_estimates, message, error_bits);
        lik[j] = m_likelihoods;

        double curr_progress = float((j+1))/float((decoding_trials));
        if(curr_progress >= progress+0.1)
        {
            progress += 0.1;
            cout<<progress*100<<"%"<<endl;
        }
    }

    cout<<"Ideal likelihoods: "<<endl;
    print_stats(lik);
}

vector<double> PolarCode::get_likelihoods() const
{
    return m_likelihoods;
}

vector<double> PolarCode::get_bits_success() const
{
    return m_bits_success;
}

vector<int> PolarCode::get_active_bits() const
{
    return m_active_bits;
}

vector<double> PolarCode::compute_operation(const vector<double> &a, const vector<double> &b, bool parity)
{
    vector<double> result;

    // XOR
    if(parity == 0)
        result = {a[0]*b[0] + a[1]*b[1], a[0]*b[1] + a[1]*b[0]};
    // EQU
    else
        result = {a[0]*b[0], a[1]*b[1]};

    double norm = result[0] + result[1];
    result[0] /= norm;
    result[1] /= norm;

    result.push_back(a[2]+b[2]+log(norm));

    return result;
}

vector<double> PolarCode::compute_block(int j , int index, DIRECTION dir)
{
    m_iterations++;

    if(j == m_m)
        return m_decoder_estimates_h[index][j].b;

    bool parity;

    if(dir == DIRECTION::LEFT)
    {
        int block_size = pow(2, m_m-j);
        parity = (index % block_size) >= (block_size / 2);
    }

    else if(dir == DIRECTION::RIGHT)
    {
        int block_size = pow(2, m_m-j+1);
        parity = (index % block_size) >= (block_size / 2);
    }

    else
        parity = dir != DIRECTION::BOTTOM;

    // XOR
    if(parity == 0)
    {
        if(dir == DIRECTION::LEFT)
        {
            int i1 = ((int) index / (int)pow(2, m_m-j))*pow(2, m_m-1-j) + index % (int)pow(2, m_m-j);

            vector<double>  a = compute_block(j+1, index, DIRECTION::LEFT);
            vector<double>  b = compute_block(j, i1, DIRECTION::TOP);

            m_decoder_estimates_h[index][j].b_init = true;
            return m_decoder_estimates_h[index][j].b = compute_operation(a, b, parity);
        }
    }

    // EQU
    else
    {
        if(dir == DIRECTION::TOP)
        {
            int i4 = ((int) index / (int)pow(2, m_m-j-1))*pow(2, m_m-j) + index % (int)pow(2, m_m-j-1) + pow(2, m_m-j-1);

            vector<double> a = m_decoder_estimates_h[i4][j].f;
            vector<double> b = compute_block(j+1, i4, DIRECTION::LEFT);

            // v(index,j).b
            return compute_operation(a, b, parity);
        }

        else if(dir == DIRECTION::LEFT)
        {
            int i3 = ((int)(index - pow(2, m_m-j-1)) / (int)pow(2, m_m-j))*pow(2, m_m-j-1) + (int)(index - pow(2, m_m-j-1)) % (int)pow(2, m_m-j);
            int i2 = ((int) i3 / (int)pow(2, m_m-j-1))*pow(2, m_m-j) + i3 % (int)pow(2, m_m-j-1);

            if(m_decoder_estimates_h[i2][j+1].b_init == false)
                m_decoder_estimates_h[i2][j+1].b  = compute_block(j+1, i2, DIRECTION::LEFT);

            vector<double> b = m_decoder_estimates_h[i2][j+1].b;
            vector<double> c = m_decoder_estimates_h[i2][j].f;

            // v(i3, j).f
            vector<double> d = compute_operation(b, c, 0);

            if(m_decoder_estimates_h[index][j+1].b_init == false)
                m_decoder_estimates_h[index][j+1].b  = compute_block(j+1, index, DIRECTION::LEFT);

            vector<double> a = m_decoder_estimates_h[index][j+1].b;

            m_decoder_estimates_h[index][j].b_init = true;
            return m_decoder_estimates_h[index][j].b = compute_operation(a, d, parity);
        }
    }
}

void PolarCode::update_block(int j , int index, DIRECTION dir)
{
    m_iterations++;

    if(j == m_m)
        return;

    int block_size = pow(2, m_m-j);
    bool parity = (index % block_size) >= (block_size / 2);

    // EQU
    if(parity == 1 & dir == DIRECTION::RIGHT)
    {
        int i3 = ((int)(index - pow(2, m_m-j-1)) / (int)pow(2, m_m-j))*pow(2, m_m-j-1) + (int)(index - pow(2, m_m-j-1)) % (int)pow(2, m_m-j);
        int i2 = ((int) i3 / (int)pow(2, m_m-j-1))*pow(2, m_m-j) + i3 % (int)pow(2, m_m-j-1);

        vector<double> a = m_decoder_estimates_h[i2][j].f;
        vector<double> b = m_decoder_estimates_h[index][j].f;

        m_decoder_estimates_h[i2][j+1].f = compute_operation(a, b, 0);
        m_decoder_estimates_h[index][j+1].f = m_decoder_estimates_h[index][j].f;

        update_block(j+1, i2, DIRECTION::RIGHT);
        update_block(j+1, index, DIRECTION::RIGHT);
    }
}

void PolarCode::initialize_data(const vector<vector<double> > &encoder_estimates)
{
    for(int i = 0; i < m_n; i++)
    {
        for(int j = 0; j < m_m+1; j++)
        {
            if(j != m_m)
                m_decoder_estimates_h[i][j] = GraphFlow({UN_MESSAGE, UN_MESSAGE, false});

            else
            {
                m_decoder_estimates_h[i][m_m].f = UN_MESSAGE;
                m_decoder_estimates_h[i][m_m].b = encoder_estimates[i];
                m_decoder_estimates_h[i][m_m].b_init = true;
            }
        }
    }
}

void PolarCode::evaluate_decoder_performance(AbstractChannel &channel, vector<double> noise_levels, int frozen_trials, int decoding_trials, DECODER decoder_type, bool show_progress, vector<vector<int>> frozen_bits, vector<vector<int>> fixed_candidates, bool lower_ml_bound, bool show_l)
{
    if(show_l)
        noise_levels.resize(1);

    if(frozen_bits.size() > 0 & noise_levels.size() > frozen_bits.size())
    {
        cout<<"Frozen bits list incomplete"<<endl;
        return;
    }

    else
        frozen_bits.resize(noise_levels.size());

    if(decoder_type == PARITY_COR_SCD)
        cout<<"Evaluating corrected parity-check SCD decoding performance"<<endl<<endl;
    else if(decoder_type == PARITY_SCD)
        cout<<"Evaluating parity-check SCD decoding performance"<<endl<<endl;
    else if(decoder_type == COR_SCD)
        cout<<"Evaluating corrected SCD decoding performance"<<endl<<endl;
    else
        cout<<"Evaluating SCD decoding performance"<<endl<<endl;

    cout<<"Code: n: "<<m_n<<"  k: "<<m_k<<"  Rate: "<<float(m_k)/float(m_n)<<endl;
    cout<<"Channel: "<<channel.get_name()<<endl;
    cout<<"Decoding trials: "<<decoding_trials<<endl;

    if(!frozen_bits[0].empty())
    {
        cout<<"User defined frozen bits"<<endl;
        frozen_trials = 1;
    }
    else
        cout<<"Frozen trials: "<<frozen_trials<<endl;

    if(!fixed_candidates.empty() && (decoder_type==PARITY_COR_SCD || decoder_type==COR_SCD))
        cout<<"Nb candidates: "<<fixed_candidates[0].size()<<endl;

    cout<<"Noise levels: ";
    print_list(noise_levels);
    cout<<endl;

    auto decoding_time=0;

    double progress = 0;
    vector<vector<double>> results;
    double current_noise = channel.get_noise();

    for(int j = 0; j < noise_levels.size(); j++)
    {
        double fail_b = 0;
        double fail_m = 0;
        double fail_ml = 0;
        double iterations = 0;

        channel.set_noise(noise_levels[j]);

        if(frozen_bits[0].size() == 0)
            estimate_frozen_bits(channel, frozen_trials, false);

        else
            m_frozen_bits = frozen_bits[j];

        for(int i = 0; i < decoding_trials; i++)
        {
            auto start = chrono::high_resolution_clock::now();

            vector<bool> message = random_bool_list(m_k);
            vector<bool> message_e = encode(message);
            vector<double> message_t = channel.transmit(message_e);
            vector<vector<double>> encoder_estimate = channel.estimate(message_t);
            vector<bool> message_estimate;

            if(decoder_type == PARITY_COR_SCD)
            {
                if(fixed_candidates.size() == 0)
                    message_estimate = parity_decode_alg(encoder_estimate, message, {}, {}, {}, false, true);
                else
                    message_estimate = parity_decode_alg(encoder_estimate, message, {}, {}, fixed_candidates[j], false, true);
            }
            else if(decoder_type == PARITY_SCD)
                message_estimate = parity_decode(encoder_estimate, message, {}, {}, true);
            else if(decoder_type == COR_SCD)
            {
                if(fixed_candidates.size() == 0)
                    message_estimate = parity_decode_alg(encoder_estimate, message, {}, {}, {}, false, false);
                else
                    message_estimate = parity_decode_alg(encoder_estimate, message, {}, {}, fixed_candidates[j], false, false);
            }
            else
                message_estimate = decode(encoder_estimate, true);

            auto stop = chrono::high_resolution_clock::now();

            decoding_time += chrono::duration_cast<chrono::microseconds>(stop-start).count();

            int d = nb_differences(message, message_estimate);

            iterations += m_iterations_parity;
            fail_b += d/double(m_k);

            if(d != 0)
                fail_m += 1;

            if(lower_ml_bound && d > 0)
            {
                double l;

                if(decoder_type == SCD || decoder_type == PARITY_SCD)
                    l= m_likelihoods[m_n-1];
                else
                    l = m_candidate_likelihoods[m_candidate_likelihoods.size()-1][m_n-1];

                if(evaluate_ideal_likelihood(encoder_estimate, message)[m_n-1] <= l)
                    fail_ml+=1;
            }

            double curr_progress = float((i+1)*(j+1))/float((noise_levels.size()*decoding_trials));
            if(show_progress && curr_progress >= progress+0.1)
            {
                progress += 0.1;
                cout<<progress*100<<"%"<<endl;
            }
        }

        if(lower_ml_bound)
        {
            if(decoder_type == PARITY_COR_SCD)
                results.push_back({noise_levels[j], fail_m/double(decoding_trials), fail_b/double(decoding_trials), fail_ml/double(decoding_trials),  iterations/double(decoding_trials)});
            else
                results.push_back({noise_levels[j], fail_m/double(decoding_trials), fail_b/double(decoding_trials), fail_ml/double(decoding_trials)});
        }
        else
        {
            if(decoder_type == PARITY_COR_SCD)
                results.push_back({noise_levels[j], fail_m/double(decoding_trials), fail_b/double(decoding_trials), iterations/double(decoding_trials)});
            else
                results.push_back({noise_levels[j], fail_m/double(decoding_trials), fail_b/double(decoding_trials)});
        }
    }

    channel.set_noise(current_noise);

    if(show_progress)
        cout<<"Duration encoding/transmitting/decoding (avg): "<<decoding_time/(decoding_trials*noise_levels.size())*10e-6<<" s"<<endl<<endl;

    print_list(results);
}

void PolarCode::evaluate_parity_decoder_behavior(AbstractChannel &channel, double noise_level, int frozen_trials, int decoding_trials, vector<int> frozen_bits, bool ev_success, bool show_progress)
{
    cout<<"Evaluating parity-check decoding behavior"<<endl<<endl;
    cout<<"Code: n: "<<m_n<<"  k: "<<m_k<<"  Rate: "<<float(m_k)/float(m_n)<<endl;
    cout<<"Channel: "<<channel.get_name()<<endl;
    cout<<"Decoding trials: "<<decoding_trials<<endl;
    if(!frozen_bits.empty())
    {
        cout<<"User defined frozen bits"<<endl;
        frozen_trials = 1;
    }
    else
        cout<<"Frozen trials: "<<frozen_trials<<endl;
    cout<<"Noise level: "<<noise_level<<endl;

    if(frozen_bits.size() == 0)
        estimate_frozen_bits(channel, frozen_trials, false);

    else
        m_frozen_bits = frozen_bits;

    double current_noise = channel.get_noise();
    channel.set_noise(noise_level);

    vector<double> lik_ideal;
    vector<double> lik_scd;
    double progress;
    vector<vector<double>> lik;
    vector<bool> message_fail;
    vector<bool> message_fail_enc;
    vector<vector<double>> message_e_fail;
    vector<double> lik_final;
    vector<int> pos_fail_a;

    do{
        do {
            lik = vector<vector<double>>(decoding_trials, vector<double>(m_n, 0));
            progress = 0;
            vector<int> pos_fail;
            pos_fail_a.clear();

            do
            {
                message_fail = random_bool_list(m_k);
                message_fail_enc = encode(message_fail);
                vector<double> message_t = channel.transmit(message_fail_enc);
                message_e_fail = channel.estimate(message_t);
                vector<bool> message_dec = decode(message_e_fail, true);
                pos_fail = pos_difference(message_fail, message_dec);
            }    while(pos_fail.size() == 0);



            int f = 0;
            for(int i = 0; i < m_n; i++)
            {
                if(find(m_frozen_bits.begin(), m_frozen_bits.end(), i) != m_frozen_bits.end())
                    f++;

                else if(find(pos_fail.begin(), pos_fail.end(), i-f) != pos_fail.end())
                    pos_fail_a.push_back(i);
            }

            cout<<"Failure position: "<<endl;
            print_list(pos_fail_a);

            cout<<"Ideal - Likelihoods: "<<endl;
            lik_ideal = evaluate_ideal_likelihood(message_e_fail, message_fail);
            print_list(lik_ideal);

            decode(message_e_fail, true);
            cout<<"SCD (fail) - Likelihoods: "<<endl;
            lik_scd = m_likelihoods;
            print_list(m_likelihoods);
        }while(!(lik_ideal[m_n-1] > lik_scd[m_n-1]));


        progress = 0;
        lik = vector<vector<double>>(decoding_trials, vector<double>(m_n, 0));

        for(int i = 0; i < decoding_trials; i++)
        {
            vector<bool> message_estimate = parity_decode(message_e_fail, message_fail);

            lik[i] = m_likelihoods;

            double curr_progress = float(i+1)/float(decoding_trials);
            if(show_progress & curr_progress >= progress+0.1)
            {
                progress += 0.1;
                cout<<"Evaluating parity-check SCD behavior: "<<progress*100<<"%"<<endl;
            }
        }
        cout<<"Reference parity-check SCD - Likelihoods: "<<endl;
        lik_final = print_stats(lik)[m_n-1];
    }while(!(lik_final[0]+lik_final[2] < lik_ideal[m_n-1]));

    cout<<"Difference: "<<lik_ideal[m_n-1]-lik_final[0]+lik_final[2]<<endl;

    channel.set_noise(current_noise);
}

void PolarCode::evaluate_parity_decoder_test(AbstractChannel &channel, double noise_level, int frozen_trials, int decoding_trials, bool compute_ideal, bool show_progress, vector<int> frozen_bits, bool parity_check)
{
    cout<<"Evaluating parity-check decoding test"<<endl<<endl;
    cout<<"Code: n: "<<m_n<<"  k: "<<m_k<<"  Rate: "<<float(m_k)/float(m_n)<<endl;
    cout<<"Channel: "<<channel.get_name()<<endl;
    cout<<"Decoding trials: "<<decoding_trials<<endl;

    if(frozen_bits.size() == 0)
        cout<<"Frozen trials: "<<frozen_trials<<endl;
    else
    {
        cout<<"User defined frozen bits"<<endl;
        frozen_trials = 1;
        m_frozen_bits = frozen_bits;
    }
    cout<<"Noise level: "<<noise_level<<endl<<endl;

    channel.set_noise(noise_level);

    double current_noise = channel.get_noise();
    vector<vector<double>> lik(decoding_trials, vector<double>(m_n, 0));
    vector<int> error_track(m_n);

    if(frozen_bits.size() == 0)
        estimate_frozen_bits(channel, frozen_trials, false);

    vector<int> pos_fail_single;
    vector<int> pos_fail_final;
    double ideal;
    double reference;
    double final;


    vector<int> pos_fail_a;
    bool keep=false;

    do
    {
        vector<int> pos_fail;

        vector<bool> message_fail;
        vector<bool> message_fail_enc;
        vector<vector<double>> message_fail_est;

        while(pos_fail.size() == 0)
        {
            message_fail = random_bool_list(m_k);
            message_fail_enc = encode(message_fail);
            vector<double> message_t = channel.transmit(message_fail_enc);
            message_fail_est = channel.estimate(message_t);
            vector<bool> message_estimate = decode(message_fail_est, true);

            pos_fail = pos_difference(message_fail, message_estimate);
        }


        cout<<endl<<"SCD likelihoods:"<<endl;
        print_list(m_likelihoods);


        pos_fail_a.clear();
        int f = 0;
        for(int i = 0; i < m_n; i++)
        {
            if(find(m_frozen_bits.begin(), m_frozen_bits.end(), i) != m_frozen_bits.end())
                f++;

            else if(find(pos_fail.begin(), pos_fail.end(), i-f) != pos_fail.end())
                pos_fail_a.push_back(i);
        }

        cout<<endl<<"Failure position ("<<pos_fail.size()<<"): "<<endl;
        print_list(pos_fail_a);

        vector<bool> random_bits = random_bool_list(m_n);

        cout<<endl<<"Ideal likelihoods:"<<endl;
        print_list(evaluate_ideal_likelihood(message_fail_est, message_fail));
        ideal = m_likelihoods[m_n-1];

        parity_decode(message_fail_est, message_fail, random_bits, {}, parity_check);
        cout<<endl<<"Reference likelihood:"<<endl;
        print_list(m_likelihoods);
        reference = m_likelihoods[m_n-1];
        final = reference;


        vector<bool>message_parity_dec =  parity_decode_alg(message_fail_est, message_fail, {}, random_bits, {}, false, parity_check);
        pos_fail_final = pos_difference(message_fail, message_parity_dec);

        if(m_candidates.size() != 0)
        {
            cout<<endl<<"Final candidates:"<<endl;
            print_list(m_candidates);
            cout<<endl<<"Final likelihoods:"<<endl;
            print_list(m_candidate_likelihoods);
            final = m_candidate_likelihoods[m_candidate_likelihoods.size()-1][m_n-1];
        }

        message_parity_dec = parity_decode(message_fail_est, message_fail, random_bits, {pos_fail_a[0]}, parity_check);
        pos_fail_single = pos_difference(message_fail, message_parity_dec);

        cout<<endl<<"Error ratio (final): "<<float(pos_fail_final.size())/float(pos_fail.size())<<endl;
        cout<<"Iterations: "<<m_iterations_parity<<endl;
        cout<<endl<<"Error ratio (single): "<<float(pos_fail_single.size())/float(pos_fail.size())<<endl;
        cout<<endl<<"Ideal / Final / Reference: "<<ideal<<":"<<final<<":"<<reference<<endl;

        if((ideal <= reference | ideal < final))
            cout<<"Impossible to decode"<<endl;
        else if(pos_fail_final.size() != 0)
        {
            cout<<"Missed decoding opportunity"<<endl;
            keep=true;
        }
        else
        {
            keep=true;
        }
    }
    while(!((pos_fail_single.size() != 0) && (pos_fail_final.size() != 0) && ideal > final));

    channel.set_noise(current_noise);
}

void PolarCode::evaluate_parity_decoder_single_error(AbstractChannel &channel, double noise_level, int frozen_trials, int decoding_trials, bool compute_ideal, bool show_progress, vector<int> frozen_bits)
{
    cout<<"Evaluating parity-check decoding test"<<endl<<endl;
    cout<<"Code: n: "<<m_n<<"  k: "<<m_k<<"  Rate: "<<float(m_k)/float(m_n)<<endl;
    cout<<"Channel: "<<channel.get_name()<<endl;
    cout<<"Decoding trials: "<<decoding_trials<<endl;

    if(frozen_bits.size() == 0)
        cout<<"Frozen trials: "<<frozen_trials<<endl;
    else
    {
        cout<<"User defined frozen bits"<<endl;
        frozen_trials = 1;
        m_frozen_bits = frozen_bits;
    }
    cout<<"Noise level: "<<noise_level<<endl<<endl;

    channel.set_noise(noise_level);

    double current_noise = channel.get_noise();
    vector<vector<double>> lik(decoding_trials, vector<double>(m_n, 0));
    vector<int> error_track(m_n);

    if(frozen_bits.size() == 0)
        estimate_frozen_bits(channel, frozen_trials, false);

    vector<int> pos_fail_single;
    vector<int> pos_fail_final;
    double ideal;
    double reference;
    double final;

    for(int k=0; k < 1e4; k++)
    {
        vector<int> pos_fail_a;
        bool keep=false;

        do
        {
            vector<int> pos_fail;

            vector<bool> message_fail;
            vector<bool> message_fail_enc;
            vector<vector<double>> message_fail_est;

            while(pos_fail.size() == 0)
            {
                message_fail = random_bool_list(m_k);
                message_fail_enc = encode(message_fail);
                vector<double> message_t = channel.transmit(message_fail_enc);
                message_fail_est = channel.estimate(message_t);
                vector<bool> message_estimate = decode(message_fail_est, true);

                pos_fail = pos_difference(message_fail, message_estimate);
            }

            if(show_progress)
            {
                cout<<endl<<"SCD Failure:"<<endl;
                print_list(m_likelihoods);
            }

            pos_fail_a.clear();
            int f = 0;
            for(int i = 0; i < m_n; i++)
            {
                if(find(m_frozen_bits.begin(), m_frozen_bits.end(), i) != m_frozen_bits.end())
                    f++;

                else if(find(pos_fail.begin(), pos_fail.end(), i-f) != pos_fail.end())
                    pos_fail_a.push_back(i);
            }

            cout<<endl<<"Failure position ("<<pos_fail.size()<<"): "<<endl;
            print_list(pos_fail_a);

            vector<bool> random_bits = random_bool_list(m_n);

            cout<<endl<<"Ideal likelihoods:"<<endl;
            print_list(evaluate_ideal_likelihood(message_fail_est, message_fail));
            ideal = m_likelihoods[m_n-1];

            parity_decode(message_fail_est, message_fail, random_bits, {}, true);
            cout<<endl<<"Reference likelihood:"<<endl;
            print_list(m_likelihoods);
            reference = m_likelihoods[m_n-1];
            final = reference;


            vector<bool>message_parity_dec =  parity_decode_alg(message_fail_est, message_fail, {}, random_bits, {}, {}, true);
            pos_fail_final = pos_difference(message_fail, message_parity_dec);

            if(m_candidates.size() != 0)
            {
                cout<<endl<<"Final candidates:"<<endl;
                print_list(m_candidates);
                cout<<endl<<"Final likelihoods:"<<endl;
                print_list(m_candidate_likelihoods);
                final = m_candidate_likelihoods[m_candidate_likelihoods.size()-1][m_n-1];
            }

            message_parity_dec = parity_decode(message_fail_est, message_fail, random_bits, {pos_fail_a[0]});
            pos_fail_single = pos_difference(message_fail, message_parity_dec);

            cout<<endl<<"Error ratio (final): "<<float(pos_fail_final.size())/float(pos_fail.size())<<endl;
            cout<<"Iterations: "<<m_iterations_parity<<endl;
            cout<<endl<<"Error ratio (single): "<<float(pos_fail_single.size())/float(pos_fail.size())<<endl;
            cout<<endl<<"Ideal / Final / Reference: "<<ideal<<":"<<final<<":"<<reference<<endl;

            if((ideal <= reference | ideal < final))
                cout<<"Impossible to decode"<<endl;
            else if(pos_fail_final.size() != 0)
            {
                cout<<"Missed decoding opportunity"<<endl;
                keep=true;
            }
            else
            {
                keep=true;
            }
        }
        while(!((pos_fail_single.size() == 0) && keep));

        if(keep)
        {
            error_track[pos_fail_a[0]] = error_track[pos_fail_a[0]]+1;
            cout<<"Count:"<<k+1<<endl;
        }
    }

    cout<<"Error counts:"<<endl;
    print_list(error_track);

    channel.set_noise(current_noise);
}

void PolarCode::compare_frozen_likelihood(AbstractChannel &channel, double noise_level, int frozen_trials, int decoding_trials)
{
    vector<vector<double>> frozen_likelihoods;
    double current_noise = channel.get_noise();
    channel.set_noise(noise_level);

    vector<double> likelihoods(m_n, 0);
    vector<double> ratio_l(m_n, 0);
    vector<double> bits_success = vector<double>(m_n, 0);
    int k_temp = m_k;
    m_k = m_n;
    double r = 1;
    double l = 1;

    for(int j = 0; j < frozen_trials; j++)
    {
        vector<bool> message = random_bool_list(m_k);
        vector<bool> message_encod = encode(message);
        vector<double> message_trans = channel.transmit(message_encod);
        vector<vector<double>> encoder_estimates = channel.estimate(message_trans);

        initialize_data(encoder_estimates);

        for(int i = 0; i < m_n; i++)
        {
            int i_r = reverse_bits(i, m_m);

            vector<double> marginal = compute_block(0, i_r, DIRECTION::LEFT);

            if(message[i] == 0)
            {
                bits_success[i] += marginal[0];
                l = exp(marginal[2])*marginal[0];
                m_decoder_estimates_h[i_r][0].f = ZERO_MESSAGE;
            }

            else
            {
                bits_success[i] += marginal[1];
                l = exp(marginal[2])*marginal[1];
                m_decoder_estimates_h[i_r][0].f = ONE_MESSAGE;
            }

            likelihoods[i] += l;
            ratio_l[i] += l/r;
            r = l;

            update_block(0, i_r, DIRECTION::RIGHT);
        }
    }

    m_k = k_temp;
    m_bits_success = bits_success;

    channel.set_noise(current_noise);

    for(int i = 1; i < m_n; i++)
        frozen_likelihoods.push_back({(double)i, 2*m_bits_success[i]/frozen_trials, ratio_l[i]/frozen_trials});

    print_list(frozen_likelihoods);

}

