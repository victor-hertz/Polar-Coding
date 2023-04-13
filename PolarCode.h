#ifndef POLARCODE_H
#define POLARCODE_H

#include "AbstractPolarCode.h"

struct GraphFlow{
    vector<double> f={0.5,0.5, 0};
    vector<double> b={0.5,0.5, 0};
    bool b_init = false;
};

enum DECODER {SCD=0, PARITY_COR_SCD=1, PARITY_SCD=2, COR_SCD=3};
enum DIRECTION { LEFT = 0, RIGHT, TOP, BOTTOM };

class PolarCode: public AbstractPolarCode
{
public:
    PolarCode(int m, int k, vector<int> frozen_bits = vector<int>(), bool debug = false);

    vector<bool> decode(const vector<vector<double> > &encoder_estimates);
    vector<bool> decode(const vector<vector<double> > &encoder_estimates, bool compute_l);
    vector<bool> parity_decode(const vector<vector<double>> &message_est, const vector<bool> &message, vector<bool> random_bits={},  vector<int> bit_flips={}, bool parity_check=true);
    vector<bool> parity_decode_alg(const vector<vector<double>> &message_est, const vector<bool> &message, const vector<vector<double>> &success_lik, vector<bool> random_bits={},  vector<int> fixed_candidates={},bool alg_avg = false, bool parity_check=true, bool debug=false);

    vector<int> estimate_frozen_bits(const AbstractChannel &channel, int trials, bool enable_log = false, bool show_success=false);

    vector<double> get_likelihoods() const;
    vector<double> get_bits_success() const;
    vector<int> get_active_bits() const;

    vector<double> evaluate_ideal_likelihood(const vector<vector<double> > &encoder_estimates, vector<bool> message, vector<int> error_bits = {});
    void evaluate_ideal_likelihoods(const AbstractChannel &channel, int decoding_trials, int frozen_trials, vector<int> error_bits = {}, vector<int> frozen_bits = {});
    void evaluate_decoder_performance(AbstractChannel &channel, vector<double> noise_levels, int frozen_trials, int decoding_trials, DECODER decoder_type, bool show_progress = 0, vector<vector<int>> frozen_bits = {}, vector<vector<int> > fixed_candidates={},  bool lower_ml_bound = false, bool show_l = false);
    void evaluate_parity_decoder_behavior(AbstractChannel &channel, double noise_level, int frozen_trials, int decoding_trials, vector<int> frozen_bits = {}, bool success_ev = true, bool show_progress = 0);
    void evaluate_parity_decoder_behavior2(AbstractChannel &channel, double noise_level, int frozen_trials, int decoding_trials, vector<int> frozen_bits = {}, bool success_ev = true, bool show_progress = 0);
    void evaluate_parity_decoder_test(AbstractChannel &channel, double noise_level, int frozen_trials, int decoding_trials, bool compute_ideal, bool show_progress = 0, vector<int> frozen_bits = {}, bool parity_check=true);
    void evaluate_parity_decoder_single_error(AbstractChannel &channel, double noise_level, int frozen_trials, int decoding_trials, bool compute_ideal, bool show_progress = 0, vector<int> frozen_bits = {});
    void compare_frozen_likelihood(AbstractChannel &channel, double noise_level, int frozen_trials, int decoding_trials);

private:
    vector<double> compute_operation(const vector<double> &a, const vector<double> &b, bool parity);

    vector<double> compute_block(int j, int index, DIRECTION dir);
    void update_block(int j, int index, DIRECTION dir);

    void initialize_data(const vector<vector<double> > &encoder_estimates);

    vector<int> m_active_bits;
    int m_iterations;
    int m_iterations_parity;

    vector<vector<double>> m_candidate_likelihoods;
    vector<int> m_candidates;

    vector<vector<GraphFlow>> m_decoder_estimates_h;
    vector<double> m_likelihoods;

    vector<double> m_bits_success;
    double m_time;
};

#endif // POLARCODE_H
