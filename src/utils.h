template<class Type>
Type dpois_robust(Type x, Type log_lambda, int give_log = 0)
{
    Type log_dpois = x * log_lambda - exp(log_lambda) - lfactorial(x);
    return ( give_log ? log_dpois : exp(log_dpois) );
}

template<class Type>
Type rnbinom_robust(Type log_mu, Type log_size)
{
    Type log_prob = log_size - logspace_add(log_mu, log_size);
    return rnbinom(exp(log_size), exp(log_prob));
    // usage: rnbinom(size, prob)
}

// https://github.com/kaskr/adcomp/issues/59
template<class Type>
bool is_NA_real_(Type x)
{
    return R_IsNA(asDouble(x));
}

template<class Type>
bool is_finite(Type x)
{
    return R_finite(asDouble(x));
}

int nchar(int i)
{
    if (i > 0)
    {
        return 1 + (int) log10((double) i);
    }
    else if (i < 0)
    {
        return 1 + nchar(-i);
    }
    else
    {
        return 1;
    }
}

template<class Type>
vector<Type> logspace_diff_1(vector<Type> log_x)
{
    vector<Type> log_diff_x(log_x.size() - 1);
    for (int i = 0; i < log_x.size() - 1; i++)
    {
        log_diff_x(i) = logspace_sub(log_x(i+1), log_x(i));
    }
    return log_diff_x;
}

template<class Type>
vector<Type> logspace_diff_n(vector<Type> log_x, vector<int> len)
{
    vector<Type> log_diff_x(log_x.size() - len.size());
    for (int s = 0, i = 0; s < len.size(); s++) // loop over segments
    {
        vector<Type> log_x_segment = log_x.segment(i + s, len(s));
        log_diff_x.segment(i, len(s) - 1) = logspace_diff_1(log_x_segment); 
        i += len(s) - 1; // increment reference index
    }
    return log_diff_x;
}

template<class Type>
vector<Type> logspace_cumsum_1(vector<Type> log_x)
{
    vector<Type> log_cumsum_x(log_x.size());
    log_cumsum_x(0) = log_x(0);
    for (int i = 1; i < log_x.size(); i++)
    {
        log_cumsum_x(i) = logspace_add(log_cumsum_x(i-1), log_x(i));
    }
    return log_cumsum_x;
}

template<class Type>
vector<Type> logspace_cumsum_n(vector<Type> log_x, vector<int> len)
{
    vector<Type> log_cumsum_x(log_x.size());
    for (int s = 0, i = 0; s < len.size(); s++) // loop over segments
    {
        vector<Type> log_x_segment = log_x.segment(i, len(s));
        log_cumsum_x.segment(i, len(s)) = logspace_cumsum_1(log_x_segment); 
        i += len(s); // increment reference index
    }
    return log_cumsum_x;
}

template<class Type>
vector<Type> logspace_sub_n_1(vector<Type> log_x, Type log_a)
{
    vector<Type> log_x_minus_a(log_x.size());
    for (int i = 0; i < log_x.size(); i++)
    {
        log_x_minus_a = logspace_sub(log_x, log_a);
    }
    return log_x_minus_a;
}
