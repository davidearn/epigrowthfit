/* Similar to R function `is.na` */
template<class Type>
bool is_na(Type x)
{
    return R_IsNA(asDouble(x));
}

/* Similar to R function `is.finite` */
template<class Type>
bool is_finite(Type x)
{
    return R_finite(asDouble(x));
}

/* Decrease all elements of an integer vector `x` by an integer `k` */
vector<int> decrement(vector<int> x, int k = 1)
{
    for (int i = 0; i < x.size(); i++)
    {
        x(i) -= k;
    }
    return x;
}

/* Get the number of characters in an integer */
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

/* Check if an integer vector has at least one non-negative element */
bool any_geq_zero(vector<int> x)
{
    int n = x.size();
    if (n == 0)
    {
        return false;
    }
    bool yes = false;
    for (int i = 0; i < n; i++)
    {
        if (x(i) >= 0)
	{
	    yes = true;
	    break;
	}
    }
    return yes;
}

/* Poisson density with robust parametrization */
template<class Type>
Type dpois_robust(Type x, Type log_lambda, int give_log = 0)
{
    Type log_dpois = x * log_lambda - exp(log_lambda) - lfactorial(x);
    return ( give_log ? log_dpois : exp(log_dpois) );
}

/* Negative binomial sampling with robust parametrization */
template<class Type>
Type rnbinom_robust(Type log_mu, Type log_size)
{
    Type log_prob = log_size - logspace_add(log_mu, log_size);
    return rnbinom(exp(log_size), exp(log_prob));
    /* usage: rnbinom(size, prob) */
}

/* Compute `log(diff(x))` given `log(x)` */
template<class Type>
vector<Type> logspace_diff(vector<Type> log_x)
{
    int n = log_x.size() - 1;
    vector<Type> log_diff_x(n);
    for (int i = 0; i < n; i++)
    {
        log_diff_x(i) = logspace_sub(log_x(i+1), log_x(i));
    }
    return log_diff_x;
}

/* Compute `c(log(diff(x1)), ..., log(diff(xn)))`
   given `c(log(x1), ..., log(xn)) and `c(length(x1), ..., length(xn))`
*/
template<class Type>
vector<Type> logspace_diff(vector<Type> log_x, vector<int> len)
{
    int N = len.size();
    int n;
    vector<Type> log_diff_x(log_x.size() - N);
    for (int s = 0, i = 0; s < N; s++)
    {
        n = len(s);
        vector<Type> log_x_segment = log_x.segment(i + s, n);
        log_diff_x.segment(i, n - 1) = logspace_diff(log_x_segment);
	i += n - 1;
    }
    return log_diff_x;
}

/* Compute `log(cumsum(x))` given `log(x)` */
template<class Type>
vector<Type> logspace_cumsum(vector<Type> log_x)
{
    int n = log_x.size();
    vector<Type> log_cumsum_x(n);
    log_cumsum_x(0) = log_x(0);
    for (int i = 1; i < n; i++)
    {
        log_cumsum_x(i) = logspace_add(log_cumsum_x(i-1), log_x(i));
    }
    return log_cumsum_x;
}

/* Compute `c(log(cumsum(x1)), ..., log(cumsum(xn)))`
   given `c(log(x1), ..., log(xn)) and `c(length(x1), ..., length(xn))`
*/
template<class Type>
vector<Type> logspace_cumsum(vector<Type> log_x, vector<int> len)
{
    int N = len.size();
    int n;
    vector<Type> log_cumsum_x(log_x.size());
    for (int s = 0, i = 0; s < N; s++)
    {
        n = len(s);
        vector<Type> log_x_segment = log_x.segment(i, n);
        log_cumsum_x.segment(i, n) = logspace_cumsum(log_x_segment);
	i += n;
    }
    return log_cumsum_x;
}

/* Compute `log(x + a)` given vector `log(x)` and scalar `log(a)` */
template<class Type>
vector<Type> logspace_add(vector<Type> log_x, Type log_a)
{
    int n = log_x.size();
    vector<Type> log_x_plus_a(n);
    for (int i = 0; i < n; i++)
    {
        log_x_plus_a(i) = logspace_add(log_x(i), log_a);
    }
    return log_x_plus_a;
}

/* Compute `log(x + a)` given vectors `log(x)` and `log(a)` of equal length */
template<class Type>
vector<Type> logspace_add(vector<Type> log_x, vector<Type> log_a)
{
    int n = log_x.size();
    vector<Type> log_x_plus_a(n);
    for (int i = 0; i < n; i++)
    {
        log_x_plus_a(i) = logspace_add(log_x(i), log_a(i));
    }
    return log_x_plus_a;
}

/* Compute `log(x - a)` given vector `log(x)` and scalar `log(a)` */
template<class Type>
vector<Type> logspace_sub(vector<Type> log_x, Type log_a)
{
    int n = log_x.size();
    vector<Type> log_x_minus_a(n);
    for (int i = 0; i < n; i++)
    {
        log_x_minus_a(i) = logspace_sub(log_x(i), log_a);
    }
    return log_x_minus_a;
}

/* Compute `log(x - a)` given vectors `log(x)` and `log(a)` of equal length */
template<class Type>
vector<Type> logspace_sub(vector<Type> log_x, vector<Type> log_a)
{
    int n = log_x.size();
    vector<Type> log_x_minus_a(n);
    for (int i = 0; i < n; i++)
    {
        log_x_minus_a(i) = logspace_sub(log_x(i), log_a(i));
    }
    return log_x_minus_a;
}

/* Recycle a vector `x` to length `len` starting from index `from` */
template<class Type>
vector<Type> rep_len_from(vector<Type> x, int len, int from = 0)
{
    int n = x.size();
    if (from > 0)
    {
	vector<Type> x1 = x.segment(0, from);
	vector<Type> x2 = x.segment(from, n - from);
	x << x2,x1;
    }
    if (len < n)
    {
        return x.segment(0, len);
    }
    if (len == n)
    {
        return x;
    }
    
    int len0 = n * ((int) len / n);
    int len_minus_len0 = len - len0;

    vector<Type> x_rep(len);
    for (int i = 0; i < len0; i += n)
    {
	x_rep.segment(i, n) = x;
    }
    if (len_minus_len0 > 0)
    {
        x_rep.segment(len0, len_minus_len0) = x.segment(0, len_minus_len0);
    }
    return x_rep;
}

template<class Type>
bool is_nll_term_ok(Type nll_term, double tol = 1.0e+09)
{
    return is_finite(nll_term) && asDouble(nll_term) < tol;
}
