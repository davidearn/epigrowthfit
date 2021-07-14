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
    int n = x.size();
    for (int i = 0; i < n; ++i)
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
    if (i < 0)
    {
        return 1 + nchar(-i);
    }
    return 1;
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
    for (int i = 0; i < n; ++i)
    {
        log_x(i) = logspace_sub(log_x(i+1), log_x(i));
    }
    log_x.conservativeResize(n);
    return log_x;
}

/* Compute `c(log(diff(x1)), ..., log(diff(xn)))`
   given `c(log(x1), ..., log(xn)) and `c(length(x1), ..., length(xn))`
*/
template<class Type>
vector<Type> logspace_diff(vector<Type> log_x, vector<int> len)
{
    int N = len.size();
    int n;
    for (int s = 0, i = 0; s < N; ++s)
    {
        n = len(s);
        log_x.segment(i, n - 1) = logspace_diff((vector<Type>) log_x.segment(i + s, n));
	i += n - 1;
    }
    log_x.conservativeResize(log_x.size() - N);
    return log_x;
}

/* Compute `log(cumsum(x))` given `log(x)` */
template<class Type>
vector<Type> logspace_cumsum(vector<Type> log_x)
{
    int n = log_x.size();
    for (int i = 1; i < n; ++i)
    {
        log_x(i) = logspace_add(log_x(i-1), log_x(i));
    }
    return log_x;
}

/* Compute `c(log(cumsum(x1)), ..., log(cumsum(xn)))`
   given `c(log(x1), ..., log(xn)) and `c(length(x1), ..., length(xn))`
*/
template<class Type>
vector<Type> logspace_cumsum(vector<Type> log_x, vector<int> len)
{
    int N = len.size();
    int n;
    for (int s = 0, i = 0; s < N; ++s)
    {
        n = len(s);
        log_x.segment(i, n) = logspace_cumsum((vector<Type>) log_x.segment(i, n));
	i += n;
    }
    return log_x;
}

/* Recycle a vector `x` to length `len` starting from index `from` */
template<class Type>
vector<Type> rep_len_from(vector<Type> x, int len, int from = 0)
{
    int n = x.size();
    if (from > 0)
    {
	vector<Type> x1 = x.head(from);
	vector<Type> x2 = x.tail(n - from);
	x << x2,x1;
    }
    if (len == n)
    {
        return x;
    }
    if (len < n)
    {
        return x.head(len);
    }
    int len0 = n * (len / n);
    vector<Type> res(len);
    for (int i = 0; i < len0; i += n)
    {
        res.segment(i, n) = x;
    }
    if (len0 < len)
    {
        res.tail(len - len0) = x.head(len - len0);
    }
    return res;
}

template<class Type>
bool is_nll_term_ok(Type nll_term, double tol = 1.0e+09)
{
    return is_finite(nll_term) && asDouble(nll_term) < tol;
}

template<class Type>
Eigen::SparseMatrix<Type> as_sparse_matrix(vector<Type> x, vector<int> index, vector<int> nnz)
{
    int n = x.size();
    Eigen::SparseMatrix<Type> res(x.size(), nnz.size());
    res.reserve(nnz);
    for (int i = 0; i < n; ++i)
    {
        res.insert(i, index(i)) = x(i);
    }
    return res;
}
