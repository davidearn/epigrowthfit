namespace egf
{

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

/* Compute `log(diff(x))` given `log(x)` */
template<class Type>
void logspace_diff(vector<Type> &log_x)
{
    int n = log_x.size() - 1;
    for (int i = 0; i < n; ++i)
    {
        log_x(i) = logspace_sub(log_x(i+1), log_x(i));
    }
    log_x.conservativeResize(n);
}

/* Log multivariate gamma function, requires `x > 0.5 * (n - 1)` */
template<class Type>
Type mvlgamma(Type x, int n)
{
    Type res = lgamma(x);
    if (n == 1)
    {
        return res;
    }
    for (int i = 1; i < n; ++i)
    {
        res += lgamma(x - Type(0.5 * i));
    }
    res += Type(0.25 * n * (n - 1)) * log(Type(M_PI));
    return res;
}

/* Compute `log(diag(L %*% t(L)))` for n-by-n unit diagonal matrix `L`
   given elements of `lower.tri(L)` in row-major order 
*/
template<class Type>
vector<Type> log_diag_LLT(const vector<Type> &x)
{
    int len = x.size();
    int n = 0.5 * (1.0 + sqrt(1.0 + 8.0 * len));
    vector<Type> res(n);
    res(0) = Type(0.0);
    for (int i = 1, k = 0; i < n; ++i)
    { /* loop over rows of L */
        /* Compute log inner product of row with itself */
        res(i) = logspace_add(Type(0.0), Type(2.0) * log(x(k)));
	++k;
	for (int j = 1; j < i; ++j, ++k)
	{
	    res(i) = logspace_add(res(i), Type(2.0) * log(x(k)));
	}
    } /* loop over rows of L */
    return res;
}

} // namespace egf
