namespace egf
{

/* Lewandowski-Kurowicka-Joe (LKJ) density function (non-normalized),
   where `x` parametrizes an n-by-n SPD correlation matrix `X`. 
   `x` should contain the `n*(n-1)/2` lower triangular elements 
   of the unit diagonal Cholesky factor `L` of `X` in row-major order.
*/
template<class Type>
Type dlkj(const vector<Type> &x, Type eta, int give_log = 0)
{
    int len = x.size();
    if (len == 0)
    {
        return ( give_log ? Type(0.0) : Type(1.0) );
    }
    Type log_det_X = -log_diag_LLT(x).sum();
    Type log_res = (eta - Type(1.0)) * log_det_X;
    return ( give_log ? log_res : exp(log_res) );
}

/* Wishart density function, where `x` and `scale` 
   parametrize corresponding n-by-n SPD covariance matrices `X` and `S`.
   `x` should contain the `n` log standard deviations associated
   with `X` followed by the `n*(n-1)/2` lower triangular elements 
   of the unit diagonal Cholesky factor `L` of `X` in row-major order
   (ditto for `scale` and `S`).
*/
template<class Type>
Type dwishart(const vector<Type> &x, Type df, const vector<Type> &scale, int give_log = 0)
{
    int len = x.size();
    int n = 0.5 * (-1.0 + sqrt(1.0 + 8.0 * len));

    vector<Type> log_diag_LLT_X = log_diag_LLT((vector<Type>) x.tail(len - n));
    vector<Type> log_diag_LLT_S = log_diag_LLT((vector<Type>) scale.tail(len - n));

    Type log_det_X = Type(2.0) * x.head(n).sum() - log_diag_LLT_X.sum();
    Type log_det_S = Type(2.0) * scale.head(n).sum() - log_diag_LLT_S.sum();

    /* Remains to compute `tr(invS * X)` ... */
    
    matrix<Type> L_X(n, n);
    L_X.setIdentity();
    for (int i = 0, k = n; i < n; ++i)
    {
        for (int j = 0; j < i; ++j, ++k)
	{
	    L_X(i, j) = scale(k);
	}
    }

    matrix<Type> L_S(n, n);
    L_S.setIdentity();
    for (int i = 0, k = n; i < n; ++i)
    {
        for (int j = 0; j < i; ++j, ++k)
	{
	    L_S(i, j) = x(k);
	}
    }
    
    Type log_det_L_S; /* gets 0, hopefully */
    matrix<Type> invL_S = atomic::matinvpd(L_S, log_det_L_S);
    
    matrix<Type> A = (invL_S.transpose() * invL_S).array() * (L_X * L_X.transpose()).array();
    vector<Type> log_diag_D = x.head(n) - scale.head(n) - Type(0.5) * (log_diag_LLT_X - log_diag_LLT_S);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
	{
	    A(i, j) *= exp(log_diag_D(i) + log_diag_D(j));
	}
    }

    Type log_res = -Type(0.5) *
        (
	 df * log_det_S +
	 (-df + Type(n + 1)) * log_det_X +
	 Type(n) * df * Type(M_LN2) +
	 Type(2.0) * mvlgamma(Type(0.5) * df, n) +
	 A.sum() /* this term is `tr(invS * X)` */
	);
    return ( give_log ? log_res : exp(log_res) ); 
}

/* Inverse Wishart density, where `x` and `scale` 
   parametrize corresponding n-by-n SPD covariance matrices `X` and `S`.
   `x` should contain the `n` log standard deviations associated
   with `X` followed by the `n*(n-1)/2` lower triangular elements 
   of the unit diagonal Cholesky factor `L` of `X` in row-major order
   (ditto for `scale` and `S`).
*/
template<class Type>
Type dinvwishart(const vector<Type> &x, Type df, const vector<Type> &scale, int give_log = 0)
{
    int len = x.size();
    int n = 0.5 * (-1.0 + sqrt(1.0 + 8.0 * len));

    vector<Type> log_diag_LLT_X = log_diag_LLT((vector<Type>) x.tail(len - n));
    vector<Type> log_diag_LLT_S = log_diag_LLT((vector<Type>) scale.tail(len - n));

    Type log_det_X = Type(2.0) * x.head(n).sum() - log_diag_LLT_X.sum();
    Type log_det_S = Type(2.0) * scale.head(n).sum() - log_diag_LLT_S.sum();

    /* Remains to compute `tr(S * invX)` ... */
    
    matrix<Type> L_X(n, n);
    L_X.setIdentity();
    for (int i = 0, k = n; i < n; ++i)
    {
        for (int j = 0; j < i; ++j, ++k)
	{
	    L_X(i, j) = scale(k);
	}
    }

    matrix<Type> L_S(n, n);
    L_S.setIdentity();
    for (int i = 0, k = n; i < n; ++i)
    {
        for (int j = 0; j < i; ++j, ++k)
	{
	    L_S(i, j) = x(k);
	}
    }
    
    Type log_det_L_X; /* gets 0, hopefully */
    matrix<Type> invL_X = atomic::matinvpd(L_X, log_det_L_X);
    
    matrix<Type> A = (L_S * L_S.transpose()).array() * (invL_X.transpose() * invL_X).array();
    vector<Type> log_diag_D = scale.head(n) - x.head(n) - Type(0.5) * (log_diag_LLT_S - log_diag_LLT_X);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
	{
	    A(i, j) *= exp(log_diag_D(i) + log_diag_D(j));
	}
    }

    Type log_res = -Type(0.5) *
        (
	 df * log_det_S +
	 (df + Type(n + 1)) * log_det_X +
	 Type(n) * df * Type(M_LN2) -
	 Type(2.0) * mvlgamma(Type(0.5) * df, n) -
	 A.sum() /* this term is `tr(S * invX)` */
	);
    return ( give_log ? log_res : exp(log_res) ); 
}

/* Poisson density function with robust parametrization */
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
}

} // namespace egf
