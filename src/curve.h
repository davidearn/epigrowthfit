double log_exponential(vector<double> x)
{
    // log(c(t))
    // = log(c0 * exp(r * t))
    // = log(c0) + r * t
    double t = x(0);
    double log_r = x(1);
    double log_c0 = x(2);
    return log_c0 + exp(log_r) * t;
}


TMB_ATOMIC_VECTOR_FUNCTION(
    // ATOMIC_NAME
    log_exponential
    ,
    // OUTPUT_DIM
    1,
    // ATOMIC_DOUBLE
    ty(0) = log_exponential(tx);
    ,
    // ATOMIC_REVERSE
    Type f = ty[0];
    Type df = 1. / (exp(W) * (1. + W)); // Derivative
    px[0] = DW * py[0];                 // Reverse mode chain rule
)



template<class Type>
Type log_exponential(Type t, Type log_r, Type log_c0)
{
    // log(c(t))
    // = log(c0 * exp(r * t))
    // = log(c0) + r * t
    return log_c0 + exp(log_r) * t;
}

template<class Type>
Type log_subexponential(Type t1, Type log_alpha, Type log_c0, Type logit_p)
{
    // log(c(t))
    // = log((c0^(1 - p) + (1 - p) * alpha * t)^(1 / (1 - p)))
    // = log(c0^(1 - p) + (1 - p) * alpha * t) / (1 - p)
    Type one_minus_p = Type(1) / (Type(1) + exp(logit_p));
    return logspace_add(one_minus_p * log_c0, log(one_minus_p) + log_alpha + log(t)) / one_minus_p;
}

template<class Type>
Type log_gompertz(Type t, Type log_alpha, Type log_c0, Type log_K)
{
    // log(c(t))
    // = log(K * (c0 / K)^(exp(-alpha * t)))
    // = log(K) + exp(-alpha * t) * (log(c0) - log(K))
    return log_K + exp(-exp(log_alpha) * t) * (log_c0 - log_K);
}

template<class Type>
Type log_logistic(Type t, Type log_r, Type log_tinfl, Type log_K)
{
    // log(c(t))
    // = log(K / (1 + exp(-r * (t - tinfl))))
    // = log(K) - log(1 + exp(-r * (t - tinfl)))
    return log_K - logspace_add(Type(0), -exp(log_r) * (t - exp(log_tinfl)));
}

template<class Type>
Type log_richards(Type t, Type log_r, Type log_tinfl, Type log_K, Type log_a)
{
    // log(c(t))
    // = log(K / (1 + a * exp(-r * a * (t - tinfl)))^(1 / a))
    // = log(K) - log(1 + a * exp(-r * a * (t - tinfl))) / a
    Type a = exp(log_a);
    return log_K - logspace_add(Type(0), log_a - exp(log_r) * a * (t - exp(log_tinfl))) / a;
}

template<class Type>
vector<Type> eval_log_cum_inc(vector<Type> t,
			      matrix<Type> Y,
			      int curve_flag,
			      bool excess,
			      int j_log_r,
			      int j_log_alpha,
			      int j_log_c0,
			      int j_log_tinfl,
			      int j_log_K,
			      int j_logit_p,
			      int j_log_a,
			      int j_log_b)
{
    vector<Type> log_cum_inc(t.size());
    for (int i = 0; i < t.size(); i++)
    {
        switch (curve_flag)
	{
	case exponential:
	    log_cum_inc(i) = eval_exponential(t(i), Y(i, j_log_r), Y(i, j_log_c0));
	    break;
	case subexponential:
	    log_cum_inc(i) = eval_subexponential(t(i), Y(i, j_log_alpha), Y(i, j_log_c0), Y(i, j_logit_p));
	    break;
	case gompertz:
	    log_cum_inc(i) = eval_gompertz(t(i), Y(i, j_log_alpha), Y(i, j_log_c0), Y(i, j_log_K));
	    break;
	case logistic:
	    log_cum_inc(i) = eval_logistic(t(i), Y(i, j_log_r), Y(i, j_log_tinfl), Y(i, j_log_K));
	    break;
	case richards:
	    log_cum_inc(i) = eval_richards(t(i), Y(i, j_log_r), Y(i, j_log_tinfl), Y(i, j_log_K), Y(i, j_log_a));
	    break;
	}
	if (excess)
	{
	    log_cum_inc(i) = logspace_add(Y(i, j_log_b) + log(t(i)), log_cum_inc(i));
	}
    }
    return log_cum_inc;
}

template<class Type>
vector<Type> eval_log_int_inc(vector<Type> log_cum_inc, vector<int> wl)
{
    vector<Type> log_int_inc(log_cum_inc.size() - wl.size());
    for (int i1 = 0, i2 = 1, k = 0; k < wl.size(); k++)
    {
        for (int j = 0; j < wl(k) - 1; j++)
	{
	    log_int_inc(i1+j) = logspace_sub(log_cum_inc(i2+j), log_cum_inc(i2+j-1));
	}
	i1 += wl(k) - 1;
	i2 += wl(k);
    }
    return log_int_inc;
}

template<class Type>
vector<Type> eval_log_rt(vector<Type> log_cum_inc,
		         matrix<Type> Y,
		         int curve_flag,
			 int j_log_r,
			 int j_log_alpha,
			 int j_log_K,
			 int j_logit_p,
			 int j_log_a)
{
    vector<Type> log_rt(log_cum_inc.size());
    Type one_minus_p;
    for (int i = 0; i < log_cum_inc.size(); i++)
    {
        switch (curve_flag)
	{
	case exponential:
	    // log(c'(t) / c(t)) = log(r)
	    log_rt(i) = Y(i, j_log_r);
	    break;
	case subexponential:
	    // log(c'(t) / c(t)) = log(alpha) - (1 - p) * log(c(t))
	    one_minus_p = Type(1) / (Type(1) + exp(Y(i, j_logit_p)));
	    log_rt(i) = Y(i, j_log_alpha) - one_minus_p * log_cum_inc(i);
	    break;
	case gompertz:
	    // log(c'(t) / c(t)) = log(alpha) + log(log(K) - log(c(t)))
	    log_rt(i) = Y(i, j_log_alpha) + log(Y(i, j_log_K) - log_cum_inc(i));
	    break;
	case logistic:
	    // log(c'(t) / c(t)) = log(r) + log(1 - c(t) / K)
	    log_rt(i) = Y(i, j_log_r) + logspace_sub(Type(0), log_cum_inc(i) - Y(i, j_log_K));
	    break;
	case richards:
	    // log(c'(t) / c(t)) = log(r) + log(1 - (c(t) / K)^a)
	    log_rt(i) = Y(i, j_log_r) + logspace_sub(Type(0), exp(Y(i, j_log_a)) * (log_cum_inc(i) - Y(i, j_log_K)));
	    break;
	}
    }
    return log_rt;
}
