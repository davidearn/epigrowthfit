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
			      vector<int> slen,
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
    for (int s = 0, k = 0; s < slen.size(); s++) // loop over segments
    {
        for (int i = 0; i < slen(s) + 1; i++) // loop over within-segment index
	{
	    switch (curve_flag)
	    {
	    case exponential:
	        log_cum_inc(k + i) = log_exponential(t(k + i), Y(s, j_log_r), Y(s, j_log_c0));
		break;
	    case subexponential:
	        log_cum_inc(k + i) = log_subexponential(t(k + i), Y(s, j_log_alpha), Y(s, j_log_c0), Y(s, j_logit_p));
		break;
	    case gompertz:
	        log_cum_inc(k + i) = log_gompertz(t(k + i), Y(s, j_log_alpha), Y(s, j_log_c0), Y(s, j_log_K));
		break;
	    case logistic:
	        log_cum_inc(k + i) = log_logistic(t(k + i), Y(s, j_log_r), Y(s, j_log_tinfl), Y(s, j_log_K));
		break;
	    case richards:
	        log_cum_inc(k + i) = log_richards(t(k + i), Y(s, j_log_r), Y(s, j_log_tinfl), Y(s, j_log_K), Y(s, j_log_a));
		break;
	    }
	    if (excess)
	    {
	        log_cum_inc(i) = logspace_add(Y(i, j_log_b) + log(t(i)), log_cum_inc(i));
	    }
	}
	k += slen(s); // increment `t` index
    }
    return log_cum_inc;
}

template<class Type>
vector<Type> eval_log_int_inc(vector<Type> log_cum_inc,
			      vector<int> slen,
			      matrix<Type> Y,
			      vector<int> dow0,
			      bool weekday,
			      int j_log_w1,
			      int j_log_w2,
			      int j_log_w3,
			      int j_log_w4,
			      int j_log_w5,
			      int j_log_w6)
{
    vector<Type> log_int_inc(log_cum_inc.size() - slen.size());
    for (int s = 0, k1 = 0, k2 = 0; s < slen.size(); s++) // loop over segments
    {
        for (int i = 0; i < slen(s); i++) // loop over within-segment index
	{
	    log_int_inc(k1+i) = logspace_sub(log_cum_inc(k2+i+1), log_cum_inc(k2+i));
	}
	k1 += slen(s); // increment `log_int_inc` index
	k2 += slen(s) + 1; // increment `log_cum_inc` index
    }

    if (weekday)
    {
        // NB: Below assumes that observations are daily, so we rely on
        //     R back-end to enforce this when `weekday = TRUE`
      
        vector<Type> log_w(7);
	for (int s = 0, k = 0; s < slen.size(); s++) // loop over segments
	{
	    log_w << Type(0),
	             Y(s, j_log_w1), Y(s, j_log_w2), Y(s, j_log_w3),
	             Y(s, j_log_w4), Y(s, j_log_w5), Y(s, j_log_w6);

	    for (int i = 0, d = dow0(s); i < slen(s); i++, d++) // loop over within-segment index
	    {
	        log_int_inc(k+i) += log_w(d % 7);
	    }
	    k += slen(s); // increment `log_int_inc` index
	}
    }
    return log_int_inc;
}

template<class Type>
vector<Type> log_rt(vector<Type> t,
		    vector<Type> log_cum_inc,
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
