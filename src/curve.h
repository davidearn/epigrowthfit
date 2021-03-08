template<class Type>
Type eval_log_exponential(Type t, Type log_r, Type log_c0)
{
    // log(c(t))
    // = log(c0 * exp(r * t))
    // = log(c0) + r * t
    return log_c0 + exp(log_r) * t;
}

template<class Type>
Type eval_log_subexponential(Type t, Type log_alpha, Type log_c0, Type logit_p)
{
    // log(c(t))
    // = log((c0^(1 - p) + (1 - p) * alpha * t)^(1 / (1 - p)))
    // = log(c0^(1 - p) + (1 - p) * alpha * t) / (1 - p)
    Type one_minus_p = Type(1) / (Type(1) + exp(logit_p));
    return logspace_add(one_minus_p * log_c0, log(one_minus_p) + log_alpha + log(t)) / one_minus_p;
}

template<class Type>
Type eval_log_gompertz(Type t, Type log_alpha, Type log_c0, Type log_K)
{
    // log(c(t))
    // = log(K * (c0 / K)^(exp(-alpha * t)))
    // = log(K) + exp(-alpha * t) * (log(c0) - log(K))
    return log_K + exp(-exp(log_alpha) * t) * (log_c0 - log_K);
}

template<class Type>
Type eval_log_logistic(Type t, Type log_r, Type log_tinfl, Type log_K)
{
    // log(c(t))
    // = log(K / (1 + exp(-r * (t - tinfl))))
    // = log(K) - log(1 + exp(-r * (t - tinfl)))
    return log_K - logspace_add(Type(0), -exp(log_r) * (t - exp(log_tinfl)));
}

template<class Type>
Type eval_log_richards(Type t, Type log_r, Type log_tinfl, Type log_K, Type log_a)
{
    // log(c(t))
    // = log(K / (1 + a * exp(-r * a * (t - tinfl)))^(1 / a))
    // = log(K) - log(1 + a * exp(-r * a * (t - tinfl))) / a
    Type a = exp(log_a);
    return log_K - logspace_add(Type(0), log_a - exp(log_r) * a * (t - exp(log_tinfl))) / a;
}

template<class Type>
vector<Type> eval_log_curve(vector<Type> t,
			    vector<int> t_seg_len,
			    int curve_flag,
			    bool excess,
			    matrix<Type> Y,
			    int j_log_r,
			    int j_log_alpha,
			    int j_log_c0,
			    int j_log_tinfl,
			    int j_log_K,
			    int j_logit_p,
			    int j_log_a,
			    int j_log_b)
{
    vector<Type> log_curve(t.size());
    for (int s = 0, i = 0; s < t_seg_len.size(); s++) // loop over segments
    {
        for (int k = 0; k < t_seg_len(s); k++) // loop over within-segment index
	{
	    switch (curve_flag)
	    {
	    case exponential:
	        log_curve(i+k) = eval_log_exponential(t(i+k), Y(s, j_log_r), Y(s, j_log_c0));
		break;
	    case subexponential:
	        log_curve(i+k) = eval_log_subexponential(t(i+k), Y(s, j_log_alpha), Y(s, j_log_c0), Y(s, j_logit_p));
		break;
	    case gompertz:
	        log_curve(i+k) = eval_log_gompertz(t(i+k), Y(s, j_log_alpha), Y(s, j_log_c0), Y(s, j_log_K));
		break;
	    case logistic:
	        log_curve(i+k) = eval_log_logistic(t(i+k), Y(s, j_log_r), Y(s, j_log_tinfl), Y(s, j_log_K));
		break;
	    case richards:
	        log_curve(i+k) = eval_log_richards(t(i+k), Y(s, j_log_r), Y(s, j_log_tinfl), Y(s, j_log_K), Y(s, j_log_a));
		break;
	    }
	    if (excess)
	    {
	        log_curve(i+k) = logspace_add(Y(s, j_log_b) + log(t(i+k)), log_curve(i+k));
	    }
	}
	i += t_seg_len(s); // increment reference index
    }
    return log_curve;
}

template<class Type>
vector<Type> eval_log_cases(vector<Type> log_curve,
			    vector<int> t_seg_len,
			    bool weekday,
			    // for weekday=true:
			    vector<int> dow,
			    matrix<Type> Y,
			    int j_log_w1,
			    int j_log_w2,
			    int j_log_w3,
			    int j_log_w4,
			    int j_log_w5,
			    int j_log_w6)
{
    vector<Type> log_cases = logspace_diff_n(log_curve, t_seg_len);

    if (weekday)
    {
        // NB: Below assumes that observations are daily. We rely
        //     on R back-end to enforce this when weekday=true.

	vector<Type> log_w(7);
	for (int s = 0, i = 0; s < t_seg_len.size(); s++) // loop over segments
	{
	    log_w << Type(0),
		     Y(s, j_log_w1), Y(s, j_log_w2), Y(s, j_log_w3),
		     Y(s, j_log_w4), Y(s, j_log_w5), Y(s, j_log_w6);

	    for (int k = 0, d = dow(s); k < t_seg_len(s) - 1; k++, d++) // loop over within-segment index
	    {
		log_cases(i+k) += log_w(d % 7);
	    }
	    i += t_seg_len(s) - 1; // increment reference index
	}
    }
    return log_cases;
}

template<class Type>
vector<Type> eval_log_rt(vector<Type> t,
			 vector<Type> log_curve,
			 vector<Type> log_cases,
			 vector<int> t_seg_len,
			 int curve_flag,
			 bool excess,
			 bool weekday,
			 matrix<Type> Y,
			 int j_log_r,
			 int j_log_alpha,
			 int j_log_K,
			 int j_logit_p,
			 int j_log_a,
			 int j_log_b)
{
    if (weekday)
    {
        // Local linear regression on log cases
        vector<Type> x(7);
	vector<Type> y(7);
	vector<Type> ybar(7);
	for (int i = 0; i < 7; i++)
	{
	    x(i) = Type(i - 3);
	}

	// In each segment, we lose 1 element due to differencing
	// and 6 elements due to insufficient data at edges
        vector<Type> log_rt(t.size() - 7 * t_seg_len.size());
	for (int s = 0, i1 = 0, i2 = 0; s < t_seg_len.size(); s++) // loop over segments
	{
	    for (int k = 0; k < t_seg_len(s) - 7; k++) // loop over within-segment index
	    {
	        y = log_cases.segment(i2+k, 7);
		ybar.fill(y.sum() / Type(7));
		log_rt(i1+k) = log((x * (y - ybar)).sum()) - log(Type(28));
	    }
	    i1 += t_seg_len(s) - 7;
	    i2 += t_seg_len(s) - 1;
	}
	return log_rt;
    }

    vector<Type> log_rt(t.size());
    Type one_minus_p;
    for (int s = 0, i = 0; s < t_seg_len.size(); s++)
    {
	for (int k = 0; k < t_seg_len(s); k++)
	{
	    if (excess)
	    {
		log_curve(i+k) = logspace_sub(log_curve(i+k), Y(s, j_log_b) + log(t(i+k)));
	    }
	    switch (curve_flag)
	    {
	    case exponential:
		// log(c'(t) / c(t)) = log(r)
		log_rt(i+k) = Y(s, j_log_r);
		break;
	    case subexponential:
		// log(c'(t) / c(t)) = log(alpha) - (1 - p) * log(c(t))
		one_minus_p = Type(1) / (Type(1) + exp(Y(s, j_logit_p)));
		log_rt(i+k) = Y(s, j_log_alpha) - one_minus_p * log_curve(i+k);
		break;
	    case gompertz:
		// log(c'(t) / c(t)) = log(alpha) + log(log(K) - log(c(t)))
		log_rt(i+k) = Y(s, j_log_alpha) + log(Y(s, j_log_K) - log_curve(i+k));
		break;
	    case logistic:
		// log(c'(t) / c(t)) = log(r) + log(1 - c(t) / K)
		log_rt(i+k) = Y(s, j_log_r) + logspace_sub(Type(0), log_curve(i+k) - Y(s, j_log_K));
		break;
	    case richards:
		// log(c'(t) / c(t)) = log(r) + log(1 - (c(t) / K)^a)
		log_rt(i+k) = Y(s, j_log_r) + logspace_sub(Type(0), exp(Y(s, j_log_a)) * (log_curve(i+k) - Y(s, j_log_K)));
		break;
	    }
	}
    }
    return log_rt;
}
