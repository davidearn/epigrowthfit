template<class Type>
Type eval_log_exponential(Type t, Type log_r, Type log_c0)
{
    /* log(c(t))
       = log(c0 * exp(r * t))
       = log(c0) + r * t
    */
    return log_c0 + exp(log_r) * t;
}

template<class Type>
Type eval_log_subexponential(Type t, Type log_alpha, Type log_c0, Type logit_p)
{
    /* log(c(t))
       = log((c0^(1 - p) + (1 - p) * alpha * t)^(1 / (1 - p)))
       = log(c0^(1 - p) + (1 - p) * alpha * t) / (1 - p)
    */
    Type log_one_minus_p = -logspace_add(Type(0), logit_p);
    Type one_minus_p = exp(log_one_minus_p);
    return logspace_add(one_minus_p * log_c0, log_one_minus_p + log_alpha + log(t)) / one_minus_p;
}

template<class Type>
Type eval_log_gompertz(Type t, Type log_alpha, Type log_c0, Type log_K)
{
    /* log(c(t))
       = log(K * (c0 / K)^(exp(-alpha * t)))
       = log(K) + exp(-alpha * t) * (log(c0) - log(K))
    */
    return log_K + exp(-exp(log_alpha) * t) * (log_c0 - log_K);
}

template<class Type>
Type eval_log_logistic(Type t, Type log_r, Type log_tinfl, Type log_K)
{
    /* log(c(t))
       = log(K / (1 + exp(-r * (t - tinfl))))
       = log(K) - log(1 + exp(-r * (t - tinfl)))
    */
    return log_K - logspace_add(Type(0), -exp(log_r) * (t - exp(log_tinfl)));
}

template<class Type>
Type eval_log_richards(Type t, Type log_r, Type log_tinfl, Type log_K, Type log_a)
{
    /* log(c(t))
       = log(K / (1 + a * exp(-r * a * (t - tinfl)))^(1 / a))
       = log(K) - log(1 + a * exp(-r * a * (t - tinfl))) / a
    */
    Type a = exp(log_a);
    return log_K - logspace_add(Type(0), log_a - exp(log_r) * a * (t - exp(log_tinfl))) / a;
}

template<class Type>
vector<Type> eval_log_curve(vector<Type> t,
			    vector<int> t_seg_len,
			    int curve_flag,
			    bool do_excess,
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
    int N = t_seg_len.size();
    int n;
    
    vector<Type> log_curve(t.size());
    for (int s = 0, i = 0; s < N; s++) /* loop over segments */
    {
        n = t_seg_len(s);
        for (int k = 0; k < n; k++) /* loop over within-segment index */
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
	    if (do_excess && t(i+k) > 0)
	    {
	        log_curve(i+k) = logspace_add(Y(s, j_log_b) + log(t(i+k)), log_curve(i+k));
	    }
	}
	i += n;
    }
    return log_curve;
}

template<class Type>
vector<Type> get_day_of_week_offsets(vector<int> log_cases_seg_len,
				     vector<int> day1,
				     matrix<Type> Y,
				     int j_log_w1,
				     int j_log_w2,
				     int j_log_w3,
				     int j_log_w4,
				     int j_log_w5,
				     int j_log_w6)
{
    int N = log_cases_seg_len.size();
    int n;
    vector<Type> log_w(7);

    vector<Type> offsets(log_cases_seg_len.sum());
    for (int s = 0, i = 0; s < N; s++) /* loop over segments */
    {
        n = log_cases_seg_len(s);
	log_w << Type(0),
		 Y(s, j_log_w1), Y(s, j_log_w2), Y(s, j_log_w3),
		 Y(s, j_log_w4), Y(s, j_log_w5), Y(s, j_log_w6);
	offsets.segment(i, n) = rep_len_from(log_w, n, day1(s));
	i += n;
    }
    return offsets;
}

template<class Type>
vector<Type> eval_log_rt_approx(vector<Type> log_cases,
				vector<int> log_cases_seg_len)
{
    int N = log_cases_seg_len.size();
    int n;
    vector<Type> t(7);
    for (int i = 0; i < 7; i++)
    {
	t(i) = Type(i - 3);
    }
    vector<Type> x(7);
    vector<Type> xbar(7);
    Type log_28 = log(28.0);
    
    /* In each segment, 6 elements are lost due to edge effects */
    vector<Type> log_rt(log_cases.size() - 6 * N);
    for (int s = 0, i1 = 0, i2 = 0; s < N; s++) /* loop over segments */
    {
        n = log_cases_seg_len(s);
        for (int k = 0; k < n - 6; k++) /* loop over within-segment index */
	{
	    x = log_cases.segment(i2+k, 7);
	    xbar.fill(x.sum() / Type(7));
	    log_rt(i1+k) = log((t * (x - xbar)).sum()) - log_28;
	}
	i1 += n - 6;
	i2 += n;
    }
    return log_rt;
}		

template<class Type>
vector<Type> eval_log_rt_exact(vector<Type> t,
			       vector<Type> log_curve,
			       vector<int> t_seg_len,
			       int curve_flag,
			       bool do_excess,
			       matrix<Type> Y,
			       int j_log_r,
			       int j_log_alpha,
			       int j_log_K,
			       int j_logit_p,
			       int j_log_a,
			       int j_log_b)
{
    int N = t_seg_len.size();
    int n;
    Type one_minus_p;
    
    vector<Type> log_rt(t.size());
    for (int s = 0, i = 0; s < N; s++)
    {
        n = t_seg_len(s);
	for (int k = 0; k < n; k++)
	{
	    if (do_excess)
	    {
		log_curve(i+k) = logspace_sub(log_curve(i+k), Y(s, j_log_b) + log(t(i+k)));
	    }
	    switch (curve_flag)
	    {
	    case exponential:
	        /*  log(c'(t) / c(t)) = log(r)  */
		log_rt(i+k) = Y(s, j_log_r);
		break;
	    case subexponential:
	        /*  log(c'(t) / c(t)) = log(alpha) - (1 - p) * log(c(t))  */
	        one_minus_p = exp(-logspace_add(Type(0), Y(s, j_logit_p)));
		log_rt(i+k) = Y(s, j_log_alpha) - one_minus_p * log_curve(i+k);
		break;
	    case gompertz:
	        /*  log(c'(t) / c(t)) = log(alpha) + log(log(K) - log(c(t)))  */
		log_rt(i+k) = Y(s, j_log_alpha) + log(Y(s, j_log_K) - log_curve(i+k));
		break;
	    case logistic:
	        /*  log(c'(t) / c(t)) = log(r) + log(1 - c(t) / K)  */
		log_rt(i+k) = Y(s, j_log_r) + logspace_sub(Type(0), log_curve(i+k) - Y(s, j_log_K));
		break;
	    case richards:
	        /*  log(c'(t) / c(t)) = log(r) + log(1 - (c(t) / K)^a)  */
		log_rt(i+k) = Y(s, j_log_r) + logspace_sub(Type(0), exp(Y(s, j_log_a)) * (log_curve(i+k) - Y(s, j_log_K)));
		break;
	    }
	}
	i += n;
    }
    return log_rt;
}
