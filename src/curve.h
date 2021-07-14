template<class Type>
vector<Type> eval_log_exponential(vector<Type> time, Type log_r, Type log_c0)
{
    /* log(c(t))
       = log(c0 * exp(r * t))
       = log(c0) + r * t
    */

    int n = time.size();
    Type r = exp(log_r);

    vector<Type> res(n);
    for (int i = 0; i < n; i++)
    {
        res(i) = log_c0 + r * time(i);
    }
    return res;
}

template<class Type>
vector<Type> eval_log_subexponential(vector<Type> time, Type log_alpha, Type log_c0, Type logit_p)
{
    /* log(c(t))
       = log((c0^(1 - p) + (1 - p) * alpha * t)^(1 / (1 - p)))
       = log(c0^(1 - p) + (1 - p) * alpha * t) / (1 - p)
    */

    int n = time.size();
    Type log_one_minus_p = -logspace_add(Type(0), logit_p);
    Type one_minus_p = exp(log_one_minus_p);

    vector<Type> res(n);
    for (int i = 0; i < n; i++)
    {
        res(i) = logspace_add(one_minus_p * log_c0, log_one_minus_p + log_alpha + log(time(i))) / one_minus_p;
    }
    return res;
}

template<class Type>
vector<Type> eval_log_gompertz(vector<Type> time, Type log_alpha, Type log_c0, Type log_K)
{
    /* log(c(t))
       = log(K * (c0 / K)^(exp(-alpha * t)))
       = log(K) + exp(-alpha * t) * (log(c0) - log(K))
    */

    int n = time.size();
    Type alpha = exp(log_alpha);

    vector<Type> res(n);
    for (int i = 0; i < n; i++)
    {
        res(i) = log_K + exp(-alpha * time(i)) * (log_c0 - log_K);
    }
    return res;
}

template<class Type>
vector<Type> eval_log_logistic(vector<Type> time, Type log_r, Type log_tinfl, Type log_K)
{
    /* log(c(t))
       = log(K / (1 + exp(-r * (t - tinfl))))
       = log(K) - log(1 + exp(-r * (t - tinfl)))
    */

    int n = time.size();
    Type r = exp(log_r);
    Type tinfl = exp(log_tinfl);

    vector<Type> res(n);
    for (int i = 0; i < n; i++)
    {
        res(i) = log_K - logspace_add(Type(0), -r * (time(i) - tinfl));
    }
    return res;
}

template<class Type>
vector<Type> eval_log_richards(vector<Type> time, Type log_r, Type log_tinfl, Type log_K, Type log_a)
{
    /* log(c(t))
       = log(K / (1 + a * exp(-r * a * (t - tinfl)))^(1 / a))
       = log(K) - log(1 + a * exp(-r * a * (t - tinfl))) / a
    */
  
    int n = time.size();
    Type r = exp(log_r);
    Type tinfl = exp(log_tinfl);
    Type a = exp(log_a);

    vector<Type> res(n);
    for (int i = 0; i < n; i++)
    {
        res(i) = log_K - logspace_add(Type(0), log_a - r * a * (time(i) - tinfl)) / a;
    }
    return res;
}

template<class Type>
vector<Type> eval_log_curve(vector<Type> time,
			    vector<int> time_seg_len,
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
    int N = time_seg_len.size();
    int n;
    
    vector<Type> log_curve(time.size());
    for (int s = 0, i = 0; s < N; s++)
    {
        n = time_seg_len(s);
	vector<Type> time0 = time.segment(i, n);
	switch (curve_flag)
	{
	case exponential:
	    log_curve.segment(i, n) = eval_log_exponential(time0, Y(s, j_log_r), Y(s, j_log_c0));
	    break;
	case subexponential:
	    log_curve.segment(i, n) = eval_log_subexponential(time0, Y(s, j_log_alpha), Y(s, j_log_c0), Y(s, j_logit_p));
	    break;
	case gompertz:
	    log_curve.segment(i, n) = eval_log_gompertz(time0, Y(s, j_log_alpha), Y(s, j_log_c0), Y(s, j_log_K));
	    break;
	case logistic:
	    log_curve.segment(i, n) = eval_log_logistic(time0, Y(s, j_log_r), Y(s, j_log_tinfl), Y(s, j_log_K));
	    break;
	case richards:
	    log_curve.segment(i, n) = eval_log_richards(time0, Y(s, j_log_r), Y(s, j_log_tinfl), Y(s, j_log_K), Y(s, j_log_a));
	    break;
	}
	if (do_excess)
	{
	    for (int k = 0; k < n; k++)
	    {
	        if (asDouble(time(i+k)) > 0.0)
		{
		    log_curve(i+k) = logspace_add(log_curve(i+k), Y(s, j_log_b) + log(time(i+k)));
		}
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
    for (int s = 0, i = 0; s < N; s++)
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
    vector<Type> time(7);
    for (int i = 0; i < 7; i++)
    {
	time(i) = Type(i - 3);
    }
    vector<Type> x(7);
    vector<Type> xbar(7);
    Type log_28 = log(28.0);
    
    /* In each segment, 6 elements are lost due to edge effects */
    vector<Type> log_rt(log_cases.size() - 6 * N);
    for (int s = 0, i1 = 0, i2 = 0; s < N; s++)
    {
        n = log_cases_seg_len(s);
        for (int k = 0; k < n - 6; k++)
	{
	    x = log_cases.segment(i2+k, 7);
	    xbar.fill(x.sum() / Type(7.0));
	    log_rt(i1+k) = log((time * (x - xbar)).sum()) - log_28;
	}
	i1 += n - 6;
	i2 += n;
    }
    return log_rt;
}		

template<class Type>
vector<Type> eval_log_rt_exact(vector<Type> time,
			       vector<Type> log_curve,
			       vector<int> time_seg_len,
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
    int N = time_seg_len.size();
    int n;
    Type one_minus_p;
    
    vector<Type> log_rt(time.size());
    for (int s = 0, i = 0; s < N; s++)
    {
        n = time_seg_len(s);
	for (int k = 0; k < n; k++)
	{
	    if (do_excess && asDouble(time(i+k)) > 0.0)
	    {   
		log_curve(i+k) = logspace_sub(log_curve(i+k), Y(s, j_log_b) + log(time(i+k)));
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
