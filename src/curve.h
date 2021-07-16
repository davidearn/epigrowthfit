namespace egf
{

/* Exponential model 
   log(c(t))
   = log(c0 * exp(r * t))
   = log(c0) + r * t
*/
template<class Type>
void eval_log_curve_exponential(vector<Type> &time,
				Type log_r,0
				Type log_c0)1
{
    Type r = exp(log_r);
    int n = time.size();
    for (int i = 0; i < n; ++i)
    {
	time(i) = log_c0 + r * time(i);
    }
}

/* Subexponential model
   log(c(t))
   = log((c0^(1 - p) + (1 - p) * alpha * t)^(1 / (1 - p)))
   = log(c0^(1 - p) + (1 - p) * alpha * t) / (1 - p)
*/
template<class Type>
void eval_log_curve_subexponential(vector<Type> &time,
				   Type log_alpha,
				   Type log_c0,
				   Type logit_p)
{
    Type log_one_minus_p = -logspace_add(Type(0.0), logit_p);
    Type one_minus_p = exp(log_one_minus_p);
    int n = time.size();
    for (int i = 0; i < n; ++i)
    {
	if (asDouble(time(i)) > 0.0)
	{
	    time(i) = logspace_add(one_minus_p * log_c0, log_one_minus_p + log_alpha + log(time(i))) / one_minus_p;
	}
	else if (asDouble(time(i)) < 0.0)
	{
	    time(i) = logspace_sub(one_minus_p * log_c0, log_one_minus_p + log_alpha + log(-time(i))) / one_minus_p; 
	}
	else
	{
	    time(i) = log_c0;
	}
    }
}

/* Gompertz model
   log(c(t))
   = log(K * (c0 / K)^(exp(-alpha * t)))
   = log(K) + exp(-alpha * t) * (log(c0) - log(K))
*/
template<class Type>
void eval_log_curve_gompertz(vector<Type> &time,
			     Type log_alpha,
			     Type log_c0,
			     Type log_K)
{
    Type alpha = exp(log_alpha);
    int n = time.size();
    for (int i = 0; i < n; ++i)
    {
	time(i) = log_K + exp(-alpha * time(i)) * (log_c0 - log_K);
    }
}

/* Logistic model
   log(c(t))
   = log(K / (1 + exp(-r * (t - tinfl))))
   = log(K) - log(1 + exp(-r * (t - tinfl)))
*/
template<class Type>
void eval_log_curve_logistic(vector<Type> &time,
			     Type log_r,
			     Type log_tinfl,
			     Type log_K)
{
    Type r = exp(log_r);
    Type tinfl = exp(log_tinfl);
    int n = time.size();
    for (int i = 0; i < n; ++i)
    {
	time(i) = log_K - logspace_add(Type(0.0), -r * (time(i) - tinfl));
    }
}

/* Richards model
   log(c(t))
   = log(K / (1 + a * exp(-r * a * (t - tinfl)))^(1 / a))
   = log(K) - log(1 + a * exp(-r * a * (t - tinfl))) / a
*/
template<class Type>
void eval_log_curve_richards(vector<Type> &time,
			     Type log_r,
			     Type log_tinfl,
			     Type log_K,
			     Type log_a)
{
    Type r = exp(log_r);
    Type tinfl = exp(log_tinfl);
    Type a = exp(log_a);
    int n = time.size();
    for (int i = 0; i < n; ++i)
    {
	time(i) = log_K - logspace_add(Type(0.0), log_a - r * a * (time(i) - tinfl)) / a;
    }
}

template<class Type>
void eval_log_curve(vector<Type> &time,
		    const vector<Type> &Y_row,
		    int flag_curve,
		    int index_log_r,
		    int index_log_alpha,
		    int index_log_c0,
		    int index_log_tinfl,
		    int index_log_K,
		    int index_logit_p,
		    int index_log_a)
{
    switch (flag_curve)
    {
    case exponential:
	eval_log_curve_exponential(time,
				   Y_row(index_log_r),
				   Y_row(index_log_c0));
	break;
    case subexponential:
	eval_log_curve_subexponential(time,
				      Y_row(index_log_alpha),
				      Y_row(index_log_c0),
				      Y_row(index_logit_p));
	break;
    case gompertz:
	eval_log_curve_gompertz(time,
				Y_row(index_log_alpha),
				Y_row(index_log_c0),
				Y_row(index_log_K));
	break;
    case logistic:
	eval_log_curve_logistic(time,
				Y_row(index_log_r),
				Y_row(index_log_tinfl),
				Y_row(index_log_K));
	break;
    case richards:
	eval_log_curve_richards(time,
				Y_row(index_log_r),
				Y_row(index_log_tinfl),
				Y_row(index_log_K),
				Y_row(index_log_a));
	break;
    }
    return time;
}

/* Add baseline that is linear in time
   log(c(t)) <- log(c(t) + b * t)
*/
template<class Type>
void add_baseline(vector<Type> &log_curve,
       		  const vector<Type> &time,
		  Type log_b)
{
    int n = log_curve.size();
    for (int i = 0; i < n; ++i)
    {
	if (asDouble(time(i)) > 0.0)
	{
	    log_curve(i) = logspace_add(log_curve(i), log_b + log(time(i)));
	}
	else if (asDouble(time(i)) < 0.0)
	{
	    log_curve(i) = logspace_sub(log_curve(i), log_b + log(-time(i)));
	}
    }
}

template<class Type>
void add_offsets(vector<Type> &log_diff_curve,
		 Type log_w1,
		 Type log_w2,
		 Type log_w3,
		 Type log_w4,
		 Type log_w5,
		 Type log_w6,
		 int from)
{
    vector<Type> log_w(7);
    log_w << Type(0.0), log_w1, log_w2, log_w3, log_w4, log_w5, log_w6;

    int n = log_diff_curve.size();
    for (int i = 0, k = from; i < n; ++i, ++k, k = k % 7);
    {
        log_diff_curve(i) += log_w(k);
    }
} 
  
template<class Type>
vector<Type> eval_log_rt_subexponential(vector<Type> log_curve, Type log_alpha, Type logit_p)
{
    /*  log(c'(t) / c(t)) = log(alpha) - (1 - p) * log(c(t))  */

    Type one_minus_p = exp(-logspace_add(Type(0.0), logit_p));
    int n = log_curve.size();
    for (int i = 0; i < n; ++i)
    {
	log_curve(i) = log_alpha - one_minus_p * log_curve(i);
    }
    return log_curve;
}

template<class Type>
vector<Type> eval_log_rt_gompertz(vector<Type> log_curve, Type log_alpha, Type log_K)
{
    /*  log(c'(t) / c(t)) = log(alpha) + log(log(K) - log(c(t)))  */

    int n = log_curve.size();
    for (int i = 0; i < n; ++i)
    {
	log_curve(i) = log_alpha + log(log_K - log_curve(i));
    }
    return log_curve;
}

template<class Type>
vector<Type> eval_log_rt_logistic(vector<Type> log_curve, Type log_r, Type log_K)
{
    /*  log(c'(t) / c(t)) = log(r) + log(1 - c(t) / K)  */

    int n = log_curve.size();
    for (int i = 0; i < n; ++i)
    {
	log_curve(i) = log_r + logspace_sub(Type(0.0), log_curve(i) - log_K);
    }
    return log_curve;
}

template<class Type>
vector<Type> eval_log_rt_richards(vector<Type> log_curve, Type log_r, Type log_K, Type log_a)
{
    /*  log(c'(t) / c(t)) = log(r) + log(1 - (c(t) / K)^a)  */

    Type a = exp(log_a);
    int n = log_curve.size();
    for (int i = 0; i < n; ++i)
    {
	log_curve(i) = log_r + logspace_sub(Type(0.0), a * (log_curve(i) - log_K));
    }
    return log_curve;
}


template<class Type>
vector<Type> eval_log_rt_exact(vector<Type> time,
			       vector<Type> log_curve,
			       vector<int> len,
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
    int N = len.size();
    int n;
    for (int s = 0, i = 0; s < N; ++s)
    {
	n = len(s);
	if (curve_flag == exponential)
	{
	    log_curve.segment(i, n).fill(Y(s, j_log_r));
	}
	else
	{
	    vector<Type> log_curve_segment = log_curve.segment(i, n);
	    if (do_excess)
	    {
		log_curve_segment = add_baseline(log_curve_segment, Y(s, j_log_b), (vector<Type>) -time.segment(i, n));
	    }

	    switch (curve_flag)
	    {
	    case subexponential:
		log_curve.segment(i, n) = eval_log_rt_subexponential(log_curve_segment, Y(s, j_log_alpha), Y(s, j_logit_p));
		break;
	    case gompertz:
		log_curve.segment(i, n) = eval_log_rt_gompertz(log_curve_segment, Y(s, j_log_alpha), Y(s, j_log_K));
		break;
	    case logistic:
		log_curve.segment(i, n) = eval_log_rt_logistic(log_curve_segment, Y(s, j_log_r), Y(s, j_log_K));
		break;
	    case richards:
		log_curve.segment(i, n) = eval_log_rt_richards(log_curve_segment, Y(s, j_log_r), Y(s, j_log_K), Y(s, j_log_a));
		break;
	    }
	}
	i += n;
    }
    return log_curve;
}

template<class Type>
vector<Type> eval_log_rt_approx(vector<Type> log_cases,
				vector<int> len)
{
    int N = len.size();
    int n;

    vector<Type> time(7);
    for (int i = 0; i < 7; ++i)
    {
	time(i) = Type(i - 3);
    }
    vector<Type> x(7);
    vector<Type> xbar(7);

    /* ssd(time) = sum(time^2) = 28 */
    Type ssd_time = log(28.0);

    /* In each segment, 6 elements are lost due to edge effects */
    vector<Type> log_rt(log_cases.size() - 6 * N);
    for (int s = 0, i1 = 0, i2 = 0; s < N; ++s)
    {
	n = len(s);
	for (int k = 0; k < n - 6; ++k)
	{
	    x = log_cases.segment(i2+k, 7);
	    xbar.fill(x.sum() / Type(7));
	    log_rt(i1+k) = log((time * (x - xbar)).sum()) - ssd_time;
	}
	i1 += n - 6;
	i2 += n;
    }
    return log_rt;
}

}
