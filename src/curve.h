namespace egf
{

/* Exponential model
   log(c(t))
   = log(c0 * exp(r * t))
   = log(c0) + r * t
*/
template<class Type>
void eval_log_curve_exponential(
	vector<Type> &time,
	Type log_r,
	Type log_c0)
{
	Type r = exp(log_r);
	int n = time.size();
	for (int i = 0; i < n; ++i)
	{
		time(i) = log_c0 + r * time(i);
	}
	return;
}

/* Subexponential model
   log(c(t))
   = log(c0 * (1 + (1 - p) * alpha * t / c0^(1 - p))^(1 / (1 - p)))
   = log(c0) + log(1 + (1 - p) * alpha * t / c0^(1 - p)) / (1 - p)
*/
template<class Type>
void eval_log_curve_subexponential(
	vector<Type> &time,
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
			time(i) = log_c0 + logspace_add(Type(0.0), log_one_minus_p + log_alpha + log(time(i)) - one_minus_p * log_c0) / one_minus_p;
		}
		else if (asDouble(time(i)) < 0.0)
		{
			time(i) = log_c0 + logspace_sub(Type(0.0), log_one_minus_p + log_alpha + log(-time(i)) - one_minus_p * log_c0) / one_minus_p;
		}
		else
		{
			time(i) = log_c0;
		}
	}
	return;
}

/* Gompertz model
   log(c(t))
   = log(K * exp(-exp(-alpha * (t - tinfl))))
   = log(K) - exp(-alpha * (t - tinfl))
*/
template<class Type>
void eval_log_curve_gompertz(
	vector<Type> &time,
	Type log_alpha,
	Type log_tinfl,
	Type log_K)
{
	Type alpha = exp(log_alpha);
	Type tinfl = exp(log_tinfl);
	int n = time.size();
	for (int i = 0; i < n; ++i)
	{
		time(i) = log_K - exp(-alpha * (time(i) - tinfl));
	}
	return;
}

/* Logistic model
   log(c(t))
   = log(K / (1 + exp(-r * (t - tinfl))))
   = log(K) - log(1 + exp(-r * (t - tinfl)))
*/
template<class Type>
void eval_log_curve_logistic(
	vector<Type> &time,
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
	return;
}

/* Richards model
   log(c(t))
   = log(K / (1 + a * exp(-a * r * (t - tinfl)))^(1 / a))
   = log(K) - log(1 + a * exp(-a * r * (t - tinfl))) / a
*/
template<class Type>
void eval_log_curve_richards(
	vector<Type> &time,
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
		time(i) = log_K - logspace_add(Type(0.0), log_a - a * r * (time(i) - tinfl)) / a;
	}
	return;
}

template<class Type>
void eval_log_curve(
	vector<Type> &time,
	const vector<Type> &Y_row,
	const indices_t<Type> &indices,
	int flag)
{
	switch (flag)
	{
	case exponential:
		eval_log_curve_exponential(
			time,
			Y_row(indices.log_r),
			Y_row(indices.log_c0));
		break;
	case subexponential:
		eval_log_curve_subexponential(
			time,
			Y_row(indices.log_alpha),
			Y_row(indices.log_c0),
			Y_row(indices.logit_p));
		break;
	case gompertz:
		eval_log_curve_gompertz(
			time,
			Y_row(indices.log_alpha),
			Y_row(indices.log_tinfl),
			Y_row(indices.log_K));
		break;
	case logistic:
		eval_log_curve_logistic(
			time,
			Y_row(indices.log_r),
			Y_row(indices.log_tinfl),
			Y_row(indices.log_K));
		break;
	case richards:
		eval_log_curve_richards(
			time,
			Y_row(indices.log_r),
			Y_row(indices.log_tinfl),
			Y_row(indices.log_K),
			Y_row(indices.log_a));
		break;
	}
	return;
}

/* Add baseline that is linear in time
   log(c(t)) <- log(c(t) + b * t)
*/
template<class Type>
void logspace_add_baseline(
	vector<Type> &log_curve,
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
	return;
}

template<class Type>
void logspace_add_offsets(
	vector<Type> &log_diff_curve,
	Type log_w1,
	Type log_w2,
	Type log_w3,
	Type log_w4,
	Type log_w5,
	Type log_w6,
	int from = 0)
{
	vector<Type> log_w(7);
	log_w << Type(0.0), log_w1, log_w2, log_w3, log_w4, log_w5, log_w6;

	int n = log_diff_curve.size();
	for (int i = 0, k = from; i < n; ++i, ++k, k = k % 7)
	{
		log_diff_curve(i) += log_w(k);
	}
	return;
}

/* Subexponential model
   log(c'(t) / c(t)) = log(alpha) - (1 - p) * log(c(t))
*/
template<class Type>
void eval_log_rt_subexponential(
	vector<Type> &log_curve,
	Type log_alpha,
	Type logit_p)
{
	Type one_minus_p = exp(-logspace_add(Type(0.0), logit_p));
	int n = log_curve.size();
	for (int i = 0; i < n; ++i)
	{
		log_curve(i) = log_alpha - one_minus_p * log_curve(i);
	}
	return;
}

/* Gompertz model
   log(c'(t) / c(t)) = log(alpha) + log(log(K) - log(c(t)))
*/
template<class Type>
void eval_log_rt_gompertz(
	vector<Type> &log_curve,
	Type log_alpha,
	Type log_K)
{
	int n = log_curve.size();
	for (int i = 0; i < n; ++i)
	{
		log_curve(i) = log_alpha + log(log_K - log_curve(i));
	}
	return;
}

/* Logistic model
   log(c'(t) / c(t)) = log(r) + log(1 - c(t) / K)
*/
template<class Type>
void eval_log_rt_logistic(
	vector<Type> &log_curve,
	Type log_r,
	Type log_K)
{
	int n = log_curve.size();
	for (int i = 0; i < n; ++i)
	{
		log_curve(i) = log_r + logspace_sub(Type(0.0), log_curve(i) - log_K);
	}
	return;
}

/* Richards model
   log(c'(t) / c(t)) = log(r) + log(1 - (c(t) / K)^a)
*/
template<class Type>
void eval_log_rt_richards(
	vector<Type> &log_curve,
	Type log_r,
	Type log_K,
	Type log_a)
{
	Type a = exp(log_a);
	int n = log_curve.size();
	for (int i = 0; i < n; ++i)
	{
		log_curve(i) = log_r + logspace_sub(Type(0.0), a * (log_curve(i) - log_K));
	}
	return;
}


template<class Type>
void eval_log_rt_exact(
	vector<Type> &log_curve,
	const vector<Type> &Y_row,
	const indices_t<Type> &indices,
	int flag)
{
	switch (flag)
	{
	case exponential:
		log_curve.fill(Y_row(indices.log_r));
		break;
	case subexponential:
		eval_log_rt_subexponential(
			log_curve,
			Y_row(indices.log_alpha),
			Y_row(indices.logit_p));
		break;
	case gompertz:
		eval_log_rt_gompertz(
			log_curve,
			Y_row(indices.log_alpha),
			Y_row(indices.log_K));
		break;
	case logistic:
		eval_log_rt_logistic(
			log_curve,
			Y_row(indices.log_r),
			Y_row(indices.log_K));
		break;
	case richards:
		eval_log_rt_richards(
			log_curve,
			Y_row(indices.log_r),
			Y_row(indices.log_K),
			Y_row(indices.log_a));
		break;
	}
	return;
}

} /* namespace egf */
