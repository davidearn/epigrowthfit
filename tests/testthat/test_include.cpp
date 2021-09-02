#include <TMB.hpp>
#include "enums.h"
#include "../../src/structs.h"
#include "../../src/utils.h"
#include "../../src/distributions.h"
#include "../../src/enums.h"
#include "../../src/curve.h"

template<class Type>
Type objective_function<Type>::operator() ()
{
    DATA_INTEGER(flag_test);
  
    switch (flag_test)
    {
    /* structs.h */
    case list_of_vectors_t:
    {
        DATA_STRUCT(x, egf::list_of_vectors_t);
	vector< vector<Type> > res = x;
	REPORT(res);
	break;
    }
    /* utils.h */
    case is_NA_real_:
    {
        DATA_VECTOR(x);
	int n = x.size();
	vector<int> res(n);
	for (int i = 0; i < n; ++i)
	{
	    res(i) = egf::is_NA_real_(x(i));
	}
	REPORT(res);
        break;
    }
    case is_finite:
    {
        DATA_VECTOR(x);
	int n = x.size();
	vector<int> res(n);
	for (int i = 0; i < n; ++i)
	{
	    res(i) = egf::is_finite(x(i));
	}
	REPORT(res);
        break;
    }
    case logspace_diff:
    {
        DATA_VECTOR(log_x);
	egf::logspace_diff(log_x);
	vector<Type> res = log_x;
	REPORT(res);
        break;
    }
    /* distributions.h */
    case mvlgamma:
    {
        DATA_VECTOR(x);
	DATA_INTEGER(p);
	int n = x.size();
	vector<Type> res(n);
	for (int i = 0; i < n; ++i)
	{
	    res(i) = egf::mvlgamma(x(i), p);
	}
	REPORT(res);
        break;
    }
    case dlkj:
    {
        DATA_VECTOR(x);
	DATA_SCALAR(eta);
	DATA_INTEGER(give_log);
	Type res = egf::dlkj(x, eta, give_log);
	REPORT(res);
        break;
    }
    case dwishart:
    {
        DATA_VECTOR(x);
	DATA_SCALAR(df);
	DATA_VECTOR(scale)
	DATA_INTEGER(give_log);
	Type res = egf::dwishart(x, df, scale, give_log);
	REPORT(res);
        break;
    }
    case dinvwishart:
    {
        DATA_VECTOR(x);
	DATA_SCALAR(df);
	DATA_VECTOR(scale)
	DATA_INTEGER(give_log);
	Type res = egf::dinvwishart(x, df, scale, give_log);
	REPORT(res);
        break;
    }
    case dpois_robust:
    {
        DATA_VECTOR(x);
        DATA_VECTOR(log_lambda);
        DATA_INTEGER(give_log);
	int n = x.size();
	vector<Type> res(n);
	for (int i = 0; i < n; ++i)
	{
	    res(i) = egf::dpois_robust(x(i), log_lambda(i), give_log);
	}
	REPORT(res);
        break;
    }
    case rnbinom_robust:
    {
        DATA_SCALAR(log_mu);
        DATA_SCALAR(log_size);
	DATA_INTEGER(n);
	vector<Type> res(n);
	for (int i = 0; i < n; ++i)
	{
	    res(i) = egf::rnbinom_robust(log_mu, log_size);
	}
	REPORT(res);
        break;
    }
    /* curve.h */
    case eval_log_curve_exponential:
    {
        DATA_VECTOR(time);
	DATA_SCALAR(log_r);
	DATA_SCALAR(log_c0);
	egf::eval_log_curve_exponential(time, log_r, log_c0);
	vector<Type> res = time;
	REPORT(res);
        break;
    }
    case eval_log_curve_subexponential:
    {
        DATA_VECTOR(time);
	DATA_SCALAR(log_alpha);
	DATA_SCALAR(log_c0);
	DATA_SCALAR(logit_p);
	egf::eval_log_curve_subexponential(time, log_alpha, log_c0, logit_p);
	vector<Type> res = time;
	REPORT(res);
        break;
    }
    case eval_log_curve_gompertz:
    {
        DATA_VECTOR(time);
	DATA_SCALAR(log_alpha);
	DATA_SCALAR(log_tinfl);
	DATA_SCALAR(log_K);
	egf::eval_log_curve_gompertz(time, log_alpha, log_tinfl, log_K);
	vector<Type> res = time;
	REPORT(res);
        break;
    }
    case eval_log_curve_logistic:
    {
        DATA_VECTOR(time);
	DATA_SCALAR(log_r);
	DATA_SCALAR(log_tinfl);
	DATA_SCALAR(log_K);
	egf::eval_log_curve_logistic(time, log_r, log_tinfl, log_K);
	vector<Type> res = time;
	REPORT(res);
        break;
    }
    case eval_log_curve_richards:
    {
        DATA_VECTOR(time);
	DATA_SCALAR(log_r);
	DATA_SCALAR(log_tinfl);
	DATA_SCALAR(log_K);
	DATA_SCALAR(log_a);
	egf::eval_log_curve_richards(time, log_r, log_tinfl, log_K, log_a);
	vector<Type> res = time;
	REPORT(res);
        break;
    }
    case logspace_add_baseline:
    {
        DATA_VECTOR(log_curve);
	DATA_VECTOR(time);
	DATA_SCALAR(log_b);
	egf::logspace_add_baseline(log_curve, time, log_b);
	vector<Type> res = log_curve;
	REPORT(res);
        break;
    }
    case logspace_add_offsets:
    {
        DATA_VECTOR(log_diff_curve);
        DATA_VECTOR(log_w);
        DATA_INTEGER(from);
	egf::logspace_add_offsets(log_diff_curve,
				  log_w(1), log_w(2), log_w(3),
				  log_w(4), log_w(5), log_w(6),
				  from);
	vector<Type> res = log_diff_curve;
	REPORT(res);
        break;
    }
    case eval_log_rt_subexponential:
    {
        DATA_VECTOR(log_curve);
        DATA_SCALAR(log_alpha);
	DATA_SCALAR(logit_p);
	egf::eval_log_rt_subexponential(log_curve, log_alpha, logit_p);
	vector<Type> res = log_curve;
	REPORT(res);
        break;
    }
    case eval_log_rt_gompertz:
    {
        DATA_VECTOR(log_curve);
        DATA_SCALAR(log_alpha);
	DATA_SCALAR(log_K);
	egf::eval_log_rt_gompertz(log_curve, log_alpha, log_K);
	vector<Type> res = log_curve;
	REPORT(res);
        break;
    }
    case eval_log_rt_logistic:
    {
        DATA_VECTOR(log_curve);
        DATA_SCALAR(log_r);
	DATA_SCALAR(log_K);
	egf::eval_log_rt_logistic(log_curve, log_r, log_K);
	vector<Type> res = log_curve;
	REPORT(res);
        break;
    }
    case eval_log_rt_richards:
    {
        DATA_VECTOR(log_curve);
        DATA_SCALAR(log_r);
	DATA_SCALAR(log_K);
	DATA_SCALAR(log_a);
	egf::eval_log_rt_richards(log_curve, log_r, log_K, log_a);
	vector<Type> res = log_curve;
	REPORT(res);
        break;
    }
    }
    
    return Type(0.0);
}
