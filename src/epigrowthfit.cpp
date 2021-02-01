#define TMB_LIB_INIT R_init_epigrowthfit
#include <TMB.hpp>

enum curve
{
    exponential,
    subexponential,
    gompertz,
    logistic,
    richards
};

enum distr
{
    pois,
    nbinom
};

template<class Type>
Type eval_exponential(Type t, Type log_r, Type log_c0)
{
    // log(c(t))
    // = log(c0 * exp(r * t))
    // = log(c0) + r * t
    return log_c0 + exp(log_r) * t;
}

template<class Type>
Type eval_subexponential(Type t, Type log_alpha, Type log_c0, Type logit_p)
{
    // log(c(t))
    // = log((c0^(1 - p) + (1 - p) * alpha * t)^(1 / (1 - p)))
    // = log(c0^(1 - p) + (1 - p) * alpha * t) / (1 - p)
    Type one_minus_p = Type(1) / (Type(1) + exp(logit_p));
    return logspace_add(one_minus_p * log_c0, log(one_minus_p) + log_alpha + log(t)) / one_minus_p;
}

template<class Type>
Type eval_gompertz(Type t, Type log_alpha, Type log_c0, Type log_K)
{
    // log(c(t))
    // = log(K * (c0 / K)^(exp(-alpha * t)))
    // = log(K) + exp(-alpha * t) * (log(c0) - log(K))
    return log_K + exp(-exp(log_alpha) * t) * (log_c0 - log_K);
}

template<class Type>
Type eval_logistic(Type t, Type log_r, Type log_tinfl, Type log_K)
{
    // log(c(t))
    // = log(K / (1 + exp(-r * (t - tinfl))))
    // = log(K) - log(1 + exp(-r * (t - tinfl)))
    return log_K - logspace_add(Type(0), -exp(log_r) * (t - exp(log_tinfl)));
}

template<class Type>
Type eval_richards(Type t, Type log_r, Type log_tinfl, Type log_K, Type log_a)
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

template<class Type>
Type dpois_robust(Type x, Type log_lambda, int give_log = 0)
{
    Type log_dpois = x * log_lambda - exp(log_lambda) - lfactorial(x);
    return ( give_log ? log_dpois : exp(log_dpois) );
}

template<class Type>
Type rnbinom_robust(Type log_mu, Type log_disp)
{
    Type log_a = log_disp - logspace_add(log_mu, log_disp);
    return rnbinom(exp(log_disp), exp(log_a));
}

// https://github.com/kaskr/adcomp/issues/59
template<class Type>
bool isNA_real_(Type x)
{
    return R_IsNA(asDouble(x));
}

template<class Type>
matrix<Type> prune_dupl_rows(matrix<Type> Y, vector<int> wl)
{
    matrix<Type> Y_short(wl.size(), Y.cols());
    for (int i = 0, k = 0; k < wl.size(); k++)
    {
        Y_short.row(k) = Y.row(i);
	i += wl(k);
    }
    return Y_short;
}

// https://kaskr.github.io/adcomp/Introduction.html
template<class Type>
Type objective_function<Type>::operator() ()
{
    // Set up ==================================================================

    // Flags
    DATA_INTEGER(curve_flag);    // curve     (enum curve)
    DATA_INTEGER(distr_flag);    // distr     (enum distr)
    DATA_INTEGER(excess_flag);   // excess    (1=yes,      0=no)
    DATA_INTEGER(sparse_X_flag); // X format  (1=sparse,   0=dense)
    DATA_INTEGER(predict_flag);  // predict   (1=predict,  0=fit)
    bool excess   = (excess_flag   == 1);
    bool sparse_X = (sparse_X_flag == 1);
    bool predict  = (predict_flag  == 1);
    
    // Data
    // time series
    DATA_VECTOR(t); // length=N
    DATA_VECTOR(x); // length=N
    // window lengths
    DATA_IVECTOR(wl); // length=nts
    int nw = wl.size();
    // factor nlevels
    DATA_IVECTOR(fnl); // length=nFE=np
    DATA_IVECTOR(rnl); // length=nRE
    // design matrices
    DATA_SPARSE_MATRIX(Xs); // nrow=N,  ncol=sum(nlevels(FE i))
    DATA_MATRIX(Xd);
    DATA_SPARSE_MATRIX(Z);  // nrow=N, ncol=sum(nlevels(RE i))
    // random effect indicators
    DATA_IMATRIX(rid); // nrow=nRE,  ncol=np,  {0,1}
    bool anyRE = (rid.rows() > 0);
    int np = rid.cols();
    vector<int> rnp = rid.rowwise().sum();
    // parameter indices
    DATA_INTEGER(j_log_r);
    DATA_INTEGER(j_log_alpha);
    DATA_INTEGER(j_log_c0);
    DATA_INTEGER(j_log_tinfl);
    DATA_INTEGER(j_log_K);
    DATA_INTEGER(j_logit_p);
    DATA_INTEGER(j_log_a);
    DATA_INTEGER(j_log_b);
    DATA_INTEGER(j_log_nbdisp);
    
    
    // Parameters
    // fixed effects
    PARAMETER_VECTOR(beta); // length=sum(nlevels(FE i))
    // log sd random effects
    PARAMETER_VECTOR(log_sd_b); // length=sum(np(RE i))
    // random effects
    PARAMETER_VECTOR(b); // length=sum(nlevels(RE i)*np(RE i)) 

    
    // Compute parameter values at every time point ============================
    // NOTE: parameters are constant within but not across fitting windows

    // (N x np) matrix of parameter values
    matrix<Type> Y(t.size(), np);
    
    // Fixed effects component:
    // matrix multiply
    // X * beta2
    // = (N x sum(nlevels(FE i))) * (sum(nlevels(FE i)) x np)
    // = (N x np)
    // need to construct block matrix `beta2` from vector `beta`

    Eigen::SparseMatrix<Type> beta2(fnl.sum(), np);
    beta2.reserve(fnl);
    
    for (int j = 0, i = 0; j < np; j++) // loop over parameters
    {
	for (int l = 0; l < fnl(j); l++) // loop over levels
	{
	    beta2.insert(i + l, j) = beta(i + l);
	}
	i += fnl(j);
    }
    
    if (sparse_X)
    {
        Y = Xs * beta2;
    }
    else
    {
        Y = Xd * beta2;
    }
    matrix<Type> Y_sim = Y; // preserve fixed effects component for SIMULATE
    
    // Random effects component:
    // matrix multiply
    // Z * b2
    // = (N x sum(nlevels(RE i))) * (sum(nlevels(RE i)) x np)
    // = (N x np)
    // need to construct block matrix `b2` from vector `b`
    // will require some setting up...
    
    matrix<Type> b2(rnl.sum(), np);
    
    // A list of (nlevels(RE i) x np(RE i)) matrices gathering related elements
    // of vector `b`: rows of each matrix follow a zero-mean, unit-variance
    // multivariate normal distribution with a common covariance matrix
    vector< matrix<Type> > re_list(rid.rows());

    // A list of s.d. vectors corresponding elementwise to `re_list`
    vector< vector<Type> > sd_list(rid.rows());

    // A list of correlation matrices corresponding elementwise to `re_list`
    vector< matrix<Type> > cor_list(rid.rows());
	
    if (anyRE)
    {

        // Form random effect matrices -----------------------------------------

	for (int r = 0, k = 0; r < re_list.size(); r++) // loop over RE
	{
	    matrix<Type> re_list_el(rnl(r), rnp(r));
	    for (int j = 0; j < rnp(r); j++) // loop over parameters
	    {
		re_list_el.col(j) = b.segment(k, rnl(r));
		k += rnl(r);
	    }
	    re_list(r) = re_list_el;
	}
	REPORT(re_list);


	// Form scale vectors --------------------------------------------------

	for (int r = 0, k = 0; r < sd_list.size(); r++)
	{
	    vector<Type> sd_list_el = exp(log_sd_b.segment(k, rnp(r)));
	    sd_list(r) = sd_list_el;
	    k += rnp(r);
	}
	REPORT(sd_list);

	
	// Form correlation matrices -------------------------------------------

	// Correlation matrix in univariate case
	matrix<Type> cor_list_el(1, 1);
	cor_list_el(0, 0) = Type(1);

	for (int r = 0; r < cor_list.size(); r++) // loop over RE
	{
	    // UNSTRUCTURED_CORR() does not tolerate `ltri.size() == 0`
	    if (rnp(r) == 1)
	    {
		cor_list(r) = cor_list_el;
		continue;
	    }
	    else
	    {
	        vector<Type> ltri(rnp(r) * (rnp(r) - 1) / 2);
		cor_list(r) = density::UNSTRUCTURED_CORR(ltri).cov();
	    }
	}
	REPORT(cor_list);


	// Construct `b2` ------------------------------------------------------

        b2.fill(Type(0));
	for (int r = 0, i = 0; r < rid.rows(); r++) // loop over RE
	{
	    for (int j = 0, k = 0; j < np; j++) // loop over parameters
	    {
		if (rid(r, j) == 1)
		{
		    b2.block(i, j, rnl(r), 1) = sd_list(r)(k) * re_list(r).col(k);
		    k++;
		}
	    }
	    i += rnl(r);
	}
	Y += Z * b2;
    }


    // Report parameter values in each fitting window ==========================
    // NOTE: parameters are constant within but not across fitting windows

    matrix<Type> Y_short = prune_dupl_rows(Y, wl);
    vector<Type> Y_short_as_vector = Y_short.vec();
    ADREPORT(Y_short_as_vector);
    

    // Compute likelihood ======================================================

    // Log cumulative incidence
    vector<Type> log_cum_inc = eval_log_cum_inc(t, Y, curve_flag, excess,
						j_log_r, j_log_alpha, j_log_c0, j_log_tinfl, j_log_K, j_logit_p, j_log_a, j_log_b);
    // Log interval incidence
    vector<Type> log_int_inc = eval_log_int_inc(log_cum_inc, wl);

    // Negative log likelihood
    Type nll = Type(0);
    Type log_var_minus_mu;
    
    for (int i1 = 0, i2 = 1, k = 0; k < wl.size(); k++) // loop over time series
    {
        for (int j = 0; j < wl(k) - 1; j++) // loop over time points
	{
	    if (!isNA_real_(x(i2+j)))
	    {
	        switch (distr_flag)
		{
		case pois:
		    nll -= dpois_robust(x(i2+j), log_int_inc(i1+j), true);
		    // usage: dpois_robust(x, log_lambda, give_log)
		    break;
		case nbinom:
		    log_var_minus_mu = Type(2) * log_int_inc(i1+j) - Y(i2+j, j_log_nbdisp);
		    nll -= dnbinom_robust(x(i2+j), log_int_inc(i1+j), log_var_minus_mu, true);
		    // usage: dnbinom_robust(x, log_mu, log_var_minus_mu, give_log)
		    break;
		}
	    }
	}
	i1 += wl(k) - 1;
	i2 += wl(k);
    }
    if (anyRE)
    {
    	for (int r = 0; r < re_list.size(); r++) // loop over RE
    	{
    	    for (int l = 0; l < rnl(r); l++) // loop over levels
    	    {
    	        nll += density::MVNORM(cor_list(r))(re_list(r).row(l));
    	    }
    	}
    }


    // Simulate incidence ======================================================

    SIMULATE
    {
        if (anyRE)
	{
	    matrix<Type> b2_sim(rnl.sum(), np);
	    b2_sim.fill(Type(0));
	    
	    for (int r = 0, i = 0; r < rid.rows(); r++) // loop over RE
	    {
		vector<Type> z(sd_list(r).size());
		for (int l = 0; l < rnl(r); l++) // loop over levels
		{
		    z = sd_list(r) * density::MVNORM(cor_list(r)).simulate();
		    for (int j = 0, k = 0; j < np; j++) // loop over parameters
		    {
			if (rid(r, j) == 1)
			{
			    b2_sim(i + l, j) = z(k);
			    k++;
			}
		    }
		}
		i += rnl(r);
	    }
	    Y_sim += Z * b2_sim;
	}
	// matrix<Type> Y_sim_short = prune_dupl_rows(Y_sim, wl);
	// REPORT(Y_sim_short);

        log_cum_inc = eval_log_cum_inc(t, Y_sim, curve_flag, excess,
				       j_log_r, j_log_alpha, j_log_c0, j_log_tinfl, j_log_K, j_logit_p, j_log_a, j_log_b);
	log_int_inc = eval_log_int_inc(log_cum_inc, wl);

	for (int i1 = 0, i2 = 1, k = 0; k < wl.size(); k++) // loop over time series
	{
	    for (int j = 0; j < wl(k) - 1; j++) // loop over time points
	    {
	        switch(distr_flag)
		{
		case pois:
		    x(i2+j) = rpois(exp(log_int_inc(i1+j)));
		    // usage: rpois(mu)
		    break;
		case nbinom:
		    x(i2+j) = rnbinom_robust(log_int_inc(i1+j), Y_sim(i2+j, j_log_nbdisp));
		    // usage: rnbinom_robust(log_mu, log_disp)
		    break;
		}
	    }
	    i1 += wl(k) - 1;
	    i2 += wl(k);
	}
	REPORT(x);
    }
    

    // Predict incidence =======================================================
    // NOTE: enforcing on R side that there is only one time series here
    
    if (predict)
    {
        DATA_VECTOR(t_new); // length=N_new
	DATA_SPARSE_MATRIX(Xs_new); // nrow=N_new,  ncol=sum(nlevels(FE i))
    	DATA_MATRIX(Xd_new); 
	DATA_SPARSE_MATRIX(Z_new); // nrow=N_new,  ncol=sum(nlevels(RE i))
	DATA_IVECTOR(predict_lci_lii_lrt_flag); // predict (log_cum_inc, log_int_inc, log_rt)  (1=do,  0=don't)
	DATA_INTEGER(se_flag); // report  (1=ADREPORT,  0=REPORT)

	matrix<Type> Y_new(t_new.size(), np);
	if (sparse_X)
	{
	    Y_new = Xs_new * beta2;
	}
	else
	{
	    Y_new = Xd_new * beta2;
	}
        if (anyRE)
    	{
	    Y_new += Z_new * b2;
    	}

	vector<Type> log_cum_inc_new = eval_log_cum_inc(t_new, Y_new, curve_flag, excess,
							j_log_r, j_log_alpha, j_log_c0, j_log_tinfl, j_log_K, j_logit_p, j_log_a, j_log_b);
	
	if (predict_lci_lii_lrt_flag(0) == 1)
	{
	    vector<int> wl_new(1);
	    wl_new(0) = t_new.size();
	    vector<Type> log_int_inc_new = eval_log_int_inc(log_cum_inc_new, wl_new);
	    if (se_flag == 1)
	    {
	        ADREPORT(log_int_inc_new);
	    }
	    else
	    {
	        REPORT(log_int_inc_new);
	    }
	}
	if (predict_lci_lii_lrt_flag(1) == 1)
	{
	    if (se_flag == 1)
	    {
	        ADREPORT(log_cum_inc_new);
	    }
	    else
	    {
	        REPORT(log_cum_inc_new);
	    }
	}
	if (predict_lci_lii_lrt_flag(2) == 1)
	{
	    vector<Type> log_rt_new = eval_log_rt(log_cum_inc_new, Y_new, curve_flag,
						  j_log_r, j_log_alpha, j_log_K, j_logit_p, j_log_a);
	    if (se_flag == 1)
	    {
	        ADREPORT(log_rt_new);
	    }
	    else
	    {
	        REPORT(log_rt_new);
	    }
	}
    }

    return nll;
}
