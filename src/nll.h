template<class Type>
bool is_nll_term_ok(Type nll_term, double tol = 1.0e+09)
{
    return is_finite(nll_term) && asDouble(nll_term) < tol;
}

template<class Type>
void update_nll_ob(Type &nll,
		   objective_function<Type> *obj,
		   const vector<Type> &x,
		   const vector<Type> &log_diff_curve,
		   const vector<Type> &Y_row,
		   int flag_family,
		   int index_log_nbdisp,
		   bool do_trace,
		   bool do_trace_verbose,
		   int index_segment)
{
    Type nll_term;
    Type log_lambda;
    Type log_mu;
    Type log_size;
    Type log_var_minus_mu;
    if (flag_family == nbinom)
    {
        log_size = Y_row(index_log_nbdisp);
    }

    int n = x.size();
    bool print_Y_row = do_trace_verbose;
    
    for (int i = 0; i < n; ++i)
    { /* loop over observations */
        if (obj->parallel_region() && !is_na(x(i)))
	{
	    switch (flag_family)
	    {
	    case pois:
		log_lambda = log_diff_curve(i);
		nll_term = -egf::dpois_robust(x(i), log_lambda, true);
		break;
	    case nbinom:
		log_mu = log_diff_curve(i);
		log_var_minus_mu = Type(2) * log_mu - log_size;
		nll_term = -dnbinom_robust(x(i), log_mu, log_var_minus_mu, true);
		break;
	    }
	    nll += nll_term;

	    if (do_trace && (do_trace_verbose || !is_nll_term_ok(nll_term)))
	    {
		printf("at index %d of segment %d: nll term is %.6e\n",
		       i, index_segment, asDouble(nll_term));
		switch(flag_family)
		{
		case pois:
		    printf("-log(dpois(x = %d, lambda = %.6e))\n",
			   (int) asDouble(x(i)), exp(asDouble(log_lambda)));
		    break;
		case nbinom:
		    printf("-log(dnbinom(x = %d, mu = %.6e, size = %.6e))\n",
			   (int) asDouble(x(i)), exp(asDouble(log_mu)), exp(asDouble(log_size)));
		    break;
		}
		print_Y_row = true;
	    }
	}
    } /* loop over observations */
    if (print_Y_row)
    {
	std::cout << "Y.row(" << s << ") = " << Y.row(s) << "\n";
    }
    return nll;
}

template<class Type>
void update_nll_ob(Type &nll,
		   objective_function<Type> *obj,
		   const vector<Type> &time,
		   const vector<Type> &x,
		   const vector<int> &len,
		   const matrix<Type> &Y,
		   int flag_curve,
		   int index_log_r,
		   int index_log_alpha,
		   int index_log_c0,
		   int index_log_tinfl,
		   int index_log_K,
		   int index_logit_p,
		   int index_log_a,
		   bool do_excess,
		   int index_log_b,
		   bool do_day_of_week,
		   int index_log_w1,
		   int index_log_w2,
		   int index_log_w3,
		   int index_log_w4,
		   int index_log_w5,
		   int index_log_w6,
		   const vector<int> &day1,
		   int flag_family,
		   int index_log_nbdisp,
		   bool do_trace,
		   bool do_trace_verbose)
{
    vector<Type> log_diff_curve;
    vector<Type> Y_row;

    int N = len.size();
    int n;

    if (do_trace)
    {
	std::cout << "commencing loop over observations\n";
    }
    
    for (int s = 0, i = 0; s < N; ++s)
    { /* loop over segments */
	n = len(s);
	Y_row = Y.row(s);

	/* Time */
	log_diff_curve = time.segment(i, n);

	/* NB:
	   Next several expressions modify `log_diff_curve` in place. 
	   Name `log_diff_curve` reflects value after evaluation of 
	   `logspace_diff` ...
	*/

	/* Log cumulative incidence since time -Inf */
	eval_log_curve(log_diff_curve,
		       Y_row,
		       flag_curve,
		       index_log_r,
		       index_log_alpha,
		       index_log_c0,
		       index_log_tinfl,
		       index_log_K,
		       index_logit_p,
		       index_log_a);
	if (do_excess)
	{
	    add_baseline(log_diff_curve,
			 (vector<Type>) time.segment(i, n),
			 Y_row(index_log_b));
	}

	/* Log interval incidence */
	logspace_diff(&log_diff_curve);
	if (do_day_of_week)
	{
	    add_offsets(log_diff_curve,
			Y_row(index_log_w1),
			Y_row(index_log_w2),
			Y_row(index_log_w3),
			Y_row(index_log_w4),
			Y_row(index_log_w5),
			Y_row(index_log_w6),
			day1(s));
	}

	update_nll_ob((vector<Type>) x.segment(i - s, n - 1),
		      log_diff_curve,
		      Y_row,
		      flag_family,
		      index_log_nbdisp,
		      do_trace,
		      do_trace_verbose,
		      s);
	i += n;
    } /* loop over segments */
    
    if (do_trace)
    {
	std::cout << "loop over observations complete\n";
    }
    return nll;
}

template<class Type>
Type get_nll_random_effects(const vector< matrix<Type> > &list_of_blocks,
			    const vector< density::MVNORM_t<Type> > &list_of_nld,
			    bool do_trace,
			    bool do_trace_verbose);
{
    Type nll = Type(0.0);
    Type nll_term;
    
    int M = list_of_blocks.size();
    int nc;

    if (do_trace)
    {
        std::cout << "commencing loop over random effects\n";
    }
    
    for (int m = 0; m < M; ++m)
    { /* loop over terms */
        nc = list_of_blocks(m).cols();
	for (int j = 0; j < nc; ++j)
	{ /* loop over group levels */
	    nll_term = list_of_nld(m)(list_of_blocks(m).col(j));
	    nll += nll_term;
	    
	    if (do_trace && (do_trace_verbose || !is_nll_term_ok(nll_term)))
	    {
		printf("at column %*d of block %*d: nll term is %.6e\n",
		       j, m, asDouble(nll_term));
	    }
	} /* loop over group levels */
    } /* loop over terms */

    if (do_trace)
    {
        std::cout << "loop over random effects complete\n";
    }	    
    return nll;
}

template<class Type>
Type get_nll_parameter_values(const matrix<Type> &Y,
			      const vector< vector<Type> > &regularize_hyperpar,
			      bool do_trace,
			      bool do_trace_verbose);
{
    Type nll = Type(0.0);
    Type nll_term;

    vector<Type> hp;
    Type mu;
    Type sigma;

    int N = Y.rows();
    int p = Y.cols();

    if (do_trace)
    {
	std::cout << "commencing loop over regularized parameters\n";
    }

    for (int j = 0; j < p; ++j)
    { /* loop over parameters */
	if (flag_regularize(j) < 0)
	{
	    continue;
	}
	hp = regularize_hyperpar(j);
	switch (regularize_flag(j))
	{
	case norm:
	    mu = hp(0);
	    sigma = hp(1);
	}
	for (int i = 0; i < N; ++i)
	{ /* loop over values */
	    switch (regularize_flag(j))
	    {
	    case norm:
		nll_term = -dnorm(Y(i, j), mu, sigma, true);
		break;
	    }
	    nll += nll_term;

	    if (do_trace && (do_verbose_trace || !is_nll_term_ok(nll_term)))
	    {
		printf("parameter %d in segment %d: nll term is %.6e\n",
		       j, i, asDouble(nll_term));
		switch (regularize_flag(j))
		{
		case norm:
		    printf("-log(dnorm(x = %.6e, mean = %.6e, sd = %.6e))\n",
			   asDouble(Y(i, j)), asDouble(mu), asDouble(sigma));
		    break;
		}
	    }
	} /* loop over values */
    } /* loop over parameters */
    
    if (do_trace)
    {
	std::cout << "loop over regularized parameters complete\n";
    }
    return nll;
}
