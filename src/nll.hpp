template<class Type>
bool is_nll_term_ok(Type nll_term, double tol = 1.0e+09)
{
    return egf::is_finite(nll_term) && asDouble(nll_term) < tol;
}

template<class Type>
void add_nll_ob(Type &nll,
		objective_function<Type> *obj,
		const vector<Type> &time,
		const vector<int> &time_seg_len,
		vector<Type> &x,
		const matrix<Type> &Y,
		const egf::indices_t<Type> &indices,
		const egf::flags_t<Type> &flags,
		const vector<int> &day1)
{
    int N = time_seg_len.size();
    int n;
    vector<Type> Y_row(Y.cols());
    vector<Type> log_diff_curve;
    
    for (int s = 0, i = 0; s < N; ++s)
    { /* loop over segments */
	n = time_seg_len(s);
	Y_row = Y.row(s);

	/* Time */
	log_diff_curve = time.segment(i, n);

	/* NB:
	   Next several expressions modify `log_diff_curve` in place. 
	   Name `log_diff_curve` reflects value after later evaluation 
	   of `logspace_diff` ...
	*/

	/* Log cumulative incidence since time -Inf */
	egf::eval_log_curve(log_diff_curve, Y_row, indices, flags.flag_curve);
	if (flags.do_excess)
	{
	    egf::add_baseline(log_diff_curve,
			      (vector<Type>) time.segment(i, n),
			      Y_row(indices.index_log_b));
	}

	/* Log interval incidence */
	egf::logspace_diff(&log_diff_curve);
	if (flags.do_day_of_week)
	{
	    egf::add_offsets(log_diff_curve,
			     Y_row(indices.index_log_w1),
			     Y_row(indices.index_log_w2),
			     Y_row(indices.index_log_w3),
			     Y_row(indices.index_log_w4),
			     Y_row(indices.index_log_w5),
			     Y_row(indices.index_log_w6),
			     day1(s));
	}

	add_nll_ob(nll, obj, x, log_diff_curve, Y_row, indices, flags,
		   i - s, s);
	i += n;
    } /* loop over segments */
}


template<class Type>
void add_nll_ob(Type &nll,
		objective_function<Type> *obj,
		vector<Type> &x,
		const vector<Type> &log_diff_curve,
		const vector<Type> &Y_row,
		const egf::indices_t<Type> &indices,
		const egf::flags_t<Type> &flags,
		int ix,
		int sx)
{
    Type nll_term;
    int n = log_diff_curve.size();
    bool print_Y_row = flags.do_trace_verbose;

    Type log_lambda;
    Type log_mu;
    Type log_size;
    Type log_var_minus_mu;
    if (flags.flag_family == nbinom)
    {
        log_size = Y_row(indices.index_log_nbdisp);
    }

    for (int i = ix, k = 0; k < n; ++i, ++k)
    { /* loop over observations */
        if (obj->parallel_region() && !egf::is_na(x(i)))
	{
	    switch (flags.flag_family)
	    {
	    case pois:
		log_lambda = log_diff_curve(k);
		break;
	    case nbinom:
		log_mu = log_diff_curve(k);
		log_var_minus_mu = Type(2) * log_mu - log_size;
		break;
	    }
	    if (flags.do_simulate)
	    {
	        switch(flags.flag_family)
		{
		case pois:
		    x(i) = rpois(exp(log_lambda));
		    break;
		case nbinom:
		    x(i) = egf::rnbinom_robust(log_mu, log_size);
		    break;
		}
		continue;
	    }
	    else
	    {
	        switch(flags.flag_family)
		{
		case pois:
		    nll_term = -egf::dpois_robust(x(i), log_lambda, true);
		    break;
		case nbinom:
		    nll_term = -dnbinom_robust(x(i), log_mu, log_var_minus_mu, true);
		    break;
		}
	    }
	    nll += nll_term;

	    if (flags.do_trace && (flags.do_trace_verbose || !is_nll_term_ok(nll_term)))
	    {
		printf("at index %d of segment %d: nll term is %.6e\n",
		       k, sx, asDouble(nll_term));
		switch(flags.flag_family)
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
	std::cout << "Y.row(" << s << ") = " << Y_row << "\n";
    }
}


template<class Type>
void add_nll_re(Type &nll,
		objective_function<Type> *obj,
		const vector< matrix<Type> > &list_of_blocks,
		const vector< density::MVNORM_t<Type> > &list_of_nld,
		const egf::flags_t<Type> &flags)
{
    Type nll_term;
    int M = list_of_blocks.size();
    int nc;
    
    for (int m = 0; m < M; ++m)
    { /* loop over terms */
        nc = list_of_blocks(m).cols();
        for (int j = 0; j < nc; ++j)
	{ /* loop over group levels */
	    if (obj->parallel_region())
	    {
	        nll_term = list_of_nld(m)(list_of_blocks(m).col(j));
		nll += nll_term;

		if (flags.do_trace && (flags.do_trace_verbose || !is_nll_term_ok(nll_term)))
		{
		    printf("at column %*d of block %*d: nll term is %.6e\n",
			   j, m, asDouble(nll_term));
		}
	    }
	} /* loop over group levels */
    } /* loop over terms */
}

template<class Type>
void add_nll_pv(Type &nll,
		objective_function<Type> *obj,
		const matrix<Type> &Y,
		const vector< vector<Type> > &hyperparameters,
		const egf::flags_t<Type> &flags)
{
    Type nll_term;
    int nr = Y.rows();
    int nc = Y.cols();

    vector<Type> hp;
    Type mu;
    Type sigma;

    for (int j = 0; j < nc; ++j)
    { /* loop over parameters */
	if (flags.flag_regularize(j) >= 0)
	{
	    hp = hyperparameters(j);
	    switch (flags.flag_regularize(j))
	    {
	    case norm:
		mu = hp(0);
		sigma = hp(1);
	    }
	    
	    for (int i = 0; i < nr; ++i)
	    { /* loop over values */
	        if (obj->parallel_region())
		{
		    switch (flag.flag_regularize(j))
		    {
		    case norm:
			nll_term = -dnorm(Y(i, j), mu, sigma, true);
			break;
		    }
		    nll += nll_term;

		    if (flags.do_trace && (flags.do_trace_verbose || !is_nll_term_ok(nll_term)))
		    {
			printf("parameter %d in segment %d: nll term is %.6e\n",
			       j, i, asDouble(nll_term));
			switch (flags.flag_regularize(j))
			{
			case norm:
			    printf("-log(dnorm(x = %.6e, mean = %.6e, sd = %.6e))\n",
				   asDouble(Y(i, j)), asDouble(mu), asDouble(sigma));
			    break;
			}
		    }
		}
	    } /* loop over values */
	}
    } /* loop over parameters */
}
