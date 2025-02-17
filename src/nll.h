template<class Type>
bool is_nll_term_ok(Type nll_term, double tol = 1.0e+9)
{
	return egf::is_finite(nll_term) && asDouble(nll_term) < tol;
}

template<class Type>
Type nll_obs(objective_function<Type> *obj,
             const vector<Type> &time,
             const vector<int> &time_seg_len,
             vector<Type> &x,
             const matrix<Type> &Y,
             const egf::indices_t<Type> &indices,
             const egf::flags_t<Type> &flags,
             const vector<int> &day1)
{
	Type res = Type(0.0);

	int N = time_seg_len.size();
	int n;
	vector<Type> Y_row;
	vector<Type> predict;

	for (int s = 0, i = 0; s < N; ++s)
	{ /* loop over segments */
		n = time_seg_len(s);
		Y_row = Y.row(s);
		/* t */
		predict = time.segment(i, n);
		/* t <- log(c(t)) */
		egf::eval_log_curve(predict, Y_row, indices, flags.curve);
		if (flags.do_excess)
		{
			/* log(c(t)) <- log(b * t + c(t)) */
			egf::logspace_add_baseline(predict,
			                           (vector<Type>) time.segment(i, n),
			                           Y_row(indices.log_b));
		}
		/* log(c(t)) <- log(diff(c(t))) */
		egf::logspace_diff(predict);
		if (flags.do_day_of_week)
		{
			/* log(diff(c(t))) <- log(diff(c(t)) * w(t[-n], t[-1])) */
			egf::logspace_add_offsets(predict,
			                          Y_row(indices.log_w1),
			                          Y_row(indices.log_w2),
			                          Y_row(indices.log_w3),
			                          Y_row(indices.log_w4),
			                          Y_row(indices.log_w5),
			                          Y_row(indices.log_w6),
			                          day1(s));
		}
		res += nll_obs(obj, x, predict, Y_row, indices, flags, i - s, s);
		i += n;
	} /* loop over segments */

	return res;
}


template<class Type>
Type nll_obs(objective_function<Type> *obj,
             vector<Type> &x,
             const vector<Type> &log_diff_curve,
             const vector<Type> &Y_row,
             const egf::indices_t<Type> &indices,
             const egf::flags_t<Type> &flags,
             int ix,
             int sx)
{
	Type res = Type(0.0);
	Type nll_term = Type(0.0);

	int n = log_diff_curve.size();
	bool print_Y_row = flags.do_trace_verbose;

	Type log_lambda = Type(0.0);
	Type log_mu = Type(0.0);
	Type log_size = Type(0.0);
	Type log_var_minus_mu = Type(0.0);
	if (flags.family == nbinom)
	{
		log_size = Y_row(indices.log_disp);
	}

	for (int i = ix, k = 0; k < n; ++i, ++k)
	{ /* loop over observations */
		if (obj->parallel_region() && !egf::is_NA_real_(x(i)))
		{
			switch (flags.family)
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
				switch(flags.family)
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
			switch(flags.family)
			{
			case pois:
				nll_term = -egf::dpois_robust(x(i), log_lambda, true);
				break;
			case nbinom:
				nll_term = -dnbinom_robust(x(i), log_mu, log_var_minus_mu, true);
				break;
			}
			res += nll_term;

			if (flags.do_trace &&
			    (flags.do_trace_verbose || !is_nll_term_ok(nll_term)))
			{
				Rprintf("at index %d of segment %d: nll term is %.6e\n",
				        k, sx, asDouble(nll_term));
				switch(flags.family)
				{
				case pois:
					Rprintf("-dpois(x = %d, lambda = %.6e, give_log = true)\n",
					        (int) asDouble(x(i)), exp(asDouble(log_lambda)));
					break;
				case nbinom:
					Rprintf("-dnbinom(x = %d, mu = %.6e, size = %.6e, give_log = true)\n",
					        (int) asDouble(x(i)), exp(asDouble(log_mu)), exp(asDouble(log_size)));
					break;
				}
				print_Y_row = true;
			}
		}
	} /* loop over observations */

	if (print_Y_row)
	{
		Rcout << "Y.row(" << sx << ") =\n" << Y_row << "\n";
	}

	return res;
}


template<class Type>
Type nll_ran(objective_function<Type> *obj,
             const vector< matrix<Type> > &list_of_blocks,
             vector< density::MVNORM_t<Type> > list_of_nld,
             const egf::flags_t<Type> &flags)
{
	Type res = Type(0.0);
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
				res += nll_term;

				if (flags.do_trace &&
				    (flags.do_trace_verbose || !is_nll_term_ok(nll_term)))
				{
					Rprintf("at column %d of block %d: nll term is %.6e\n",
					        j, m, asDouble(nll_term));
				}
			}
		} /* loop over group levels */
	} /* loop over terms */

	return res;
}

template<class Type>
Type nll_top(objective_function<Type> *obj,
             const matrix<Type> &Y,
             const vector< vector<Type> > &hyperparameters_top,
             const egf::flags_t<Type> &flags)
{
	Type res = Type(0.0);
	Type nll_term = Type(0.0);

	int nr = Y.rows();
	int nc = Y.cols();

	vector<Type> hp;
	Type mu = Type(0.0);
	Type sigma = Type(0.0);

	for (int j = 0; j < nc; ++j)
	{ /* loop over parameters */
		if (flags.regularize_top(j) >= 0)
		{
			hp = hyperparameters_top(j);
			switch (flags.regularize_top(j))
			{
			case norm:
				mu = hp(0);
				sigma = hp(1);
				break;
			}

			for (int i = 0; i < nr; ++i)
			{ /* loop over values */
				if (obj->parallel_region())
				{
					switch (flags.regularize_top(j))
					{
					case norm:
						nll_term = -dnorm(Y(i, j), mu, sigma, true);
						break;
					}
					res += nll_term;

					if (flags.do_trace &&
					    (flags.do_trace_verbose || !is_nll_term_ok(nll_term)))
					{
						Rprintf("parameter %d in segment %d: nll term is %.6e\n",
						        j, i, asDouble(nll_term));
						switch (flags.regularize_top(j))
						{
						case norm:
							Rprintf("-dnorm(x = %.6e, mean = %.6e, sd = %.6e, give_log = true)\n",
							        asDouble(Y(i, j)), asDouble(mu), asDouble(sigma));
							break;
						}
					}
				}
			} /* loop over values */
		}
	} /* loop over parameters */

	return res;
}

template<class Type>
Type nll_bot(objective_function<Type> *obj,
             const vector<Type> &beta,
             const vector<Type> &theta,
             const vector< vector<Type> > &list_of_sd,
             const vector< vector<Type> > &list_of_chol,
             const vector< vector<Type> > &hyperparameters_bottom,
             const egf::flags_t<Type> &flags)
{
	Type res = Type(0.0);
	Type nll_term = Type(0.0);
	vector<Type> hp;

	Type mu = Type(0.0);
	Type sigma = Type(0.0);

	int i = 0;
	for (int j = 0; j < beta.size(); ++i, ++j)
	{ /* loop over 'beta' elements */
		if (obj->parallel_region() && flags.regularize_bottom(i) >= 0)
		{
			hp = hyperparameters_bottom(i);
			switch (flags.regularize_bottom(i))
			{
			case norm:
				mu = hp(0);
				sigma = hp(1);
				nll_term = -dnorm(beta(j), mu, sigma, true);
				break;
			}
			res += nll_term;

			if (flags.do_trace &&
			    (flags.do_trace_verbose || !is_nll_term_ok(nll_term)))
			{
				Rprintf("beta element %d: nll term is %.6e\n",
				        j, asDouble(nll_term));
				switch (flags.regularize_bottom(i))
				{
				case norm:
					Rprintf("-dnorm(x = %.6e, mean = %.6e, sd = %.6e, give_log = true)\n",
					        asDouble(beta(j)), asDouble(mu), asDouble(sigma));
					break;
				}
			}
		}
	} /* loop over 'beta' elements */

	if (!flags.do_random_effects)
	{
		return res;
	}

	for (int j = 0; j < theta.size(); ++i, ++j)
	{ /* loop over 'theta' elements */
		if (obj->parallel_region() && flags.regularize_bottom(i) >= 0)
		{
			hp = hyperparameters_bottom(i);
			switch (flags.regularize_bottom(i))
			{
			case norm:
				mu = hp(0);
				sigma = hp(1);
				nll_term = -dnorm(theta(j), mu, sigma, true);
				break;
			}
			res += nll_term;

			if (flags.do_trace &&
			    (flags.do_trace_verbose || !is_nll_term_ok(nll_term)))
			{
				Rprintf("theta element %d: nll term is %.6e\n",
				        j, asDouble(nll_term));
				switch (flags.regularize_bottom(i))
				{
				case norm:
					Rprintf("-dnorm(x = %.6e, mean = %.6e, sd = %.6e, give_log = true)\n",
					        asDouble(theta(j)), asDouble(mu), asDouble(sigma));
					break;
				}
			}
		}
	} /* loop over 'theta' elements */

	vector<Type> x;
	Type eta = Type(0.0);
	Type df = Type(0.0);
	vector<Type> scale;

	for (int j = 0; j < list_of_chol.size(); ++i, ++j)
	{ /* loop over correlation/covariance matrices */
		if (obj->parallel_region() && flags.regularize_bottom(i) >= 0)
		{
			hp = hyperparameters_bottom(i);
			switch (flags.regularize_bottom(i))
			{
			case lkj:
				x = list_of_chol(j);
				eta = hp(0);
				nll_term = -egf::dlkj(x, eta, true);
				break;
			default:
				x.resize(list_of_sd(j).size() + list_of_chol(j).size());
				x << log(list_of_sd(j)),list_of_chol(j);
				df = hp(0);
				scale = hp.tail(x.size());
				switch (flags.regularize_bottom(i))
				{
				case wishart:
					nll_term = -egf::dwishart(x, df, scale, true);
					break;
				case invwishart:
					nll_term = -egf::dinvwishart(x, df, scale, true);
					break;
				}
				break;
			}
			res += nll_term;

			if (flags.do_trace &&
			    (flags.do_trace_verbose || !is_nll_term_ok(nll_term)))
			{
				Rprintf("correlation/covariance matrix %d: nll term is %.6e\n",
				        j, asDouble(nll_term));
				switch (flags.regularize_bottom(i))
				{
				case lkj:
					Rcout << "-dlkj(x = " << x << ", eta = " << eta << ", give_log = true)\n";
					break;
				case wishart:
					Rcout << "-dwishart(x = " << x << ", df = " << df << ", scale = " << scale << ", give_log = true)\n";
					break;
				case invwishart:
					Rcout << "-dinvwishart(x = " << x << ", df = " << df << ", scale = " << scale << ", give_log = true)\n";
					break;
				}
			}
		}
	} /* loop over correlation/covariance matrices */

	return res;
}
