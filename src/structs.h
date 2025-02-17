namespace egf
{

int get_list_integer(SEXP list, const char *name)
{
	SEXP el = getListElement(list, name, &isNumericScalar);
	return (int) REAL(el)[0];
}

template<class Type>
struct indices_t
{
	int log_r;
	int log_alpha;
	int log_c0;
	int log_tinfl;
	int log_K;
	int logit_p;
	int log_a;
	int log_b;
	int log_disp;
	int log_w1;
	int log_w2;
	int log_w3;
	int log_w4;
	int log_w5;
	int log_w6;

	indices_t(SEXP x)
	{
		log_r     = get_list_integer(x, "log_r");
		log_alpha = get_list_integer(x, "log_alpha");
		log_c0    = get_list_integer(x, "log_c0");
		log_tinfl = get_list_integer(x, "log_tinfl");
		log_K     = get_list_integer(x, "log_K");
		logit_p   = get_list_integer(x, "logit_p");
		log_a     = get_list_integer(x, "log_a");
		log_b     = get_list_integer(x, "log_b");
		log_disp  = get_list_integer(x, "log_disp");
		log_w1    = get_list_integer(x, "log_w1");
		log_w2    = get_list_integer(x, "log_w2");
		log_w3    = get_list_integer(x, "log_w3");
		log_w4    = get_list_integer(x, "log_w4");
		log_w5    = get_list_integer(x, "log_w5");
		log_w6    = get_list_integer(x, "log_w6");
	}
};

template<class Type>
struct flags_t
{
	/* Model of cumulative incidence (enum.h) */
	int curve;
	/* Baseline term in model of cumulative incidence (1=yes, 0=no) */
	int excess;
	/* Model of observation error (enum.h) */
	int family;
	/* Day of week effects (1=yes, 0=no) */
	int day_of_week;
	/* Priors on top and bottom level parameters (enum.h, -1=no prior) */
	vector<int> regularize_top;
	vector<int> regularize_bottom;
	/* Trace
	   0=nothing
	   1=degenerate nll terms
	   2=all nll terms
	*/
	int trace;
	/* X format (1=sparse, 0=dense) */
	int sparse_X;
	/* predict (1=yes, 0=no) */
	int predict;

	bool do_excess;
	bool do_day_of_week;
	bool do_regularize_top;
	bool do_regularize_bottom;
	bool do_trace;
	bool do_trace_verbose;
	bool do_sparse_X;
	bool do_random_effects;
	bool do_predict;
	bool do_simulate;

	flags_t(SEXP x)
	{
		curve       = get_list_integer(x, "curve");
		excess      = get_list_integer(x, "excess");
		family      = get_list_integer(x, "family");
		day_of_week = get_list_integer(x, "day_of_week");
		trace       = get_list_integer(x, "trace");
		sparse_X    = get_list_integer(x, "sparse_X");
		predict     = get_list_integer(x, "predict");

		do_excess        = (excess      == 1);
		do_day_of_week   = (day_of_week == 1);
		do_trace         = (trace       >= 1);
		do_trace_verbose = (trace       >= 2);
		do_sparse_X      = (sparse_X    == 1);
		do_predict       = (predict     == 1);

		regularize_top    = asVector<int>(getListElement(x, "regularize_top"   , &Rf_isNumeric));
		regularize_bottom = asVector<int>(getListElement(x, "regularize_bottom", &Rf_isNumeric));

		do_regularize_top = false;
		for (int i = 0; !do_regularize_top    && i < regularize_top   .size(); ++i)
		{
			do_regularize_top = regularize_top(i) >= 0;
		}

		do_regularize_bottom = false;
		for (int i = 0; !do_regularize_bottom && i < regularize_bottom.size(); ++i)
		{
			do_regularize_bottom = regularize_bottom(i) >= 0;
		}
	}
};

template<class Type>
struct list_of_vectors_t : vector< vector<Type> >
{
	list_of_vectors_t(SEXP x)
	{
		(*this).resize(LENGTH(x));
		for (int i = 0; i < LENGTH(x); ++i)
		{
			SEXP v = VECTOR_ELT(x, i);
			(*this)(i) = asVector<Type>(v);
		}
	}
};

} /* namespace egf */
