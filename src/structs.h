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
    int index_log_r;
    int index_log_alpha;
    int index_log_c0;
    int index_log_tinfl;
    int index_log_K;
    int index_logit_p;
    int index_log_a;
    int index_log_b;
    int index_log_nbdisp;
    int index_log_w1;
    int index_log_w2;
    int index_log_w3;
    int index_log_w4;
    int index_log_w5;
    int index_log_w6;
  
    indices_t(SEXP x)
    {
        index_log_r      = get_list_integer(x, "index_log_r");
	index_log_alpha  = get_list_integer(x, "index_log_alpha");
	index_log_c0     = get_list_integer(x, "index_log_c0");
	index_log_tinfl  = get_list_integer(x, "index_log_tinfl");
	index_log_K      = get_list_integer(x, "index_log_K");
	index_logit_p    = get_list_integer(x, "index_logit_p");
	index_log_a      = get_list_integer(x, "index_log_a");
	index_log_b      = get_list_integer(x, "index_log_b");
	index_log_nbdisp = get_list_integer(x, "index_log_nbdisp");
	index_log_w1     = get_list_integer(x, "index_log_w1");
	index_log_w2     = get_list_integer(x, "index_log_w2");
	index_log_w3     = get_list_integer(x, "index_log_w3");
	index_log_w4     = get_list_integer(x, "index_log_w4");
	index_log_w5     = get_list_integer(x, "index_log_w5");
	index_log_w6     = get_list_integer(x, "index_log_w6");
    }
};

template<class Type>
struct flags_t
{
    /* Model of cumulative incidence (enums.hpp) */
    int flag_curve;
    /* Baseline term in model of cumulative incidence (1=yes, 0=no) */
    int flag_excess;
    /* Model of observation error (enums.hpp) */
    int flag_family;
    /* Day of week effects (1=yes, 0=no) */
    int flag_day_of_week;
    /* Priors on nonlinear model parameters (enums.hpp, -1=no prior) */
    vector<int> flag_regularize;
    /* Trace
       0=nothing
       1=degenerate nll terms
       2=all nll terms
    */
    int flag_trace;
    /* X format (1=sparse, 0=dense) */
    int flag_sparse_X;
    /* predict (1=yes, 0=no) */
    int flag_predict;

    bool do_excess;
    bool do_day_of_week;
    bool do_regularize;
    bool do_trace;
    bool do_trace_verbose;
    bool do_sparse_X;
    bool do_random_effects;
    bool do_predict;
    bool do_simulate;
  
    flags_t(SEXP x)
    {
        flag_curve       = get_list_integer(x, "flag_curve");
	flag_excess      = get_list_integer(x, "flag_excess");
	flag_family      = get_list_integer(x, "flag_family");
	flag_day_of_week = get_list_integer(x, "flag_day_of_week");
	flag_regularize  = asVector<int>(getListElement(x, "flag_regularize", &Rf_isNumeric));
	flag_trace       = get_list_integer(x, "flag_trace");
	flag_sparse_X    = get_list_integer(x, "flag_sparse_X");
	flag_predict     = get_list_integer(x, "flag_predict");

	do_excess        = (flag_excess         == 1);
	do_day_of_week   = (flag_day_of_week    == 1);
	do_regularize    = false;
	do_trace         = (flag_trace          >= 1);
	do_trace_verbose = (flag_trace          >= 2);
	do_sparse_X      = (flag_sparse_X       == 1);
	do_predict       = (flag_predict        == 1);
	for (int i = 0; !do_regularize && i < flag_regularize.size(); ++i)
	{
	    do_regularize = flag_regularize(i) >= 0;
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

} // namespace egf
