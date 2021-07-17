namespace egf
{

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
    int log_nbdisp;
    int log_w1;
    int log_w2;
    int log_w3;
    int log_w4;
    int log_w5;
    int log_w6;
  
    indices_t(SEXP x)
    {
        index_log_r      = getListInteger(x, "index_log_r",      -1);
	index_log_alpha  = getListInteger(x, "index_log_alpha",  -1);
	index_log_c0     = getListInteger(x, "index_log_c0",     -1);
	index_log_tinfl  = getListInteger(x, "index_log_tinfl",  -1);
	index_log_K      = getListInteger(x, "index_log_K",      -1);
	index_logit_p    = getListInteger(x, "index_logit_p",    -1);
	index_log_a      = getListInteger(x, "index_log_a",      -1);
	index_log_b      = getListInteger(x, "index_log_b",      -1);
	index_log_nbdisp = getListInteger(x, "index_log_nbdisp", -1);
	index_log_w1     = getListInteger(x, "index_log_w1",     -1);
	index_log_w2     = getListInteger(x, "index_log_w2",     -1);
	index_log_w3     = getListInteger(x, "index_log_w3",     -1);
	index_log_w4     = getListInteger(x, "index_log_w4",     -1);
	index_log_w5     = getListInteger(x, "index_log_w5",     -1);
	index_log_w6     = getListInteger(x, "index_log_w6",     -1);
    }
};

template<class Type>
struct flags_t
{
    /* Model of cumulative incidence (enum.h) */
    int flag_curve;
    /* Baseline term in model of cumulative incidence (1=yes, 0=no) */
    int flag_excess;
    /* Model of observation error (enum.h) */
    int flag_family;
    /* Day of week effects (1=yes, 0=no) */
    int flag_day_of_week;
    /* Priors on nonlinear model parameters (enum.h, -1=no prior) */
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
    bool do_predict;
    bool do_predict_lci;
    bool do_predict_lii;
    bool do_predict_lrt;
    bool do_simulate;
    bool do_random_effects;
  
    flags_t(SEXP x)
    {
        flag_curve       = getListInteger(x, "flag_curve",       -1);
	flag_excess      = getListInteger(x, "flag_excess",      -1);
	flag_family      = getListInteger(x, "flag_family",      -1);
	flag_day_of_week = getListInteger(x, "flag_day_of_week", -1);
	flag_regularize  = asVector<int>(getListElement(x, "flag_regularize", &Rf_isNumeric));
	flag_trace       = getListInteger(x, "flag_trace",       -1);
	flag_sparse_X    = getListInteger(x, "flag_sparse_X",    -1);
	flag_predict     = getListInteger(x, "flag_predict",     -1);

	do_excess        = (excess      == 1);
	do_day_of_week   = (day_of_week == 1);
	do_regularize    = false;
	do_trace         = (trace       >= 1);
	do_trace_verbose = (trace       >= 2);
	do_sparse_X      = (sparse_X    == 1);
	do_predict       = (predict     == 1);
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
