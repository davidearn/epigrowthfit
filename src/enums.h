enum curve
{
    exponential,
    subexponential,
    gompertz,
    logistic,
    richards
};

enum family
{
    pois,
    nbinom
};

enum prior
{
    norm,
    lkj,
    wishart,
    invwishart
};

enum test
{
    list_of_vectors_t,
    is_na,
    is_finite,
    logspace_diff,
    mvlgamma,
    log_diag_LLT,
    dlkj,
    dwishart,
    dinvwishart,
    dpois_robust,
    rnbinom_robust,
    eval_log_curve_exponential,
    eval_log_curve_subexponential,
    eval_log_curve_gompertz,
    eval_log_curve_logistic,
    eval_log_curve_richards,
    logspace_add_baseline,
    logspace_add_offsets,
    eval_log_rt_subexponential,
    eval_log_rt_gompertz,
    eval_log_rt_logistic,
    eval_log_rt_richards
};
