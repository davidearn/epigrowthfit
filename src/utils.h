namespace egf
{

/* Return TRUE if x = NA_real_ and FALSE otherwise */
template<class Type>
bool is_NA_real_(Type x)
{
	return R_IsNA(asDouble(x));
}

/* Similar to R function is.finite */
template<class Type>
bool is_finite(Type x)
{
	return R_finite(asDouble(x));
}

/* Compute log(diff(x)) given log(x) */
template<class Type>
void logspace_diff(vector<Type> &log_x)
{
	int n = log_x.size() - 1;
	for (int i = 0; i < n; ++i)
	{
		log_x(i) = logspace_sub(log_x(i + 1), log_x(i));
	}
	log_x.conservativeResize(n);
	return;
}

} /* namespace egf */
