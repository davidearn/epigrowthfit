What do we want to keep? Throw away?

## definitely keep

- `getInits` machinery
- ability to try multiple response distributions (dnbinom, Poisson) and epidemic curves (exponential, logistic, or Richards) 
     - integer flags in the data (or `enum` values) - specify which distribution you're going to use (e.g. https://github.com/glmmTMB/glmmTMB/blob/master/glmmTMB/src/glmmTMB.cpp )
	 - pass unused arguments as `NA` and test for them internally (BMB will make an example!)
	       - think about an interface that will work for all models; maybe we always pass `x0`, `thalf`, *and* `K`, but at least one of them is always `NA`?
     - if parameterization is such that all models are nested/special cases then you can use the `map=` argument in `MakeADFun` to fix some parameters.

## throw away

- S4 machinery
- custom optimizer (hopefully not necessary! combination of (1) reparameterization (2) using TMB (3) maybe??? regularization should make it unnecessary to fight with convergence as much ...

## maybe add

- generalized Richards curve?




