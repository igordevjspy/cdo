cdo
===

Matlab code of a simplified computation of expected payoffs and default probability of CDO
This code reproduce some results of the following paper: 

The Economics of Structured Finance - J. Coval, J. Jurek, E. Stafford. 


============
The code is in Matlab and contains 3 different files: 
- Coprnd: where the copulas are computed. Only the Gaussian is available, the others can easily be implemented. 
- Monte:  MC simulations to compute expected payoffs and default proba
- Main: for non-noisy and noisy simulations on different parameters: default proba of the asset, recovery rate, rho. 


============
