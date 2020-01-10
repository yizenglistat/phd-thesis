#********************************************#
#
#
#
#
#
#
#
#
#
#
#********************************************#
#     0. Required R packages                 #
#********************************************#
library(Matrix)
library(MASS)
library(splines)
library(splines2)
library(nloptr)
library(lsei)
library(arm)
library(Rsolnp)
library(pspline)
#********************************************#
#     1. Help functions                      #
#********************************************#
source('./loading/help_funs.r')
#********************************************#
#     2. B-splines                           #
#********************************************#
source('./loading/basis_splines.r')
#********************************************#
#     3. B-splines                           #
#********************************************#
source('./loading/basis_splines_deriv.r')
#********************************************#
#     4. Link functions                      #
#********************************************#
source('./loading/basis_splines_link_funs.r')
#********************************************#
#     5. True Beta and delta                 #
#********************************************#
source('./loading/parameters_true.r')
#********************************************#
#     6. True links functions                #
#********************************************#
source('./loading/link_funs_true.r')
#********************************************#
#     7. Simulation data                     #
#********************************************#
source('./loading/simulation_data.r')
#********************************************#
#     8. Conditional expectation functions   #
#********************************************#
source('./loading/cond_expectation_funs.r')
#********************************************#
#     9. Expected log-likelihood function    #
#********************************************#
source('./loading/expected_loglike.r')
#********************************************#
#     10. Get initials for GEM algorithm     #
#********************************************#
source('./loading/get_initials.r')
#********************************************#
#     11. Optimization alpha function        #
#********************************************#
source('./loading/optim_alpha.r')
#********************************************#
#     12. Optimization delta function        #
#********************************************#
source('./loading/optim_delta.r')
#********************************************#
#     13. Optimization beta function         #
#********************************************#
source('./loading/optim_beta.r')
#********************************************#
#     14. Estimate MLE by GEM algorithm      #
#********************************************#
source('./loading/estimate_mle.r')
#********************************************#
#     15. End                                #
#********************************************#
source('./loading/nn_trans.r')