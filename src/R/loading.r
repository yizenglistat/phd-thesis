# load pacakge
library(MASS)
library(Matrix)
library(splines2)
library(lsei)
library(crayon)

# load customized functions
source('./R/help_funs.r')
source('./R/basis_splines.r')
source('./R/basis_splines_deriv.r')
source('./R/basis_splines_link_funs.r')
source('./R/parameters_true.r')
source('./R/link_funs_true.r')
source('./R/simulate.r')
source('./R/cond_expectation_funs.r')
source('./R/expected_loglike.r')



source('./R/get_initials.r')
source('./R/mle.r')
source('./R/update_alpha.r')
source('./R/update_beta.r')
source('./R/update_delta.r')