# load pacakge
require(MASS)
require(Matrix)
require(splines2)
require(lsei)
require(lattice)
require(crayon)
require(GeneralizedHyperbolic)
# install.packages('splines2')
# install.packages('lsei')
# install.packages('lattice')
# install.packages('crayon')


# load customized functions
source('./R/help_funs.r')
source('./R/parameters_true.r')
source('./R/link_funs_true.r')
source('./R/simulate.r')
source('./R/cond_expectation_funs.r')
source('./R/expected_loglike.r')

source('./R/mon_splines.r')

source('./R/mle.r')
source('./R/get_initials.r')
source('./R/update_alpha.r')
source('./R/update_beta.r')
source('./R/update_delta.r')

source('./R/output_header.r')
source('./R/output_body.r')
source('./R/output_tail.r')
source('./R/output_figure.r')
