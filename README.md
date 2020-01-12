## Part of My Ph.D. Dissertation at University of South Carolina

### Introduction
Medical researchers are often interested in modeling the disease infection status of individuals to identify important risk factors and to estimate subject-specific risk probabilities. In many cases, pooling specimens (e.g., blood, urine, swabs, etc.) through group testing offers a novel approach to significantly reduce the number of tests, the time expended, and the overall costs. This has led to the adoption of group testing in a number of infectious disease applications

### Libraries
```r
require(MASS)
require(Matrix)
require(splines2)
require(lsei)
require(crayon)
```

### Usage
```
.
+-- R/
+-- output/
|   +-- application/
|       +-- figures/
|       +-- output_(seednumber).txt
|   +-- simulation/
|       +-- figures/
|       +-- output_(seednumber).txt
+-- run.r
+-- README.md
```
Open `run.r` in R or RStudio and it will generate outputs.

### Example
An illustrative example is provided. With the default simulation setting and seed number `set.seed(1452)`, we could obtain results in the following table. Its convergence is much faster classical EM algorithm and will converge to global minimal (cost) more precisely.

```r
-----------------------------------------------------------
setting
-----------------------------------------------------------
N                 : 5000
ord               : 4
niknots           : 8
Se                : 0.95 0.95
Sp                : 0.95 0.95
true beta         : 2.0000 -1.0000 -3.0000 4.0000 0.0000
true delta        : 0.3000
-----------------------------------------------------------
accelerated EM algorithm
-----------------------------------------------------------
last three costs  : 4965.91
EM cost desc      : NA%
nest EM cost diff : 168.4178
est beta          : -0.0015 3e-04 -0.0058 -0.0034 -0.0064
est delta         : 0.4322
-----------------------------------------------------------
last three costs  : 4965.91 3617.643
EM cost desc      : 27.15%
nest EM cost diff : 301.4924
est beta          : 0.4444 -0.2909 -0.6975 1.0093 -0.4708
est delta         : 0.2225
-----------------------------------------------------------
last three costs  : 4965.91 3617.643 1821.714
EM cost desc      : 49.64%
nest EM cost diff : 196.152
est beta          : 2.2047 -1.2505 -3.2984 4.9309 -0.0743
est delta         : 0.1996
-----------------------------------------------------------
last three costs  : 3617.643 1821.714 932.084
EM cost desc      : 48.83%
nest EM cost diff : 3.6008
est beta          : 1.9164 -1.2128 -2.5559 4.0992 -0.1078
est delta         : 0.2496
-----------------------------------------------------------
last three costs  : 1821.714 932.084 789.81
EM cost desc      : 15.26%
nest EM cost diff : 1.2037
est beta          : 1.9604 -1.1346 -2.6969 4.1075 -0.1073
est delta         : 0.2856
-----------------------------------------------------------
last three costs  : 932.084 789.81 752.957
EM cost desc      : 4.67%
nest EM cost diff : 0.5154
est beta          : 2.0333 -1.0656 -2.7998 4.1702 -0.1039
est delta         : 0.3096
-----------------------------------------------------------
last three costs  : 932.084 789.81 752.957
EM cost desc      : 4.67%
nest EM cost diff : 0.0975
est beta          : 2.0311 -1.1237 -2.8358 4.1702 -0.1322
est delta         : 0.3202
-----------------------------------------------------------
last three costs  : 789.81 752.957 727.22
EM cost desc      : 3.42%
nest EM cost diff : 0.0182
est beta          : 2.0613 -1.1146 -2.8844 4.2143 -0.1063
est delta         : 0.3264
-----------------------------------------------------------
last three costs  : 752.957 727.22 697.625
EM cost desc      : 4.07%
nest EM cost diff : 0.0279
est beta          : 2.0738 -1.165 -2.9152 4.2199 -0.1303
est delta         : 0.3267
-----------------------------------------------------------
last three costs  : 727.22 697.625 698.298
EM cost desc      : -0.1%
nest EM cost diff : 0.004
est beta          : 2.0885 -1.1556 -2.9391 4.2216 -0.138
est delta         : 0.327
-----------------------------------------------------------
last four costs   : 697.625 698.298 693.966
EM cost desc      : 0.62%
nest EM cost diff : 0.004
est beta          : 2.0885 -1.1556 -2.9391 4.2216 -0.138
est delta         : 0.327
-----------------------------------------------------------
```

### Conclusions
To be continue. 