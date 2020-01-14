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
|       +-- output_(seednumber).md
|       +-- output_(seednumber).csv
|   +-- simulation/
|       +-- figures/
|       +-- output_(seednumber).md
|       +-- output_(seednumber).csv
+-- run.r
+-- README.md
```
Open `run.r` in R or RStudio and it will generate outputs.

### Example
An illustrative example is provided. With the default simulation setting and seed number `set.seed(1452)`, we could obtain results in the following table. Its convergence is much faster classical EM algorithm and will converge to global minimal (cost) more precisely.

```r
-----------------------------------------------------------
                          Setting
-----------------------------------------------------------
N                 : 5000
ord               : 5
niknots           : 10
Se                : 0.95 0.95
Sp                : 0.95 0.95
true beta         : 2.0000 -1.0000 -3.0000 4.0000 0.0000
true delta        : 0.3000
-----------------------------------------------------------
                 Accelerated EM Algorithm
```

### Conclusions
To be continue. 