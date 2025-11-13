# hSSALT
hSSALT Rpackage: Analyze a simple heterogeneous SSALT model with exponential lifetime distributions.

## Installation
```r
install.packages("hSSALT")
```

## Example
```r
library(hSSALT)
data(hSSALTdata) 
mle <- MLEhSSALT(data = hSSALTdata$data, 
                 n = 35, 
                 censoring = 1, 
                 tau = c(8, 20), 
                 theta21 = 1, 
                 theta22 = 8, 
                 p = 0.4) 
```

## Author
Yao Lu (Author, Maintainer) — yao.lu1@rwth-aachen.de  

## License
GPL-3 © Yao Lu
