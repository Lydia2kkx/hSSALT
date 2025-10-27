library(coin)
set.seed(123)

d <- data.frame(
  y = rnorm(20),
  g = gl(2, 10)
)

tst <- wilcox_test(y ~ g, data = d, distribution = approximate(B = 5000))
qperm(tst, c(0.025, 0.975))

########################################################

library(perm)
set.seed(1)

x <- rnorm(20)
g <- gl(2, 10)

res <- permTS(x ~ g,
              alternative = "two.sided",
              method = "exact.mc",
              control = permControl(nmc = 10000, setSEED=F))

res$p.value