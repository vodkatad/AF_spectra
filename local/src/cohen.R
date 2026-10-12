
library(pwr)
n <- 18+17
u <- 2
v <- n - u -1
pwr.f2.test(u = u, v = v, f2 = NULL, sig.level = 0.05, power = 0.8)

n <- 10+8
u <- 2
v <- n - u -1
pwr.f2.test(u = u, v = v, f2 = NULL, sig.level = 0.05, power = 0.8)