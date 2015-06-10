## ---- fig.width=7, fig.height=5------------------------------------------
library(HyperCube)

# The package includes the data set canadian.earnings.
# The age is considered as factor in the data set.
canadian.earnings$age <- factor(canadian.earnings$age)

# Plot the data
plot(as.numeric(as.character(canadian.earnings$age)), canadian.earnings$log.income,
     xlab = "age", ylab = "log(income)")

# The number of ages in the data set, p, in Example 1 in Beran (2014).
p <- length(unique(canadian.earnings[,1]))

# D_5 as in equation (3.10) in Beran (2014)
D <- diffMatrix(p, 5)

# The parametor nu in equation (3.11) in Beran (2014)
nu <- c(0, 10^c(2,5,8,11))

# Plotting Hypercube Estimator fits for varying nu
lcolor <- 1:5
for(k in 1:5) {
  
  # The matrix W in equation (3.11) in Beran (2014)
  W <- nu[k] * t(D) %*% D
  
  # Convert W to V, as described in (1.6) in Beran (2014)
  V <- plsW2V(W)
  
  # Hyperpercube Estimator Fit
  hcmod <- hypercube( log.income ~ age -1, data=canadian.earnings, V)
  
  # Plot the fits
  lines(as.numeric(levels(canadian.earnings$age)), 
        hcmod$coefficients, col = lcolor[k])
  legend("topleft", cex = 0.8, 
         legend = c("0", "10^2", "10^5", "10^8", "10^11"), 
         lty = rep(1,5), col=1:5)
}

## ------------------------------------------------------------------------
library(HyperCube)

# The package includes the data set litter.
# The formula specifying the two-way layout is "weight ~ mother:infant -1".
# hypercubeOptimization computes the optimal d
hcmodopt <- hypercubeOptimization( weight ~ mother:infant -1, data = litter)

# The optimal d
# Same result as stated in Example 2 in Beran (2014) 
hcmodopt$projcoef

# Compare the estimated risk
summary(hcmodopt$est)

