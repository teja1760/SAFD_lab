y = rnorm(2000,0,sqrt(5))
beta = 0.5*sum(y^2)+0.5
alpha = (length(y)+5)/2
cat("MAP estimate : ")
cat(beta/(alpha+1),"\n")
cat("Bayesian estimate : ")
cat(beta/(alpha-1))