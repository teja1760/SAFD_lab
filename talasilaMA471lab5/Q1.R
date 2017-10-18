library(coda)

N = 100
y = rnorm(N, 0, sqrt(5))

a = (N+5)/2
b = (1+sum(y^2))/2

Theta = rgamma(N, shape = a, scale = 1/b)
Theta = 1/Theta
Theta = sort(Theta)

alpha = 1-0.95

Min_interval_length = 9999999

theta_a = 0
theta_b = 0

for(j in 1:(N - as.integer(N*(1-alpha))))
{
  temp = Theta[j+as.integer(N*(1-alpha))] - Theta[j]
  if(temp < Min_interval_length)
  {
    Min_interval_length = temp
    theta_a = Theta[j]
    theta_b = Theta[j+as.integer(N*(1-alpha))]
  }
}

cat("95% Credible Interval for sigma^2 :\n")
cat("Confidence Interval = 
    [ ",theta_a, ", ", theta_b, " ]\n\n")