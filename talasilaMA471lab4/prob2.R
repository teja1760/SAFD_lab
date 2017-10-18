data = read.table('d-csp0108.txt',header=T)
attach(data)
C_log = log(C+1)
SP_log = log(SP+1)

h = function(p)
{
    log(p/(1-p))
}

h_prime = function(p)
{
    1/(p*(1-p))
}

h_inv = function(p)
{
    exp(p)/(1 + exp(p))
}

conf_interval = function(X,dist="normal")
{
    mu = mean(X)
    sigma = sd(X)
    n = length(X)
    if(dist=="normal")
    {
        return(c(mu-1.96*sigma/sqrt(n),mu+1.96*sigma/sqrt(n)))
    }
    else if(dist=="bernoulli")
    {
        p_hat = mu
        return(c(p_hat-1.96*sqrt(p_hat*(1-p_hat)/n),p_hat+1.96*sqrt(p_hat*(1-p_hat)/n)))
    }
    else if(dist=="normalized_bernoulli")
    {
        p_hat = mu
        if(p_hat==0)
        	return(c(0,0))
        L = h(p_hat) - (1.96*(h_prime(p_hat))*sqrt(p_hat*(1-p_hat)/n))
        U = h(p_hat) + (1.96*(h_prime(p_hat))*sqrt(p_hat*(1-p_hat)/n))
        return(c(h_inv(L),h_inv(U)))
    }
}

interval = conf_interval(C)
cat("Confidence Interval of C: [",interval[1],",",interval[2],"]\n")
mu_c = mean(C)
sigma_c = sd(C)
n = length(C)
N = 1000
count = 0
for(i in 1:N)
{
    sample = rnorm(n,mu_c,sigma_c)
    interval = conf_interval(sample)
    if(interval[1]<=mu_c && mu_c<=interval[2])
        count = count + 1
}

coverage_prob_C = count/N
cat("Coverage Probability of C = ",coverage_prob_C)

p = 0.1
sample_size = c(20,50,100,1000)
for(size in sample_size)
{
    count = 0
    false_count = 0
    for(i in 1:1000)
    {
        sample = rbinom(size,1,p)
        interval = conf_interval(sample,dist="bernoulli")
        if(interval[1]<0 || interval[2]>1)
            false_count = false_count + 1
        if(interval[1]<=p && p<=interval[2])
            count = count + 1
    }
    cat("\nFor sample size = ",size)
    cat("\nCoverage Probability = ",count/1000)
    cat("\nNo. of intervals outside parameter space = ",false_count)
}

for(size in sample_size)
{
    count = 0
    false_count = 0
    for(i in 1:1000)
    {
        sample = rbinom(size,1,p)
        interval = conf_interval(sample,dist="normalized_bernoulli")
        if(interval[1]<0 || interval[2]>1)
            false_count = false_count + 1
        if(interval[1]<=p && p<=interval[2])
            count = count + 1
    }
    cat("\nFor sample size = ",size)
    cat("\nCoverage Probability = ",count/1000)
    cat("\nNo. of intervals outside parameter space = ",false_count)
}