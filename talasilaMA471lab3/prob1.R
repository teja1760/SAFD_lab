rmixnormal <- function(n,p,mu1,s1,mu2,s2)
{
	res = vector(,n)
	for(i in 1:n)
	{
		u = runif(1)
		if(u<p)
			res[i] = rnorm(1,mu1,s1)
		else
			res[i] = rnorm(1,mu2,s2)
	}
	return(res)
}
dmixnorm <- function(x,p,mu1,s1,mu2,s2)
{
	return(p*dnorm(x,mu1,s1) + (1-p)*dnorm(x,mu2,s2))
}
dmynorm <- function(X,mu,s)
{
	res = vector(,length(X))
	for(i in 1:length(X))
	{
		res[i] = dnorm(X[i],mu,s)
	}
	return(res)
}
EMmixnorm <- function(X,maxiter=1000)
{
	n = length(X)
	p = runif(1)
	mu1 = rnorm(1,0,1)
	mu2 = rnorm(1,0,1)
	s1 = rexp(1,1)
	s2 = rexp(1,1)
	m = vector(,n)
	for(t in 1:maxiter)
	{
		m = (p*dmynorm(X,mu1,s1))/(p*dmynorm(X,mu1,s1) + (1-p)*dmynorm(X,mu2,s2))
		m[is.na(m)] = 0.5
		
		p = mean(m)
		mu1 = (m%*%X)/sum(m)
		mu2 = ((1-m)%*%X)/sum(1-m)
		s1 = (m%*%((X-mu1)^2))/sum(m)
		s1 = sqrt(s1)
		s2 = ((1-m)%*%((X-mu2)^2))/sum(1-m)
		s2 = sqrt(s2)
	}
	return(c(p,mu1,s1,mu2,s2))
}
X = rmixnormal(200,0.4,0,1,0,5)
v = EMmixnorm(X,1000)
params = c("p","mu1","s1","mu2","s2")
cat("Estimates\n")
for(i in 1:5)
{
	cat(params[i]," = ",v[i],"\n")
}
p = v[1]
mu1 = v[2]
s1 = v[3]
mu2 = v[4]
s2 = v[5]
L = sum(log(dmixnorm(X,0.4,0,1,0,5)))
cat("Log loglihood : ",L)