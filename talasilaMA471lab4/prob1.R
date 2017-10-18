A = read.table("d-csp0108.txt",header = T)
A[,c(2,3)] = log(1+A[,c(2,3)])
A = as.data.frame(A)
colMeans(A[,c(2,3)])
m = colMeans(A[,c(2,3)])
s = c(sd(A[,2]),sd(A[,3]))

count = 0
for(i in 1:1000)
{
	Z = rnorm(1000,m[1],s[1])
	mu = mean(Z)
	sigma = sd(Z)
	L = mu-1.96*sigma/sqrt(length(Z))
	U = mu+1.96*sigma/sqrt(length(Z))
	if(L < m[1] && m[1] < U)
		count = count + 1
}
cat("Coverage Probability : ",count/1000)