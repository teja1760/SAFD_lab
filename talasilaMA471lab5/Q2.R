# data
data = read.table(file = "d-csp0108.txt", header = TRUE);

C_rtn = data[,2]
SP_rtn = data[,3]

# Log returns calculated
log_C = log(1+C_rtn)
log_SP = log(1+SP_rtn)

alpha = 1-0.95

z = qnorm(1-alpha/2, mean=0, sd=1)

# --------------------------- PART (a) -----------------------

cat("95% Confidence Inervals Assuming 
    data follows normal\n\n")

# ---------------For log returns for Citi Group----------------
# (i) Based on first 50 samples

cat("For mean of daily log return 
    of CitiGroup Stock:\n\n")

sample1_50 = log_C[1:50]
x1_50 = mean(sample1_50)
s1_50 = sd(sample1_50)

a1_50 = x1_50 - s1_50*z/sqrt(50);
b1_50 = x1_50 + s1_50*z/sqrt(50);
  
cat("Based on first 50 sample points:\n")
cat("Confidence Interval: 
    [ ", a1_50, ", ",b1_50," ]\n")

# (ii) Based on Complete sample

X1 = mean(log_C)
s1 = sd(log_C)

a1 = X1 - s1*z/sqrt(length(log_C))
b1 = X1 + s1*z/sqrt(length(log_C))

cat("Based on the Complete sample points:\n")
cat("Confidence Interval: 
    [ ", a1, ", ",b1," ]\n")


# ----------------For log returns for S&P-----------------
# (i) Based on first 50 samples

cat("For mean of daily log 
    return of S&P Stock:\n\n")

sample2_50 = log_SP[1:50]
x2_50 = mean(sample2_50)
s2_50 = sd(sample2_50)

a2_50 = x2_50 - s2_50*z/sqrt(50);
b2_50 = x2_50 + s2_50*z/sqrt(50);

cat("Based on first 50 sample points:\n")
cat("Confidence Interval: 
    [ ", a2_50, ", ",b2_50," ]\n")

# (ii) Based on Complete sample

X2 = mean(log_SP)
s2 = sd(log_SP)

a2 = X2 - s2*z/sqrt(length(log_SP))
b2 = X2 + s2*z/sqrt(length(log_SP))

cat("Based on the Complete sample points:\n")
cat("Confidence Interval: 
    [ ", a2, ", ",b2," ]\n")


# ---------------------- PART (b) ------------------------


# ================= PERCENTILE BOOTSTRAP =================

GenReSample <- function(X)
{
  n = length(X)
  K = sample.int(n, n, replace = TRUE)
  Y = X[K]
  return(Y)
}

cat("95% Percentile Bootstrap Confidence Inervals\n\n")

# --------------For log returns for Citi Group---------------
# (i) Based on first 50 samples

cat("For mean of daily log return 
    of CitiGroup Stock:\n\n")

sample1_50 = log_C[1:50]
n = length(sample1_50)
x1_50 = mean(sample1_50)
s1_50 = sd(sample1_50)

Y1_50 = GenReSample(sample1_50)
Y1_50 = sort(Y1_50)

ay1_50 = Y1_50[as.integer(n*alpha/2)+1]
by1_50 = Y1_50[as.integer(n*(1-alpha/2))+1]

cat("Based on first 50 sample points:\n")
cat("Confidence Interval: 
    [ ", ay1_50, ", ",by1_50," ]\n")


# (ii) Based on Complete sample

n = length(log_C)
X1 = mean(log_C)
s1 = sd(log_C)

Y1 = GenReSample(log_C)
Y1 = sort(Y1)

ay1 = Y1[as.integer(n*alpha/2)+1]
by1 = Y1[as.integer(n*(1-alpha/2))+1]

cat("Based on Complete sample points:\n")
cat("Confidence Interval: 
    [ ", ay1, ", ",by1," ]\n")

# ---------------- For log returns for S&P-----------------
# (i) Based on first 50 samples

cat("For mean of daily log 
    return of S&P Stock:\n\n")

sample2_50 = log_SP[1:50]
n = length(sample2_50)
x2_50 = mean(sample2_50)
s2_50 = sd(sample2_50)

Y2_50 = GenReSample(sample2_50)
Y2_50 = sort(Y2_50)

ay2_50 = Y2_50[as.integer(n*alpha/2)+1]
by2_50 = Y2_50[as.integer(n*(1-alpha/2))+1]

cat("Based on first 50 sample points:\n")
cat("Confidence Interval: 
    [ ", ay2_50, ", ",by2_50," ]\n")


# (ii) Based on Complete sample

n = length(log_SP)
X2 = mean(log_SP)
s2 = sd(log_SP)

Y2 = GenReSample(log_SP)
Y2 = sort(Y2)

ay2 = Y2[as.integer(n*alpha/2)+1]
by2 = Y2[as.integer(n*(1-alpha/2))+1]

cat("Based on Complete sample points:\n")
cat("Confidence Interval: 
    [ ", ay2, ", ",by2," ]\n")


# ==================== BOOT-t =====================

GenerateBootstrapSamples <- function(X, B)
{
  n = length(X)
  
  muX = vector(,B)
  sdX = vector(,B)
  for(i in 1:B)
  {
    K = sample.int(n, n, replace = TRUE)
    xB = X[K]
    muX[i] = mean(xB)
    sdX[i] = sd(xB)
  }
  
  params <- list("mu"=muX, "sd"=sdX)
  return(params)
}


cat("95% Boot-t Confidence Inervals\n\n")

# --------------For log returns for Citi Group--------------
# (i) Based on first 50 samples

cat("For mean of daily log return
    of CitiGroup Stock:\n\n")

sample1_50 = log_C[1:50]
n = length(sample1_50)
x1_50 = mean(sample1_50)
s1_50 = sd(sample1_50)

B = 100
params = GenerateBootstrapSamples(sample1_50, B)
muxb1_50 = params$mu
sdxb1_50 = params$sd

tb1_50 = (x1_50 - muxb1_50)/(sdxb1_50/sqrt(length(n)))
tb1_50 = sort(tb1_50)

L = as.integer(B*alpha/2)+1
U = as.integer(B*(1-alpha/2))+1

ab1_50 = x1_50 + s1_50*tb1_50[L]/sqrt(n) 
bb1_50 = x1_50 + s1_50*tb1_50[U]/sqrt(n)

cat("Based on first 50 sample points:\n")
cat("Confidence Interval: [ ", ab1_50, ", ",bb1_50," ]\n")


# (ii) Based on Complete sample

n = length(log_C)
X1 = mean(log_C)
s1 = sd(log_C)

B = 100
params = GenerateBootstrapSamples(log_C, B)
muXb1 = params$mu
sdXb1 = params$sd

tb1 = (X1 - muXb1)/(sdXb1/sqrt(length(n)))
tb1 = sort(tb1)

L = as.integer(B*alpha/2)+1
U = as.integer(B*(1-alpha/2))+1

ab1 = X1 + s1*tb1[L]/sqrt(n) 
bb1 = X1 + s1*tb1[U]/sqrt(n)

cat("Based on Complete sample points:\n")
cat("Confidence Interval: [ ", ab1, ", ",bb1," ]\n")


# ----------------For log returns for S&P-----------------
# (i) Based on first 50 samples

cat("For mean of daily log return of S&P Stock:\n\n")

sample2_50 = log_SP[1:50]
n = length(sample2_50)
x2_50 = mean(sample2_50)
s2_50 = sd(sample2_50)

B = 100
params = GenerateBootstrapSamples(sample2_50, B)
muxb2_50 = params$mu
sdxb2_50 = params$sd

tb2_50 = (x2_50 - muxb2_50)/(sdxb2_50/sqrt(length(n)))
tb2_50 = sort(tb2_50)

L = as.integer(B*alpha/2)+1
U = as.integer(B*(1-alpha/2))+1

ab2_50 = x2_50 + s2_50*tb2_50[L]/sqrt(n) 
bb2_50 = x2_50 + s2_50*tb2_50[U]/sqrt(n)

cat("Based on first 50 sample points:\n")
cat("Confidence Interval: [ ", ab2_50, ", ",bb2_50," ]\n")


# (ii) Based on Complete sample

n = length(log_SP)
X2 = mean(log_SP)
s2 = sd(log_SP)

B = 100
params = GenerateBootstrapSamples(log_SP, B)
muXb2 = params$mu
sdXb2 = params$sd

tb2 = (X2 - muXb2)/(sdXb2/sqrt(length(n)))
tb2 = sort(tb2)

L = as.integer(B*alpha/2)+1
U = as.integer(B*(1-alpha/2))+1

ab2 = X2 + s2*tb2[L]/sqrt(n) 
bb2 = X2 + s2*tb2[U]/sqrt(n)

cat("Based on Complete sample points:\n")
cat("Confidence Interval: [ ", ab2, ", ",bb2," ]\n")
