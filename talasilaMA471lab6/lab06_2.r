#############Ques 2#######################
temp<-read.table("d-csp0108.txt",header=TRUE)
conf<-95     ########95% confidence interval


cat("\n\n\n\t\t\t\tQuestion 2\n\n")


C_ret=array(dim=c(dim(temp)[1],2))
n=dim(temp)[1]
C_ret[,1]=temp[,2]
C_ret[,2]=log(1+temp[,2])

C_mean=mean(C_ret[,2])
C_sd=sd(C_ret[,2])
C_skew=mean(((C_ret[,2]-C_mean)/C_sd)^3)
C_kurt=mean(((C_ret[,2]-C_mean)/C_sd)^4)

z=qnorm((1-(conf/100))/2)
C_low=C_skew+(z*sqrt(6/n))
C_high=C_skew-(z*sqrt(6/n))
cat(paste("Skewness of C data\t\t",C_skew,"\n",conf,"% Skewness interval of C\t [",C_low,",",C_high,"]\n"))
if((C_low<=0)&&(C_high>=0))
  cat(paste("\nSkewness measure of the log-returns for C data is zero.\n\n"))
if((C_low>=0)||(C_high<=0))
  cat(paste("\nSkewness measure of the log-returns for C data is not zero.\n\n"))
C_low=C_kurt+(z*sqrt(24/n))-3
C_high=C_kurt-(z*sqrt(24/n))-3
cat(paste("Excess Kurtosis of C data\t",C_kurt-3,"\n",conf,"% Kurtosis interval of C\t [",C_low,",",C_high,"]\n\n"))
if((C_low<=0)&&(C_high>=0))
  cat(paste("\nExcess Kurtosis measure of the log-returns for C data is zero.\n\n"))
if((C_low>=0)||(C_high<=0))
  cat(paste("\nExcess Kurtosis measure of the log-returns for C data is not zero.\n\n"))



SP_ret=array(dim=c(dim(temp)[1],2))
n=dim(temp)[1]
SP_ret[,1]=temp[,3]
SP_ret[,2]=log(1+temp[,3])

SP_mean=mean(SP_ret[,2])
SP_sd=sd(SP_ret[,2])
SP_skew=mean(((SP_ret[,2]-SP_mean)/SP_sd)^3)
SP_kurt=mean(((C_ret[,2]-C_mean)/C_sd)^4)

z=qnorm((1-(conf/100))/2)
SP_low=SP_skew+(z*sqrt(6/n))
SP_high=SP_skew-(z*sqrt(6/n))
cat(paste("Skewness of SP data\t\t",SP_skew,"\n",conf,"% Skewness interval of SP\t [",SP_low,",",SP_high,"]\n"))
if((SP_low<=0)&&(SP_high>=0))
  cat(paste("\nSkewness measure of the log-returns for SP data is zero.\n\n"))
if((SP_low>=0)||(SP_high<=0))
  cat(paste("\nSkewness measure of the log-returns for SP data is not zero.\n\n"))
SP_low=SP_kurt+(z*sqrt(24/n))-3
SP_high=SP_kurt-(z*sqrt(24/n))-3
cat(paste("Excess Kurtosis of SP data\t",SP_kurt-3,"\n",conf,"% Kurtosis interval of SP\t [",SP_low,",",SP_high,"]\n"))
if((SP_low<=0)&&(SP_high>=0))
  cat(paste("\nExcess Kurtosis measure of the log-returns for SP data is zero.\n\n"))
if((SP_low>=0)||(SP_high<=0))
  cat(paste("\nExcess Kurtosis measure of the log-returns for SP data is not zero.\n\n"))