#############Ques 1#######################
cat("\t\t\t\tQuestion 1\n\n")
N_ar<-array(c(10,20,100,500,1000))
size_ar<-array(c(20,50,100,500,1000))
#N_ar<-array(c(10,100,1000,5000,10000))
#size_ar<-array(c(20,50,100,500,1000))   ######### given in question
p_seed<-.4                         ######### just a random value
conf_seed<-95                      ######### 95% confidence interval
z_up<-qchisq(conf_seed/100,1)     ##########  quantile point  



coveragep_ar<-array(dim=c(length(size_ar),2))
p_ar<-seq(0,1,by = .01)
for(k_v in 1:length(N_ar)){
N<-N_ar[k_v]
for(j in 1:length(size_ar)){
	counterp<-0
	for(i in 1:N){
		sample_ar<-rbinom(size_ar[j],1,p_seed)
		samp_mean<-mean(sample_ar)
		samp_sum<-sum(sample_ar)
		Likelihood_mle<-(log(samp_mean)*samp_sum)+(log(1-samp_mean)*(size_ar[j]-samp_sum))
		Likelihood_p<-(log(p_ar)*samp_sum)+(log(1-p_ar)*(size_ar[j]-samp_sum))
		p_final<-2*(Likelihood_mle-Likelihood_p)<z_up
		ci_low<-min(p_ar[p_final])
		ci_high<-max(p_ar[p_final])
		#cat(paste("The confidence interval is\t\t[",ci_low,",",ci_high,"]\n"))
		if(ci_high!=Inf && ci_low!=Inf && ci_low!=-Inf)
			if((p_seed<=ci_high)&&(p_seed>=ci_low))
				counterp<-counterp+1
		
	}
	coveragep_ar[j,1]<-size_ar[j]
	coveragep_ar[j,2]<-counterp/N
}
cat(paste("\n_____________________________________________________________\n"))
cat(paste("\nN value is ",N,"\n"))
#print(coveragep_ar)
cat(paste("The theoretical coverage probability is ",conf_seed/100,"\n\n"))

for(i in 1:length(size_ar)){
	cat(paste("Size=",coveragep_ar[i,1],",\nThe calculated coverage probability is \t\t\t",coveragep_ar[i,2],"\n\n"))
}

}
