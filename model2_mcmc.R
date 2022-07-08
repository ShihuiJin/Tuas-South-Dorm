#a combination of all models
#incorporate different treatments in the model
ID=t(read.csv("ID.csv"))
onset_c=t(read.csv('onset_c.csv'))
cinfected=t(read.csv("cinfected.csv")) #for symptom data
cinfected_ID=ID[cinfected] #in the order of cinfected
location=read.csv("location_outcomes.csv")
data=read.csv('2_data.csv')
missing_matrix_c=as.matrix(read.csv('missing_matrix_c.csv'))

#dose allocation information by place in location file
allocation=as.matrix(read.csv('allocation.csv'))
allocation_l=rep(0,nrow(data))
for(i in 1 :nrow(data))
{
  if (sum(allocation[,1]==location$ID[i])==0){
    allocation_l[i]=0
    next
  }
  t=allocation[which(allocation[,1]==location$ID[i]),2]
  if (t=='A'){
    allocation_l[i]=1
  }else if (t=='B'){
    allocation_l[i]=2
  }else if (t=='C'){
    allocation_l[i]=3
  }else if (t=='D'){
    allocation_l[i]=4
  }else{
    allocation_l[i]=5
    }
}
#missing data
{
  allocation_l[15]=4
  allocation_l[2194]=3
  allocation_l[2606]=3
  allocation_l[2618]=3
  allocation_l[3104]=allocation_l[3112]=5
  allocation_l[4200]=3
}

#log_posterior functions
source("model2_functions.R")
#mh for proposing new b1, b2, b3
mh=function(current,old, data)
{
  reject=FALSE
  if(current$b1<0|current$b2<0|current$b3<0|current$a2<0|current$a3<0|current$theta<0)reject=TRUE
  if(sum(current$t<=0)>0)reject=TRUE
  if (mut_c(exp(current$a1-(current$a2)^2), current$a1, current$a2, current$a3)>=1) reject=TRUE
  if(!reject)
  {
    current$m1=tau_c(current)
    current=logposterior(current, data)
    logaccprob=current$logposterior[1]-old$logposterior[1]
    lu=-rexp(1)
    if(lu>logaccprob)reject=TRUE
  }
  if(reject)current=old
  current
}

#mh for proposing new infection/removal time
mh_inf=function(data,old_data, k)
{
  reject=FALSE
  if(data$t_inf[k]<=data$t_inf[first_inf[(data$block[k]-3)/2]]) reject=TRUE
  if(max(sort(data$t_inf[which(data$block==data$block[k])])[2:sum(data$block==data$block[k]&data$t_inf<=end)]-sort(data$t_inf[which(data$block==data$block[k])])[2:sum(data$block==data$block[k]&data$t_inf<=end)-1])>=10) reject=TRUE
  if(data$S0[k]==0){
    if(is.na(data$T1[k])==0&data$S1[k]==1){
      if(data$t_inf[k]>data$T1[k]-139)reject=TRUE
    }else if(is.na(data$T1[k])==1&data$S1[k]==1){
      if(data$t_inf[k]>end)reject=TRUE
    }
  }else if(data$S0[k]==1){
    if(data$t_inf[k]>data$T0[k]-139) reject=TRUE
  }
  if(data$t_inf[k]<=start)reject=TRUE
  if(issym[k]>0){
    if (data$t_inf[k]+1>onset_c[issym[k]]|data$t_inf[k]+20<=onset_c[issym[k]]) reject=TRUE
    }
  if(!reject)
  {
    data=logposterior_inf(data, old_data, k)
    logaccprob=data$delta[1]
    if(is.finite(logaccprob)==0){
      reject=TRUE
      }else{
        lu=-rexp(1)
        if(lu>logaccprob)reject=TRUE
      }
  }
  if(reject)data=old_data
  data
}


library(MASS)
first_inf=rep(0,5) #fix the logposterior of these as 0
for (i in 1:5)
{
  first_inf[i]=which(data$t_inf==min(data$t_inf[which(data$block==i*2+3)]))
}
issym=rep(0, nrow(data)) #symtomatic infections with value >0, place in onset_c
onset_sym=rep(0,nrow(data))
for(i in 1:nrow(data))
{
  if (length(which(cinfected_ID==data$ID[i]))>0){
    issym[i]=which(cinfected_ID==data$ID[i])
    onset_sym[i]=onset_c[issym[i]]
  }
}
start=-28; end=60
current=list(a1=2.4, a2=0.76, a3=1.87, theta=0.105, m1=1, b1=rnorm(1,7.5e-5, 0.5e-5), b2=rnorm(1,1.5e-4, 2e-5), b3=rnorm(1,0.012, 2e-3))
current$m1=tau_c(current)
current$t=rep(1,5)
current$inf_order=c(start, sort(unique(c(data$t_inf[which(data$t_inf>start&data$t_inf<=end)],data$t_rem[which(data$t_rem>start&data$t_rem<=end)]))))
var=diag(c(0.01,0.004,0.01,0.001,5e-6,1e-5,2e-4,rep(0.01,length(current$t)-1)),6+length(current$t),6+length(current$t))/100
var[1:7,1:7]=cov(store[1:MCMCiterations,c(1,2,3,4,6,7,8)])/100
MCMCiterations=10000
acceptance_rate=0
s_inf=lapply(1:5, function(k) matrix(0,nrow=MCMCiterations, ncol=end-start+1) )
logpos=matrix(0,nrow=MCMCiterations,ncol=15)
store=data.frame(a1=rep(0,MCMCiterations), a2=rep(0,MCMCiterations), a3=rep(0,MCMCiterations), theta=rep(0,MCMCiterations), m1=rep(0,MCMCiterations), b1=rep(0,MCMCiterations), b2=rep(0,MCMCiterations), b3=rep(0,MCMCiterations),t2=rep(0,MCMCiterations),t3=rep(0,MCMCiterations),t4=rep(0,MCMCiterations),t5=rep(0,MCMCiterations))
mean_incub=rep(0,MCMCiterations)
for(iteration in 1:MCMCiterations)
{
  if(iteration%%1e3==0)print(paste0('iteration: ',iteration, ';  time: ', Sys.time()))
  old=logposterior(current, data)
  logpos[iteration,]=old$logposterior
  new=mvrnorm(1,c(as.numeric(current[c(1,2,3,4,6,7,8)]),current$t[-1]),var)
  current[c(1,2,3,4,6,7,8)]=new[1:7]
  current$t[-1]=new[8:(6+length(current$t))]
  current=mh(current,old, data)
  if(current$logposterior[1]!=old$logposterior[1]) acceptance_rate=acceptance_rate+1
  #variables for calculating the log-posterior of augmented t_inf&t_rem
  inf_order=current$inf_order[-1]
  lambda_i=current$lambda_i
  data$lambda_inf=rep(1, nrow(data)) #if infected before -10, lambda_inf=1
  for (i in which(data$t_inf>start&data$t_inf<=end)) #one-time use, for calculating lambda_inf
  {
    data$lambda_inf[i]=lambda_i[i, which(inf_order==data$t_inf[i])]
  }
  for(k in sample(seq(1,nrow(data)),500)) #change infection/removal time of the infection group only
  {
    if(k%in%first_inf)next
    old_data=data
    data$t_inf[k]=rnorm(1,data$t_inf[k],0.5)
    data$t_rem[k]=data$t_inf[k]+10
    data=mh_inf(data, old_data, k)
  }
  store[iteration,1]=current$a1
  store[iteration,2]=current$a2
  store[iteration,3]=current$a3
  store[iteration,4]=current$theta
  store[iteration,5]=current$m1
  store[iteration,6]=current$b1
  store[iteration,7]=current$b2
  store[iteration,8]=current$b3
  store[iteration,9:(7+length(current$t))]=current$t[-1]
  mean_incub[iteration]=mean(data$t_inf[which(issym>0)]-onset_sym[which(issym>0)])
  for(k in 1:5)
  {
    s_inf[[k]][iteration,]=as.vector(do.call("cbind", lapply(start:end, function(t) sum(data$t_inf[data$block==k*2+3]<=t)-sum(data$t_inf[data$block==k*2+3]<=t-1) ))) #new infections in the time interval (t-1,t]
  }
}

r=sample(seq(1,1e4),1)
lapply(1:5, function(k) write.csv(s_inf[[k]], paste0(r,"_inf_result_",k,".csv"), row.names = FALSE) )
write.csv(logpos, paste0(r,"_logpos_result.csv"), row.names = FALSE)
write.csv(store, paste0(r,"_store_result.csv"), row.names = FALSE)
