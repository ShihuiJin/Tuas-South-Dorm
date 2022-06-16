#include log_normal sero-conversion rate in the likelihood calculation, LN(2.75,0.62^2)
#allow some seronegative individuals to be infected

#adjacency matrix
distance=matrix(0, nrow=nrow(data), ncol=nrow(data))
for (i in 1:nrow(data))
{
  for (j in 1:nrow(data))
  {
    if (data$block[i]==data$block[j]) 
    {
      distance[i,j]=1
      if (data$floor[i]==data$floor[j]){
        distance[i,j]=2
        if (data$unit[i]==data$unit[j]) distance[i,j]=3
      }
    }
  }
}

repre=rep(0, 5*9*14)
for(i in 1:5)
{
  for(j in 1:9)
  {
    for (k in 1:14)
    {
      repre[(i-1)*9*14+(j-1)*14+k]=which(data$floor==j&data$unit==(i*14+14+k))[1]
      if (is.na(repre[(i-1)*9*14+(j-1)*14+k])==1) {
        repre[(i-1)*9*14+(j-1)*14+k]=0
      }
    }
  }
}

matching=rep(0,nrow(data))
for(i in 1:nrow(data))
{
  matching[i]=(data$block[i]-5)*9*7+(data$floor[i]-1)*14+data$unit[i]-data$block[i]*7+7
}

distance_repre=matrix(0, nrow=length(repre), ncol=length(repre))
for (i in 1:length(repre))
{
  for (j in 1:length(repre))
  {
    if (repre[i]>0&repre[j]>0){
      if (data$block[repre[i]]==data$block[repre[j]]) 
      {
        distance_repre[i,j]=1
        if (data$floor[repre[i]]==data$floor[repre[j]]){
          distance_repre[i,j]=2
          if (data$unit[repre[i]]==data$unit[repre[j]]) distance_repre[i,j]=3
        }
      }
    }
  }
}

Rcpp::sourceCpp("model1_functions_rcpp.cpp")

#lambda matrix
lambda_matrix=function(current, data){
  lambda_i=matrix(0, nrow=nrow(data), ncol=length(current$inf_order))
  lambda_small=matrix(0,nrow=length(repre), ncol=length(current$inf_order))
  #lambda at the beginning,before the first event takes place
  for (i in 1:length(repre))
  {
    if(repre[i]==0) next
    l_unit=sum(data$t_inf<=start&data$t_rem>start&distance[repre[i],]==3)
    l_floor=sum(data$t_inf<=start&data$t_rem>start&distance[repre[i],]==2)
    l_block=sum(data$t_inf<=start&data$t_rem>start&distance[repre[i],]==1)
    lambda_small[i,1]=(l_unit*(current$d1+current$d2+1)+l_floor*(current$d1+current$d2)+l_block*current$d1)*current$b3
  }
  
  for (t in 2:ncol(lambda_small))
  {
    lambda_small[,t]=lambda_small[,t-1]
    if (length(which(data$t_inf==current$inf_order[t]))>0)
    {
      ii=matching[which(data$t_inf==current$inf_order[t])]
      lambda_small[which(distance_repre[ii,]==3),t]=lambda_small[which(distance_repre[ii,]==3),t]+(current$d1+current$d2+1)*current$b3
      lambda_small[which(distance_repre[ii,]==2),t]=lambda_small[which(distance_repre[ii,]==2),t]+(current$d1+current$d2)*current$b3
      lambda_small[which(distance_repre[ii,]==1),t]=lambda_small[which(distance_repre[ii,]==1),t]+current$d1*current$b3
    }
    if (length(which(data$t_rem==current$inf_order[t]))>0)
    {
      ii=matching[which(data$t_rem==current$inf_order[t])]
      lambda_small[which(distance_repre[ii,]==3),t]=lambda_small[which(distance_repre[ii,]==3),t]-(current$d1+current$d2+1)*current$b3
      lambda_small[which(distance_repre[ii,]==2),t]=lambda_small[which(distance_repre[ii,]==2),t]-(current$d1+current$d2)*current$b3
      lambda_small[which(distance_repre[ii,]==1),t]=lambda_small[which(distance_repre[ii,]==1),t]-current$d1*current$b3
    }
  }
  lambda_i=copy_c(lambda_i, lambda_small, matching, nrow(lambda_i), ncol(lambda_i))
}


logposterior=function(current, data)
{
  inf_order=c(start, sort(unique(c(data$t_inf[which(data$t_inf>start&data$t_inf<=end)],data$t_rem[which(data$t_rem>start&data$t_rem<=end)]))))
  current$inf_order=inf_order
  deltat=current$inf_order[-1]-current$inf_order[-length(current$inf_order)]
  lambda_i=lambda_matrix(current,data)
  current$lambda_i=lambda_i
  inf_place=rep(0,nrow(data))
  for (i in 1:nrow(data)){
    if(sum(current$inf_order==data$t_inf[i])>0) inf_place[i]=which(current$inf_order==data$t_inf[i])
  }
  end_place=rep(0,end+1)
  for (i in 1:(end+1))
  {
    end_place[i]=sum(current$inf_order<=i+138)
  }
  current$logposterior = logposterior_c(current, missing_matrix_c, onset_c,rep(0,15),first_inf, issym, data,end_place,inf_place,inf_order, deltat, lambda_i,nrow(data), start, end)
  current 
}


#updated logposterior (based on lambda_i(t_inf) calculation)
likelihood_contr=function(data, k){
  inf_order=c(start,sort(c(data$t_inf[which(data$t_inf>start&data$t_inf<=data$t_inf[k])],data$t_rem[which(data$t_rem>start&data$t_rem<=data$t_inf[k])])))
  inf_list=data.frame(place=rep(0, length(inf_order)-1), time=rep(0, length(inf_order)-1), inf=rep(0, length(inf_order)-1))
  inf_list$place=c(which(data$t_inf>start&data$t_inf<=data$t_inf[k]),which(data$t_rem>start&data$t_rem<=data$t_inf[k]))
  inf_list$time=c(data$t_inf[which(data$t_inf>start&data$t_inf<=data$t_inf[k])],data$t_rem[which(data$t_rem>start&data$t_rem<=data$t_inf[k])])
  inf_list$inf[1:sum(data$t_inf>start&data$t_inf<=data$t_inf[k])]=rep(1, sum(data$t_inf>start&data$t_inf<=data$t_inf[k]))          
  inf_list=inf_list[order(inf_list$time),]
  lambda_new=rep(0, length(inf_order)-1)
  lambda_new[1]=lambda_i[k,1]
  if (length(inf_order)==2) return(-lambda_new[1]*(inf_order[2]-start))
  lambda_new=lambda_new_cal(lambda_new, distance, inf_list, current$d1*current$b3, current$d2*current$b3, current$b3, k, length(lambda_new))
  deltat=inf_order[-1]-inf_order[-length(inf_order)]
  loglikelihood_c(lambda_new, deltat, length(deltat))
}
likelihood_contr1=function(data, k, end_time){
  inf_order=c(start,sort(c(data$t_inf[which(data$t_inf>start&data$t_inf<=end_time)],data$t_rem[which(data$t_rem>start&data$t_rem<=end_time)])))
  inf_list=data.frame(place=rep(0, length(inf_order)-1), time=rep(0, length(inf_order)-1), inf=rep(0, length(inf_order)-1))
  inf_list$place=c(which(data$t_inf>start&data$t_inf<=end_time),which(data$t_rem>start&data$t_rem<=end_time))
  inf_list$time=c(data$t_inf[which(data$t_inf>start&data$t_inf<=end_time)],data$t_rem[which(data$t_rem>start&data$t_rem<=end_time)])
  inf_list$inf[1:sum(data$t_inf>start&data$t_inf<=end_time)]=rep(1, sum(data$t_inf>start&data$t_inf<=end_time))          
  inf_list=inf_list[order(inf_list$time),]
  lambda_new=rep(0, length(inf_order))
  lambda_new[1]=lambda_i[k,1]
  lambda_new=lambda_new_cal(lambda_new, distance, inf_list, current$d1*current$b3, current$d2*current$b3, current$b3, k, length(lambda_new))
  deltat=inf_order[-1]-inf_order[-length(inf_order)]
  loglikelihood_c(lambda_new[-length(inf_order)], deltat, length(deltat))-lambda_new[length(inf_order)]*(end_time-inf_order[length(inf_order)])
}

#first_inf=c(277,597,2280,2504,3883) #fix the logposterior of these as 0
#include prior for infection times of symptomatic infections
logposterior_inf=function(data, old_data, k){
  if(data$t_inf[k]>end&old_data$t_inf[k]>end){
    data$delta=0
    return(data)
  }
  if(issym[k]>0){
    delta=incubation_c(current, data$t_inf[k], old_data$t_inf[k], missing_matrix_c, onset_c, issym[k]-1)
  }else{
    delta=0
  }
  if(data$S0[k]==1&is.na(location$T0[k])==0) delta=delta+plnorm(location$T0[k]-139-data$t_inf[k], 2.75,0.62, log=TRUE)-plnorm(location$T0[k]-139-old_data$t_inf[k], 2.75,0.62, log=TRUE)
  if(data$S0[k]==0&is.na(data$T0[k])==0) delta=delta+log(1-plnorm(location$T0[k]-139-data$t_inf[k], 2.75,0.62))-log(1-plnorm(location$T0[k]-139-old_data$t_inf[k], 2.75,0.62))
  if(data$S0[k]==0){
    if(data$S1[k]==1&is.na(location$T1[k])==0) delta=delta+plnorm(location$T1[k]-139-data$t_inf[k], 2.75,0.62, log=TRUE)-plnorm(location$T1[k]-139-old_data$t_inf[k], 2.75,0.62, log=TRUE)
    if(data$S1[k]==0&is.na(location$T1[k])==0) delta=delta+log(1-plnorm(location$T1[k]-139-data$t_inf[k], 2.75,0.62))-log(1-plnorm(location$T1[k]-139-old_data$t_inf[k], 2.75,0.62))
  }
  inf_new=min(data$t_inf[k], end)
  inf_old=min(data$t_inf[k], end)
  rem_new=min(data$t_rem[k], end)
  rem_old=min(old_data$t_rem[k],end)
  if(data$t_inf[k]>old_data$t_inf[k]){
    for(i in which(distance[,k]==3&data$t_inf>=inf_old))
    {
      if(i%in%first_inf) next
      if(i==k){
        l_unit=sum(data$t_inf<inf_new&data$t_rem>=inf_new&distance[,k]==3)
        l_floor=sum(data$t_inf<inf_new&data$t_rem>=inf_new&distance[,k]==2)
        l_block=sum(data$t_inf<inf_new&data$t_rem>=inf_new&distance[,k]==1)
        data$lambda_inf[k]=(current$d1*l_block+(current$d1+current$d2)*l_floor+(current$d1+current$d2+1)*l_unit)*current$b3
        if (data$lambda_inf[i]<=0){
          delta=-999999
          data$delta=delta
          return(data)
        }else{
          if (data$t_inf[k]<=end){
            delta=delta-log(old_data$lambda_inf[i])+log(data$lambda_inf[i])+likelihood_contr(data,k)-likelihood_contr(old_data,k)
          }else{
            data$lambda_inf[i]=1
            delta=delta-log(old_data$lambda_inf[i])+likelihood_contr1(data,k,end)-likelihood_contr(old_data,k)
          }
        }
      }else if (data$S1[i]==2&data$S0[i]==0&data$t_inf[i]>end){
        if (data$T0[i]-139>rem_new) next
        if (data$T0[i]-139>inf_old) delta=delta+current$b3*(current$d1+current$d2+1)*(min(data$T0[i]-139,inf_new)-inf_old)
        if (data$T0[i]-139>rem_old) delta=delta-current$b3*(current$d1+current$d2+1)*((data$T0[i]-139)-rem_old)
      }else if(data$S1[i]==0&data$t_inf[i]>end){
        if (rem_new<data$T1[i]-139) next
        if (data$T1[i]-139>inf_old) delta=delta+current$b3*(current$d1+current$d2+1)*(min(data$T1[i]-139,inf_new)-inf_old)
        if (data$T1[i]-139>rem_old) delta=delta-current$b3*(current$d1+current$d2+1)*((data$T1[i]-139)-rem_old)
      }else if(data$t_inf[i]<=inf_new){
          data$lambda_inf[i]=old_data$lambda_inf[i]-current$b3*(current$d1+current$d2+1)
          if (data$lambda_inf[i]<=0){
            delta=-999999
            data$delta=delta
            return(data)
          }else{
            delta=delta-log(old_data$lambda_inf[i])+log(data$lambda_inf[i])+current$b3*(current$d1+current$d2+1)*(data$t_inf[i]-inf_old)
          }
      }else if(data$t_inf[i]>rem_old&data$t_inf[i]<=rem_new){
        data$lambda_inf[i]=old_data$lambda_inf[i]+current$b3*(current$d1+current$d2+1)
        if (data$lambda_inf[i]<=0){
          delta=-999999
          data$delta=delta
          return(data)
        }else{
          delta=delta-log(old_data$lambda_inf[i])+log(data$lambda_inf[i])+current$b3*(current$d1+current$d2+1)*(rem_new-data$t_inf[i])
        }
      }else if (data$t_inf[i]>inf_new&data$t_inf[i]<=rem_old){
        delta=delta+current$b3*(current$d1+current$d2+1)*(inf_new-inf_old)
      }
    }
    for(i in which(distance[,k]==2&data$t_inf>=inf_old))
    {
      if(i%in%first_inf) next
      if (data$S1[i]==2&data$S0[i]==0&data$t_inf[i]>end){
        if (data$T0[i]-139>rem_new) next
        if (data$T0[i]-139>inf_old) delta=delta+current$b3*(current$d1+current$d2)*(min(data$T0[i]-139,inf_new)-inf_old)
        if (data$T0[i]-139>rem_old) delta=delta-current$b3*(current$d1+current$d2)*((data$T0[i]-139)-rem_old)
      }else if(data$S1[i]==0&data$t_inf[i]>end){
        if (data$T1[i]-139>rem_new) next
        if (data$T1[i]-139>inf_old) delta=delta+current$b3*(current$d1+current$d2)*(min(data$T1[i]-139,inf_new)-inf_old)
        if (data$T1[i]-139>rem_old) delta=delta-current$b3*(current$d1+current$d2)*((data$T1[i]-139)-rem_old)        
      }else if (data$t_inf[i]<=inf_new){
        delta=delta-log(old_data$lambda_inf[i])
        data$lambda_inf[i]=old_data$lambda_inf[i]-current$b3*(current$d1+current$d2)
        if (data$lambda_inf[i]<=0){
          delta=-999999
          data$delta=delta
          return(data)
        }else{
          delta=delta+log(data$lambda_inf[i])+current$b3*(current$d1+current$d2)*(data$t_inf[i]-inf_old)
        }
      }else if(data$t_inf[i]>rem_old&data$t_inf[i]<=rem_new){
        data$lambda_inf[i]=old_data$lambda_inf[i]+current$b3*(current$d1+current$d2)
        if (data$lambda_inf[i]<=0){
          delta=-999999
          data$delta=delta
          return(data)
        }else{
          delta=delta-log(old_data$lambda_inf[i])+log(data$lambda_inf[i])+current$b3*(current$d1+current$d2)*(rem_new-data$t_inf[i])
        }
      }else if (data$t_inf[i]>inf_new&data$t_inf[i]<=rem_old){
        delta=delta+current$b3*(current$d1+current$d2)*(inf_new-inf_old)
      }
    }
    for(i in which(distance[,k]==1&data$t_inf>inf_old))
    {
      if(i%in%first_inf) next
      if (data$S1[i]==2&data$S0[i]==0&data$t_inf[i]>end){
        if (data$T0[i]-139>rem_new) next
        if (data$T0[i]-139>inf_old) delta=delta+current$b3*current$d1*(min(data$T0[i]-139,inf_new)-inf_old)
        if (data$T0[i]-139>rem_old) delta=delta-current$b3*current$d1*((data$T0[i]-139)-rem_old)
      }else if(data$S1[i]==0&data$t_inf[i]>end){
        if (data$T1[i]-139>rem_new) next
        if (data$T1[i]-139>inf_old) delta=delta+current$b3*current$d1*(min(data$T1[i]-139,inf_new)-inf_old)
        if (data$T1[i]-139>rem_old) delta=delta-current$b3*current$d1*((data$T1[i]-139)-rem_old)        
      }else if(data$t_inf[i]<=inf_new){
        delta=delta-log(old_data$lambda_inf[i])
        data$lambda_inf[i]=old_data$lambda_inf[i]-current$b3*current$d1
        if (data$lambda_inf[i]<=0){
          delta=-999999
          data$delta=delta
          return(data)
        }else{
          delta=delta+log(data$lambda_inf[i])+current$b3*current$d1*(data$t_inf[i]-inf_old)
        }
      }else if(data$t_inf[i]>rem_old&data$t_inf[i]<=rem_new){
        data$lambda_inf[i]=old_data$lambda_inf[i]+current$b3*current$d1
        if (data$lambda_inf[i]<=0){
          delta=-999999
          data$delta=delta
          return(data)
        }else{
          delta=delta-log(old_data$lambda_inf[i])+log(data$lambda_inf[i])+current$b3*current$d1*(rem_new-data$t_inf[i])
        }
      }else if (data$t_inf[i]>inf_new&data$t_inf[i]<=rem_old){
        delta=delta+current$b3*current$d1*(inf_new-inf_old)
      }
    }
  }
  if(data$t_inf[k]<old_data$t_inf[k]){
    for(i in which(distance[,k]==3&data$t_inf>=inf_new))
    {
      if(i%in%first_inf) next
      if(i==k){
          l_unit=sum(data$t_inf<inf_new&data$t_rem>=inf_new&distance[,k]==3)
          l_floor=sum(data$t_inf<inf_new&data$t_rem>=inf_new&distance[,k]==2)
          l_block=sum(data$t_inf<inf_new&data$t_rem>=inf_new&distance[,k]==1)
          data$lambda_inf[k]=(current$d1*l_block+(current$d1+current$d2)*l_floor+(current$d1+current$d2+1)*l_unit)*current$b3
          if (data$lambda_inf[i]<=0){
            delta=-999999
            data$delta=delta
            return(data)
          }else{
            if (old_data$t_inf[k]<=end){
              delta=delta-log(old_data$lambda_inf[i])+log(data$lambda_inf[i])+likelihood_contr(data,k)-likelihood_contr(old_data,k)
            }else{
              delta=delta+log(data$lambda_inf[i])+likelihood_contr(data,k)-likelihood_contr1(old_data,k, end)
            }
          }
      }else if (data$S1[i]==2&data$S0[i]==0&data$t_inf[i]>end){
        if (data$T0[i]-139>rem_old) next
        if (data$T0[i]-139>inf_new) delta=delta-current$b3*(current$d1+current$d2+1)*(min(data$T0[i]-139,inf_old)-inf_new)
        if (data$T0[i]-139>rem_new) delta=delta+current$b3*(current$d1+current$d2+1)*((data$T0[i]-139)-rem_new)
      }else if(data$S1[i]==0&data$t_inf[i]>end){
        if (rem_old<data$T1[i]-139) next
        if (data$T1[i]-139>inf_new) delta=delta-current$b3*(current$d1+current$d2+1)*(min(data$T1[i]-139,inf_old)-inf_new)
        if (data$T1[i]-139>rem_new) delta=delta+current$b3*(current$d1+current$d2+1)*((data$T1[i]-139)-rem_new)
      }else if(data$t_inf[i]<=inf_old){
          data$lambda_inf[i]=old_data$lambda_inf[i]+current$b3*(current$d1+current$d2+1)
          if (data$lambda_inf[i]<=0){
            delta=-999999
            data$delta=delta
            return(data)
          }else{
            delta=delta-log(old_data$lambda_inf[i])+log(data$lambda_inf[i])-current$b3*(current$d1+current$d2+1)*(data$t_inf[i]-inf_new)
          }
        
      }else if(data$t_inf[i]>rem_new&data$t_inf[i]<=rem_old){
        data$lambda_inf[i]=old_data$lambda_inf[i]-current$b3*(current$d1+current$d2+1)
        if (data$lambda_inf[i]<=0){
          delta=-999999
          data$delta=delta
          return(data)
        }else{
          delta=delta-log(old_data$lambda_inf[i])+log(data$lambda_inf[i])-current$b3*(current$d1+current$d2+1)*(rem_old-data$t_inf[i])
        }
      }else if (data$t_inf[i]>inf_old&data$t_inf[i]<=rem_new){
        delta=delta-current$b3*(current$d1+current$d2+1)*(inf_old-inf_new)
      }
    }
    for(i in which(distance[,k]==2&data$t_inf>inf_new))
    {
      if(i%in%first_inf) next
      if (data$S1[i]==2&data$S0[i]==0&data$t_inf[i]>end){
        if (data$T0[i]-139>rem_old) next
        if (data$T0[i]-139>inf_new) delta=delta-current$b3*(current$d1+current$d2)*(min(data$T0[i]-139,inf_old)-inf_new)
        if (data$T0[i]-139>rem_new) delta=delta+current$b3*(current$d1+current$d2)*((data$T0[i]-139)-rem_new)
      }else if(data$S1[i]==0&data$t_inf[i]>end){
        if (rem_old<data$T1[i]-139) next
        if (data$T1[i]-139>inf_new) delta=delta-current$b3*(current$d1+current$d2)*(min(data$T1[i]-139,inf_old)-inf_new)
        if (data$T1[i]-139>rem_new) delta=delta+current$b3*(current$d1+current$d2)*((data$T1[i]-139)-rem_new)
      }else if(data$t_inf[i]<=inf_old){
        data$lambda_inf[i]=old_data$lambda_inf[i]+current$b3*(current$d1+current$d2)
        if (data$lambda_inf[i]<=0){
          delta=-999999
          data$delta=delta
          return(data)
        }else{
          delta=delta-log(old_data$lambda_inf[i])+log(data$lambda_inf[i])-current$b3*(current$d1+current$d2)*(data$t_inf[i]-inf_new)
        }
      }else if(data$t_inf[i]>rem_new&data$t_inf[i]<=rem_old){
        data$lambda_inf[i]=old_data$lambda_inf[i]-current$b3*(current$d1+current$d2)
        if (data$lambda_inf[i]<=0){
          delta=-999999
          data$delta=delta
          return(data)
        }else{
          delta=delta-log(old_data$lambda_inf[i])+log(data$lambda_inf[i])-current$b3*(current$d1+current$d2)*(rem_old-data$t_inf[i])
        }
      }else if (data$t_inf[i]>inf_old&data$t_inf[i]<=rem_new){
        delta=delta-current$b3*(current$d1+current$d2)*(inf_old-inf_new)
      }
    }
    for(i in which(distance[,k]==1&data$t_inf>inf_new))
    {
      if(i%in%first_inf) next
      if (data$S1[i]==2&data$S0[i]==0&data$t_inf[i]>end){
        if (data$T0[i]-139>rem_old) next
        if (data$T0[i]-139>inf_new) delta=delta-current$b3*current$d1*(min(data$T0[i]-139,inf_old)-inf_new)
        if (data$T0[i]-139>rem_new) delta=delta+current$b3*current$d1*((data$T0[i]-139)-rem_new)
      }else if(data$S1[i]==0&data$t_inf[i]>end){
        if (rem_old<data$T1[i]-139) next
        if (data$T1[i]-139>inf_new) delta=delta-current$b3*current$d1*(min(data$T1[i]-139,inf_old)-inf_new)
        if (data$T1[i]-139>rem_new) delta=delta+current$b3*current$d1*((data$T1[i]-139)-rem_new)
      }else if(data$t_inf[i]<=inf_old){
        data$lambda_inf[i]=old_data$lambda_inf[i]+current$b3*current$d1
        if (data$lambda_inf[i]<=0){
          delta=-999999
          data$delta=delta
          return(data)
        }else{
          delta=delta-log(old_data$lambda_inf[i])+log(data$lambda_inf[i])-current$b3*current$d1*(data$t_inf[i]-inf_new)
        }
      }else if(data$t_inf[i]>rem_new&data$t_inf[i]<=rem_old){
        data$lambda_inf[i]=old_data$lambda_inf[i]-current$b3*current$d1
        if (data$lambda_inf[i]<=0){
          delta=-999999
          data$delta=delta
          return(data)
        }else{
          delta=delta-log(old_data$lambda_inf[i])+log(data$lambda_inf[i])-current$b3*current$d1*(rem_old-data$t_inf[i])
        }
      }else if (data$t_inf[i]>inf_old&data$t_inf[i]<=rem_new){
        delta=delta-current$b3*current$d1*(inf_old-inf_new)
      }
    }
  }
  data$delta=delta
  data
}
