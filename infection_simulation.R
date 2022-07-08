#simulate infection  by changing block/floor numbers
#randomly assign 4257 persons to different blocks, floors and units

infection_simulation=function(nblock, nfloor,b1,b2,b3,simulation=1e3){
  start=0; end=180
  fake=data.frame(no=seq(1,4257), block=rep(0,4257), floor=rep(0,4257), unit=rep(0,4257), inf=rep(end+10, 4257), rem=rep(end+10, 4257))
  for (i in 1:nrow(fake))
  {
    fake$block[i]=sample(seq(1,nblock),1)
    fake$floor[i]=sample(seq(1,nfloor),1)
    fake$unit[i]=sample(seq(1,14),1)
  }
  distance=matrix(0, nrow=nrow(fake), ncol=nrow(fake))
  for (i in 1:nrow(fake))
  {
    for (j in 1:nrow(fake))
    {
      if (fake$block[i]==fake$block[j]) 
      {
        distance[i,j]=1
        if (fake$floor[i]==fake$floor[j]){
          distance[i,j]=2
          if (fake$unit[i]==fake$unit[j]) distance[i,j]=3
        }
      }
    }
  }
  repre=rep(0, nfloor*14*nblock)
  for(i in 1:nblock)
  {
    for(j in 1:nfloor)
    {
      for (k in 1:14)
      {
        repre[(i-1)*14*nfloor+(j-1)*14+k]=which(fake$block==i&fake$floor==j&fake$unit==k)[1]
        if (is.na(repre[(i-1)*14*nfloor+(j-1)*14+k])==1) {
          repre[(i-1)*14*nfloor+(j-1)*14+k]=0
        }
      }
    }
  }
  matching=rep(0,nrow(fake))
  for(i in 1:nrow(fake))
  {
    matching[i]=(fake$block[i]-1)*14*nfloor+(fake$floor[i]-1)*14+fake$unit[i]
  }
  distance_repre=matrix(0, nrow=length(repre), ncol=length(repre))
  for (i in 1:length(repre))
  {
    for (j in 1:length(repre))
    {
      if (repre[i]>0&repre[j]>0){
        if (fake$block[repre[i]]==fake$block[repre[j]]) 
        {
          distance_repre[i,j]=1
          if (fake$floor[repre[i]]==fake$floor[repre[j]]){
            distance_repre[i,j]=2
            if (fake$unit[repre[i]]==fake$unit[repre[j]]) distance_repre[i,j]=3
          }
        }
      }
    }
  }
  sero1=rep(0,simulation) 
  #simulate the first 50 infections, let them happen 5 days before the start
  inf=50
  initial_inf=rep(floor(inf/nblock),nblock)
  if (inf-floor(inf/nblock)*nblock>0){
    initial_inf[sample(seq(1,nblock),inf-floor(inf/nblock)*nblock)]=initial_inf[sample(seq(1,nblock),inf-floor(inf/nblock)*nblock)]+1
  }
  for (iteration in 1:simulation)
  {
    fake$t_inf=rep(end+10,nrow(fake))
    fake$t_rem=rep(end+10,nrow(fake))
    for(i in 1:nblock)
    {
      s=sample(which(fake$block==i),initial_inf[i])
      fake$t_inf[s]=runif(initial_inf[i],-10,0)
      fake$t_rem[s]=fake$t_inf[s]+10
    }
    event_time=c(fake$t_inf[which(fake$t_inf>start&fake$t_inf<end)], fake$t_rem[which(fake$t_rem>start&fake$t_rem<end)])
    candidate=which(fake$t_inf>end)
    #calculate lambda
    lambda=rep(0, length(repre))
    #initial value
    for (i in 1:length(repre))
    {
      l_unit=sum(fake$t_inf<=start&fake$t_rem>start&distance[repre[i],]==3)
      l_floor=sum(fake$t_inf<=start&fake$t_rem>start&distance[repre[i],]==2)
      l_block=sum(fake$t_inf<=start&fake$t_rem>start&distance[repre[i],]==1)
      lambda[i]=l_unit*(b1+b2+b3)+l_floor*(b1+b2)+l_block*b1
    }
    t=max(fake$t_inf[which(fake$t_inf<end)])
    while(t<end&sum(lambda>0))
    {
      if(length(candidate)<=0) break
      t_new=rep(t, length(candidate))
      for (i in 1:length(candidate))
      {
        m=matching[candidate[i]]
        if(lambda[m]<=0){
          t_new[i]=9999
          lambda[m]=0
        }else{
          t_new[i]=t+rexp(1,lambda[m])
        } 
      }
      if(min(t_new, event_time)>end) break
      if(min(t_new)>=min(event_time)){
        t=min(event_time)
        if(sum(fake$t_inf==min(event_time))>0){
          j=matching[which(fake$t_inf==min(event_time))]
          lambda[j]=lambda[j]+(b1+b2+b3)
          lambda[which(distance_repre[j,]==2)]=lambda[which(distance_repre[j,]==2)]+(b1+b2)
          lambda[which(distance_repre[j,]==1)]=lambda[which(distance_repre[j,]==1)]+b1
        }else{
          j=matching[which(fake$t_rem==min(event_time))]
          lambda[j]=lambda[j]-(b1+b2+b3)
          lambda[which(distance_repre[j,]==2)]=lambda[which(distance_repre[j,]==2)]-(b1+b2)
          lambda[which(distance_repre[j,]==1)]=lambda[which(distance_repre[j,]==1)]-b1
        }
        event_time=event_time[-which.min(event_time)]
      }else{
        t=min(t_new)
        i=which.min(t_new)
        fake$t_inf[candidate[i]]=t
        fake$t_rem[candidate[i]]=t+10
        event_time=c(event_time, t+10)
        j=matching[candidate[i]]
        lambda[j]=lambda[j]+(b1+b2+b3)
        lambda[which(distance_repre[j,]==2)]=lambda[which(distance_repre[j,]==2)]+(b1+b2)
        lambda[which(distance_repre[j,]==1)]=lambda[which(distance_repre[j,]==1)]+b1
        candidate=candidate[-i]
      }
    }
    sero1[iteration]=sum(fake$t_inf<=end)
  }
  c(mean(sero1), quantile(sero1,c(0.5,0.025,0.975,0.25,0.75)))
}

block_simulation=function(b){
  b1=b[1];b2=b[2];b3=b[3]
  out=c()
  nfloor=9
  for(nblock in 4:15)
  {
    print(paste0('nblock = ',nblock, ';  nfloor = ', nfloor, ";  Time = ", Sys.time()) )
    out=rbind(out, infection_simulation(nblock, nfloor,b1,b2,b3, 100))
  }
  write.csv(out, paste0('simulation_', nblock,'_', nfloor, '_', i, '.csv' ), row.names = FALSE)
  out
}

floor_simulation=function(b){
  b1=b[1];b2=b[2];b3=b[3]
  out=c()
  nblock=5
  for(nfloor in 8:19)
  {
    print(paste0('nblock = ',nblock, ';  nfloor = ', nfloor, ";  Time = ", Sys.time()) )
    out=rbind(out, infection_simulation(nblock, nfloor,b1,b2,b3,100))
    write.csv(out, paste0('simulation_', nblock,'_', nfloor, '_', i, '.csv' ), row.names = FALSE)
  }
  out
}

#100 different posterior draws
{
  library(data.table)
  r=sample(1:nrow(store), 100)
  b=cbind(store[r,'d1']*store[r,'b3'],store[r,'d2']*store[r,'b3'] ,store[r,'b3'])
  colMeans(b)
}

#posterior mean
b=matrix(c(7.6e-5, b2=0.00016, b3=0.012),1, 3) #posterior mean


stat=c()
for(i in rev(50:nrow(b))){
  print(paste0('simulation = ', i, ';  Time = ', Sys.time()) )
  out=block_simulation(b[i,])
  #out=floor_simulation(b[i,])
  write.csv(out, paste0('simulation_summary_', i, '.csv' ), row.names = FALSE)
  stat=rbind(stat,out[,1])
}
colMeans(stat)/4257
write.csv(stat,'summary_stat.csv',row.names = FALSE)



#Figure 3
{
  png('new_infections.png',height=8,width=16,units='cm', res=300, pointsize=7)
  pushViewport(plotViewport(c(1,1,1,1)))
  pushViewport(viewport(layout=grid.layout(nrow=1,ncol=2))) 
  pushViewport(viewport(layout.pos.col=2,layout.pos.row=1)) 
  #change block plot
  sero=c(77,59,30,13,5.6,3.7,2.9,2.5,2.2,1.9,1.9,1.8)
  nblock=seq(4,15)
  pushViewport(plotViewport(c(3.5,4.5,3.5,1), xscale=c(0,length(sero)+1),yscale=c(0,85)))
  grid.rect()
  grid.text('(b)', x=0.05, y=0.95)
  for(i in 1:length(sero))
  {
    if(i==2){
      grid.rect(x=i, y=0, height=sero[i], width=0.6, just = c("center","bottom"), default.units='native', gp=gpar(fill=scales::alpha('navyblue',0.8), col=scales::alpha('navyblue',0.1)))
    }else{
      grid.rect(x=i, y=0, height=sero[i], width=0.6, just = c("center","bottom"), default.units='native', gp=gpar(fill=scales::alpha('navyblue',0.5), col=scales::alpha('navyblue',0.1)))
    }
  }
  grid.yaxis()
  grid.xaxis(at=seq(1,length(sero)), label=nblock)
  grid.xaxis(at=seq(1,length(sero)), label=round(4257/9/nblock),main='FALSE')
  grid.text('Final attack rates (%)',x=unit(-3,'lines'),rot=90)
  grid.text('Number of blocks (number of floors = 9)',y=unit(-2.8,'lines'))
  grid.text('Number of men per floor',y=unit(1,'npc')+unit(2.8,'lines'))
  popViewport()
  popViewport()
  pushViewport(viewport(layout.pos.col=1,layout.pos.row=1)) 
  #change floor plot
  sero=c(66,59,52,45,39,32,28,23,19,17,15,13)
  nfloor=seq(8,19)
  pushViewport(plotViewport(c(3.5,4,3.5,1.5), xscale=c(0,length(sero)+1),yscale=c(0,85)))
  grid.rect()
  grid.text('(a)', x=0.05, y=0.95)
  for(i in 1:length(sero))
  {
    if(i==2){
      grid.rect(x=i, y=0, height=sero[i], width=0.6, just = c("center","bottom"), default.units='native', gp=gpar(fill=scales::alpha('indianred',0.8), col=scales::alpha('indianred',0.1)))
    }else{
      grid.rect(x=i, y=0, height=sero[i], width=0.6, just = c("center","bottom"), default.units='native', gp=gpar(fill=scales::alpha('indianred',0.5), col=scales::alpha('indianred',0.1)))
    }
  }
  grid.yaxis()
  grid.xaxis(at=seq(1,length(sero)), label=nfloor)
  grid.xaxis(at=seq(1,length(sero)), label=round(4257/5/nfloor),main='FALSE')
  grid.text('Final attack rates (%)',x=unit(-3,'lines'),rot=90)
  grid.text('Number of floors (number of blocks = 5)',y=unit(-2.8,'lines'))
  grid.text('Number of men per floor',y=unit(1,'npc')+unit(2.8,'lines'))
  popViewport()
  popViewport()
  popViewport()
  dev.off()
}

