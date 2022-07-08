#for paper
#loading the basic data
{
  library(data.table)
  library(grid)
  library(scales)
  data=as.data.frame(fread('data.csv'))
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
  start=-28; end=60
  
  stat_cal=function(m){
    c(mean(m),quantile(m, c(0.025,0.975)))
  }
}

#Figure 1
#infection harzard
{
  data$t_inf=s_inf[sample(1:nrow(s_inf),1), ]
  k=sample(seq(1,nrow(data)),1)
  start0=-14
  inf=matrix(0, 3, end-start0+1)
  para=c(7.6e-5, 1.4e-4, 0.012)
  for(i in start0:end)
  {
    for(j in 1:3)
    {
      if(j==1){
        inf[j,i-start0+1]=length(which(data$t_inf<=i&data$t_rem>i&distance[k,]==j))*para[j]
      }else{
        inf[j,i-start0+1]=inf[j-1,i-start0+1]+length(which(data$t_inf<=i&data$t_rem>i&distance[k,]==j))*para[j]
      }
    }
  }
  text=c('Block','Level','Room')
  y_place=c(max(inf[1,])*0.4, max(inf[2,])*0.8, max(inf[3,])*0.5)
  png('hazard.png',height=8,width=8.5,units='cm', res=300, pointsize=10) 
  pushViewport(plotViewport(c(4,5,1,1), xscale=c(start0,end),yscale=c(0,max(inf)*1.05)))
  grid.rect()
  grid.text('(b)', x=0.08, y=0.95)
  for (i in 1:3)
  {
    grid.polygon(x=c(seq(start0,end),rev(c(start0,end))), y=c(inf[i,],0,0), default.units='native', gp = gpar(fill = scales::alpha('darkorange',0.8-i*0.2),col=scales::alpha('tomato',0.1) , alpha = 1))
  }
  for (i in 1:3)
  {
    grid.text(text[i], x=which.max(inf[i,])+start0, y=y_place[i], default.units='native', gp=gpar(col='navyblue', fontsize=10))
  }
  grid.yaxis(at=seq(0,max(inf),0.02))
  grid.xaxis(at=seq(ceiling(start0/14)*14,end,14))
  grid.text('Example total',x=unit(-4.5,'lines'),rot=90)
  grid.text('infection hazard (/d)',x=unit(-3.5,'lines'),rot=90)
  grid.text('Time of study (d)',y=unit(-3,'lines'))
  popViewport()
  dev.off()
}


#seroconversion
{
  png('serocoversion.png',height=8,width=8,units='cm', res=300, pointsize=10) 
  pushViewport(plotViewport(c(4,4,1,1), xscale=c(0,35),yscale=c(0,100)))
  grid.rect()
  grid.text('(d)', x=0.08, y=0.95)
  grid.lines(seq(0,35,0.5), plnorm(seq(0,35,0.5), 2.75,0.62)*100, default.units='native', gp=gpar(col=scales::alpha('firebrick3',0.6), lwd=2, cex=3))
  grid.yaxis()
  grid.xaxis(at=seq(0,35,7))
  grid.text('Fraction seroconverted (%)',x=unit(-3,'lines'),rot=90)
  grid.text('Time since infection (d)',y=unit(-3,'lines'))
  popViewport()
  dev.off()
}


#incubation distribution
{
  png('incubation_distribution.png',height=8,width=8,units='cm', res=300, pointsize=10) 
  pushViewport(plotViewport(c(4,4,1,1), xscale=c(0,14),yscale=c(0,0.4)))
  grid.rect()
  grid.text('(c)', x=0.08, y=0.95)
  grid.lines(seq(0,14,0.1), dlnorm(seq(0,14,0.1), 1,sqrt(0.3)), default.units='native', gp=gpar(col=scales::alpha('navyblue',0.6), lwd=2, cex=3))
  grid.yaxis()
  grid.xaxis(at=seq(0,14,7))
  grid.text('Probability density',x=unit(-3,'lines'),rot=90)
  grid.text('Time of onset after infection (d)',y=unit(-3,'lines'))
  popViewport()
  dev.off()
}


#infection data
{
  s_inf=as.matrix(fread(paste0('~/OneDrive - National University of Singapore/dorm/0_v2/2_1e4/',8419,'_inf_result.csv')))
  s_inf=rbind(s_inf,as.matrix(fread(paste0('~/OneDrive - National University of Singapore/dorm/0_v2/1e4/',7955,'_inf_result.csv'))))
  s_inf=rbind(s_inf,as.matrix(fread(paste0('~/OneDrive - National University of Singapore/dorm/0_v2/1e4/',9473,'_inf_result.csv'))))
}



#Figure 2 statistics
{
  infection_t=function(n,t){
    as.vector(do.call('sum', lapply(1:5, function(k) sum(s_inf[[k]][n,1:(t-start+1)]))))
  }
  infectious_t=function(n,t){
    as.matrix(do.call('sum', lapply(1:5, function(k) sum(s_inf[[k]][n,(t-start-8):(t-start+1)]))))
  }
  infection_all=function(n){
    as.matrix(do.call('cbind', lapply(start:end, function(t) infection_t(n,t) )))
  }
  infectious_all=function(n){
    as.matrix(do.call('cbind', lapply(start:end, function(t) infectious_t(n,t) )))
  }
  infection_t_k=function(n,k){
    as.matrix(do.call('cbind', lapply(start:end, function(t) sum(s_inf[[(k-3)/2]][n,1:(t-start+1)]))))
  }
  infectious_t_k=function(n,k){
    as.matrix(do.call('cbind', lapply(start:end, function(t) sum(s_inf[[(k-3)/2]][n,(t-start-8):(t-start+1)]))))
  }
  sero_k=function(k){
    temp=as.matrix(do.call('rbind', lapply(1:nrow(s_inf[[1]]), function(n) infection_t_k(n,k))))
    as.matrix(do.call("cbind", lapply(1:ncol(temp), function(i) stat_cal(temp[,i]))))
  }
  inf_k=function(k){
    temp=as.matrix(do.call('rbind', lapply(1:nrow(s_inf[[1]]), function(n) infectious_t_k(n,k))))
    as.matrix(do.call("cbind", lapply(1:ncol(temp), function(i) stat_cal(temp[,i]))))
  }
  conv_k=function(k){
    as.matrix(do.call('rbind', lapply(1:nrow(s_inf[[1]]), function(n) sum(s_inf[[k]][n,(1-start):(end-start+1)]) )))
  }
  
}
#overall stat
temp=as.matrix(do.call('rbind', lapply(1:nrow(s_inf[[1]]), function(n) infection_all(n))))
sero=as.matrix(do.call("cbind", lapply(1:ncol(temp), function(i) stat_cal(temp[,i]))))
temp=as.matrix(do.call('rbind', lapply(1:nrow(s_inf[[1]]), function(n) infectious_all(n))))
inf=as.matrix(do.call("cbind", lapply(1:ncol(temp), function(i) stat_cal(temp[,i]))))
remove(temp)
#block stat
block_count=rep(0,5)
for (i in 1:5)
{
  block_count[i]=sum(data$block==2*i+3)
}
num_sero=lapply(1:5*2+3, function(k) sero_k(k))
num_inf=lapply(1:5*2+3, function(k) inf_k(k))

#seroconvergence 
stat_cal(as.vector(do.call('sum', lapply(1:5, function(k) conv_k(k) )))/nrow(data))
#for each block
as.matrix(do.call('cbind', lapply(1:5, function(k) stat_cal(conv_k(k))/block_count[k] )))
#sero-prevalence at a certain time point
as.matrix(do.call('rbind', lapply(1:5, function(k) num_sero[[k]][,29]/block_count[k])))
#maximum infectious day
as.matrix(do.call('rbind', lapply(1:5, function(i) c(which.max(num_inf[[i]][1,])-29,round(num_inf[[i]][,which.max(num_inf[[i]][1,])]*100/block_count[i],4)) )))

#plot
{
  colour = list(navy = '#00306b',navy2 = '#406390',navy3 = '#7e96b5',navy4 = '#c0cbda',red = '#e01f33',red2 = '#ea5766',red3 = '#f18f99',red4 = '#f8c7cb',blue = '#0072bc',blue2 = '#3f96cd',blue3 = '#80b9de',blue4 = '#c0dcef',orange = '#ee6f01',orange2 = '#f2943f',orange3 = '#f7b780',orange4 = '#fcdac0',sky = '#10a9e0',sky2 = '#4bbee8',sky3 = '#88d4f0',sky4 = '#c3eaf8',pink = '#970e53',pink2 = '#d31876',pink3 = '#ff42a1',pink4 = '#ff95ca',gray = '#536778',gray2 = '#8195a7',gray3 = '#d8dee7')
  clr=c(colour$navy, colour$red, colour$blue, colour$orange, colour$pink,'black')
  legend=c('A','B','C','D','E')
  png('Figure2.png',height=8,width=16,units='cm', res=300, pointsize=8) 
  pushViewport(plotViewport(c(1,1,1,1)))
  pushViewport(viewport(layout=grid.layout(nrow=1,ncol=2))) 
  pushViewport(viewport(layout.pos.col=1,layout.pos.row=1)) 
  pushViewport(plotViewport(c(4,3,1,4), xscale=c(start,end),yscale=c(0,90)))
  grid.rect()
  grid.text('(a)',x=0.05,y=0.95,just='left')
  grid.lines(seq(start,end), sero[1,]/nrow(data)*100, default.units='native', gp=gpar(col=scales::alpha(clr[6],1), cex=0.5))
  grid.polygon(x=c(seq(start,end),rev(seq(start,end))), y=c(sero[2,]/nrow(data)*100,rev(sero[3,]/nrow(data)*100)), default.units='native', gp = gpar(fill = scales::alpha(clr[6],0.5),col=scales::alpha(clr[6],0.1) , alpha = 1))
  grid.rect(x=0.05, y=unit(1,'npc')-unit(3,'lines'), width=0.05, height=0.005, default.units="npc", gp=gpar(fill=scales::alpha(clr[6],0.8), col=scales::alpha(clr[6],0.01)))
  grid.rect(x=0.05, y=unit(1,'npc')-unit(3,'lines'), width=0.05, height=0.02, default.units="npc", gp=gpar(fill=scales::alpha(clr[6],0.5), col=scales::alpha(clr[6],0.01)))
  grid.text('All', x=0.1, y=unit(1,'npc')-unit(3,'lines'), just='left')
  for (i in 1:5)
  {
    grid.lines(start:end, num_sero[[i]][1,]/block_count[i]*100, default.units='native', gp=gpar(col=scales::alpha(clr[i],0.25), cex=0.5))
    grid.polygon(x=c(seq(start,end),rev(seq(start,end))), y=c(num_sero[[i]][2,]/block_count[i]*100,rev(num_sero[[i]][3,]/block_count[i]*100)), default.units='native', gp = gpar(fill = scales::alpha(clr[i],0.1),col=scales::alpha(clr[i],0.01) , alpha = 1))
    #grid.circle(x=0.05, y=unit(1,'npc')-unit(i+3,'lines'), r=0.02, default.units="npc", gp=gpar(fill=clr[i], col=clr[i]))
    grid.rect(x=0.05, y=unit(1,'npc')-unit(i+3,'lines'), width=0.05, height=0.005, default.units="npc", gp=gpar(fill = scales::alpha(clr[i],0.8),col=scales::alpha(clr[i],0.01)))
    grid.rect(x=0.05, y=unit(1,'npc')-unit(i+3,'lines'), width=0.05, height=0.02, default.units="npc", gp=gpar(fill = scales::alpha(clr[i],0.5),col=scales::alpha(clr[i],0.01)))
    grid.text(paste0('Block ',legend[i]), x=0.1, y=unit(1,'npc')-unit(i+3,'lines'), just='left')
    
  }
  grid.yaxis(at=seq(0,3500/4257*100,500/4257*100), label=seq(0, 3500, 500), main=FALSE)
  grid.yaxis(at=seq(0,90,10))
  grid.xaxis(at=seq(ceiling(start/14)*14,end,14),label=seq(-4,end/7,2))
  grid.text('Prevalence (%)',x=unit(-2.8,'lines'),rot=90)
  grid.text('Prevalence (n)',x=unit(1,'npc')-unit(-3.2,'lines'),rot=90)
  grid.text('Time of study (weeks)',y=unit(-3,'lines'))
  popViewport()
  popViewport()
  pushViewport(viewport(layout.pos.col=2,layout.pos.row=1)) 
  pushViewport(plotViewport(c(4,3.5,1,3.5), xscale=c(start,end),yscale=c(0,33)))
  grid.rect()
  grid.text('(b)',x=0.05,y=0.95,just='left')
  grid.lines(seq(start,end), inf[1,]/nrow(data)*100, default.units='native', gp=gpar(col=scales::alpha(clr[6],1), cex=0.5))
  grid.polygon(x=c(seq(start,end),rev(seq(start,end))), y=c(inf[2,]/nrow(data)*100,rev(inf[3,]/nrow(data)*100)), default.units='native', gp = gpar(fill = scales::alpha(clr[6],0.5),col=scales::alpha(clr[6],0.1) , alpha = 1))
  grid.rect(x=0.05, y=unit(1,'npc')-unit(3,'lines'), width=0.05, height=0.005, default.units="npc", gp=gpar(fill=scales::alpha(clr[6],0.8), col=scales::alpha(clr[6],0.01)))
  grid.rect(x=0.05, y=unit(1,'npc')-unit(3,'lines'), width=0.05, height=0.02, default.units="npc", gp=gpar(fill=scales::alpha(clr[6],0.5), col=scales::alpha(clr[6],0.01)))
  grid.text('All', x=0.1, y=unit(1,'npc')-unit(3,'lines'), just='left')
  for (i in 1:5)
  {
    grid.lines(start:end, num_inf[[i]][1,]/block_count[i]*100, default.units='native', gp=gpar(col=scales::alpha(clr[i],0.25), cex=0.5))
    grid.polygon(x=c(seq(start,end),rev(seq(start,end))), y=c(num_inf[[i]][2,]/block_count[i]*100,rev(num_inf[[i]][3,]/block_count[i]*100)), default.units='native', gp = gpar(fill = scales::alpha(clr[i],0.1),col=scales::alpha(clr[i],0.01) , alpha = 1))
    grid.rect(x=0.05, y=unit(1,'npc')-unit(i+3,'lines'), width=0.05, height=0.005, default.units="npc", gp=gpar(fill = scales::alpha(clr[i],0.8),col=scales::alpha(clr[i],0.01)))
    grid.rect(x=0.05, y=unit(1,'npc')-unit(i+3,'lines'), width=0.05, height=0.02, default.units="npc", gp=gpar(fill = scales::alpha(clr[i],0.5),col=scales::alpha(clr[i],0.01)))
    grid.text(paste0('Block ',legend[i]), x=0.1, y=unit(1,'npc')-unit(i+3,'lines'), just='left')
  }
  grid.yaxis(at=seq(0,1250/4257*100,250/4257*100), label=seq(0, 1250, 250), main=FALSE)
  grid.yaxis(at=seq(0,30,5))
  grid.xaxis(at=seq(ceiling(start/14)*14,end,14),label=seq(-4,end/7,2))
  grid.text('Number infectious (n)',x=unit(1,'npc')-unit(-3.2,'lines'),rot=90)
  grid.text('Proportion infectious (%)',x=unit(-2.5,'lines'),rot=90)
  grid.text('Time of study (weeks)',y=unit(-3,'lines'))
  popViewport()
  popViewport()
  popViewport()
  dev.off()
}



#trace-plot for model without diverse treatment effects
{
  MCMCiterations=1e4
  text=c(expression(a[1]), expression(a[2]),expression(a[3]),expression(theta),'incublation length',expression(beta[block]),expression(beta[floor]),expression(beta[unit]))
  ylab1=c(2.1,0.65,1.5,0.08,3,5e-5,0,0.009)
  ylab2=c(2.5,0.95,2.1,0.16,6,1e-4,4e-4,0.014)
  clr=c("indianred3", "darkorange", "thistle4","darkseagreen", "purple", "slateblue3",'navyblue','lightpink3', 'mediumpurple4','tan2','palevioletred4')
  png('MCMC_traceplot_1.png',height=16,width=32,units='cm', res=300, pointsize=10) 
  pushViewport(plotViewport(c(1,1,1,1)))
  pushViewport(viewport(layout=grid.layout(nrow=2,ncol=4)))
  j=0
  k_total=floor(nrow(store)/MCMCiterations)
  for(i in c(1:4,6:8))
  {
    j=j+1
    pushViewport(viewport(layout.pos.col=j-(ceiling(j/4)-1)*4,layout.pos.row=ceiling(j/4))) 
    pushViewport(plotViewport(c(4,4,1,1), xscale=c(0,MCMCiterations),yscale=c(ylab1[i],ylab2[i])))
    grid.rect()
    for(k in 1:k_total)
    {
      idx=(k-1)*MCMCiterations+1:MCMCiterations
      grid.lines(seq(1,MCMCiterations), store[idx,i], default.units='native', gp=gpar(col=scales::alpha(clr[i],0.5), cex=3))
    }
    grid.lines(c(1,MCMCiterations), c(mean(store[,i]),mean(store[,i])), default.units='native', gp=gpar(col=scales::alpha('tomato',0.6), cex=3))
    grid.rect(x=0.05, y=unit(1,'npc')-unit(1,'lines'), width=0.03, height=0.025, default.units="npc", gp=gpar(fill=clr[i], col=scales::alpha(clr[i],0.6), fontsize=8))  
    grid.text(text[i], x=0.09, y=unit(1,'npc')-unit(1,'lines'), just='left')
    grid.yaxis()
    grid.xaxis()
    popViewport()
    popViewport()
  }
  popViewport()
  dev.off()
}


#trace-plot for model with diverse treatment effects
{
  MCMCiterations=1e4
  text=c(expression(a[1]), expression(a[2]),expression(a[3]),expression(theta),'incublation length',expression(beta[block]),expression(beta[floor]),expression(beta[unit]), expression(gamma[B]),expression(gamma[C]), expression(gamma[D]),expression(gamma[E]))
  ylab1=c(2.2,0.7,1.6,0.08,3,5e-5,0,0.009,rep(0.7,4))
  ylab2=c(2.5,0.9,2.0,0.16,6,1e-4,4e-4,0.015,rep(1.3,4))
  clr=c("indianred3", "darkorange", "thistle4","darkseagreen", "purple", "slateblue3",'navyblue','lightpink3', 'mediumpurple4','tan2','palevioletred4','lightsalmon3','powderblue','navajowhite3','yellowgreen')
  png('MCMC_traceplot_2.png',height=24,width=32,units='cm', res=300, pointsize=10) 
  pushViewport(plotViewport(c(1,1,1,1)))
  pushViewport(viewport(layout=grid.layout(nrow=3,ncol=4)))
  j=0
  k_total=floor(nrow(store)/MCMCiterations)
  for(i in c(1:4,6:12))
  {
    j=j+1
    if(j<=7){
      pushViewport(viewport(layout.pos.col=j-(ceiling(j/4)-1)*4,layout.pos.row=ceiling(j/4)))
    }else{
      pushViewport(viewport(layout.pos.col=i-(ceiling(i/4)-1)*4,layout.pos.row=ceiling(i/4))) 
    }
    pushViewport(plotViewport(c(4,4,1,1), xscale=c(0,MCMCiterations),yscale=c(ylab1[i],ylab2[i])))
    grid.rect()
    for(k in 1:k_total)
    {
      idx=(k-1)*MCMCiterations+1:MCMCiterations
      grid.lines(seq(1,MCMCiterations), store[idx,i], default.units='native', gp=gpar(col=scales::alpha(clr[i],0.5), cex=3))
    }
    grid.lines(c(1,MCMCiterations), c(mean(store[,i]),mean(store[,i])), default.units='native', gp=gpar(col=scales::alpha('tomato',0.6), cex=3))
    
    grid.rect(x=0.05, y=unit(1,'npc')-unit(1,'lines'), width=0.035, height=0.025, default.units="npc", gp=gpar(fill=clr[i], col=scales::alpha(clr[i],0.6), fontsize=8))  
    grid.text(text[i], x=0.09, y=unit(1,'npc')-unit(1,'lines'), just='left')
    grid.yaxis()
    grid.xaxis()
    popViewport()
    popViewport()
  }
  popViewport()
  dev.off()
}


#convergence diagnostic test
library(coda)
for(i in 1:ncol(store))
{
  print(heidel.diag(store[,i]))
}


#trace-plot for cumulative infections
{
  t=0:7*7
  temp=as.matrix(do.call('rbind', lapply(1:nrow(s_inf[[1]]), function(n) infection_all(n))))
  s=as.matrix(do.call("cbind", lapply(t-start+1, function(i) stat_cal(temp[,i]))))
  library(viridis)
  clr=viridis_pal(option='D')(12)[-1]
  png('infection_traceplot.png',height=8,width=9,units='cm', res=300, pointsize=10) 
  pushViewport(plotViewport(c(4,4,1,4), xscale=c(0,MCMCiterations),yscale=c(0,3000)))
  grid.rect()
  for(i in c(1:6,8))
  {
    k_total=floor(nrow(s_inf[[1]])/MCMCiterations)
    for(k in 1:k_total){
      idx=(k-1)*MCMCiterations+1:MCMCiterations  
      grid.lines(seq(1,MCMCiterations), s[idx,i], default.units='native', gp=gpar(col=scales::alpha(clr[i],0.6), cex=3))
    } 
    grid.text(paste0('t = ',t[i]), x=unit(1,'npc')+unit(0.5,'lines'), y=unit(mean(s[,i]),'native'), just='left')
  }
  grid.yaxis(at=0:6*500)
  grid.xaxis()
  grid.text('cumulative infections', x=-unit(3.5,'lines'), rot=90)
  grid.text('Time (d)', x=unit(1,'npc')+unit(0.5,'lines'), y=unit(1,'npc')+unit(0.5,'lines'), just='left')
  grid.text('iteration', y=-unit(2.5,'lines'))
  popViewport()
  dev.off()
}







