#simulation functions-------------
#main
sim=function(n.inf){
  #simulate infections
  n.inf_1=n.inf*3
  t.inf=runif(n.inf_1, -29.99,0)
  incub.inf=rlnorm(n.inf_1, incub_paras[1], incub_paras[2])
  inf.inf=rlnorm(n.inf_1, inf_paras[1], inf_paras[2])
  infectious.end=t.inf+incub.inf+inf.inf
  inf.travel=which(infectious.end>0)[1:n.inf]
  #n.inf=length(inf.travel)
  
  t.inf=t.inf[inf.travel]; incub.inf=incub.inf[inf.travel]
  inf.inf=inf.inf[inf.travel]; infectious.end=infectious.end[inf.travel]
  test.res=function(test.t){
    test.res1=predict(pcr1,newdata=data.frame(days_pcr=test.t), type='response')
    test.res1[(test.t<=0)|(test.t>inf.inf)]=0
    test.res1[is.na(test.res1)]=0
    test.res2=predict(pcr2,newdata=data.frame(days_pcr=test.t), type='response')
    test.res2[(test.t<=0)|(test.t>inf.inf)]=0
    test.res2[is.na(test.res2)]=0
    test.res=test.res1*symptom+test.res2*(1-symptom)
    test.res[test.res>1]=1
    rbinom(n.inf,1,test.res)
  }
  test1.t=-3-t.inf-incub.inf
  test2.t=0-t.inf-incub.inf
  test3.t=7-t.inf-incub.inf
  
  #simulate test
  symptom=rbinom(n.inf,1,p.sym)
  test1.res=test.res(test1.t)
  test2.res=test.res(test2.t)
  test3.res=test.res(test3.t)
  
  #output
  missed=c(
    sum(infectious.end>0), #S1
    sum(infectious.end>3 & test2.res==0), #S2
    sum(infectious.end>0 & test1.res==0), #S3
    sum(infectious.end>3 & test1.res==0 & test2.res==0), #S4
    sum(infectious.end>10 & test1.res==0 & test3.res==0), #S5
    sum(infectious.end>7), #S6
    sum(infectious.end>14), #S7
    sum(infectious.end>21), #S8
    sum(infectious.end>28) #S9
  )
  
  recovered=c(
    sum(infectious.end<=0), #S1
    sum(infectious.end>0 & infectious.end<=3 & test2.res==0), #S2
    sum(infectious.end<=0), #S3
    sum(infectious.end>0 & infectious.end<=3 & test1.res==0 & test2.res==0),
    sum(infectious.end>0 & infectious.end<=10 & test1.res==0 & test3.res==0), #S5
    sum(infectious.end<=7), #S6
    sum(infectious.end<=14), #S7
    sum(infectious.end<=21), #S8
    sum(infectious.end<=28) #S9
  )
  
  pre.dep=c(
    rep(0,2), #S1&2
    rep(sum(test1.res==1),3), #S3-5
    rep(0,4) #S6-10
  )
  
  on.arr=c(
    0, sum(test2.res==1&infectious.end>0),0,
    sum(test1.res==0&test2.res==1&infectious.end>0),rep(0,5)
  )
  
  post.qua=c(rep(0,4),sum(test1.res==0&test3.res==1&infectious.end>7),rep(0,4))
  
  out=rbind(recovered,pre.dep,on.arr,post.qua,missed)
  out
}

#supplementary analysis: quarantine scenarios (pre-departure, on-arrival, both)
sim_add=function(n.inf){
  #simulate infections
  n.inf_1=n.inf*3
  t.inf=runif(n.inf_1, -29.99,0)
  incub.inf=rlnorm(n.inf_1, incub_paras[1], incub_paras[2])
  inf.inf=rlnorm(n.inf_1, inf_paras[1], inf_paras[2])
  infectious.end=t.inf+incub.inf+inf.inf
  inf.travel=which(infectious.end>0)[1:n.inf]
  #n.inf=length(inf.travel)
  
  t.inf=t.inf[inf.travel]; incub.inf=incub.inf[inf.travel]
  inf.inf=inf.inf[inf.travel]; infectious.end=infectious.end[inf.travel]
  test.res=function(test.t){
    test.res1=predict(pcr1,newdata=data.frame(days_pcr=test.t), type='response')
    test.res1[(test.t<=0)|(test.t>inf.inf)]=0
    test.res1[is.na(test.res1)]=0
    test.res2=predict(pcr2,newdata=data.frame(days_pcr=test.t), type='response')
    test.res2[(test.t<=0)|(test.t>inf.inf)]=0
    test.res2[is.na(test.res2)]=0
    test.res=test.res1*symptom+test.res2*(1-symptom)
    test.res[test.res>1]=1
    rbinom(n.inf,1,test.res)
  }
  test1.t=-3-t.inf-incub.inf
  test2.t=0-t.inf-incub.inf
  
  #simulate test
  symptom=rbinom(n.inf,1,p.sym)
  test1.res=test.res(test1.t)
  test2.res=test.res(test2.t)
  
  #output
  n.quar=4
  t.quar=(1:n.quar)*7
  t.test=3
  test3.res=lapply(1:n.quar,function(t){
    test3.t=t.quar[t]-t.inf-incub.inf
    test.res(test3.t)
  })
  
  
  missed=c(
    #pre-departure
    unlist(lapply(1:n.quar, function(i){
      sum(infectious.end>t.quar[i] & test1.res==0)
    })),
    #post-quarantine
    unlist(lapply(1:n.quar, function(i){
      sum(infectious.end>(t.quar[i]+t.test) & test3.res[[i]]==0)
    })),
    #both
    unlist(lapply(1:n.quar, function(i){
      sum(infectious.end>(t.quar[i]+t.test) & test1.res==0 & test3.res[[i]]==0)
    }))
  )
  
  recovered=c(
    unlist(lapply(1:n.quar, function(i){
      sum(infectious.end>0 & infectious.end<=t.quar[i] & test1.res==0)
    })),
    unlist(lapply(1:n.quar, function(i){
      sum(infectious.end>0 & infectious.end<=(t.quar[i]+t.test) & test3.res[[i]]==0)
    })),
    unlist(lapply(1:n.quar, function(i){
      sum(infectious.end>0 & infectious.end<=(t.quar[i]+t.test) & test1.res==0 & test3.res[[i]]==0)
    }))
  )
  
  pre.dep=c(
    rep(sum(test1.res==1),n.quar), 
    rep(0,n.quar),
    rep(sum(test1.res==1),n.quar)
  )
  
  on.arr=rep(0,n.quar*3)
  
  post.qua=c(rep(0,4),
             unlist(lapply(1:n.quar, function(i){
               sum(test3.res[[i]]==1)
             })),
             unlist(lapply(1:n.quar, function(i){
               sum(test1.res==0&test3.res[[i]]==1)
             })) )
  
  
  out=rbind(recovered,pre.dep,on.arr,post.qua,missed)
  out
}


#sensitivity analysis: infectious before symptom onset
sim_preinf=function(n.inf, preinf=F){
  #simulate infections
  n.inf_1=n.inf*2
  t.inf=runif(n.inf_1, -29.99,0)
  incub.inf=rlnorm(n.inf_1, incub_paras[1], incub_paras[2])
  inf.inf=rlnorm(n.inf_1, inf_paras[1], inf_paras[2])
  infectious.end=t.inf+incub.inf+inf.inf
  inf.travel=which(infectious.end>0)[1:n.inf]
  #n.inf=length(inf.travel)
  
  t.inf=t.inf[inf.travel]; incub.inf=incub.inf[inf.travel]
  inf.inf=inf.inf[inf.travel]; infectious.end=infectious.end[inf.travel]
  t0=rep(0,n.inf)
  if(preinf){
    t0=runif(n.inf, -4,0)*rbinom(n.inf,1,0.8)*(incub.inf>4)
  }
  
  test.res=function(test.t){
    test.res1=predict(pcr1,newdata=data.frame(days_pcr=test.t), type='response')
    test.res1[(test.t<=0)|(test.t>inf.inf)]=0
    test.res2=predict(pcr2,newdata=data.frame(days_pcr=test.t), type='response')
    test.res3=rep(0,n.inf)
    test.res3[(test.t<=0)&(test.t>t0)]=test.res2[(test.t<=0)&(test.t>t0)]
    test.res2[(test.t<=0)|(test.t>inf.inf)]=0
    test.res=test.res1*symptom+test.res2*(1-symptom)+test.res3
    rbinom(n.inf,1,test.res)
  }
  test1.t=-3-t.inf-incub.inf
  test2.t=0-t.inf-incub.inf
  test3.t=7-t.inf-incub.inf
  
  #simulate test
  symptom=rbinom(n.inf,1,p.sym)
  test1.res=test.res(test1.t)
  test2.res=test.res(test2.t)
  test3.res=test.res(test3.t)
  
  #output
  missed=c(
    sum(infectious.end>0), #S1
    sum(infectious.end>3 & test2.res==0), #S2
    sum(infectious.end>0 & test1.res==0), #S3
    sum(infectious.end>3 & test1.res==0 & test2.res==0), #S4
    sum(infectious.end>10 & test1.res==0 & test3.res==0), #S5
    sum(infectious.end>7), #S6
    sum(infectious.end>14), #S7
    sum(infectious.end>21), #S8
    sum(infectious.end>28) #S9
  )
  
  recovered=c(
    sum(infectious.end<=0), #S1
    sum(infectious.end>0 & infectious.end<=3 & test2.res==0), #S2
    sum(infectious.end<=0), #S3
    sum(infectious.end>0 & infectious.end<=3 & test1.res==0 & test2.res==0),
    sum(infectious.end>0 & infectious.end<=10 & test1.res==0 & test3.res==0), #S5
    sum(infectious.end<=7), #S6
    sum(infectious.end<=14), #S7
    sum(infectious.end<=21), #S8
    sum(infectious.end<=28) #S9
  )
  
  pre.dep=c(
    rep(0,2), #S1&2
    rep(sum(test1.res==1),3), #S3-5
    rep(0,4) #S6-10
  )
  
  on.arr=c(
    0, sum(test2.res==1&infectious.end>0),0,
    sum(test1.res==0&test2.res==1&infectious.end>0),rep(0,5)
  )
  
  post.qua=c(rep(0,4),sum(test1.res==0&test3.res==1&infectious.end>7),rep(0,4))
  
  out=rbind(recovered,pre.dep,on.arr,post.qua,missed)
  out
}

#sensitivity analysis: identify asymptomatic cases
sim_asym=function(n.inf){
  #simulate infections
  n.inf_1=n.inf*3
  t.inf=runif(n.inf_1, -29.99,0)
  incub.inf=rlnorm(n.inf_1, incub_paras[1], incub_paras[2])
  inf.inf=rlnorm(n.inf_1, inf_paras[1], inf_paras[2])
  infectious.end=t.inf+incub.inf+inf.inf
  inf.travel=which(infectious.end>0)[1:n.inf]
  #n.inf=length(inf.travel)
  
  t.inf=t.inf[inf.travel]; incub.inf=incub.inf[inf.travel]
  inf.inf=inf.inf[inf.travel]; infectious.end=infectious.end[inf.travel]
  test.res=function(test.t){
    test.res1=predict(pcr1,newdata=data.frame(days_pcr=test.t), type='response')
    test.res1[(test.t<=0)|(test.t>inf.inf)]=0
    test.res1[is.na(test.res1)]=0
    test.res2=predict(pcr2,newdata=data.frame(days_pcr=test.t), type='response')
    test.res2[(test.t<=0)|(test.t>inf.inf)]=0
    test.res2[is.na(test.res2)]=0
    test.res=test.res1*symptom+test.res2*(1-symptom)
    test.res[test.res>1]=1
    rbinom(n.inf,1,test.res)
  }
  test1.t=-3-t.inf-incub.inf
  test2.t=0-t.inf-incub.inf
  test3.t=7-t.inf-incub.inf
  
  #simulate test
  symptom=rbinom(n.inf,1,p.sym)
  #asymptomatic infections
  asymptom=rep(0,n.inf); asymptom[symptom==0]=rbinom(sum(symptom==0),1,0.1/(1-p.sym))
  test1.res=test.res(test1.t)
  test2.res=test.res(test2.t)
  test3.res=test.res(test3.t)
  
  #output
  missed=c(
    sum(infectious.end>0 & asymptom==0), #S1
    sum(infectious.end>3 & test2.res==0 & asymptom==0), #S2
    sum(infectious.end>0 & test1.res==0 & asymptom==0), #S3
    sum(infectious.end>3 & test1.res==0 & test2.res==0 & asymptom==0), #S4
    sum(infectious.end>10 & test1.res==0 & test3.res==0 & asymptom==0), #S5
    sum(infectious.end>7 & asymptom==0), #S6
    sum(infectious.end>14 & asymptom==0), #S7
    sum(infectious.end>21 & asymptom==0), #S8
    sum(infectious.end>28 & asymptom==0) #S9
  )
  
  asymptom=c(
    sum(infectious.end>0 & asymptom==1), #S1
    sum(infectious.end>3 & test2.res==0 & asymptom==1), #S2
    sum(infectious.end>0 & test1.res==0 & asymptom==1), #S3
    sum(infectious.end>3 & test1.res==0 & test2.res==0 & asymptom==1), #S4
    sum(infectious.end>10 & test1.res==0 & test3.res==0 & asymptom==1), #S5
    sum(infectious.end>7 & asymptom==1), #S6
    sum(infectious.end>14 & asymptom==1), #S7
    sum(infectious.end>21 & asymptom==1), #S8
    sum(infectious.end>28 & asymptom==1) #S9
  )
  
  recovered=c(
    sum(infectious.end<=0), #S1
    sum(infectious.end>0 & infectious.end<=3 & test2.res==0), #S2
    sum(infectious.end<=0), #S3
    sum(infectious.end>0 & infectious.end<=3 & test1.res==0 & test2.res==0),
    sum(infectious.end>0 & infectious.end<=10 & test1.res==0 & test3.res==0), #S5
    sum(infectious.end<=7), #S6
    sum(infectious.end<=14), #S7
    sum(infectious.end<=21), #S8
    sum(infectious.end<=28) #S9
  )
  
  pre.dep=c(
    rep(0,2), #S1&2
    rep(sum(test1.res==1),3), #S3-5
    rep(0,4) #S6-10
  )
  
  on.arr=c(
    0, sum(test2.res==1&infectious.end>0),0,
    sum(test1.res==0&test2.res==1&infectious.end>0),rep(0,5)
  )
  
  post.qua=c(rep(0,4),sum(test1.res==0&test3.res==1&infectious.end>7),rep(0,4))
  
  out=rbind(recovered,pre.dep,on.arr,post.qua,asymptom,missed)
  out
}

#sensitivity analysis: home isolation
sim_home=function(n.inf, n.det=0.9){
  #simulate infections
  n.inf_1=n.inf*3
  t.inf=runif(n.inf_1, -29.99,0)
  incub.inf=rlnorm(n.inf_1, incub_paras[1], incub_paras[2])
  inf.inf=rlnorm(n.inf_1, inf_paras[1], inf_paras[2])
  infectious.end=t.inf+incub.inf+inf.inf
  inf.travel=which(infectious.end>0)[1:n.inf]
  #n.inf=length(inf.travel)
  
  t.inf=t.inf[inf.travel]; incub.inf=incub.inf[inf.travel]
  inf.inf=inf.inf[inf.travel]; infectious.end=infectious.end[inf.travel]
  
  
  #output
  detect=rbinom(n.inf, 1,n.det)
  p.home=seq(0,0.9,0.3)
  n.quar=4;t.quar=1:n.quar*7
  missed=lapply(p.home, function(p){
    home=rep(0,n.inf)
    home[detect==F]=rbinom(sum(detect==F),1,p)
    lapply(t.quar, function(t){
      sum(infectious.end>t & (detect==T|(home==T&detect==F)))+
        sum(infectious.end>0 & home==F & detect==F)
    })%>%unlist
  })%>%do.call('rbind',.)%>%as.matrix()
  missed
}

#summary statistics--------------
#count
sum_num=function(res){
  n.vars=nrow(res[[1]])
  lapply(1:n.vars, function(i){ #five summary statistics
    count=as.matrix(do.call('rbind',lapply(1:n.sim, function(m){
      res[[m]][i,]
    })))
    as.matrix(do.call('rbind',lapply(1:ncol(count), function(s){
      dat=count[,s]%>%na.omit()
      c(mean(dat), quantile(dat,c(0.5,0.025,0.975,0.25,0.75)))
    })))
  })
}

#proportion
sum_prop=function(res){
  n.vars=nrow(res[[1]])
  lapply(1:n.vars, function(i){ #five summary statistics
    count=as.matrix(do.call('rbind',lapply(1:n.sim, function(m){
      res[[m]][i,]
    })))
    as.matrix(do.call('rbind',lapply(1:ncol(count), function(s){
      dat=count[,s]
      c(mean(dat), quantile(dat,c(0.5,0.025,0.975,0.25,0.75)))/n.inf*100
    })))
  })
}




