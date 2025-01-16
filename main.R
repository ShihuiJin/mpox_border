library(grid)
library(dplyr)
library(viridis)
library(parallel)

#step 1: set parameter values----------
source('set_paras.R')

#step 2: load simulation functions-------------
load('paras.RData') #parameter: PCR, incubation, infectious
source('functions.R') #simulation and summary functions

#step 3: simulations-----------
n.cores=detectCores()
n.sim=1e4

#analysis 1: proportion of missed and detected cases
n.inf=1e4; n.vars=5
res=mclapply(1:n.sim, function(m) sim(n.inf),mc.cores = n.cores) #for figure3, all infections
sum_prop(res)

#analysis 2: number of missed cases per 100,000 travellers
n=1e6
prevs=c(0.01,0.05,0.1)/1000
res=lapply(prevs,function(prev){
  lapply(1:n.sim, function(m){
    n.inf=rbinom(1,n,prev)
    if(n.inf==0){
      rep(0,5)
    }else{
      sim(n.inf)
    }
  })
})
lapply(seq_along(prevs), function(j){
  sum_num(res[[j]])[[5]]/10
})

#supplementary: additional quarantine scenarios
n.s=12
res_add=mclapply(1:n.sim, function(m) sim_add(n.inf),mc.cores = n.cores)
sum_prop(res_add)

