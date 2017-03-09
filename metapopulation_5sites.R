#with ceiling-type density dependence and stochastic variation in annual survival rate 

library(MASS)


######### set up site effects and calculate probability of hatching a nest during the season #######
### functions to calculate probability of hatching at least one nest for exclosed/unexclosed ####
### this set up is a bit cumbersome but likely best way to build reactive decision tool; set these up once and then r can recalculate only when needed ####

h.un.calc<-function(alpha.a,alpha.p,alpha.f, eta.a,eta.p,r2,r3,m){
  lin.f<-exp(alpha.f)
  lin.a<-exp(alpha.a+eta.a)
  lin.p<-exp(alpha.p+eta.p)
  survp<-1/(lin.f+lin.a+lin.p+1)   #daily survival probability
  
  #period fate probabilities for each nest
  pred.p<-((lin.p/(lin.f+lin.a+lin.p+1))/(1-survp)*(1-survp^34))  #predation unexclosed
  a.p.un<-((lin.a/(lin.f+lin.a+lin.p+1))/(1-survp)*(1-survp^34))  #abandonment unexclosed
  flood.p<-((lin.f/(lin.f+lin.a+lin.p+1))/(1-survp)*(1-survp^34))  #flooding probability
  o.un<-flood.p+pred.p
  hatch.un<-survp^34  #probability of period survival unexclosed nest
  h.un<-hatch.un+r2*(o.un*hatch.un+a.p.un*(1-m)*hatch.un+o.un*(o.un*r3+a.p.un*(1-m)*r3*hatch.un)+
                       a.p.un*(1-m)*(o.un*hatch.un+a.p.un*(1-m)*r3*hatch.un))
  
  
}

h.ex.calc<-function(alpha.a,alpha.p,alpha.f, beta.a.ex, beta.p.ex, eta.a,eta.p,r2,r3,m){
  
  lin.f.ex<-exp(alpha.f)
  lin.a.ex<-exp(alpha.a + eta.a + beta.a.ex)
  lin.p.ex<-exp(alpha.p + eta.p + beta.p.ex)
  survp.ex<-1/(lin.f.ex+lin.a.ex+lin.p.ex+1)
  pred.p.ex<-((lin.p.ex/(lin.f.ex+lin.a.ex+lin.p.ex+1))/(1-survp.ex)*(1-survp.ex^34))
  a.p.ex<-((lin.a.ex/(lin.f.ex+lin.a.ex+lin.p.ex+1))/(1-survp.ex)*(1-survp.ex^34))
  flood.p.ex<-((lin.f.ex/(lin.f.ex+lin.a.ex+lin.p.ex+1))/(1-survp.ex)*(1-survp.ex^34))
  
  hatch.ex<-survp.ex^34   #probability of period survival exclosed nest
  o.ex<-flood.p.ex+pred.p.ex 
  h.ex<-hatch.ex+r2*(o.ex*hatch.ex+a.p.ex*(1-m)*hatch.ex+o.ex*(o.ex*r3+a.p.ex*(1-m)*r3*hatch.ex)+
                       a.p.ex*(1-m)*(o.ex*hatch.ex+a.p.ex*(1-m)*r3*hatch.ex))
}

aban.ex.calc<-function(alpha.a,alpha.p,alpha.f, beta.a.ex, beta.p.ex, eta.a,eta.p){
  
  lin.f.ex<-exp(alpha.f)
  lin.a.ex<-exp(alpha.a + eta.a + beta.a.ex)
  lin.p.ex<-exp(alpha.p + eta.p + beta.p.ex)
  survp.ex<-1/(lin.f.ex+lin.a.ex+lin.p.ex+1)
  a.p.ex<-((lin.a.ex/(lin.f.ex+lin.a.ex+lin.p.ex+1))/(1-survp.ex)*(1-survp.ex^34))
}
aban.un.calc<-function(alpha.a, alpha.p, alpha.f, eta.a, eta.p){
  
  lin.f<-exp(alpha.f)
  lin.a<-exp(alpha.a+eta.a)
  lin.p<-exp(alpha.p+eta.p)
  survp<-1/(lin.f+lin.a+lin.p+1)
  a.p.un<-((lin.a/(lin.f+lin.a+lin.p+1))/(1-survp)*(1-survp^34)) 
}


###################################### 

##set up:

#permutation matrix

E11<-matrix(0, nrow=2, ncol=5)
E12<-matrix(0, nrow=2, ncol=5)
E13<-matrix(0, nrow=2, ncol=5)
E14<-matrix(0, nrow=2, ncol=5)
E15<-matrix(0, nrow=2, ncol=5)
E21<-matrix(0, nrow=2, ncol=5)
E22<-matrix(0, nrow=2, ncol=5)
E23<-matrix(0, nrow=2, ncol=5)
E24<-matrix(0, nrow=2, ncol=5)
E25<-matrix(0, nrow=2, ncol=5)

E11[1,1]<-1
E12[1,2]<-1
E13[1,3]<-1
E14[1,4]<-1
E15[1,5]<-1
E21[2,1]<-1
E22[2,2]<-1
E23[2,3]<-1
E24[2,4]<-1
E25[2,5]<-1

P<-E11%x%t(E11) + E12%x%t(E12) + E13%x%t(E13) + E14%x%t(E14) + E15%x%t(E15) +
  E21%x%t(E21) + E22%x%t(E22) + E23%x%t(E23) + E24%x%t(E24) + E25%x%t(E25)


# site-independent variables  ##current values from run of "basic_approach_revise"; replace with final values after model finishes

r2<-0.7
r3<-0.7
m<-0.7   ###### should change to separate mortalities associated with exclosed-abandoned and unexclosed-abandoned ######
yt=0.99  
ys=0.68  
Phij=0.52  
Phia=0.74
f<-0.4
f.CV=0.06  #already on log-odds scale
CV=0.1
alpha.a<- -7.356
alpha.p<- -3.99
alpha.f<- -5.35
beta.a.ex<- 1.5
beta.p.ex<- -2.1
E<-0.94
rr.ad<-0.95  #adult site fidelity
rr.sy<-0.2  #second-year site fidelity

##baseline mortality risk: abandonment at unexclosed nests
a.p.un<-aban.un.calc(alpha.a, alpha.p, alpha.f, eta.a = 0, eta.p = 0)
m.a.un = yt*(a.p.un*m+a.p.un*(1-m)*r2*(a.p.un*m+a.p.un*(1-m)*r3*a.p.un*m))
m.s.un = ys*(a.p.un*m+a.p.un*(1-m)*r2*(a.p.un*m+a.p.un*(1-m)*r3*a.p.un*m))

###site-specific variables

#site 1 (=average sites) (=37/46 sites = 80% of population)
eta.p1 <- 0
eta.a1 <- 0

h.ex.1 <- h.ex.calc(alpha.a, alpha.p,alpha.f, beta.a.ex, beta.p.ex, eta.a1,eta.p1,r2,r3,m)
h.un.1 <- h.un.calc(alpha.a, alpha.p, alpha.f, eta.a1, eta.p1, r2, r3, m)
n.ex.1 <- 1
a.p.ex.1 <- aban.ex.calc(alpha.a,alpha.p,alpha.f, beta.a.ex, beta.p.ex, eta.a1, eta.p1)
a.p.un.1 <- aban.un.calc(alpha.a, alpha.p, alpha.f, eta.a1,eta.p1)

Nsy.site1.0<-0.8*2000/2*1.4   #start close to stable stage distribution - is this why there is a small jump followed by plateau or decline
Nad.site1.0<-0.8*2000/2
K.site1<-0.8*3000

#site 2 = high aban ave pred  (=4/46 = 8% of population)
eta.p2<-0
eta.a2<-1.5
h.ex.2 <- h.ex.calc(alpha.a, alpha.p,alpha.f, beta.a.ex, beta.p.ex, eta.a=eta.a2,eta.p=eta.p2,r2,r3,m)
h.un.2<-h.un.calc(alpha.a, alpha.p, alpha.f, eta.a=eta.a2, eta.p=eta.p2, r2, r3, m)
n.ex.2<-0.5
a.p.ex.2 <- aban.ex.calc(alpha.a,alpha.p,alpha.f, beta.a.ex, beta.p.ex, eta.a=eta.a2,eta.p=eta.a2)
a.p.un.2 <- aban.un.calc(alpha.a, alpha.p, alpha.f, eta.a=eta.a2, eta.p=eta.p2)

Nsy.site2.0<-0.08*2000/2
Nad.site2.0<-0.08*2000/2
K.site2<-0.08*3000

#site 3 = high pred ave aban  (=6% of population)
eta.p3<-0.5
eta.a3<-0
h.ex.3 <- h.ex.calc(alpha.a, alpha.p,alpha.f, beta.a.ex, beta.p.ex, eta.a=eta.a3,eta.p=eta.p3,r2,r3,m)
h.un.3<-h.un.calc(alpha.a, alpha.p, alpha.f, eta.a=eta.a3, eta.p=eta.p3, r2, r3, m)
n.ex.3<-1
a.p.ex.3 <- aban.ex.calc(alpha.a,alpha.p,alpha.f, beta.a.ex, beta.p.ex, eta.a=eta.a3,eta.p=eta.a3)
a.p.un.3 <- aban.un.calc(alpha.a, alpha.p, alpha.f, eta.a=eta.a3, eta.p=eta.p3)

Nsy.site3.0<-0.06*2000/2
Nad.site3.0<-0.06*2000/2
K.site3<-0.06*3000

#site 4 = ave pred low aban  (=4% of population)
eta.p4<-0
eta.a4<- -1.2
h.ex.4 <- h.ex.calc(alpha.a, alpha.p,alpha.f, beta.a.ex, beta.p.ex, eta.a=eta.a4,eta.p=eta.p4,r2,r3,m)
h.un.4<-h.un.calc(alpha.a, alpha.p, alpha.f, eta.a=eta.a4, eta.p=eta.p4, r2, r3, m)
n.ex.4<-1
a.p.ex.4 <- aban.ex.calc(alpha.a,alpha.p,alpha.f, beta.a.ex, beta.p.ex, eta.a=eta.a4,eta.p=eta.a4)
a.p.un.4 <- aban.un.calc(alpha.a, alpha.p, alpha.f, eta.a=eta.a4, eta.p=eta.p4)

Nsy.site4.0<-0.04*2000/2*1.4
Nad.site4.0<-0.04*2000/2
K.site4<-0.04*3000

#site 5 = focal site  
eta.p5<- -1.5
eta.a5<- -3
h.ex.5 <- h.ex.calc(alpha.a, alpha.p,alpha.f, beta.a.ex, beta.p.ex, eta.a=eta.a5,eta.p=eta.p5,r2,r3,m)
h.un.5<-h.un.calc(alpha.a, alpha.p, alpha.f, eta.a=eta.a5, eta.p=eta.p5, r2, r3, m)
n.ex.5<-1
a.p.ex.5 <- aban.ex.calc(alpha.a,alpha.p,alpha.f, beta.a.ex, beta.p.ex, eta.a=eta.a5,eta.p=eta.a5)
a.p.un.5 <- aban.un.calc(alpha.a, alpha.p, alpha.f, eta.a=eta.a5, eta.p=eta.p5)

Nsy.site5.0<-20
Nad.site5.0<-20
K.site5<-60

##relative site sizes to use for caluclating dispersal rates

P12<-(K.site2/sum(K.site2+K.site3+K.site4+K.site5))  #when dispersing from site 1 to 2, probability based on size of site 2 relative to all site choices (2:5)
P13<-(K.site3/sum(K.site2+K.site3+K.site4+K.site5))
P14<-(K.site4/sum(K.site2+K.site3+K.site4+K.site5))
P15<-(K.site5/sum(K.site2+K.site3+K.site4+K.site5))
P21<-(K.site1/sum(K.site1+K.site3+K.site4+K.site5))
P23<-(K.site3/sum(K.site1+K.site3+K.site4+K.site5))
P24<-(K.site4/sum(K.site1+K.site3+K.site4+K.site5))
P25<-(K.site5/sum(K.site1+K.site3+K.site4+K.site5))
P31<-(K.site1/sum(K.site1+K.site2+K.site4+K.site5))
P32<-(K.site2/sum(K.site1+K.site2+K.site4+K.site5))
P34<-(K.site4/sum(K.site1+K.site2+K.site4+K.site5))
P35<-(K.site5/sum(K.site1+K.site2+K.site4+K.site5))
P41<-(K.site1/sum(K.site1+K.site2+K.site3+K.site5))
P42<-(K.site2/sum(K.site1+K.site2+K.site3+K.site5))
P43<-(K.site3/sum(K.site1+K.site2+K.site3+K.site5))
P45<-(K.site4/sum(K.site1+K.site2+K.site3+K.site5))
P51<-(K.site1/sum(K.site1+K.site2+K.site3+K.site4))
P52<-(K.site2/sum(K.site1+K.site2+K.site3+K.site5))
P53<-(K.site3/sum(K.site1+K.site2+K.site3+K.site5))
P54<-(K.site4/sum(K.site1+K.site2+K.site3+K.site5))


#calculate SD for vital rates 
sdPhij=CV*Phij
sdPhia=CV*Phia

#sd on logit scale
sdPhijlogit=sdPhij/(Phij*(1-Phij))
sdPhialogit=sdPhia/(Phia*(1-Phia))

vcovSurv=matrix(c(sdPhijlogit^2,sdPhijlogit*sdPhialogit, sdPhijlogit*sdPhialogit, sdPhialogit^2),2,2,byrow=TRUE)

sdf=f.CV*f   #value of f.CV already on log-odds scale


###########initialize stuff
n.iter=100
n.year=100
Nmatrix<-matrix(0,nrow=n.iter,ncol=n.year)
Nmatrix1<-matrix(0,nrow=n.iter,ncol=n.year)
Nmatrix2<-matrix(0,nrow=n.iter,ncol=n.year)
Nmatrix3<-matrix(0,nrow=n.iter,ncol=n.year)
Nmatrix4<-matrix(0,nrow=n.iter,ncol=n.year)
Nmatrix5<-matrix(0,nrow=n.iter,ncol=n.year)
Extinct<-0
lambda<-matrix(0,nrow=n.iter,ncol=n.year)
lambda.site<-matrix(0,nrow=n.iter,ncol=n.year)
lambda4<-matrix(0, nrow=n.iter, ncol=n.year)


##start loop

for(iter in 1:n.iter){ 
  #set up initial population sizes (some of this can be initialized before loop)
  popmat1<-matrix(0,10,n.year)
  popmat1[1,1]<-Nsy.site1.0
  popmat1[2,1]<-Nsy.site2.0
  popmat1[3,1]<-Nsy.site3.0
  popmat1[4,1]<-Nsy.site4.0
  popmat1[5,1]<-Nsy.site5.0
  popmat1[6,1]<-Nad.site1.0
  popmat1[7,1]<-Nad.site2.0
  popmat1[8,1]<-Nad.site3.0
  popmat1[9,1]<-Nad.site4.0
  popmat1[10,1]<-Nad.site5.0
  
  
  Ntotal1<-c(popmat1[1,1]+popmat1[2,1]+popmat1[3,1]+popmat1[4,1]+popmat1[5,1]+popmat1[6,1]+popmat1[7,1]
             +popmat1[8,1]+popmat1[9,1]+popmat1[10,1],rep(0, n.year-1))
  Nsite.1<-c(popmat1[1,1]+popmat1[6,1], rep(0,n.year-1))
  Nsite.2<-c(popmat1[2,1]+popmat1[7,1], rep(0,n.year-1))
  Nsite.3<-c(popmat1[3,1]+popmat1[8,1], rep(0,n.year-1))
  Nsite.4<-c(popmat1[4,1]+popmat1[9,1], rep(0,n.year-1))
  Nsite.5<-c(popmat1[5,1]+popmat1[10,1], rep(0,n.year-1))
  
  
  for(i in 2:n.year){
    
    #generate stochastic winter survival (adult and juveniles correlated) for use at all sites
    
    survmean.w<-plogis(mvrnorm(1,c(qlogis(Phij),qlogis(Phia)),vcovSurv))
    Phij.w <- survmean.w[1]^(10/12)
    Phia.w <- survmean.w[2]^(10/12)
    survmean.b<-plogis(mvrnorm(1,c(qlogis(Phij),qlogis(Phia)),vcovSurv))
    Phijt.b<-survmean.b[1]^(2/12)
    Phiat.b<-survmean.b[2]^(2/12)
    
    ###Site 1:	
    #ceiling: sets probability of nesting to 0 for females in excess of carrying capacity
    if (Nsite.1[i-1]>K.site1){
      
      ytt1<-(popmat1[6,i-1]/Nsite.1[i-1])*(0+K.site1/Nsite.1[i-1]*yt)  #proportion of adults in population * (excess*0 + nonexcess*Fa)
      yst1<-(popmat1[1,i-1]/Nsite.1[i-1])*(0+K.site1/Nsite.1[i-1]*ys)  
      dsy21t<-0; dsy31t<-0; dsy41t<-0; dsy51t<-0  #no dispersal into site 1 if above K
      rr.sy1t<-((Nsite.1[i-1]-K.site1)*(1/0.95)+K.site1*rr.sy)/Nsite.1[i-1]  #if above K, excess birds disperse with prob = 95%
      rr.ad1t<-((Nsite.1[i-1]-K.site1)*(1/0.95)+K.site1*rr.ad)/Nsite.1[i-1]
    }#excess females 
    else{ #default values
      ytt1<-yt
      yst1<-ys
      rr.sy1t<-rr.sy
      rr.ad1t<-rr.ad
      
    }
    
    #yearly fecundity rates; 2 (=female eggs) * hatch.proportion *prob nesting * prob hatch  (=females hatched)  #demographic stochasticity: E[i] rpois(clutch.size)
    Fat1<-ytt1*(2*E*(n.ex.1*h.ex.1 + (1-n.ex.1)*h.un.1))
    Fst1<-yst1*(2*E*(n.ex.1*h.ex.1 + (1-n.ex.1)*h.un.1))
    
    #adult mortality
    m.a.ex.1=ytt1*(a.p.ex.1*m+a.p.ex.1*(1-m)*r2*(a.p.ex.1*m+a.p.ex.1*(1-m)*r3*a.p.ex.1*m))  
    #m.a.un.1=ytt1*(a.p.un*m+a.p.un*(1-m)*r2*(a.p.un*m+a.p.un*(1-m)*r3*a.p.un*m)) #removed; now a baseline, outside of loop
    m.s.ex.1=yst1*(a.p.ex.1*m+a.p.ex.1*(1-m)*r2*(a.p.ex.1*m+a.p.ex.1*(1-m)*r3*a.p.ex.1*m))  
    #m.s.un.1=yst1*(a.p.un*m+a.p.un*(1-m)*r2*(a.p.un*m+a.p.un*(1-m)*r3*a.p.un*m)) 
    
    #calculate breeding survival and stochastic fledging rates
    
    Phijt.b1<-(Phijt.b)/(1-m.s.un)*(1-m.s.ex.1*n.ex.1)
    Phiat.b1<-((Phiat.b)*(1-m.a.ex.1*n.ex.1))/(1-m.a.un)  #mean adult survival
    f.t.1<-plogis(rnorm(1,qlogis(f),sdf))
    
    #Site 2:	
    #ceiling: sets probability of nesting to 0 for females in excess of carrying capacity
    if (Nsite.2[i-1]>K.site2){
      
      ytt2<-(popmat1[7,i-1]/Nsite.2[i-1])*(0+K.site2/Nsite.2[i-1]*yt)  #proportion of adults in population * (excess*0 + nonexcess*Fa)
      yst2<-(popmat1[2,i-1]/Nsite.2[i-1])*(0+K.site2/Nsite.2[i-1]*ys)  
      
      rr.sy2t<-((Nsite.2[i-1]-K.site2)*(1/0.95)+K.site2*rr.sy)/Nsite.2[i-1]  #if above K, excess birds disperse with prob = 95%
      rr.ad2t<-((Nsite.2[i-1]-K.site2)*(1/0.95)+K.site2*rr.ad)/Nsite.2[i-1]	
      
    }#excess females 
    else{
      ytt2<-yt
      yst2<-ys
      rr.sy2t<-rr.sy
      rr.ad2t<-rr.ad
    }
    
    #yearly fecundity rates; 2 (=female eggs) * hatch.proportion *prob nesting * prob hatch  (=females hatched)  #demographic stochasticity: E[i] rpois(clutch.size)
    Fat2<-ytt2*(2*E*(n.ex.2*h.ex.2 + (1-n.ex.2)*h.un.2))
    Fst2<-yst2*(2*E*(n.ex.2*h.ex.2 + (1-n.ex.2)*h.un.2))
    
    #adult mortality
    m.a.ex.2=ytt2*(a.p.ex.2*m+a.p.ex.2*(1-m)*r2*(a.p.ex.2*m+a.p.ex.2*(1-m)*r3*a.p.ex.2*m))  
    m.s.ex.2=yst2*(a.p.ex.2*m+a.p.ex.2*(1-m)*r2*(a.p.ex.2*m+a.p.ex.2*(1-m)*r3*a.p.ex.2*m))
    
    #calculate breeding survival and stochastic fledging rates  
    
    Phijt.b2<-(Phijt.b)/(1-m.s.un)*(1-m.s.ex.2*(n.ex.2))  #mean post-fledge survival-stochastic baseline *probability of surviving abandonment 
    Phiat.b2<-(Phiat.b)/(1-m.a.un)*(1-m.a.ex.2*(n.ex.2))  #mean adult survival
    f.t.2<-plogis(rnorm(1,qlogis(f),sdf))  
    
    
    #Site 3:	
    #ceiling: sets probability of nesting to 0 for females in excess of carrying capacity
    if (Nsite.3[i-1]>K.site3){
      
      ytt3<-(popmat1[8,i-1]/Nsite.3[i-1])*(0+K.site3/Nsite.3[i-1]*yt)  #proportion of adults in population * (excess*0 + nonexcess*Fa)
      yst3<-(popmat1[3,i-1]/Nsite.3[i-1])*(0+K.site3/Nsite.3[i-1]*ys)  	
      
      
      rr.sy3t<-((Nsite.3[i-1]-K.site3)*(1/0.95)+K.site3*rr.sy)/Nsite.3[i-1]  #if above K, excess birds disperse with prob = 95%
      rr.ad3t<-((Nsite.3[i-1]-K.site3)*(1/0.95)+K.site3*rr.ad)/Nsite.3[i-1]
      
    }#excess females 
    else{
      ytt3<-yt
      yst3<-ys
      rr.sy3t<-rr.sy
      rr.ad3t<-rr.ad
      
    }
    
    #yearly fecundity rates; 2 (=female eggs) * hatch.proportion *prob nesting * prob hatch  (=females hatched)  
    Fat3<-ytt3*(2*E*(n.ex.3*h.ex.3 + (1-n.ex.3)*h.un.3))
    Fst3<-yst3*(2*E*(n.ex.3*h.ex.3 + (1-n.ex.3)*h.un.3))
    
    #adult mortality
    m.a.ex.3=ytt3*(a.p.ex.3*m+a.p.ex.3*(1-m)*r2*(a.p.ex.3*m+a.p.ex.3*(1-m)*r3*a.p.ex.3*m))  
    m.s.ex.3=yst3*(a.p.ex.3*m+a.p.ex.3*(1-m)*r2*(a.p.ex.3*m+a.p.ex.3*(1-m)*r3*a.p.ex.3*m))	
    
    #calculate breeding survival and stochastic fledging rates
    
    Phijt.b3<-(Phijt.b)/(1-m.s.un)*(1-m.s.ex.3*(n.ex.3))  #mean post-fledge survival-stochastic baseline *probability of surviving abandonment 
    Phiat.b3<-(Phiat.b)/(1-m.a.un)*(1-m.a.ex.3*(n.ex.3))  #mean adult survival
    f.t.3<-plogis(rnorm(1,qlogis(f),sdf))  
    
    #Site 4:	
    #ceiling: sets probability of nesting to 0 for females in excess of carrying capacity
    if (Nsite.4[i-1]>K.site4){
      
      ytt4<-(popmat1[9,i-1]/Nsite.4[i-1])*(0+K.site4/Nsite.4[i-1]*yt)  #proportion of adults in population * (excess*0 + nonexcess*Fa)
      yst4<-(popmat1[4,i-1]/Nsite.4[i-1])*(0+K.site4/Nsite.4[i-1]*ys)  
      
      
      rr.sy4t<-((Nsite.4[i-1]-K.site4)*(1/0.95)+K.site4*rr.sy)/Nsite.4[i-1]  #if above K, excess birds disperse with prob = 95%
      rr.ad4t<-((Nsite.4[i-1]-K.site4)*(1/0.95)+K.site4*rr.ad)/Nsite.4[i-1]	
      
    }#excess females 
    else{
      ytt4<-yt
      yst4<-ys
      rr.sy4t<-rr.sy
      rr.ad4t<-rr.ad
    }
    
    #yearly fecundity rates; 2 (=female eggs) * hatch.proportion *prob nesting * prob hatch  (=females hatched)  
    Fat4<-ytt4*(2*E*(n.ex.4*h.ex.4 + (1-n.ex.4)*h.un.4))
    Fst4<-yst4*(2*E*(n.ex.4*h.ex.4 + (1-n.ex.4)*h.un.4))
    
    #adult mortality
    m.a.ex.4=ytt4*(a.p.ex.4*m+a.p.ex.4*(1-m)*r2*(a.p.ex.4*m+a.p.ex.4*(1-m)*r3*a.p.ex.4*m))  
    m.s.ex.4=yst4*(a.p.ex.4*m+a.p.ex.4*(1-m)*r2*(a.p.ex.4*m+a.p.ex.4*(1-m)*r3*a.p.ex.4*m))
    
    #calculate breeding survival and stochastic fledging rates 
    
    Phijt.b4<-(Phijt.b)/(1-m.s.un)*(1-m.s.ex.4*(n.ex.4))  #mean post-fledge survival-stochastic baseline *probability of surviving abandonment 
    Phiat.b4<-(Phiat.b)*(1-m.a.un)*(1-m.a.ex.4*(n.ex.4))  #mean adult survival
    f.t.4<-plogis(rnorm(1,qlogis(f),sdf))  
    
    #Site 5:	
    #ceiling: sets probability of nesting to 0 for females in excess of carrying capacity
    if (Nsite.5[i-1]>K.site5){
      
      ytt5<-(popmat1[10,i-1]/Nsite.5[i-1])*(0+K.site5/Nsite.5[i-1]*yt)  #proportion of adults in population * (excess*0 + nonexcess*Fa)
      yst5<-(popmat1[5,i-1]/Nsite.5[i-1])*(0+K.site5/Nsite.5[i-1]*ys)  
      
      
      rr.sy5t<-((Nsite.5[i-1]-K.site5)*(1/0.95)+K.site5*rr.sy)/Nsite.5[i-1]  #if above K, excess birds disperse with prob = 95%
      rr.ad5t<-((Nsite.5[i-1]-K.site5)*(1/0.95)+K.site5*rr.ad)/Nsite.5[i-1]	
      
      
    }#excess females 
    else{
      ytt5<-yt
      yst5<-ys
      rr.sy5t<-rr.sy
      rr.ad5t<-rr.ad
      
    }
    
    #yearly fecundity rates; 2 (=female eggs) * hatch.proportion *prob nesting * prob hatch  (=females hatched)  
    Fat5<-ytt5*(2*E*(n.ex.5*h.ex.5 + (1-n.ex.5)*h.un.5))
    Fst5<-yst5*(2*E*(n.ex.5*h.ex.5 + (1-n.ex.5)*h.un.5))
    
    #adult mortality
    m.a.ex.5=ytt5*(a.p.ex.5*m+a.p.ex.5*(1-m)*r2*(a.p.ex.5*m+a.p.ex.5*(1-m)*r3*a.p.ex.5*m))  
    m.s.ex.5=yst5*(a.p.ex.5*m+a.p.ex.5*(1-m)*r2*(a.p.ex.5*m+a.p.ex.5*(1-m)*r3*a.p.ex.5*m))
    
    #calculate breeding survival and stochastic fledging rates 
    
    Phijt.b5<-(Phijt.b)/(1-m.s.un)*(1-m.s.ex.5*(n.ex.5))  #mean post-fledge survival-stochastic baseline *probability of surviving abandonment 
    Phiat.b5<-(Phiat.b)/(1-m.a.un)*(1-m.a.ex.5*(n.ex.5))  #mean adult survival
    f.t.5<-plogis(rnorm(1,qlogis(f),sdf))  
    
    #calculate dispersal rates of sy birds
    dsy12t<-(1-rr.sy1t)*P12*(1-(1/(1+exp(-(Nsite.2[i-1]-K.site2) ) ) ) )  #(1-return rate site 1)*relative size of site 2 * (multiplier; moves rapidly to 1 when site.2 < K; to 1 when site.2 > K)
    
    dsy13t<-(1-rr.sy1t)*P13*(1-(1/(1+exp(-(Nsite.3[i-1]-K.site3) ) ) ) )
    dsy14t<-(1-rr.sy1t)*P14*(1-(1/(1+exp(-(Nsite.4[i-1]-K.site4) ) ) ) )
    dsy15t<-(1-rr.sy1t)*P15*(1-(1/(1+exp(-(Nsite.5[i-1]-K.site5) ) ) ) )	
    dsy21t<-(1-rr.sy2t)*P21*(1-(1/(1+exp(-(Nsite.1[i-1]-K.site1) ) ) ) )
    dsy23t<-(1-rr.sy2t)*P23*(1-(1/(1+exp(-(Nsite.3[i-1]-K.site3) ) ) ) )
    dsy24t<-(1-rr.sy2t)*P24*(1-(1/(1+exp(-(Nsite.4[i-1]-K.site4) ) ) ) )
    dsy25t<-(1-rr.sy2t)*P25*(1-(1/(1+exp(-(Nsite.5[i-1]-K.site5) ) ) ) )
    dsy31t<-(1-rr.sy3t)*P31*(1-(1/(1+exp(-(Nsite.1[i-1]-K.site1) ) ) ) )
    dsy32t<-(1-rr.sy3t)*P31*(1-(1/(1+exp(-(Nsite.2[i-1]-K.site2) ) ) ) )
    dsy34t<-(1-rr.sy3t)*P34*(1-(1/(1+exp(-(Nsite.4[i-1]-K.site4) ) ) ) )
    dsy35t<-(1-rr.sy3t)*P35*(1-(1/(1+exp(-(Nsite.5[i-1]-K.site5) ) ) ) )
    dsy41t<-(1-rr.sy4t)*P41*(1-(1/(1+exp(-(Nsite.1[i-1]-K.site1) ) ) ) )
    dsy42t<-(1-rr.sy4t)*P42*(1-(1/(1+exp(-(Nsite.2[i-1]-K.site2) ) ) ) )
    dsy43t<-(1-rr.sy4t)*P43*(1-(1/(1+exp(-(Nsite.3[i-1]-K.site3) ) ) ) )
    dsy45t<-(1-rr.sy4t)*P45*(1-(1/(1+exp(-(Nsite.5[i-1]-K.site5) ) ) ) )
    dsy51t<-(1-rr.sy5t)*P51*(1-(1/(1+exp(-(Nsite.1[i-1]-K.site1) ) ) ) )
    dsy52t<-(1-rr.sy5t)*P52*(1-(1/(1+exp(-(Nsite.2[i-1]-K.site2) ) ) ) )
    dsy53t<-(1-rr.sy5t)*P53*(1-(1/(1+exp(-(Nsite.3[i-1]-K.site3) ) ) ) )
    dsy54t<-(1-rr.sy5t)*P54*(1-(1/(1+exp(-(Nsite.4[i-1]-K.site4) ) ) ) )
    
    #calculate dispersal rates of adult birds
    dad12t<-(1-rr.ad1t)*P12*(1-(1/(1+exp(-(Nsite.2[i-1]-K.site2) ) ) ) )  
    dad13t<-(1-rr.ad1t)*P13*(1-(1/(1+exp(-(Nsite.3[i-1]-K.site3) ) ) ) )
    dad14t<-(1-rr.ad1t)*P14*(1-(1/(1+exp(-(Nsite.4[i-1]-K.site4) ) ) ) )
    dad15t<-(1-rr.ad1t)*P15*(1-(1/(1+exp(-(Nsite.5[i-1]-K.site5) ) ) ) )	
    dad21t<-(1-rr.ad2t)*P21*(1-(1/(1+exp(-(Nsite.1[i-1]-K.site1) ) ) ) )
    dad23t<-(1-rr.ad2t)*P23*(1-(1/(1+exp(-(Nsite.3[i-1]-K.site3) ) ) ) )
    dad24t<-(1-rr.ad2t)*P24*(1-(1/(1+exp(-(Nsite.4[i-1]-K.site4) ) ) ) )
    dad25t<-(1-rr.ad2t)*P25*(1-(1/(1+exp(-(Nsite.5[i-1]-K.site5) ) ) ) )
    dad31t<-(1-rr.ad3t)*P31*(1-(1/(1+exp(-(Nsite.1[i-1]-K.site1) ) ) ) )
    dad32t<-(1-rr.ad3t)*P31*(1-(1/(1+exp(-(Nsite.2[i-1]-K.site2) ) ) ) )
    dad34t<-(1-rr.ad3t)*P34*(1-(1/(1+exp(-(Nsite.4[i-1]-K.site4) ) ) ) )
    dad35t<-(1-rr.ad3t)*P35*(1-(1/(1+exp(-(Nsite.5[i-1]-K.site5) ) ) ) )
    dad41t<-(1-rr.ad4t)*P41*(1-(1/(1+exp(-(Nsite.1[i-1]-K.site1) ) ) ) )
    dad42t<-(1-rr.ad4t)*P42*(1-(1/(1+exp(-(Nsite.2[i-1]-K.site2) ) ) ) )
    dad43t<-(1-rr.ad4t)*P43*(1-(1/(1+exp(-(Nsite.3[i-1]-K.site3) ) ) ) )
    dad45t<-(1-rr.ad4t)*P45*(1-(1/(1+exp(-(Nsite.5[i-1]-K.site5) ) ) ) )
    dad51t<-(1-rr.ad5t)*P51*(1-(1/(1+exp(-(Nsite.1[i-1]-K.site1) ) ) ) )
    dad52t<-(1-rr.ad5t)*P52*(1-(1/(1+exp(-(Nsite.2[i-1]-K.site2) ) ) ) )
    dad53t<-(1-rr.ad5t)*P53*(1-(1/(1+exp(-(Nsite.3[i-1]-K.site3) ) ) ) )
    dad54t<-(1-rr.ad5t)*P54*(1-(1/(1+exp(-(Nsite.4[i-1]-K.site4) ) ) ) )
    
    M.grand<-matrix(c(rr.sy1t, dsy21t, dsy31t, dsy41t, dsy51t, 0,0,0,0,0,
                      dsy12t, rr.sy2t, dsy32t, dsy42t, dsy52t, 0,0,0,0,0,
                      dsy13t, dsy23t, rr.sy3t, dsy43t, dsy53t, 0,0,0,0,0,
                      dsy14t, dsy24t, dsy34t, rr.sy4t, dsy45t, 0,0,0,0,0,
                      dsy15t, dsy25t, dsy35t, dsy45t, rr.sy5t, 0,0,0,0,0,
                      0,0,0,0,0, rr.ad1t, dad21t, dad31t, dad41t, dad51t,
                      0,0,0,0,0, dad12t, rr.ad2t, dad32t, dad42t, dad52t,
                      0,0,0,0,0, dad13t, dad23t, rr.ad3t, dad43t, dad53t,
                      0,0,0,0,0, dad14t, dad24t, dad34t, rr.ad4t, dad54t,
                      0,0,0,0,0, dad15t, dad25t, dad35t, dad45t, rr.ad5t), nrow=10, byrow=TRUE)
    
    
    B.grand<-matrix(c(Phijt.b1*Phij.w*Fst1*f.t.1, Fat1*Phiat.b1*Phia.w,0,0,0,0,0,0,0,0,
                      Phijt.b1*Phij.w*f.t.1,Phiat.b1*Phia.w,0,0,0,0,0,0,0,0,
                      0,0,Phijt.b2*Phij.w*Fst2*f.t.2, Fat2*Phiat.b2*Phia.w,0,0,0,0,0,0,
                      0,0,Phijt.b2*Phij.w*f.t.2,Phiat.b2*Phia.w,0,0,0,0,0,0,
                      0,0,0,0,Phijt.b3*Phij.w*Fst3*f.t.3, Fat3*Phiat.b3*Phia.w,0,0,0,0,
                      0,0,0,0,Phijt.b3*Phij.w*f.t.3,Phiat.b3*Phia.w,0,0,0,0,
                      0,0,0,0,0,0,Phijt.b4*Phij.w*Fst4*f.t.4, Fat4*Phiat.b4*Phia.w,0,0,
                      0,0,0,0,0,0,Phijt.b4*Phij.w*f.t.4,Phiat.b4*Phia.w,0,0,
                      0,0,0,0,0,0,0,0,Phijt.b5*Phij.w*Fst5*f.t.5, Fat5*Phiat.b5*Phia.w,
                      0,0,0,0,0,0,0,0,Phijt.b5*Phij.w*f.t.5,Phiat.b5*Phia.w
    ),
    nrow=10,byrow=TRUE)
    
    
    popmat1[1:10,i]<-P%*%B.grand%*%t(P)%*%M.grand%*%popmat1[1:10,i-1]   #where columns arranged as Nsy.site1, Nsy.site2, Nad.site1, Nad.site2
    Ntotal1[i]<-popmat1[1,i]+popmat1[2,i]+popmat1[3,i]+popmat1[4,i]+popmat1[5,i]+popmat1[6,i]+popmat1[7,i]+popmat1[8,i]+popmat1[9,i]+popmat1[10,i]
    Nsite.1[i]<-popmat1[1,i]+popmat1[6,i]
    Nsite.2[i]<-popmat1[2,i]+popmat1[7,i]
    Nsite.3[i]<-popmat1[3,i]+popmat1[8,i]
    Nsite.4[i]<-popmat1[4,i]+popmat1[9,i]
    Nsite.5[i]<-popmat1[5,i]+popmat1[10,i]
    
    
    lambda[iter,i]<-eigen(B.grand)$values[1]
    lambda4[iter,i]<-eigen(B.grand[7:8,7:8])$values[1]
    
  }#n.year
  
  Extinct[Ntotal1[100]<1]=Extinct+1
  
  
  
  Nmatrix[iter,]<-Ntotal1[]
  Nmatrix1[iter,]<-Nsite.1[]
  Nmatrix2[iter,]<-Nsite.2[]
  Nmatrix3[iter,]<-Nsite.3[]
  Nmatrix4[iter,]<-Nsite.4[]
  Nmatrix5[iter,]<-Nsite.5[]
  
  
}#n.iter

Nave=colMeans(Nmatrix)
Nsd<-apply(Nmatrix,2,sd)
Nlow<-Nave-1.96*Nsd
Nhi<-Nave+1.96*Nsd

plot(c(1:100), ylim=c(0,100+max(Nhi)), Nave, type="l", lwd=2)
lines(c(1:100), Nlow, lty=2)
lines(c(1:100), Nhi, lty=2)

Nave1<-colMeans(Nmatrix1)
Nave2<-colMeans(Nmatrix2)
Nave3<-colMeans(Nmatrix3)
Nave4<-colMeans(Nmatrix4)
Nave5<-colMeans(Nmatrix5)

plot(c(1:100), ylim=c(0,100+max(Nave1)), Nave1, type="l", lwd=2)
lines(c(1:100), Nave2, type="l", lwd=2, col="red")
lines(c(1:100), Nave3, type="l", lwd=2, col="blue")
lines(c(1:100), Nave4, type="l", lwd=2, col="green")
lines(c(1:100), Nave5, type="l", lwd=2, col="purple")

