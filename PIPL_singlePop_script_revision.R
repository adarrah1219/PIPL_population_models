#revision following some aspects of Erik Osnas' code

library(MASS)

# functions####

hatch.calc<-function(parms, eta.a=0, eta.p=0, ex = 0){
	lin.f<-exp(parms$alpha.f)
	lin.a<-exp(parms$alpha.a + eta.a + parms$beta.a.ex*ex)
	lin.p<-exp(parms$alpha.p + eta.p + parms$beta.p.ex*ex)
	survp<-1/(lin.f+lin.a+lin.p+1)   #daily survival probability

	#period fate probabilities for each nest
	pred.p<-((lin.p/(lin.f+lin.a+lin.p+1))/(1-survp)*(1-survp^34))  #predation unexclosed
	a<-((lin.a/(lin.f+lin.a+lin.p+1))/(1-survp)*(1-survp^34))  #abandonment unexclosed
	flood.p<-((lin.f/(lin.f+lin.a+lin.p+1))/(1-survp)*(1-survp^34))  #flooding probability
	o<-flood.p+pred.p
	h<-survp^34 
	m=parms$m
	r2=parms$r2
	r3=parms$r3
	season.hatch <- h + (o*r2*h + a*(1-m)*r2*h) + o*r2*(o*r3*h + a*(1-m)*r3*h) + a*(1-m)*r2*(o*r3*h+a*(1-m)*r3*h) #revised; mistake somewhere in original
	output<-list(hatch = h, season.hatch = season.hatch, a=a, p = pred.p)
	return(output)
	}


logodds=function(x){log(x/(1-x))}

#parameters to use ####

CV=0.1

parms <- list(
  Na0 = 20,
  Ns0=20,
  yt=0.99,  #probability of third-year bird nesting
  ys=0.68,  #probability of second-year bird nesting
  Phij=0.52,  
  Phia=0.74,
  f=0.4,
  m=0.78,
  E = 0.94,
  r2 = 0.7,
  r3 = 0.7,
  beta.a.ex = 1.284,#values from model "reff.siteYear" in preliminary analyses.RData
  beta.p.ex = -2.532,
  alpha.a = -7.447,
  alpha.p = -4.257,
  alpha.f = -6.040
)	
	
#standard deviations for stochasticity	
sd.parms <- list(
	vcovSurv=matrix(c((logodds(parms$Phij)*CV)^2,logodds(parms$Phij)*CV*logodds(parms$Phia)*CV, 
	           logodds(parms$Phij)*CV*logodds(parms$Phia)*CV, (logodds(parms$Phia)*CV))^2,2,2),  #had sd in here earlier, not variance
	sd.f = abs(logodds(parms$f))*0.06,
	yt = logodds(parms$yt)*CV,
	ys = logodds(parms$ys)*CV,
	r2 = abs(logodds(parms$r2))*CV,
	r3 = abs(logodds(parms$r3))*CV,
	m = abs(logodds(parms$m))*1.04,  #this might need to be changed to reflect halving the mortality risk for female-based model?
	beta.a.ex = 0.348,
	beta.p.ex = 0.272,
	alpha.a = 0.386,
	alpha.p = 0.116,
	alpha.f = 0.222,
	vcovAban = matrix(c(0.386^2, -0.092, -0.092, 0.348^2),2,2),
	vcovPred = matrix(c(0.116^2, -0.007, -0.007, 0.272^2),2,2)
)

#function to calculate lambda from draws of parameters ####

lambda.calc <- function(parms, sd.parms, eta.a=0, eta.p=0, n.ex=0, n.iter=1000){
  #define a few things
  eta.a=eta.a
  eta.p=eta.p
  lambda <- rep(NA, n.iter)
  aban.counts.ex <- rep(NA, n.iter)
  aban.counts.un <- rep(NA, n.iter)
  dead.females.ex <- rep(NA, n.iter)
  
  #run iterations
  for (i in 1:n.iter){
  #draw stochastic values
    parms.iter <- list(
      survmean=plogis(mvrnorm(1,c(qlogis(parms$Phij),qlogis(parms$Phia)),sd.parms$vcovSurv)),
      f = plogis(rnorm(1,qlogis(parms$f),sd.parms$sd.f)),
      ys = rnorm(1,logodds(parms$ys),sd.parms$ys), #probability of second-year bird nesting
     yt = parms$yt, #probability of third-year bird nesting fixed at 0.99
     r2 = rnorm(1, logodds(parms$r2), sd.parms$r2),
      r3 = rnorm(1, logodds(parms$r3), sd.parms$r3),
      m = plogis(rnorm(1,qlogis(parms$m), sd.parms$m)),
      alpha.f = rnorm(1, parms$alpha.f, sd.parms$alpha.f),
      aban.coeff = mvrnorm(1, c(parms$alpha.a, parms$beta.a.ex), sd.parms$vcovAban),
      pred.coeff = mvrnorm(1, c(parms$alpha.p, parms$beta.p.ex), sd.parms$vcovPred)
       )
    parms.iter$alpha.a <- parms.iter$aban.coeff[1] + eta.a
    parms.iter$beta.a.ex <- parms.iter$aban.coeff[2]
    parms.iter$alpha.p <- parms.iter$pred.coeff[1] + eta.p
    parms.iter$beta.p.ex <- parms.iter$pred.coeff[2]
    
  #baseline abandonment rate (unexclosed)
  a.p.un<-hatch.calc(parms=parms.iter)$a
  
	#hatch and exclosure-related abandonment probabilities
	h.un<-hatch.calc(parms=parms.iter)$season.hatch
	h.ex<-hatch.calc(parms=parms.iter, ex=1)$season.hatch
	a.p.ex<-hatch.calc(parms=parms.iter, ex=1)$a

  #fecundity
	Fa<-parms.iter$yt*(2*parms$E*(n.ex*h.ex + (1-n.ex)*h.un))
	Fs<-parms.iter$ys*(2*parms$E*(n.ex*h.ex + (1-n.ex)*h.un))

  #breeding-season mortality with abandonment
	m<-parms.iter$m; r2<-parms.iter$r2; r3<-parms.iter$r3
	m.a.un=parms.iter$yt*(a.p.un*m+a.p.un*(1-m)*r2*(a.p.un*m+a.p.un*(1-m)*r3*a.p.un*m)) 
	m.s.un=parms.iter$ys*(a.p.un*m+a.p.un*(1-m)*r2*(a.p.un*m+a.p.un*(1-m)*r3*a.p.un*m))
	m.a.ex=parms.iter$yt*(a.p.ex*m+a.p.ex*(1-m)*r2*(a.p.ex*m+a.p.ex*(1-m)*r3*a.p.ex*m))  
	m.s.ex=parms.iter$ys*(a.p.ex*m+a.p.ex*(1-m)*r2*(a.p.ex*m+a.p.ex*(1-m)*r3*a.p.ex*m))
	
	
	aban.counts.ex[i]<-a.p.ex*n.ex*sum(parms$Na0,parms$Ns0)
	aban.counts.un[i]<-a.p.un*sum(parms$Na0,parms$Ns0)  #reference point if exclosures used
	dead.females.ex[i]<- m.s.ex*parms$Ns0 + m.a.ex*parms$Na0

	#annual survival rates	
	Phij.w <- parms.iter$survmean[1]^(10/12) #winter post-fledging juvenile survival
	Phia.w <- parms.iter$survmean[2]^(10/12) #winter adult survival
	Phij.b <- parms.iter$survmean[1]^(2/12)  #second-year breeding season survival
	Phia.b <- parms.iter$survmean[2]^(2/12)  #ASY breeding survival
	Phij.b <- Phij.b/(1-m.s.un)*(1-m.s.ex*n.ex)  #add in probability of surviving exclosure-related abandonment
	Phia.b <- Phia.b/(1-m.a.un)*(1-m.a.ex*n.ex)
	Phij.ann <- Phij.b*Phij.w 
  Phia.ann <- Phia.b*Phia.w
	
  #matrix calculations
	f<-parms.iter$f
	Lesmat<-matrix(c(Phij.ann*Fs*f,Fa*Phia.ann,Phij.ann*f,Phia.ann),2,2,byrow=TRUE)  #found mistake in original; entry [2,1] had adult survival, not juv
	lambda[i]<-eigen(Lesmat)$values[1]
  } #n.iter
  output <- list(lambda=lambda, aban.counts.ex=aban.counts.ex, aban.counts.un=aban.counts.un, dead.females.ex=dead.females.ex)
  return(output)
}



###########plotting and exporting values for plotting in tool #############

n.iter<-100
n.vals <-100
eta.a.vec <- rnorm(n.vals,0,1.351)
eta.p.vec <- rnorm(n.vals,0,1.5) #variance expanded to get more extreme values for plot
eta.a.vec <- eta.a.vec[sort.list(eta.a.vec)] #sort from low to high; makes plotting later easier
eta.p.vec <- eta.p.vec[sort.list(eta.p.vec)]

a.prob.ex <- hatch.calc(parms, eta.a = eta.a.vec[], eta.p = eta.p.vec[], ex=1)$a
a.prob.un <- hatch.calc(parms, eta.a = eta.a.vec[], eta.p = eta.p.vec[], ex=0)$a
p.prob.ex <- hatch.calc(parms, eta.a = eta.a.vec[], eta.p = eta.p.vec[], ex=1)$p
p.prob.un <- hatch.calc(parms, eta.a = eta.a.vec[], eta.p = eta.p.vec[], ex=0)$p

lambda.mat.ex <- lambda.mat.un <- matrix(NA, nrow=n.vals, ncol=n.vals)

for(i in 1:n.vals){ 
  for(j in 1:n.vals){
  
  growth.ex <- lambda.calc(parms=parms, sd.parms=sd.parms, eta.a=eta.a.vec[i], n.ex=1, eta.p=eta.p.vec[j],n.iter=n.iter)
  lambda.mat.ex[i,j] <- mean(growth.ex$lambda)
  
  growth.un <- lambda.calc(parms=parms, sd.parms=sd.parms, eta.a=eta.a.vec[i], n.ex=0, eta.p=eta.p.vec[j],n.iter=n.iter)
  lambda.mat.un[i,j] <- mean(growth.un$lambda)
  }
}

#difference in lambda as a function of exclosure use, for different combinations of abandonment and predation probabilities
lambda.diff <- lambda.mat.ex  - lambda.mat.un


#plot for probabilities####
df2 <- as.data.frame(cbind(as.numeric(lambda.diff), rep(a.prob.ex,n.vals), rep(p.prob.un, each=n.vals)))
colnames(df2) <- c("lambda.diff","a.prob.ex", "p.prob.un")
write.csv(df2, "predictionFrame.csv") #export for use in tool
df2<-read.csv("predictionFrame.csv", header=T)
lambda.loess2 = loess(lambda.diff~a.prob.ex*p.prob.un, data = df2) #prediction function for creating smooth plot
values.fit2 = expand.grid(list(a.prob.ex = seq(min(a.prob.ex),max(a.prob.ex),length.out = test.vals),
                              p.prob.un = seq(min(p.prob.un),max(p.prob.un),length.out=test.vals)))
lambda.predict2 = predict(lambda.loess2, newdata = values.fit2)
values.fit2$lambda.diff <- as.numeric(lambda.predict2)

lambda.predict.mat2 <- matrix(values.fit2$lambda.diff, nrow=test.vals, ncol=test.vals)


#3D plot####
#choose x and y sequences
x <- seq(min(values.fit2$a.prob.ex),max(values.fit2$a.prob.ex),length.out=50)
y <- seq(min(values.fit2$p.prob.un),max(values.fit2$p.prob.un),length.out=50)
z <- lambda.predict.mat2[plot.vals[2:51],plot.vals[2:51]]

#make colors for plot
nrz<-length(plot.vals)-1
ncz<-length(plot.vals)-1

jet.colors <-  colorRampPalette(c("darkred","red","orange","yellow","green"))
nbcol<-64
color<-jet.colors(nbcol)
zfacet<-z[-1,-1]+z[-1,-ncz]+z[-nrz,-1]+z[-nrz,-ncz]
facetcol<-cut(zfacet,nbcol)

graph <- persp(x,y,z,
      zlab="Gain in Lambda with Exclosures", xlab="Abandonment",ylab="Predation",
      phi=20,theta=30, ticktype="detailed", col=color[facetcol])
site.pt <- c(0.2,0.6)
xval<-which.min(abs(x-site.pt[1]))
yval<-which.min(abs(y-site.pt[2]))
zval<-z[xval,yval]
add.pts<-trans3d(x[xval],y[yval],zval,pmat=graph)
points(add.pts, pch="O",col="black")

#testing area
graph <- persp(x,y,z, expand=0.9,mar=c(10,1,0,2),
               zlab="", xlab="",ylab="", border=NA,
               phi=20,theta=30, ticktype="detailed", col=color[facetcol])
x.label <- c("Abandonment")
x.label2 <- c("with Exclosures")
y.label <- c("Predation")
y.label2 <- c("without Exclosures")
x.label.pos <- trans3d(0.1, -0.1, -1.1, pmat=graph)
x.label2.pos <- trans3d(0.1,-0.18,-1.15,pmat=graph)
y.label.pos <- trans3d(0.6, 0.1, -1.05, pmat=graph)
y.label2.pos <- trans3d(0.64, 0.1, -1.1, pmat=graph)
text(x.label.pos$x,x.label.pos$y,labels=x.label, srt=-25, cex=1.2, adj=c(0,NA), xpd=T)
text(x.label2.pos$x,x.label2.pos$y,labels=x.label2,srt=-25,cex=1.2,adj=c(0,NA),xpd=T)
text(y.label.pos$x, y.label.pos$y, labels=y.label, srt=60, cex=1.2, adj=c(0,NA), xpd=T)
text(y.label2.pos$x, y.label2.pos$y, labels=y.label2, srt=60, cex=1.2, adj=c(0,NA), xpd=T)
z.label <- ("Gain in Growth Rate")
z.label2 <- ("with Exclosures")
z.label.pos <- trans3d(-0.2, 0, 0.2, pmat=graph)
z.label2.pos <- trans3d(-0.25,-0.03,0.2, pmat=graph)
text(z.label.pos$x, z.label.pos$y, labels=z.label, srt=-82, cex=1.2, adj=c(0,NA), xpd=T)
text(z.label2.pos$x, z.label2.pos$y, labels=z.label2, srt=-82, cex=1.2, adj=c(0,NA), xpd=T)

site.pt <- c(0.07,0.4)

xval<-which.min(abs(x-site.pt[1]))
yval<-which.min(abs(y-site.pt[2]))
zval<-z[xval,yval]
add.pts<-trans3d(x[xval],y[yval],zval,pmat=graph)
points(add.pts, pch=19,col="black")
site.pos <- trans3d(site.pt[1]+0.05, site.pt[2]+0.05, zval, pmat=graph)
text(site.pos,labels=c("Your Site"), adj=c(0,NA), cex=2)

#contour line

no.effect.line <- contourLines(x, y, z, nlevels=1, level=0)
no.effect <- trans3d(no.effect.line[[1]]$x, no.effect.line[[1]]$y, rep(0, length(no.effect.line[[1]]$x)), pmat=graph)
#polygon(no.effect)
lines(no.effect)

#add variance intervals
sd.site.pt <- c(0.05,0.08)
var.mat <- matrix(c(0.05^2,0,0,0.08^2), nrow=2, byrow=T)
library(car)
var.points<-ellipse(site.pt,shape=var.mat, radius=1)  
var.z <- rep(NA, length(var.points[,1]))
var.z <- predict(lambda.loess2, newdata = var.points)

var.points.trans <- trans3d(var.points[,1], var.points[,2], z=var.z+0.01, pmat=graph)
lines(var.points.trans) 

#export values to use in shiny
write.csv(z,"zvalues.csv")
write.csv(cbind(x,y), "predictors.csv")

save.image("~/Projections/single_site/3DgraphScript.RData")
###############

#threshold plot####
require(splines)

n.iter<-100
n.vals <-100
eta.a.vec <- rnorm(n.vals,0,1.351)
eta.a.vec <- eta.a.vec[sort.list(eta.a.vec)] #sort from low to high; makes plotting later easier

a.prob.ex <- hatch.calc(parms, eta.a = eta.a.vec[], eta.p = 0, ex=1)$a
a.prob.un <- hatch.calc(parms, eta.a = eta.a.vec[], eta.p = 0, ex=0)$a

lambda.mat.ex <- lambda.mat.un <- abans.ex <- abans.un <- matrix(NA, nrow=n.vals, ncol=n.vals)

for(i in 1:n.vals){ 

    growth.ex <- lambda.calc(parms=parms, sd.parms=sd.parms, eta.a=eta.a.vec[i], n.ex=1, eta.p=0,n.iter=n.iter)
    lambda.mat.ex[,i] <- mean(growth.ex$lambda)
    abans.ex[,i] <- growth.ex$aban.counts.ex
    growth.un <- lambda.calc(parms=parms, sd.parms=sd.parms, eta.a=eta.a.vec[i], n.ex=0, eta.p=0,n.iter=n.iter)
    lambda.mat.un[,i] <- mean(growth.un$lambda)
    abans.un[,i] <- growth.un$aban.counts.un
  
}

lambda.diff <- lambda.mat.ex  - lambda.mat.un

#summarize data for plotting
ref.line <- lambda.calc(parms=parms, sd.parms=sd.parms, eta.a=0, eta.p=0, 
                        n.ex=0, n.iter=n.iter) #unexclosed reference
aban.ex.means <- colMeans(abans.ex)
prob.decline <- prob.decline.ref <- rep(NA, n.vals)
prob.decline.ref <- mean(ref.line$lambda[]<1)
for (i in 1:n.vals){
  prob.decline[i] <- length(which(lambda.mat.ex[,i]<1))/n.iter 
}

par(oma=c(1,1,1,1), mar=c(6,6,2,1)) 
plot(x=NULL, y=NULL, xlab="Number of Exclosure-related Abandonments", ylab="Probability of Population Decline", 
     ylim=c(0,1), xlim=c(0, max(aban.ex.means)))
prob.curv <- smooth.spline(aban.ex.means, prob.decline)
lines(prob.curv, lwd=2)
if (length(prob.curv$y)>=n.vals) {y.coord <- prob.curv$y[1:200]} else {y.coord <- c(prob.curv$y[1:(n.vals-length(prob.curv$y))], prob.curv)}
polygon(c(aban.ex.means[1], aban.ex.means, aban.ex.means[n.vals]), c(0,y.coord,0), density=10, angle=45)

abline(h=prob.decline.ref, lty=2, lwd=2)
par(xpd=T)
legend("topright", c("Exclosed", "Unexclosed Reference"), xpd=T, lty=c(1,2), lwd=c(2,2))

diffs <- abs(prob.curv$y-ref.line$lambda)

aban.tolerance <- round(prob.curv$x[which(diffs[]==min(diffs))],0)
text(0, 1.1, labels=paste("Reassess after", aban.tolerance, "observed nest abandonments"), xpd=T, adj=c(0,NA))
