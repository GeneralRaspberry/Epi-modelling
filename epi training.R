library("epimdr")
library("rootSolve")

require(deSolve)

sirmod=function(t,y,parms){
  
  S=y[1]
  I=y[2]
  R=y[3]
  
  beta=parms["beta"]
  mu= parms["mu"]
  gamma = parms["gamma"]
  N= parms["N"]
  
  dS=mu*(N-S) - beta * S * I/N
  dI = beta * S * I/N - (mu + gamma) * I
  dR = gamma * I - mu * R
  
res = c(dS, dI, dR)

list(res)
}

times = seq(0,50*52,by=1/10)
parms = c(mu = 1/(50*52), N= 1, beta = 2, gamma = 1/2)
start = c(S = .19, I= .01, R=0.8)

out=ode(y=start, times=times, func=sirmod, parms=parms)
out=as.data.frame(out)

plot(x=out$time,y=out$S,ylab="Fraction", xlab="time", type="l")
lines(x=out$time,y=out$I,col="red")
lines(x=out$time,y=out$R,col="green")

par(mfrow=c(1,2)) #make room for side by side plot
plot(times,out$I,ylab="Fraction",xlab="Time",type="l")
plot(out$S,out$I,type="l",xlab="Susceptible",ylab="Infected")


#phase analysis

simod=function(t,y,parameters){
  S=y[1]
  I=y[2]
  beta=parameters["beta"]
  mu=parameters["mu"]
  gamma=parameters["gamma"]
  N=parameters["N"]
  
  dS=mu*(N-S) - beta * S * I/N
  dI = beta * S * I/N - (mu + gamma) * I
  res= c(dS, dI)
  list(res)
}


library("phaseR")
#Plot vector field
par(mfrow=c(1,1))
fld=flowField(simod,xlim=c(.15,.35),ylim=c(0,.01),parameters=parms,
              system="two.dim",add=FALSE,ylab="I",xlab="S")

#Add trajectory

out = as.data.frame(ode(y=c(S=.19,I=.01),times=seq(0,52*100,by=.1),func=simod,
                        parms=parms))
lines(out$S,out$I,col="red")

#Add S isocline

curve(parms["mu"]*(1/x-1)/parms["beta"],.15,.35,
       xlab="S",ylab="I",add=TRUE)

#Add I-isocline
shat=(parms["gamma"]+parms["mu"])/parms["beta"]
      lines(rep(shat,2),c(0.,0.01))

legend ("topright",legend=c("Transient","Isoclines"),
        lty=c(1,1),col=c("red","black"))

#pull values from parms vector
gamma=parms["gamma"]
beta=parms["beta"]
mu=parms["mu"]
N=parms{"N"}

#endemic equilibrium
Sstar<-(gamma+mu)/beta
Istar<-mu*(beta/(gamma+mu)-1)/beta

eql<-list(S=Sstar,I=Istar)

dS<-expression(mu*(N-S)-beta*S*I/N)
dI<-expression(beta*S*I/N-(mu+gamma)*I)

j11 = D(dS,"S")
j12 = D(dS,"I")
j21 = D(dI, "S")
j22 = D(dI, "I")


