
####################################
library("epimdr")
library("rootSolve")
library("phaseR")

library("shiny")


ui=pageWithSidebar(
  #TheTitle
  headerPanel("The SIR model"),
  #sidebar
  sidebarPanel(
    sliderInput("beta","Transmission (yr^-1)",300, 
                min=0,max=1000),
    sliderInput("infper","Infectious Period (days)",5,
                min=1,max=100),
    sliderInput("mu","birth rate", 5,
                min=0, max=100),
    sliderInput("T","Time range:",
                min=0,max=1,value = c(0,1))
  ),
  mainPanel (
    tabsetPanel(
      tabPanel("Time", plotOutput("plot1")),
      tabPanel("Phase plane", plotOutput("plot2",
                                         height = 500)),
      tabPanel("Equations",
               withMathJax(
                 helpText("Susceptible $$\\frac{dS}{dt} =
                            \\mu (N - S) - \\frac{\\beta I S} {N} $$"),
                 helpText("Infected $$\\frac{dI}{dt} = 
                            \\frac{\\beta I S} {N} - (\\mu+\\sigma) I$$"),
                 helpText("Removed $$\\frac{dR}{dt} = 
                            \\gamma I - \\mu R $$"),
                 helpText("Reproductive ratio $$R_0 =
                            \\frac{1}{\\gamma+\\mu} \\frac{\\beta N} {N}$$")
               ))
    )
  )
)
server=function(input,output){
  sirmod=function(t,x,parms){
    S=x[1]
    I=x[2]
    R=x[3]
    
    beta=parms["beta"]
    mu=parms["mu"]
    gamma=parms["gamma"]
    N=parms["N"]
    dS = mu*(N-S)-beta*S*I/N
    dI = beta * S * I / N - (mu + gamma)* I
    dR = gamma * I - mu * R
    res= c(dS,dI,dR)
    list(res)
  }
  simod=function(t,y,parameters){
    S=y[1]
    I=y[2]
    beta=parameters["beta"]
    mu=parameters["mu"]
    gamma=parameters["gamma"]
    N=parameters["N"]
    dS=mu * (N-S) - beta * S * I/N
    dI= beta * S * I/N - (mu+gamma)*I
    res = c(dS,dI)
    list(res)
  }
  
  output$plot1=renderPlot({
    times = seq(0,input$T[2], by=1/1000)
    parms=c(mu = input$mu, N= 1, beta = input$beta,
            gamma = 365/input$infper)
    start = c(S=.999,I=.001, R=0)
    R0 = round(with(as.list(parms),beta/(gamma+mu)),1)
    
    out= ode(y=start,times=times,func=sirmod,
             parms=parms)
    out=as.data.frame(out)
    
    sel=out$time>input$T[1]&out$time<input$T[2]
    plot(x=out$time[sel],y=out$S[sel],ylab="fraction",
         xlab="time",type="l",ylim=range(out[sel,-c(1,4)]))
    title(paste("R0=",R0))
    
    lines(x=out$time[sel],y=out$I[sel],col="red")
    lines(x=out$time[sel],y=out$R[sel],col="green")
    legend("right",legend=c("S","I","R"),lty=c(1,1,1),
           col=c("black","red","green"))
    
  })
  
  output$plot2 = renderPlot({
    times= seq(0,input$T[2], by=1/1000)
    parms=c(mu=input$mu,N=1,beta=input$beta,
            gamma=365/input$infper)
    start= c(S=.999,I=.001,R=0)
    R0=round(with(as.list(parms),beta/(gamma+mu)),1)
    
    out=ode(y=start[-3],times=times, func=simod,
            parms=parms)
    out=as.data.frame(out)
    
    plot(x=out$S,y=out$I,xlab="Fraction susceptible",
         ylab="Fraction infected",type="l")
    title(paste("R0=",R0))
    fld=flowField(simod, xlim=range(out$S),ylim=range(out$I),
                  parameters=parms,system="two.dim",
                  add=TRUE, ylab="I",xlab="S")
    
    abline(v=1/R0,col="green")
    
    curve(parms["mu"]*(1-x)/(parms["beta"]*x),min(out$S),
          max(out$S),add=TRUE,col="red")
    legend("topright",legend=c("I-isocline","S-isocline"),lty=c(1,1),
           col=c("green","red"))
  })
}


shinyApp(ui,server)


