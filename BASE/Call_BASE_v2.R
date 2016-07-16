  
  
  ### GRACE et al. 2015 'BASE' CODE ###
  
  # SCRIPT 1
  
  # required packages (must be downloaded and installed prior to running the model)
  library(coda)
  library(R2jags)
  library(zoo)
  
  #-------------------------------------------------------------
  
  # (A) SET LOCATION OF 'BASE' FOLDER (e.g. "C:/Desktop/Analysis")
  folder.location <- "[your directory]" 

  #-------------------------------------------------------------
  
  # (B) SET MEASUREMENT INTERVAL (SECONDS)  (must be the same for all files in a run)
  interval <- 600
  
  #-------------------------------------------------------------
 
  # (C) MCMC settings
  n.iter <- 20000   
  n.burnin <- n.iter*0.5
  

  #################################################################
  #      do not alter below here (expect to set debug = T)        #
  #################################################################
  
  start.time<-NULL; start.time<-Sys.time()
  
  wd1<- paste(folder.location, "/BASE", sep="")
  wd2<- paste(folder.location, "/BASE/input", sep="")
  wd3<- paste(folder.location, "/BASE/output", sep="")
  wd4<- paste(folder.location, "/BASE/output/validation plots", sep="")
  
  # set up output table array
  sec.per.day<-86400
  K.mean <- 0 #(reaeration.prior * interval)/sec.per.day  # 5/day 
  R.mean <- 0 #(respiration.prior * interval)/sec.per.day    # 5/day
  output.table<-NULL
  output.table<-data.frame(File=NA, GPP.mean=NA, GPP.sd=NA, ER.mean=NA, ER.sd=NA, K.mean=NA, K.sd=NA, theta.mean=NA, theta.sd=NA, p.mean=NA, p.sd=NA, 
                           R2=NA, PPP=NA, rmse=NA, rmse.relative=NA, mrl.fraction=NA, ER.K.cor=NA, convergence.check=NA, A.Rhat=NA,
                           K.Rhat=NA, theta.Rhat=NA, p.Rhat=NA, R.Rhat=NA, GPP.Rhat=NA, DIC=NA, pD=NA)
  
  
  
  # Start files looping
  setwd(wd2)
  filenames=list.files()  

  for (fname in filenames) 
    { 
  
    # read next file
  	setwd(wd2)
  	# fname <- filenames[1]
    data<-read.csv(fname, sep=",", header=T)
  
  	# define data vectors
  	num.measurements<-nrow(data)
  	tempC<-data$tempC
  	DO.meas<-data$DO.meas
  	PAR<-data$I
  	salinity<-data$salinity
  	atmo.pressure<-data$atmo.pressure
  
  	## could smooth here but better for whole dataset?
  	#library(zoo)
  	#smooth<- function(x) (rollapply(x, 5, mean, na.rm=T,align="center"))  
  	#I<-c(I[1:2],smooth(I),I[length(I)-1],I[length(I)] )
  	
  	inits <- function()
  		{
  		list(sd=0.1)
  		}
  
  	# different random seeds
  	kern=as.integer(runif(1000,min=1,max=10000))
  	iters=sample(kern,1)
  	n.chains <- 3
  	n.thin <- 10
  	data.list <- list("num.measurements","interval","tempC","DO.meas","PAR","salinity","atmo.pressure")  
  	
  	# define monitoring variables
  	params=c("A","R","K","K.day","p","theta","sd","ER","GPP", "NEP","sum.obs.resid","sum.ppa.resid","PPfit","DO.modelled", "P1", "P2")
    
  	## call jags
  	
  	# Set debug = T below to inspect each file for model convergence 
  	# (inspect the main parameters for convergence using bgr diagrams, history, density and autocorrelation)
  	setwd(wd1) 
  	metab=NULL
  	metab <- do.call(jags.parallel,
  	                 list(data=data.list, inits=inits, parameters.to.save=params, model.file = "BASE_metab_model_v2.txt",
  	                      n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin,
  	                      n.thin = n.thin, n.cluster= n.chains, DIC = TRUE,
  	                      working.directory = wd1, jags.seed = 123, digits=5))  # setwd
  	
  	# print(metab, digits=2) # to inspect results of last "metab" estimate
    
  	## diagnostic summaries
  	# Rhat (srf) test
  	srf<- metab$BUGSoutput$summary[,8]
  	Rhat.test <- NULL
  	Rhat.test <- ifelse(any(srf>1.1, na.rm=T)==TRUE,"Check convergence", "Fine")
  	
  	# autocorr test
  	metab.mcmc<-as.mcmc(metab)
  	ac.lag1 <- autocorr.diag(metab.mcmc, lags = 1)
  	auto.corr.test <- NULL
  	auto.corr.test <- ifelse(any(abs(ac.lag1)>0.1, na.rm=T)==TRUE,"Check ac", "ac OK")
  	
  	PPP <- metab$BUGSoutput$summary["PPfit","mean"] # posterior predictive p-value # updated
  	
  	DO.mod.means <- metab$BUGSoutput$summary[2:(num.measurements+1),"mean"] 
  	DO.mod.sd <- metab$BUGSoutput$summary[2:(num.measurements+1),"sd"]

  	R2 = cor(DO.mod.means,DO.meas)^2
    rmse = sqrt(sum((metab$BUGSoutput$mean$DO.modelled-DO.meas)^2)/length(DO.meas))
    post.mean.dev <- metab$BUGSoutput$mean$deviance
    pD <- metab$BUGSoutput$pD
    DIC <- metab$BUGSoutput$DIC
    
    DO.lag<-DO.meas[2:length(DO.meas)]-DO.meas[1:(length(DO.meas)-1)]
    ptpvar <- sqrt((sum((DO.lag)^2)/(length(DO.meas)-1))) # point to point variation
    rmse.relative <- rmse / ptpvar
    
    diff<-metab$BUGSoutput$mean$DO.modelled-DO.meas
  	mrl.max<-max(rle(sign(as.vector(diff)))$lengths)
  	mrl.fraction<-max(rle(sign(as.vector(diff)))$lengths)/length(DO.meas) # prop of largest run
            
  	ER.K.cor <- cor(metab$BUGSoutput$sims.list$ER,metab$BUGSoutput$sims.list$K) # plot(metab$sims.list$ER ~ metab$sims.list$K)
    
  	# append results to table and write table
  	result <- c(fname, metab$BUGSoutput$mean$GPP, metab$BUGSoutput$sd$GPP, metab$BUGSoutput$mean$ER, metab$BUGSoutput$sd$ER, metab$BUGSoutput$mean$K.day, 
  			metab$BUGSoutput$sd$K.day,  metab$BUGSoutput$mean$theta, metab$BUGSoutput$sd$theta, metab$BUGSoutput$mean$p, metab$BUGSoutput$sd$p, 
        R2, PPP, rmse, rmse.relative, mrl.fraction, ER.K.cor, Rhat.test, metab$BUGSoutput$summary["A",8] , metab$BUGSoutput$summary["K",8], 
  			metab$BUGSoutput$summary["theta",8], metab$BUGSoutput$summary["p",8], metab$BUGSoutput$summary["R",8], metab$BUGSoutput$summary["GPP",8],  DIC, pD)
  	
  	output.table<-rbind(output.table,result)
  	setwd(wd3)
  	write.csv(output.table, "BASE_results.csv")  # output file name
  
  	
  	# diagnostic multi-plot
  	setwd(wd4)
  	jpeg(file=paste(fname, ".jpg", sep=""), width=1200, height=1200, pointsize=30, quality=250)
  	
  	# fit plot
  	traceplot(metab, varname=c('A','p','R','K.day','theta'), ask=FALSE, mfrow=c(3,3), mar=c(2,2,0,8), new=FALSE)
  	
  	plot(1:num.measurements,DO.mod.means, type="l",lwd=2, ylim=c(min(DO.mod.means-DO.mod.sd)-0.5,max(DO.mod.means+DO.mod.sd)+0.5), xlab="Timestep")
  	points(1:num.measurements,DO.meas,pch=1,xlab="Timestep", col="grey60", cex=0.75)  
  	points(1:num.measurements,DO.mod.means+DO.mod.sd, type="l", lty=2)
  	points(1:num.measurements,DO.mod.means-DO.mod.sd, type="l", lty=2)
  	
  	plot(1:num.measurements,tempC,pch=1,xlab="Timestep" , typ='p')
  	plot(1:num.measurements,PAR,pch=1,xlab="Timestep" , typ='p')
  	
  	graphics.off()
  	
   	}
  
  # traceplot(metab, var=c("K.day","R","p", "A", "sd", "P1", "P2"), mfrow=c(3,3),ask=FALSE)
  
  # final csv write
  output.table<-output.table[-1,]
  setwd(wd3)
  write.csv(output.table, "BASE_results.csv", row.names=F)  # output file name
  end.time<-NULL; end.time<-Sys.time()
  start.time
  end.time
  end.time-start.time
  
  
  
  
  
