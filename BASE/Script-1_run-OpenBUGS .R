  
  
  ### GRACE et al. 2015 'BASE' CODE ###
  
  # SCRIPT 1
  
  # required packages (must be downloaded and installed prior to running the model)
  library(coda)
  library(R2OpenBUGS)
  
  #-------------------------------------------------------------
  
  # (A) SET LOCATION OF 'BASE' FOLDER (e.g. "C:/Desktop/Analysis")
  folder.location <-"[your directory]"
    
  #-------------------------------------------------------------
  
  # (B) SET MEASUREMENT INTERVAL (SECONDS)  (must be the same for all files in a run)
  interval<-300
  
  #-------------------------------------------------------------
  
  # (C) SET NUMBER OF MODEL ITERATIONS (DEFAULT 2000/1000)
  n.iter<-2000
  n.burnin<-1000
  
  #-------------------------------------------------------------
  
  
  
  #################################################################
  #      do not alter below here (expect to set debug = T)        #
  #################################################################
  
  start.time<-NULL; start.time<-Sys.time()
  
  wd1<- paste(folder.location, "/BASE", sep="")
  wd2<- paste(folder.location, "/BASE/input", sep="")
  wd3<- paste(folder.location, "/BASE/output", sep="")
  wd4<- paste(folder.location, "/BASE/output/fitting plots", sep="")
  
  # set up output table array
  sec.per.day<-86400
  K.mean <- 0 #(reaeration.prior * interval)/sec.per.day  # 5/day 
  R.mean <- 0 #(respiration.prior * interval)/sec.per.day    # 5/day
  File <- 0; GPP.mean <- 0; ER.mean <- 0; theta.mean <- 0; theta.sd <- 0; p.mean <- 0; p.sd<- 0; PPfit.mean <- 0
  GPP.sd <- 0; ER.sd <- 0; K.sd <- 0; R2<- 0; rmse<- 0; rmse.relative<- 0; mrl.fraction<- 0; ER.K.cor<-0;
  convergence.check <- 0; A.Rhat <- 0; R.Rhat <- 0; K.Rhat <- 0; theta.Rhat<- 0; p.Rhat<- 0; R.Rhat<- 0; gpp.Rhat<- 0
  DIC <- 0; pD <- 0
  
  output.table<-NULL
  output.table<-data.frame(File,GPP.mean , GPP.sd, ER.mean, ER.sd, K.mean, K.sd, theta.mean, theta.sd, 
  p.mean, p.sd, PPfit.mean, R2, rmse, rmse.relative, mrl.fraction, ER.K.cor, convergence.check, A.Rhat, R.Rhat, K.Rhat, theta.Rhat, p.Rhat, gpp.Rhat, DIC, pD)
  
  # define monitoring variables
  params=c("A","R","K","p","theta","sd","ER","gpp","sum.obs.resid","sum.ppa.resid","PPfit","DO.modelled")
  
  # Start model looping
  setwd(wd2)
  filenames=list.files()
  
  for (fname in filenames) 
  	{ 
  
    # read next file
  	setwd(wd2)
  	data<-read.csv(fname, sep=",", header=T)
  
  	# define data vectors
  	num.measurements<-nrow(data)
  	tempC<-data$tempC
  	DO.meas<-data$DO.meas
  	I<-data$I
  	salinity<-data$salinity
  	atmo.pressure<-data$atmo.pressure
  
  	inits <- function()
  		{
  		list(sd=0.1)
  		}
  
  	# different random seeds
  	kern=as.integer(runif(1000,min=1,max=10000))
  	iters=sample(kern,1)
  	
    ## TEMP
  	#library(zoo)
  	#smooth<- function(x) (rollapply(x, 5, mean, na.rm=T,align="center"))  
    #I<-c(I[1:2],smooth(I),I[length(I)-1],I[length(I)] )
    ## TEMP
    
  	data.list <- list("num.measurements","interval","tempC","DO.meas","I","salinity","atmo.pressure")  
  	
    # call OpenBUGS and store results in object "metab"
  	
  	# Set debug = T below to inspect each file for model convergence 
  	# (inspect the main parameters for convergence using bgr diagrams, history, density and autocorrelation)
  	setwd(wd1) 
  	metab=NULL
  	metab=bugs(data.list,inits,parameters.to.save=params,model.file="Script-2_metab-model.txt", n.thin=1,n.iter=n.iter,n.burnin=n.burnin, n.chains=3, 
               debug=F)  # <---------- DEBUG ARGUMENT -----------
    
  	# print(metab, digits=2) # to inspect results of last "metab" estimate
    
  	# append results to table and write table
  	R2 = cor(metab$mean$DO.modelled,DO.meas)^2
    rmse = sqrt(sum((metab$mean$DO.modelled-DO.meas)^2)/length(DO.meas))
    	
    DO.lag<-DO.meas[2:length(DO.meas)]-DO.meas[1:(length(DO.meas)-1)]
    ptpvar <- sqrt((sum((DO.lag)^2)/(length(DO.meas)-1))) # point to point variation
    rmse.relative <- rmse / ptpvar
    
    diff<-metab$mean$DO.modelled-DO.meas
  	mrl.max<-max(rle(sign(as.vector(diff)))$lengths)
  	mrl.fraction<-max(rle(sign(as.vector(diff)))$lengths)/length(DO.meas) # prop of largest run
            
  	ER.K.cor <- cor(metab$sims.list$ER,metab$sims.list$K) # plot(metab$sims.list$ER ~ metab$sims.list$K)
    
  	test<-NULL;convergence.check<-NULL
  	test<- any(c(metab$summary[1,8] , metab$summary[2,8], metab$summary[3,8], metab$summary[5,8], metab$summary[4,8], metab$summary[8,8])>1.1)
  	convergence.check<- ifelse(test==TRUE, "Check mixing", "Fine")
  
  	result <- c(fname, metab$mean$gpp, metab$sd$gpp, metab$mean$ER, metab$sd$ER, metab$mean$K*(sec.per.day/interval), 
  			metab$sd$K*(sec.per.day/interval),  metab$mean$theta, metab$sd$theta, metab$mean$p, metab$sd$p, 
        metab$mean$PPfit, R2, rmse, rmse.relative, mrl.fraction, ER.K.cor, convergence.check,
  			metab$summary[1,8] , metab$summary[2,8], metab$summary[3,8], metab$summary[5,8], metab$summary[4,8], metab$summary[8,8], 
        metab$DIC, metab$pD)
  	
  	output.table<-rbind(output.table,result)
  	setwd(wd3)
  	write.csv(output.table, "BASE_results.csv")  # output file name
  
  	# fitting plot and printing
  	setwd(wd4)
  	jpeg(file=paste(fname, ".jpg", sep=""), width=600, height=1000, quality=250, pointsize=30)
  	
    par(mfrow=c(3,1), mar=c(4,4,1,1))
    abscissa=1:num.measurements
  	plot(abscissa,DO.meas,pch=1 ,ylab="DO mg/l", col="grey50", xlab="",ylim=c(min(min(DO.meas),min(metab$mean$DO.modelled)), max(max(DO.meas),max(metab$mean$DO.modelled))), las=1, bty="l")
  	points(abscissa,metab$mean$DO.modelled, type="l", lwd=2)
  	legend("topleft", c("data","fit"), col=c("grey50","black"), lwd=c(1,2), lty=c(0,1), pch=c(1,NA),bty="n")
        
  	plot(1:num.measurements,tempC,pch=1 , typ='p', col="grey50", xlab="", las=1, bty="l")
  	plot(1:num.measurements,I,pch=1, typ='p', col="grey50",xlab="Timestep",ylab="PAR", las=1, bty="l")
    
    dev.off()
  
  	}
  
  # final csv write
  output.table<-output.table[-1,]
  setwd(wd3)
  write.csv(output.table, "BASE_results.csv", row.names=F)  # output file name
  end.time<-NULL; end.time<-Sys.time()
  start.time
  end.time
  end.time-start.time
  
  
  