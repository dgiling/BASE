  

  ### BAyesian Single-station Estimation (BASE) ###     
  ###    Grace et al. 2015 L&O: Methods         ### 

  # SCRIPT 1
  
  # required packages (must be downloaded and installed prior to running the model)
  library(coda)
  library(R2jags)
  library(zoo)
  
  #-------------------------------------------------------------
  
  # (A) SET LOCATION OF 'BASE' FOLDER (e.g. "C:/Desktop/Analysis")
  directory <- "C:/Desktop/Analysis"  # example

  #-------------------------------------------------------------
  
  # (B) SET MEASUREMENT INTERVAL (SECONDS)  (must be the same for all files in a run)
  interval <- 600
  
  #-------------------------------------------------------------
 
  # (C) MCMC settings
  n.iter <- 20000
  n.burnin <- n.iter*0.5
  
  #-------------------------------------------------------------
  
  # (D) OPTIONAL - SET SMOOTHING BEHAVIOUR
  
  # Dissolved oxygen data (fast Fourier transform)
  # Filter out this proportion of high-frequency fluctuations (0 = no smoothing)
  smooth.DO <- 0.0 # proportion

  # PAR data (moving average)
  smooth.PAR <- FALSE  # logical
  
  #-------------------------------------------------------------
  
  
  #############################################
  #      do not need to alter below here      #
  #############################################
  
  start.time<-NULL; start.time<-Sys.time()
  
  # Other functions
  smooth5 <- function(x) (rollapply(x, 5, mean, na.rm=T,align="center"))  # moving average of 5 time intervals
  
  # Data input and set up output table dataframes
  filenames<-list.files(file.path(directory,"BASE/input"))  
  
  # Set up output tables
  output.table<-NULL
  output.table<-data.frame(File=character(), Date=character(), GPP.mean=double(), GPP.sd=double(), ER.mean=double(), ER.sd=double(), 
                           K.mean=double(), K.sd=double(), theta.mean=double(), theta.sd=double(), A.mean=double(), A.sd=double(), p.mean=double(), p.sd=double(), 
                           R2=double(), PPP=double(), rmse=double(), rmse.relative=double(), mrl.fraction=double(), ER.K.cor=double(), convergence.check=double(), A.Rhat=double(),
                           K.Rhat=double(), theta.Rhat=double(), p.Rhat=double(), R.Rhat=double(), GPP.Rhat=double(), DIC=double(), pD=double(), smooth.DO=double(),
                           stringsAsFactors=FALSE)
  instant.rates<-data.frame(File=character(), Date=character(), interval=integer(), 
                            tempC=double(), I=double(), K.instant=double(), GPP.instant=double(), ER.instant=double(),
                            stringsAsFactors=FALSE)
  
  # Analyse files sequentially
  for (fname in filenames) {
    
      data<-read.csv(file.path(directory,"/BASE/input/",fname), head=T) # read next file
      seconds<-86400
      N = nrow(data)
      x = 0:(N-1)
     
      # check dates for "/" and replace with "-"
      data$Date <- gsub("/", "-", data$Date)
      
      ## Smoothing data
      plot.days <- 2       # number of consecutive days to plot for the smoothing plot
      if(smooth.DO > 0) {
        
        # fast Fourier transform smoothing - low pass filter
        DO.fft = fft(data$DO.meas)
        inx_filter = floor(N/2*(1-smooth.DO))
        filter = rep(1, N)
        filter[inx_filter:(N-inx_filter)] = 0
        DO.fft_filtered = filter * DO.fft
        data$DO.smooth <- Re( fft( DO.fft_filtered, inverse=TRUE) / N )
        noise <- data$DO.meas - data$DO.smooth

        limit <- (seconds/interval) * plot.days
        jpeg(file=file.path(directory,"BASE/output/validation plots", paste0(fname, "_DO_smoothing.jpg")), width=1200, height=600, pointsize=20, quality=400)
        plot(data$DO.meas[1:limit], typ='l', col='grey', ylab="DO concentration", xlab="Timesteps", lwd=4)
        points(noise[1:limit]+mean(data$DO.meas[1:limit]), typ='l', col="blue", lwd=2)
        points(data$DO.smooth[1:limit] , typ='l', col="red", lwd=2)
        legend("topleft", legend=c("Measured", "FFT smooth", "High-frequency noise (centered)"), col=c("grey", "red", "blue"), lwd=c(3,2,2), bty='n')
        graphics.off()
      }
  
      if(smooth.PAR == T) {
        data$I.smooth<-c(data$I[1:2],smooth5(data$I),data$I[nrow(data)-1],data$I[nrow(data)]) # moving average over 5 time intervals
        
        jpeg(file=file.path(directory,"BASE/output/validation plots", paste0(fname, "_PAR_smoothing.jpg")), width=1200, height=600, pointsize=20, quality=400)
        plot(data$I[1:limit], typ='l', col='grey', ylab="DO concentration", xlab="Timesteps", lwd=4)
        points(data$I.smooth[1:limit] , typ='l', col="red", lwd=2)
        legend("topleft", legend=c("Measured", "Smoothed"), col=c("grey", "red"), lwd=c(3,2), bty='n')
        graphics.off()
        
      }
      
      # Select dates
      data$Date <- factor(data$Date, levels = unique(data$Date))
      dates <- unique(data$Date)
      n.records <- tapply(data$Date, INDEX=data$Date, FUN=length)
      dates <- dates[n.records == (seconds/interval)] # select only dates with full days
      
      ## Analyse days sequentially
      for (d in dates) 
          { 
          data.sub <- data[data$Date == d,]
          
          # Define data vectors
        	num.measurements <- nrow(data.sub)
        	tempC <- data.sub$tempC
        	salinity <-data.sub$salinity
        	atmo.pressure <- data.sub$atmo.pressure
        	DO.meas <- if(smooth.DO > 0) data.sub$DO.smooth else data.sub$DO.meas
        	PAR     <- if(smooth.PAR == TRUE) data.sub$I.smooth else data.sub$I
        
        	# Initial values
        	# Set these to something sensible if the model is becoming stuck in a bad parameter space
        	# These values here are expressed per timestep, not per day. Divide desired initial K (/day) by the number of timesteps in a day, as shown in default below 
        	inits <- function()	{	list(K= 2 / (86400/interval) ) }

        	# Different random seeds
        	kern=as.integer(runif(1000,min=1,max=10000))
        	iters=sample(kern,1)
        	n.chains <- 3
        	n.thin <- 10
        	data.list <- list("num.measurements","interval","tempC","DO.meas","PAR","salinity","atmo.pressure")  
        	
        	# Define monitoring variables
        	params=c("A","R","K","K.day","p","theta","tau","ER","GPP", "NEP","sum.obs.resid","sum.ppa.resid","PPfit","DO.modelled")
          
        	## Call jags ##
        	
        	# Set debug = T below to inspect each file for model convergence 
        	# (inspect the main parameters for convergence using bgr diagrams, history, density and autocorrelation)
        	metabfit=NULL
        	metabfit <- do.call(jags.parallel,
        	                 list(data=data.list, inits=inits, parameters.to.save=params, model.file = "BASE_metab_model_v2.3.txt",
        	                      n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin,
        	                      n.thin = n.thin, n.cluster= n.chains, DIC = TRUE,
        	                      working.directory = file.path(directory,"/BASE"), jags.seed = 123, digits=5))
    	
        	# print(metabfit, digits=2) # to inspect results of last metabfit
          
        	## diagnostic summaries
        	# Rhat (srf) test
        	srf<- metabfit$BUGSoutput$summary[,8]
        	Rhat.test <- NULL
        	Rhat.test <- ifelse(any(srf>1.1, na.rm=T)==TRUE,"Check convergence", "Fine")
        	
            	# Check for convergence and update once if required
            	if(Rhat.test == "Check convergence") {
            	  recompile(metabfit)
            	  metabfit <- update(metabfit, n.iter=n.iter)    
            	}
        	
        	# Rhat (srf) test - second round in case metabfit is updated
        	srf<- metabfit$BUGSoutput$summary[,8]
        	Rhat.test <- NULL
        	Rhat.test <- ifelse(any(srf>1.1, na.rm=T)==TRUE,"Check convergence", "Fine")
        	##
        	
        	# autocorr test
        	metabfit.mcmc<-as.mcmc(metabfit)
        	ac.lag1 <- autocorr.diag(metabfit.mcmc, lags = 1)
        	auto.corr.test <- NULL
        	auto.corr.test <- ifelse(any(abs(ac.lag1)>0.1, na.rm=T)==TRUE,"Check ac", "ac OK")
        	
        	PPP <- metabfit$BUGSoutput$summary["PPfit","mean"] # posterior predictive p-value
        	
        	DO.mod.means <- metabfit$BUGSoutput$mean$DO.modelled
        	DO.mod.sd <- metabfit$BUGSoutput$sd$DO.modelled
      
        	R2 = cor(DO.mod.means,DO.meas)^2
          rmse = sqrt(sum((metabfit$BUGSoutput$mean$DO.modelled-DO.meas)^2)/length(DO.meas))
          post.mean.dev <- metabfit$BUGSoutput$mean$deviance
          pD <- metabfit$BUGSoutput$pD
          DIC <- metabfit$BUGSoutput$DIC
          
          DO.lag<-DO.meas[2:length(DO.meas)]-DO.meas[1:(length(DO.meas)-1)]
          ptpvar <- sqrt((sum((DO.lag)^2)/(length(DO.meas)-1))) # point to point variation
          rmse.relative <- rmse / ptpvar
          
          diff<-metabfit$BUGSoutput$mean$DO.modelled-DO.meas
        	mrl.max<-max(rle(sign(as.vector(diff)))$lengths)
        	mrl.fraction<-max(rle(sign(as.vector(diff)))$lengths)/length(DO.meas) # proportion of largest run
                  
        	ER.K.cor <- cor(metabfit$BUGSoutput$sims.list$ER,metabfit$BUGSoutput$sims.list$K) # plot(metabfit$sims.list$ER ~ metabfit$sims.list$K)
          
        	# insert results to table and write table
        	result <- c(as.character(fname), as.character(d), metabfit$BUGSoutput$mean$GPP, metabfit$BUGSoutput$sd$GPP, metabfit$BUGSoutput$mean$ER, metabfit$BUGSoutput$sd$ER, metabfit$BUGSoutput$mean$K.day, 
        			metabfit$BUGSoutput$sd$K.day,  metabfit$BUGSoutput$mean$theta, metabfit$BUGSoutput$sd$theta, metabfit$BUGSoutput$mean$A, metabfit$BUGSoutput$sd$A, metabfit$BUGSoutput$mean$p, metabfit$BUGSoutput$sd$p, 
              R2, PPP, rmse, rmse.relative, mrl.fraction, ER.K.cor, Rhat.test, metabfit$BUGSoutput$summary["A",8] , metabfit$BUGSoutput$summary["K",8], 
        			metabfit$BUGSoutput$summary["theta",8], metabfit$BUGSoutput$summary["p",8], metabfit$BUGSoutput$summary["R",8], metabfit$BUGSoutput$summary["GPP",8],  DIC, pD, smooth.DO)
        	output.table[nrow(output.table)+1,] <- result
        	write.csv(output.table, file=file.path(directory,"BASE/output/BASE_results.csv")) # output file overwritten at each iteration
        	
        	# insert results to instantaneous table and write
        	instant.result <- data.frame(File=as.character(rep(fname,seconds/interval)), Date=as.character(rep(d,seconds/interval)),interval=1:(seconds/interval),
        	                             tempC=tempC, I=PAR, 
        	                             K.instant=as.vector(metabfit$BUGSoutput$mean$K) * 1.0241^(tempC-mean(tempC)),
        	                             GPP.instant=as.vector(metabfit$BUGSoutput$mean$A) * PAR^as.vector(metabfit$BUGSoutput$mean$p),
        	                             ER.instant=as.vector(metabfit$BUGSoutput$mean$R) * as.vector(metabfit$BUGSoutput$mean$theta)^(tempC-mean(tempC)),
        	                             stringsAsFactors = FALSE)
        	instant.rates[(nrow(instant.rates)+1):(nrow(instant.rates)+(seconds/interval)),] <- instant.result
        		write.csv(instant.rates, file=file.path(directory,"BASE/output/instantaneous rates/instantaneous_rates.csv")) # output file name
        	
        	# diagnostic multi-panel plot
          jpeg(file=file.path(directory,"BASE/output/validation plots", paste0(substr(fname, 1,(nchar(fname)-4)),"_", as.character(d), ".jpg")), width=1200, height=1200, pointsize=30, quality=250)
        	
        	traceplot(metabfit, varname=c('A','p','R','K.day','theta'), ask=FALSE, mfrow=c(3,3), mar=c(2,2,0,8), new=FALSE)
        	
        	plot(1:num.measurements,DO.mod.means, type="l",lwd=2, ylim=c(min(DO.mod.means-DO.mod.sd)-0.5,max(DO.mod.means+DO.mod.sd)+0.5), xlab="Timestep")
        	points(1:num.measurements,DO.meas,pch=1,xlab="Timestep", col="grey60", cex=0.75)  
        	points(1:num.measurements,DO.mod.means+DO.mod.sd, type="l", lty=2)
        	points(1:num.measurements,DO.mod.means-DO.mod.sd, type="l", lty=2)
        	legend(x="topleft", legend=c("DO meas or smoothed", "DO modelled"), pch=c(1,NA), lty=c(NA,1), col=c("grey60", "black"), cex=0.6, bty='n')
        	
        	plot(1:num.measurements,tempC,pch=1,xlab="Timestep" , typ='p')
        	plot(1:num.measurements,PAR,pch=1,xlab="Timestep" , typ='p')
        	
        	graphics.off()
        	
         	}
  }
    
  end.time<-NULL; end.time<-Sys.time()
  end.time-start.time
  
  
  
  
  