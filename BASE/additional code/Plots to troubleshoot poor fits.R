

# run this code after a model run to plot measured DO (mg L-1 & % sat), modelled DO, PAR and Temp from the last iteration.


# sat calcs
kelvin<-0; DO.sat.pre<-0 ;DO.sat <-0; C <- 0; D<-0
S1<-0;S2<-0;S3<-0;S4<-0;sal.factor<-0;DOsalinity.corr<-0;alpha<-0;beta<-0;gamma<-0

for (i in 1:num.measurements)
	{
		kelvin[i] <- 273.15 + tempC[i]
		
		# correction for salinity
		
		S1[i] <- 157570.1 / kelvin[i]
		S2[i] <- -6.6423080E7 / (kelvin[i] * kelvin[i])
		S3[i] <- 1.2438E10 / (kelvin[i]^3)
		S4[i] <-  -8.621949E11 / (kelvin[i]^4)
		sal.factor[i] <- -1.0*salinity[i]*(0.017674-10.754/kelvin[i]+2140.7/(kelvin[i]*kelvin[i]))
		
		DOsalinity.corr[i] <- 0.0319988*exp(-135.90205+S1[i]+S2[i]+S3[i]+S4[i]+sal.factor[i])

		# correction for atmospheric pressure
		alpha[i] <- 0.009672-0.00004942*kelvin[i]+0.00000006436*(kelvin[i]^2)
		beta[i] <- exp(11.8571-3840.7/kelvin[i]-216961/(kelvin[i]^2))
		gamma[i] <- ((1-beta[i]/atmo.pressure[i]) / (1-beta[i])) * ((1-alpha[i]*atmo.pressure[i])/(1-alpha[i]))

		DO.sat[i] <- DOsalinity.corr[i]*atmo.pressure[i]*gamma[i]		

		# DO deficit or surfeit
		D[i] <- DO.sat[i]-DO.meas[i]
		
	}



DO.mod.means <- metab$mean$DO.modelled
DO.mod.sd <- metab$sd$DO.modelled


# fit plot
par(mfrow=c(2,2), mar=c(4,4,1,1))

plot(1:num.measurements, DO.mod.means, ylim=c(5.5,9.0), type="l", xlab="Timestep")
points(1:num.measurements,DO.mod.means+DO.mod.sd, type="l", lty=2)
points(1:num.measurements,DO.mod.means-DO.mod.sd, type="l", lty=2)
points(1:num.measurements,DO.meas,pch=1,xlab="Timestep" , col="grey")

plot(1:num.measurements,(DO.meas/DO.sat*100),pch=1,xlab="Timestep" , col="grey")
plot(1:num.measurements,tempC,pch=1,xlab="Timestep" , col="grey", typ='l')
plot(1:num.measurements,I,pch=1,xlab="Timestep" , typ='l')






