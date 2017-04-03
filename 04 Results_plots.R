
################################################################################
# "Evaluation of the impact of ambient temperatures on occupational injuries
#  in Spain (1994-2013)"
#   
#   ISGlobal, Center for Research in Environmental Epidemiology (CREAL) 
#   April 2017
#   
#	This code reproduces the analyses, which is based on the methodology 
#   Distributed Lag Non-Linear Models, developed by Antonio Gasparrini 
#	(http://www.ag-myresearch.com/)
#
###############################################################################

################################################################################
# PLOTS AND RESULTS
################################################################################

####################################################################
# FIGURE 1. OVERALL CUMULATIVE SUMMARY ASSOCIATION ####
# FIGURE 1A.  EXPOSURE RESPONSE-ASSOCIATION####

indlab <- predper%in%c(0,1,10,50,90,99,100)

plot(cp,ylab="Percent change (%)",xlab="Temperature percentile",axes=F,
     ylim=c(0.9,1.2),lwd=2,col="white")
ind1 <- cp$predvar<=cen
ind2 <- cp$predvar>=cen
lines(cp$predvar[ind1],cp$allRRfit[ind1],col=4,lwd=2)
lines(cp$predvar[ind2],cp$allRRfit[ind2],col=2,lwd=2)

axis(1,at=tmeancountry[indlab],labels=predper[indlab],cex.axis=0.9)
axis(2,cex.axis=0.9, at= seq(0.9,1.2, by=0.1), labels=c(-10,0,10,20))

abline(v=c(tmeancountry[cenindcountry],tmeancountry[c("1.0%","99.0%")]),
       lty=c(3,2,2))


	   
# FIGURE 1B.  LAG-RESPONSE ASSOCIATION FOR COLD####	   
plot(cplagcold,ylab="Percent change (%)",xlab="Lag, days",ylim=c(0.95,1.05),lwd=2,col="steelblue3",
     ci.arg=list(density=20,col="steelblue3"),cex.axis=0.9, yaxt='n')
axis(2,at=c(0.96,0.98,1.00,1.02,1.04), labels=c(-4,-2,0,2,4), cex.axis=0.9)

# FIGURE 1C.  LAG-RESPONSE ASSOCIATION FOR HEAT####	 
plot(cplagheat,ylab="Percent change (%)",xlab="Lag, days",ylim=c(0.95,1.05),lwd=2,col="red1",
     ci.arg=list(density=20,col="red1"),cex.axis=0.9, yaxt='n')
axis(2,at=c(0.96,0.98,1.00,1.02,1.04), labels=c(-4,-2,0,2,4), cex.axis=0.9)


################################################################################
# MAIN RESULTS
################################################################################

################################################################################


# RR AT 1st, 10th, 90TH AND 99TH VS MMP (WITH 95%CI)
predper <- c(seq(0,1,0.1),2:98,seq(99,100,0.1))

results <- matrix (NA, nrow=5, ncol=4)
results[,1] <- c("1", "10", "90", "99", "MWIP")

results[1,2] <- cp$allRRfit[predper==1];results[1,3] <- cp$allRRlow[predper==1];results[1,4] <- cp$allRRhigh[predper==1]
results[2,2] <- cp$allRRfit[predper==10];results[2,3] <- cp$allRRlow[predper==10];results[2,4] <- cp$allRRhigh[predper==10]
results[3,2] <- cp$allRRfit[predper==90];results[3,3] <- cp$allRRlow[predper==90];results[3,4] <- cp$allRRhigh[predper==90]
results[4,2] <- cp$allRRfit[predper==99];results[4,3] <- cp$allRRlow[predper==99];results[4,4] <- cp$allRRhigh[predper==99]

# Minimum occupational-injuries percentile
cenpercountry


################################################################################
#### SUPPLEMENTARY INFORMATION 
# FIGURE 1. OVERALL CUMULATIVE EXPOSURE-RESPONSE ASSOCIATIONS BY PROVINCE ####

xlab <- expression(paste("Temperature (",degree,"C)"))

pdf("figureS1.pdf",width=8,height=9)
layout(matrix(c(0,1,1,2,2,0,rep(3:8,each=2),0,9,9,10,10,0),ncol=6,byrow=T))
par(mar=c(4,3.8,3,2.4),mgp=c(2.5,1,0),las=1)

per <- t(sapply(dlist,function(x) 
  quantile(x$tempmax,c(1,2.5,10,25,50,75,90,97.5, 99)/100,na.rm=T)))


for(i in seq(length(dlist))) {
  data <- dlist[[i]]
  # 
  argvar <- list(x=data$tempmax,fun="bs",degree=2, knots=quantile(data$tempmax, 
                                                                  varper/100, na.rm=TRUE))
  
  bvar <- do.call(onebasis,argvar)
  pred <- crosspred(bvar,coef=blup3[[i]]$blup,vcov=blup3[[i]]$vcov,
                    model.link="log",by=0.1,cen=mintempcity[i])
  
  
  
  df <- data.frame(matrix(nrow=length(pred$predvar), ncol=4))
  colnames(df) <- c("temp", "RRfit", "RRlow", "RRhigh")
  df$temp <- pred$predvar
  df$RRfit <- (pred$allRRfit-1)*100
  df$RRlow <- (pred$allRRlow-1)*100
  df$RRhigh <- (pred$allRRhigh-1)*100
  
  plot(df$temp, df$RRfit, type="n", ylim=c(-55, 40), 
       ylab="Percent change (%)",xlab=xlab, lwd=2,col="white",
       cex.axis = 0.8, main=provincies_n[i], lab=c(6,5,7),
       yaxt="n", axes=F)
  polygon(c(df$temp,rev(df$temp)),c(df$RRlow,rev(df$RRhigh)),col = "grey89", border = FALSE)
  ind1 <- pred$predvar<=mintempcity[i]
  ind2 <- pred$predvar>=mintempcity[i]
  lines(df$temp[ind1],df$RRfit[ind1],col=4,lwd=2)
  lines(df$temp[ind2],df$RRfit[ind2],col=2,lwd=2)
  axis(2,at=-2:2.4*20, col.axis="black", cex.axis=0.8)
  axis(1, col.axis="black")
  abline(h=0)
  
  breaks <- c(min(data$tempmax,na.rm=T)-1,seq(pred$predvar[1],
                                              pred$predvar[length(pred$predvar)],
                                              length=30),max(data$tempmax,na.rm=T)+1)
  hist <- hist(data$tempmax,breaks=breaks,plot=F)
  hist$density <- hist$density/max(hist$density)*0.7
  prop <- max(hist$density)/max(hist$counts)
  counts <- pretty(hist$count,3)
  par(new=TRUE)
  plot(hist,ylim=c(0,max(hist$density)*2.0),axes=F,ann=F,col=grey(0.95),freq=F)
  axis(4,at=counts*prop,labels=counts,cex.axis=0.7)
  #mtext("N",4,line=-0.5,at=mean(counts*prop),cex=0.5)
  abline(v=mintempcity[i],lty=1)
  abline(v=c(per[i,c("1%","99%")]),lty=2)
}

dev.off(); 

par(opar) 











