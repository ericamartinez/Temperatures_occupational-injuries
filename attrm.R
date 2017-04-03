
################################################################################
# "Evaluation of the impact of ambient temperatures on occupational injuries
#  in Spain (1994-2013)"
#   
#   ISGlobal, Center for Research in Environmental Epidemiology (CREAL) 
#   April 2017
#   
#	This code reproduces the analyses, which is based on the methodoly 
#   Distributed Lag Non-Linear Models, developed by Antonio Gasparrini 
#	(http://www.ag-myresearch.com/)
#
###############################################################################


# COMPUTE THE ATTRIBUTABLE INJURIES FOR EACH PROVINCE, WITH EMPIRICAL CI
# ESTIMATED USING THE RE-CENTERED BASES
################################################################################

# LOAD THE FUNCTION FOR COMPUTING THE ATTRIBUTABLE RISK MEASURES
source("attrdl.R")

# CREATE THE VECTORS TO STORE THE TOTAL INJURIES (ACCOUNTING FOR MISSING)
totdeath <- rep(NA,length(provinces))
names(totdeath) <- provincies_n

# CREATE THE MATRIX TO STORE THE ATTRIBUTABLE INJURIES (simulations)
matsim <- matrix(NA,length(provinces),7,dimnames=list(provincies_n,
  c("glob","cold","heat","extreme cold","moderate cold","moderate heat","extreme heat")))

# NUMBER OF SIMULATION RUNS FOR COMPUTING EMPIRICAL CI
nsim <- 1000

# CREATE THE ARRAY TO STORE THE CI OF ATTRIBUTABLE INJURIES
arraysim <- array(NA,dim=c(length(provinces),7,nsim),dimnames=list(provincies_n,
  c("glob","cold","heat","extreme cold","moderate cold","moderate heat","extreme heat")))

################################################################################

# RUN THE LOOP
for(i in 1:length(dlist)){
  
  # PRINT
  cat(i,"")
  
  # EXTRACT THE DATA
  data <- dlist[[i]]
  
  # COMPUTE PERCENTILES 2.5 AND 97.5
  perc <- quantile(data$tempmax_compl,c(0.025,0.975))
  
  # DERIVE THE CROSS-BASIS
  # NB: CENTERING POINT DIFFERENT THAN ORIGINAL CHOICE OF 75TH
  cb <- crossbasis(data$tempmax_compl, lag=lag, argvar=list(fun="bs", degree=2, 
    knots=quantile(data$tempmax_compl, varper/100, na.rm=TRUE)), 
    arglag=arglag)
  
  # COMPUTE THE ATTRIBUTABLE INJURIES
  # NB: THE REDUCED COEFFICIENTS ARE USED HERE
  matsim[i,"glob"] <- attrdl(data$tempmax_compl,cb,data$atotal,coef=blup[[i]]$blup,
    vcov=blup[[i]]$vcov,type="an",dir="forw",cen=mintempcity[i])
  matsim[i,"cold"] <- attrdl(data$tempmax_compl,cb,data$atotal,coef=blup[[i]]$blup,
    vcov=blup[[i]]$vcov,type="an",dir="forw",cen=mintempcity[i],
    range=c(-100,mintempcity[i]))
  matsim[i,"heat"] <- attrdl(data$tempmax_compl,cb,data$atotal,coef=blup[[i]]$blup,
    vcov=blup[[i]]$vcov,type="an",dir="forw",cen=mintempcity[i],
    range=c(mintempcity[i],100))
  
  matsim[i,"extreme cold"] <- attrdl(data$tempmax_compl,cb,data$atotal,coef=blup[[i]]$blup,
                                     vcov=blup[[i]]$vcov,type="an",dir="forw",cen=mintempcity[i],
                                     range=c(-100,perc[1]))
  matsim[i,"moderate cold"] <- attrdl(data$tempmax_compl,cb,data$atotal,coef=blup[[i]]$blup,
                                     vcov=blup[[i]]$vcov,type="an",dir="forw",cen=mintempcity[i],
                                     range=c(perc[1],mintempcity[i]))
  matsim[i,"moderate heat"] <- attrdl(data$tempmax_compl,cb,data$atotal,coef=blup[[i]]$blup,
                                     vcov=blup[[i]]$vcov,type="an",dir="forw",cen=mintempcity[i],
                                     range=c(mintempcity[i],perc[2]))
  matsim[i,"extreme heat"] <- attrdl(data$tempmax_compl,cb,data$atotal,coef=blup[[i]]$blup,
                                      vcov=blup[[i]]$vcov,type="an",dir="forw",cen=mintempcity[i],
                                      range=c(perc[2],100))
  
  # COMPUTE EMPIRICAL OCCURRENCES OF THE ATTRIBUTABLE INJURIES
  # USED TO DERIVE CONFIDENCE INTERVALS
  arraysim[i,"glob",] <- attrdl(data$tempmax_compl,cb,data$atotal,coef=blup[[i]]$blup,
    vcov=blup[[i]]$vcov,type="an",dir="forw",cen=mintempcity[i],sim=T,nsim=nsim)
  arraysim[i,"cold",] <- attrdl(data$tempmax_compl,cb,data$atotal,coef=blup[[i]]$blup,
    vcov=blup[[i]]$vcov,type="an",dir="forw",cen=mintempcity[i],
    range=c(-100,mintempcity[i]),sim=T,nsim=nsim)
  arraysim[i,"heat",] <- attrdl(data$tempmax_compl,cb,data$atotal,coef=blup[[i]]$blup,
    vcov=blup[[i]]$vcov,type="an",dir="forw",cen=mintempcity[i],
    range=c(mintempcity[i],100),sim=T,nsim=nsim)
  
  arraysim[i,"extreme cold",] <- attrdl(data$tempmax_compl,cb,data$atotal,coef=blup[[i]]$blup,
                                vcov=blup[[i]]$vcov,type="an",dir="forw",cen=mintempcity[i],
                                range=c(-100,perc[1]),sim=T,nsim=nsim)
  arraysim[i,"moderate cold",] <- attrdl(data$tempmax_compl,cb,data$atotal,coef=blup[[i]]$blup,
                                        vcov=blup[[i]]$vcov,type="an",dir="forw",cen=mintempcity[i],
                                        range=c(perc[1],mintempcity[i]),sim=T,nsim=nsim)
  arraysim[i,"moderate heat",] <- attrdl(data$tempmax_compl,cb,data$atotal,coef=blup[[i]]$blup,
                                        vcov=blup[[i]]$vcov,type="an",dir="forw",cen=mintempcity[i],
                                        range=c(mintempcity[i],perc[2]),sim=T,nsim=nsim)
  arraysim[i,"extreme heat",] <- attrdl(data$tempmax_compl,cb,data$atotal,coef=blup[[i]]$blup,
                                        vcov=blup[[i]]$vcov,type="an",dir="forw",cen=mintempcity[i],
                                        range=c(perc[2],100),sim=T,nsim=nsim)
  
  
  
  # STORE THE TOTAL INJURIES (ACCOUNTING FOR MISSING)
  totdeath[i] <- sum(data$atotal,na.rm=T)
}

################################################################################
# ATTRIBUTABLE NUMBERS

# CITY-SPECIFIC
ancity <- matsim
ancitylow <- apply(arraysim,c(1,2),quantile,0.025)
ancityhigh <- apply(arraysim,c(1,2),quantile,0.975)
rownames(ancity) <- rownames(ancitylow) <- rownames(ancityhigh) <- provincies_n

# TOTAL
# NB: FIRST SUM THROUGH CITIES
antot <- colSums(matsim)
antotlow <- apply(apply(arraysim,c(2,3),sum),1,quantile,0.025)
antothigh <- apply(apply(arraysim,c(2,3),sum),1,quantile,0.975)

################################################################################
# TOTAL INJURIES

# BY COUNTRY
totdeathtot <- sum(totdeath)

################################################################################
# ATTRIBUTABLE FRACTIONS

# CITY-SPECIFIC
afcity <- ancity/totdeath*100
afcitylow <- ancitylow/totdeath*100
afcityhigh <- ancityhigh/totdeath*100

# TOTAL
aftot <- antot/totdeathtot*100
aftotlow <- antotlow/totdeathtot*100
aftothigh <- antothigh/totdeathtot*100

#





