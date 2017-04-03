
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
# FIRST-STAGE ANALYSIS: RUN THE MODEL IN EACH PROVINCE, REDUCE AND SAVE
################################################################################

################################################################################
# CREATE THE OBJECTS TO STORE THE RESULTS

# COEFFICIENTS AND COVARIANCE FOR OVERALL CUMULATIVE SUMMARY
coef <- matrix(NA,length(provinces), 3,
               dimnames= list(provincies_n))
vcov <- vector("list",length(provinces))
names(vcov) <- provincies_n


################################################################################
# RUN THE LOOP

# LOOP
time <- proc.time()[3]
for(i in seq(length(dlist))) {
  
  # PRINT
  cat(i,"")
  
  # EXTRACT THE DATA
  data <- dlist[[i]]
  
  # DEFINE THE CROSSBASIS
  
  # We define 1 internal knot placed at the 50th percentile
  cb <- crossbasis(data$tempmax_compl, lag=lag, argvar=list(fun="bs", degree=2, 
    knots=quantile(data$tempmax_compl, varper/100, na.rm=TRUE)), 
    arglag=arglag)
  
  #summary(cb)
  
  # RUN THE MODEL AND OBTAIN PREDICTIONS
  # We include month (mm) (in order to capture variability due to changes in active population every month),
  # day of the week, holidays, and the number of workers afiliatted to the SS System
  
  model <- glm(atotal ~ cb + dow + hday + phday + mm + logworkers 
    + ns(date, df=8*20), data, family=quasipoisson, na.action="na.exclude")
  
  # Predicted effects: extract from the model those parameters corresponding to cross-basis variables through functions coef and vcov
  pred <- crosspred(cb,model)
  
  # PLOT 
  plot(pred, "slices", lag=0, ylim=c(0.85,1.10), lwd=1, col=2, 
      xlab="Tempertaure", ylab="RR", main= paste("Exposure-response
       relationship for different lags. Province: ", provincies_n[i]))
  lines(pred, "slices", lag=1, col=3)
  lines(pred, "slices", lag=3, col=4)
  legend("topright", c("Lag 0", "Lags 1-2", "Lags 2-4"), col=2:4, lwd=1)
  
  plot(pred, "3d", xlab="Tempertaure", ylab="Lag", zlab="RR", 
    main= paste("Exposure-response relationship for different lags. Province: ",
      provincies_n[i]))
  
  #   
  # REDUCTION TO OVERALL CUMULATIVE
  # Sum the effects of all lags in order to eliminate one dimension of the association. 
  #Sum (acumulate) the risk during the lag period
  red <- crossreduce(cb, model)
  
  # OVERALL CUMULATIVE SUMMARY FOR THE MAIN MODEL
  coef[i,] <- coef(red)
  vcov[[i]] <- vcov(red)
  
 
}
proc.time()[3]-time

#


