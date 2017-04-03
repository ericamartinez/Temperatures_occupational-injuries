
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
# LAG-RESPONSE ASSOCIATIONS 
################################################################################

################################################################################
# CREATE THE OBJECTS TO STORE THE RESULTS
# NOW LAG-RESPONSE ASSOCIATION AT THE 1st AND 99th PERCENTILE VS THE MOIP

lag2 <- 7

# OVERALL CUMULATIVE
coeflag2 <- coeflag3 <- matrix(NA,length(provinces),8,
  dimnames= list(provincies_n))
vcovlag2 <- vcovlag3 <- vector("list",length(provinces))
names(vcovlag2) <- names(vcovlag3) <- provincies_n

################################################################################
# RUN THE MODEL FOR EACH CITY (WITH RE-CENTERED VALUES)

# LOOP
for(i in seq(length(dlist))) {
  
  # PRINT
  cat(i,"")
  
  # EXTRACT THE DATA
  data <- dlist[[i]]
  
  # DEFINE THE CROSS-BASIS (NOW CENTERED AT MINIMUM INJURIES TEMPERATURE)
  cb <- crossbasis(data$tempmax_compl, lag=lag2, 
                    argvar=list(fun="bs", degree=2, knots=quantile(data$tempmax_compl, varper/100, na.rm=TRUE)),
                    arglag=arglag)
  
  # RUN THE MODELS 
  model <- glm(atotal ~ cb + dow + hday + phday + mm + logworkers + ns(date, df=8*20), data, family=quasipoisson, na.action="na.exclude")
   
  # PREDICTIONS AND REDUCTION TO LAG-RESPONSE AT 99TH
  # NB: CENTERING NEEDED HERE, AS IT CHANGES COEF-VCOV
  cen <- quantile(data$tempmax_compl,cenpercountry/100,na.rm=T)
 
  redlag3 <- crossreduce(cb,model,"var",value=quantile(data$tempmax_compl,0.99,na.rm=T),
    cen=cen)
  coeflag3[i,] <- coef(redlag3)
  vcovlag3[[i]] <- vcov(redlag3)
  
  redlag2 <- crossreduce(cb,model,"var",value=quantile(data$tempmax_compl,0.01,na.rm=T),
    cen=cen)
  coeflag2[i,] <- coef(redlag2)
  vcovlag2[[i]] <- vcov(redlag2)
  
  
}

################################################################################
# META-ANALYSIS

# RUN THE LAG MODELS FOR COLD AND HEAT
mvlagcold <- mvmeta(coeflag2~1,vcovlag2,data=provincies_n,
  control=list(showiter=T))
summary(mvlagcold)

mvlagheat <- mvmeta(coeflag3~1,vcovlag3,data=provincies_n,
  control=list(showiter=T))
summary(mvlagheat)



################################################################################
# PREDICT THE POOLED COEFFICIENTS

mvpredlagcold <- predict(mvlagcold,datanew,vcov=T,format="list")
mvpredlagheat <- predict(mvlagheat,datanew,vcov=T,format="list")


################################################################################
# PREDICT THE POOLED LAG-RESPONSE ASSOCIATIONS

# OBTAIN THE PREDICTIONS
blag <- do.call(onebasis,c(list(x=seq(0,lag2)),attr(cb,"arglag")))

cplagheat <- crosspred(blag,coef=mvpredlagheat$fit,vcov=mvpredlagheat$vcov,
                    model.link="log",at=0:7)
cplagcold <- crosspred(blag,coef=mvpredlagcold$fit,vcov=mvpredlagcold$vcov,
                    model.link="log",at=0:7)



