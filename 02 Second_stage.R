
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


########################################################################################################
# SECOND-STAGE ANALYSIS: MULTIVARIATE META-ANALYSIS OF THE REDUCED COEF AND THEN COMPUTATION OF BLUP
########################################################################################################

# CREATE AVERAGE TEMPERATURE AND RANGE AS META-PREDICTORS
avgtmean <- sapply(dlist,function(x) mean(x$tempmax_compl,na.rm=T))
rangetmean <- sapply(dlist,function(x) diff(range(x$tempmax_compl,na.rm=T)))

################################################################################
# META-ANALYSIS
# We pool the estimated location-specific overall cumulative exposure-response associations 
# (we obtained an overall estimation - at country level)

# OVERALL CUMULATIVE SUMMARY FOR THE MAIN MODEL
mv <- mvmeta(coef~1,vcov,data=provincies_n,control=list(showiter=T))
summary(mv)


################################################################################
# OBTAIN BLUPS (best linear unbiased predictions)
# For random effects models, predictions are the sum of the average results 
# of the fixed part of the model plus the predicted random effects variances

blup <- blup(mv, vcov=T)


################################################################################
# RE-CENTERING

# GENERATE THE MATRIX FOR STORING THE RESULTS
minperccity <- mintempcity <- rep(NA,length(dlist))
names(mintempcity) <- names(minperccity) <- provincies_n

predper <- c(seq(0,1,0.1),2:98,seq(99,100,0.1))
tmeancountry <- rowMeans(sapply(dlist,function(x) quantile(jitter(x$tempmax_compl),
                                                           predper/100,na.rm=T)))

# DEFINE MINIMUM INJURIES VALUES: EXCLUDE LOW AND VERY HOT TEMPERATURE
for(i in seq(length(dlist))) {
  data <- dlist[[i]]
  predvar <- quantile(data$tempmax_compl,1:99/100,na.rm=T)
  # REDEFINE THE FUNCTION USING ALL THE ARGUMENTS (BOUNDARY KNOTS INCLUDED)
  
  # We define 1 knot at 50th percentile
  argvar <- list(x=predvar,fun="bs",
                 degree=2, knots=quantile(data$tempmax_compl, varper/100, na.rm=T),
                 Bound=range(data$tempmax_compl,na.rm=T))
  
  bvar <- do.call(onebasis,argvar)
  minperccity[i] <- (1:99)[which.min((bvar%*%blup[[i]]$blup))]
  mintempcity[i] <- quantile(data$tempmax_compl,minperccity[i]/100,na.rm=T)
}

# COUNTRY-SPECIFIC POINTS OF MINIMUM INJURIES
(minperccountry <- median(minperccity))

#


################################################################################
# PREDICT THE POOLED OVERALL CUMULATIVE ASSOCIATIONS

################################################################################

prov <- c("Alava", "Albacete", "Alicante", "Almeria", "Avila", "Badajoz", "Illes Balears", "Barcelona", "Burgos",
          "Caceres", "Cadiz", "Castellon", "Ciudad Real", "Cordoba", "A Coruna", "Cuenca", "Girona", "Granada",
          "Guadalajara", "Guipuzkoa", "Huelva", "Huesca", "Jaen", "Leon", "Lleida", "La Rioja", "Lugo", "Madrid",
          "Malaga", "Murcia", "Navarra", "Ourense", "Asturias", "Palencia", "Las Palmas", "Pontevedra", "Salamanca",
          "Santa Cruz de Tenerife", "Cantabria", "Segovia", "Sevilla", "Soria", "Tarragona", "Teruel", "Toledo",
          "Valencia", "Valladolid", "Vizcaia", "Zamora", "Zaragoza")

datanew <- data.frame(avgtmean=mean(tapply(avgtmean,prov,mean)),
                      rangetmean=mean(tapply(rangetmean,prov,mean)))

# PREDICT THE POOLED COEFFICIENTS FOR EACH MODEL
mvpred <- predict(mv,datanew,vcov=T,format="list")


# DEFINE PERCENTILES AND RELATED AVERAGE TEMPERATURES
# ADD A 'JITTER' TO MAKE PERCENTILES UNIQUE (WITH SEED)
predper <- c(seq(0,1,0.1),2:98,seq(99,100,0.1))
set.seed(13041975)
tmeancountry <- rowMeans(sapply(dlist,function(x) quantile(jitter(x$tempmax_compl),
                                                           predper/100,na.rm=T)))

# DEFINE INDICATOR FOR CENTERING PERCENTILE FROM AVERAGE ASSOCIATION
bvar <- crossbasis(tmeancountry,argvar = list(fun="bs",degree=2, knots=tmeancountry[paste(varper,".0%",sep="")]), 
  arglag = list(fun="integer"))
cenindcountry <- which.min(bvar%*%mvpred$fit)

# DEFINE CENTERING PERCENTILE FOR COUNTRY
cenpercountry <- pmin(pmax(predper[cenindcountry],10),90)


# OBTAIN THE CENTERED PREDICTIONS
cen <- tmeancountry[paste(cenpercountry,".0%",sep="")]

# PREDICT THE POOLED OVERALL CUMULATIVE ASSOCIATIONSN FOR MAIN MODEL
cp <- crosspred(bvar,coef=mvpred$fit,vcov=mvpred$vcov,model.link="log",
                 at=tmeancountry,cen=cen)






