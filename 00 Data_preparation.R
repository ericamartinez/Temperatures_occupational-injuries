
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
################################################################################

################################################################################
# PREPARE THE DATA
################################################################################

load("data.Rdata")

provincies_n <- list("Alava", "Albacete", "Alicante", "Almeria", "Avila", 
  "Badajoz", "Illes Balears", "Barcelona", "Burgos", "Caceres", "Cadiz", 
  "Castellon", "Ciudad Real", "Cordoba", "A Coruna", "Cuenca", "Girona", 
  "Granada", "Guadalajara", "Guipuzcoa", "Huelva", "Huesca", "Jaen", "Leon", 
  "Lleida", "La Rioja", "Lugo", "Madrid", "Malaga", "Murcia", "Navarra", 
  "Ourense", "Asturias", "Palencia", "Las Palmas", "Pontevedra", "Salamanca",
  "Santa Cruz de Tenerife", "Cantabria", "Segovia", "Sevilla", "Soria", 
  "Tarragona", "Teruel", "Toledo", "Valencia", "Valladolid", "Vizcaya", 
  "Zamora", "Zaragoza")


# ARRANGE THE DATA AS A LIST OF DATA SETS
provinces <- as.character(unique(tempLATS$centerprovince)) 
dlist <- lapply(provinces,function(x) tempLATS[tempLATS$centerprovince==x,]) 
# Create a list with 50 provinces 
names(dlist) <- provincies_n

#### PARAMETERS FOR THE MAIN MODEL
# 1 internal knot placed at the 50th percentile of location-specific temperature distribution
varper <- c(50)
lag <- 4
arglag <- list(fun="integer")


# LOAD THE PACKAGES
library(dlnm) ; library(mvmeta) ; library(splines) ; library(tsModel)



