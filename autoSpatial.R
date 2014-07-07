# H. Achicanoy
# CIAT, 2013

# ====================================================================== #
# autoSpatial
# ====================================================================== #

### Packages
require(geoR)
require(ape)
require(spdep)
require(geoR)
require(gstat)
require(MASS)
require(ade4)

### Indice de Moran

# Training dataset
trainDist <- as.matrix(dist(x=cbind(trainAll$lon,trainAll$lat),method='euclidean'))
trainDist_inv <- 1/trainDist
diag(trainDist_inv) <- 0

ade4_autocorr <- gearymoran(bilis=trainDist_inv,X=trainAll[,3:19],nrepet=1000,alter="two-sided")

# Testing dataset
testnAllFilter <- testnAll[!duplicated(testnAll[,1:2]),]
testDist <- as.matrix(dist(x=cbind(testnAllFilter$lon,testnAllFilter$lat),method='euclidean'))
testDist_inv <- 1/testDist
diag(testDist_inv) <- 0

ade4_autocorrTest <- gearymoran(bilis=testDist_inv,X=testnAllFilter[,3:19],nrepet=1000,alter="two-sided")

### Variogramas

testnAllFilter
dists <- dist(testnAllFilter[,1:2])
summary(as.vector(dists))
breaks = seq(0, 40, l = 11)

varList <- names(testnAllFilter)[3:17]

smv <- list()

for(var in varList){
  v1 <- variog(coords = testnAllFilter[,1:2], data = testnAllFilter[,var], breaks = breaks)
  v1.summary <- cbind(c(1:10), v1$v, v1$n)
  colnames(v1.summary) <- c("lag", "semi-variance", "pairs")
  smv[[var]] <- v1.summary
  # Sys.sleep(1) # time in seconds
  
  par(mfrow=c(9,9))
  
  plot(v1, type = "b", main = paste(var,sep=""))
  rm(v1); rm(v1.summary)
}

### Modelacion con autocorrelacion espacial

boston <- readOGR(dsn="C:/Users/haachicanoy/Downloads/examples",layer="boston")





##########################################################################
##########################################################################
# La función más importante
install.packages("biomod2")
require(biomod2)
fix(BIOMOD_FormatingData) # La función interna es BIOMOD.formated.data.PA
getAnywhere(BIOMOD.formated.data.PA)[2]
getAnywhere(.pseudo.absences.sampling)
getAnywhere(sre.pseudo.abs.selection)
getMethod(sre.pseudo.abs.selection)
getS3method(sre.pseudo.abs.selection)
getMethod("BIOMOD.formated.data.PA", "sre.pseudo.abs.selection")
##########################################################################
##########################################################################

attach(getNamespace("biomod2"))
getAnywhere(.BIOMOD_FormatingData.check.args)



################################## ==== ##################################
### MODELO LINEAL GENERALIZADO MANEJANDO AUTOCORRELACIÓN ESPACIAL
################################## ==== ##################################

#### Ajuste del modelo lineal generalizado corrigiendo la autocorrelación espacial

glm_0 <- glm(Helianthus_tuberosus~bio_1+bio_2+bio_3+bio_4+bio_5+bio_6+bio_7+bio_8+bio_9+bio_10+bio_11+bio_12+bio_13+bio_14+bio_15+bio_16+bio_17+bio_18+bio_19, data=tr_data, family=binomial(link = "logit"))

### Paso 1: excluir del modelo las variables explicativas linealmente asociadas

vifstep(tr_data[tr_data$Helianthus_tuberosus==1,4:ncol(tr_data)])
#11 variables from the 19 input variables have collinearity problem:

#  bio_5 bio_11 bio_1 bio_6 bio_7 bio_16 bio_17 bio_3 bio_12 bio_14 bio_19

#After excluding the collinear variables, the linear correlation coefficients ranges between:
#  min correlation ( bio_4 ~ bio_2 ):  -0.000540148
#max correlation ( bio_10 ~ bio_8 ):  0.6195154

#---------- VIFs of the remained variables --------
#  Variables      VIF
#1     bio_2 1.479840
#2     bio_4 3.388081
#3     bio_8 3.136583
#4     bio_9 3.698594
#5    bio_10 3.359483
#6    bio_13 2.733437
#7    bio_15 2.430054
#8    bio_18 3.067764

### Paso 2: Ajustar el modelo con las variables explicativas restantes

glm_1 <- glm(Helianthus_tuberosus~bio_2+bio_4+bio_8+bio_9+bio_10+bio_13+bio_15+bio_18, data=tr_data, family=binomial(link = "logit"))
glm_2 <- glm(Helianthus_tuberosus~auto+bio_2+bio_4+bio_8+bio_9+bio_10+bio_13+bio_15+bio_18, data=tr_data, family=binomial(link = "logit"))

### Paso 3: aplicar un procedimiento de selección de variables utilizando
### el criterio de información de Akaike

step_glm1 <- step(glm_1)
step_glm2 <- step(glm_2)

### Paso 4: Extraer los vectores propios a incluir como covariables en el modelo
### a partir del procedimiento: Moran eigenvector GLM filtering

require(spdep)
coords = tr_data[,2:3] # Extraer coordenadas
IDs <- row.names(as(tr_data, "data.frame"))
nb <- knn2nb(knearneigh(as(coords, "matrix"), k=10), row.names=IDs) # Definir el vecindario
nb.list <- nb2listw(nb)
set.seed(1234); me1 <- ME(glm_1$formula,data=tr_data,family="binomial",listw=nb.list)
set.seed(1234); me2 <- ME(glm_2$formula,data=tr_data,family="binomial",listw=nb.list)

# Extraer los vectores propios seleccionados
fitted(me1)
fitted(me2)

# Volver a ajustar el modelo, incluyendo los vectores propios seleccionados por
# Moran eigenvector GLM filtering

step_glm1
glm1_EM_all <- glm(Helianthus_tuberosus~bio_2+bio_4+bio_9+bio_10+bio_15+eigen_vec1,data=tr_data,family=binomial(link="logit"))

require(ape)

dists <- as.matrix(dist(cbind(tr_data$lon, tr_data$lat))) # Matriz de distancias
inv_dists <- 1/dists # Matriz de distancias invertidas
diag(inv_dists) <- 0
inv_dists[1:5, 1:5]

Moran.I(glm1_EM_all$residuals, inv_dists)

require(ape); require(MASS)

# Corregir la autocorrelación espacial del modelo
cSAC <- function(model,coords){
  
  require(ape)
  require(MASS)
  require(spdep)
  
  cat("Check residuals of model \n")
  
  dists <- as.matrix(dist(cbind(coords$lon,coords$lat)))
  inv.dists <- 1/dists; diag(inv.dists) <- 0; rm(dists)
  p.val <- Moran.I(model$residuals,inv.dists)
  
  if(p.val <= 0.05){ # Autocorrelación espacial
    
    cat("> Spatially autocorrelated residuals in GLM fit \n
        > Eigenvector filtering applied \n")
    
    IDs <- row.names(as(coords,"data.frame"))
    nb <- knn2nb(knearneigh(as(coords,"matrix"), k=10), row.names=IDs)
    nb.list <- nb2listw(nb)
    set.seed(1234)
    ME.filter <- ME(model$formula,data=model$data,family="binomial",listw=nb.list)
    new.model <- glm(formula=model$formula+ME.filter$vectors,
                     family=binomial(link="logit"),
                     data=model$data)
    model <- stepAIC(new.model); rm(new.model)
    
  } else {model = model}
  
  return(model)
  
}


cSAC()

pval <- Moran.I(glm_adj$residuals, inv_dists)$p.value
if(pval <= 0.05){
  cat("> Spatially autocorrelated residuals in GLM fit \n
      > Eigenvector filtering applied \n")
  coords = tr_data[,2:3] # Extraer coordenadas
  IDs <- row.names(as(tr_data, "data.frame"))
  nb <- knn2nb(knearneigh(as(coords, "matrix"), k=10), row.names=IDs) # Definir el vecindario
  nb.list <- nb2listw(nb)
  set.seed(1234); me <- ME(glm_adj$formula,data=tr_data,family="binomial",listw=nb.list)
  glm_adjM <- glm(Helianthus_tuberosus~bio_2+bio_4+bio_9+bio_10+bio_15+me$vectors,
                  data=tr_data,family=binomial(link="logit"))
  glm_adj <- stepAIC(glm_adjM,trace=0); rm(glm_adjM)
} else {
  glm_adj <- glm_adj
}

################################## ==== ##################################
### KRIGING
################################## ==== ##################################

### Kriging ordinario

## load some libraries first:
library(gstat)
## load data
d <- read.csv('F:/elev.csv')

## gstat does not like missing data, subset original data:
e <- na.omit(d)

## convert simple data frame into a spatial data frame object:
coordinates(e) <- ~ x+y

## test result with simple bubble plot:
bubble(e, zcol='elev', fill=FALSE, do.sqrt=FALSE, maxsize=2)

## create a grid onto which we will interpolate:
## first get the range in data
x.range <- as.integer(range(e@coords[,1]))
y.range <- as.integer(range(e@coords[,2]))

## now expand to a grid with 500 meter spacing:
grd <- expand.grid(x=seq(from=x.range[1], to=x.range[2], by=500), y=seq(from=y.range[1], to=y.range[2], by=500) )

## convert to SpatialPixel class
coordinates(grd) <- ~ x+y
gridded(grd) <- TRUE

## test it out:
plot(grd, cex=0.5)
points(e, pch=1, col='red', cex=0.7)
title("Interpolation Grid and Sample Points")

## make gstat object:
g <- gstat(id="elev", formula=elev ~ 1, data=e)

## the original data had a large north-south trend, check with a variogram map
plot(variogram(g, map=TRUE, cutoff=4000, width=200), threshold=10)

## another approach:
# create directional variograms at 0, 45, 90, 135 degrees from north (y-axis)
v <- variogram(g, alpha=c(0,45,90,135))

## 0 and 45 deg. look good. lets fit a linear variogram model:
## an un-bounded variogram suggests additional source of anisotropy... oh well.
v.fit <- fit.variogram(v, model=vgm(model='Lin' , anis=c(0, 0.5)))

## plot results:
plot(v, model=v.fit, as.table=TRUE)

## update the gstat object:
g <- gstat(g, id="elev", model=v.fit )

## perform ordinary kriging prediction:
p <- predict(g, model=v.fit, newdata=grd)

## visualize it:

## base graphics
par(mar=c(2,2,2,2))
image(p, col=terrain.colors(20))
contour(p, add=TRUE, drawlabels=FALSE, col='brown')
points(e, pch=4, cex=0.5)
title('OK Prediction')

## lattice graphics: thanks for R. Bivand's advice on this
## 
## alternatively plot quantiles with
## ... col.regions=terrain.colors(6), cuts=quantile(p$elev.pred) ...
##
pts <- list("sp.points", e, pch = 4, col = "black", cex=0.5)
spplot(p, zcol="elev.pred", col.regions=terrain.colors(20), cuts=19, sp.layout=list(pts), contour=TRUE, labels=FALSE, pretty=TRUE, col='brown', main='OK Prediction')

## plot the kriging variance as well
spplot(p, zcol='elev.var', col.regions=heat.colors(100), cuts=99, main='OK Variance',sp.layout=list(pts) )

######################################## / $

require(gstat)

e <- data_adj
e$prob <- glm_0$fitted.values
coordinates(e) <- ~ lon+lat
e <- e[,-c(1:20)]

bubble(e, zcol='prob', fill=FALSE, do.sqrt=FALSE, maxsize=2)

plot(myResp.ras)
points(e, pch=1, col='red', cex=0.7)

g <- gstat(id="prob", formula=prob ~ 1, data=e)
plot(variogram(g, map=TRUE, cutoff=4000, width=200), threshold=10)

v <- variogram(g, alpha=c(0,45,90,135))
v.fit <- fit.variogram(v, model=vgm(model='Lin' , anis=c(0, 0.5)))

plot(v, model=v.fit, as.table=TRUE)

g <- gstat(g, id="prob", model=v.fit )

p <- predict(g, model=v.fit, newdata=myResp.ras)

par(mar=c(2,2,2,2))
image(p, col=terrain.colors(20))
contour(p, add=TRUE, drawlabels=FALSE, col='brown')
points(e, pch=4, cex=0.5)
title('OK Prediction')







