# H. Achicanoy
# CIAT, 2013

# ====================================================================== #
# Helianthus tuberosus
# ====================================================================== #

### packages
require(raster); require(rgdal); require(biomod2); require(dismo)
require(BIOMOD); require(usdm); require(stringr)

### working directory
# setwd("/home/haachicanoy") # Home curie
setwd("/mnt/workspace_cluster_6/haachicanoy") # DAPADFS

# ======== native area
natAr.dir <- "/curie_data2/ncastaneda/cluster_variables/pruebas_aposteriori_20052013/nareas"
myResp.ras <- list.files(natAr.dir, pattern=".asc", full.names=T)
myResp.ras <- raster(myResp.ras[[3]]) # [[3]] para Helianthus tuberosus
rm(list=c("natAr.dir"))

# ======== cut native area
ext <- extent(-135,-50,20,65)
myResp.ras <- crop(myResp.ras,ext)
rm(list=c("ext"))
## Anotaciones: recortar área nativa para ubicaciones donde se cuenta
## con información

# ======== read response variable
myRespName <- "Helianthus_tuberosus"

oc.dir   <- "/curie_data2/ncastaneda/cluster_variables/pruebas_aposteriori_20052013/occurrences"
myRespXY <- read.csv(paste(oc.dir,"/",myRespName,".csv",sep=""),header=T,sep=",")
myRespXY <- myRespXY[,1:3]; rm(oc.dir)
myRespXY <- subset(myRespXY[,2:3],myRespXY$Taxon=="Helianthus_tuberosus")
dups     <- duplicated(myRespXY[])
myRespXY <- myRespXY[!dups,]; rm(dups)
myRespXY <- myRespXY[,1:2] # 178 registros

myResp   <- extract(x=myResp.ras, y=myRespXY)
myRespXY <- cbind(myRespXY, myResp)
myRespXY <- as.data.frame(myRespXY); rm(myResp)

# coordenadas de puntos dentro del área nativa
myRespXY <- myRespXY[complete.cases(myRespXY),] # 157 registros
colnames(myRespXY) <- c("lon","lat",myRespName)
rownames(myRespXY) <- 1:nrow(myRespXY)

# myRespXY[,1]: longitude
# myRespXY[,2]: latitude
# myRespXY[,3]: occurrences

# ======== read explanatory variables
envDir <- "/curie_data2/ncastaneda/cluster_variables/pruebas_aposteriori_20052013/clim/Helianthus_tuberosus-clim/"
myExpl <- paste(envDir,list.files(envDir),sep="")
myExpl <- lapply(myExpl,raster)
ext <- extent(-135,-50,20,65)
myExpl <- lapply(myExpl,FUN=function(x){crop(x,ext)})
rm(list=c("ext","envDir"))
myExpl <- raster::stack(myExpl)

# ======== BIOMOD formating data
# Pseudo-ausencias: "sre" method or "random" method
# 10.000 valores
# 10 conjuntos

set.seed(1234)
system.time({ # 42 seconds
myBiomodData <- BIOMOD_FormatingData(resp.var       = myRespXY[,3],
                                     expl.var       = myExpl,
                                     resp.xy        = myRespXY[,1:2],
                                     resp.name      = myRespName,
                                     PA.strategy    = 'sre',
                                     PA.nb.rep      = 10,
                                     PA.nb.absences = 10000,
                                     PA.sre.quant   = 0.05
                                     )
})

myBiomodData@PA # Pseudo-ausencias seleccionadas
coords <- myBiomodData@coord # Coordenadas
rownames(coords) <- 1:nrow(coords)
PA_data <- cbind(coords,myBiomodData@PA); rm(coords)

# Para las filas 1:nrow(myRespXY), seleccionar 70% train / 30% test
presences <- PA_data[1:nrow(myRespXY),]
presences <- presences[,-c(1:2)]
pr_sel <- matrix(0,ncol=10,nrow=nrow(presences),byrow=FALSE)
for(i in 1:10){
  set.seed(1234)
  sam <- sample(1:nrow(presences),size=round((nrow(presences)*0.70),0))
  pr_sel[as.numeric(rownames(presences[sam,])),i] <- 1
}
rm(list=c("i","sam")); pr_sel <- as.data.frame(pr_sel)
colnames(pr_sel) <- paste('IDtrain_PA',1:10,sep="")

# Para las filas restantes, que sean TRUE, seleccionar 70% train / 30% test
pabsences <- PA_data[-c(1:nrow(myRespXY)),]
pabsences <- pabsences[,-c(1:2)]
pabsences[pabsences==FALSE] <- NA
pa_sel <- pabsences
rownames(pabsences) <- 1:nrow(pabsences)
rownames(pa_sel) <- 1:nrow(pa_sel)
for(i in 1:10){
  rows <- which(pabsences[,i]==TRUE)
  set.seed(1234)
  sam <- sample(rows,size=round((length(rows)*0.70),0))
  pa_sel[sam,i] <- 1
  pa_sel[setdiff(rows,sam),i] <- 0
}
rm(list=c("i","sam","rows"))
colnames(pa_sel) <- paste('IDtrain_PA',1:10,sep="")

# Unir los puntos identificados como training
dt_sel <- rbind(pr_sel,pa_sel)
rm(list=c("pr_sel","pa_sel","presences","pabsences"))

PA_data <- cbind(PA_data,dt_sel)
rm(dt_sel)

PA_data[,"occ"] <- 0
PA_data[1:nrow(myRespXY),"occ"] <- 1

### Explorar funcionamiento interno de funciones en R

# fix() :::::::::::::::: Primera forma de ingresar a las funciones
# getAnywhere() :::::::: Segunda forma de ingresar a las funciones internas
# findMethods() :::::::: Tercer forma de ingresar a las funciones más complicadas

fix(BIOMOD_FormatingData)
getAnywhere()
findMethods(sre.pseudo.abs.selection)
getAnywhere(sre)

# ======== BIOMOD modeling options
# solo una variable explicativa
#system.time({ # 3 seconds
#  myBiomodOption <- BIOMOD_ModelingOptions(GLM = list(myFormula=Helianthus_tuberosus~bio_1, family='binomial'))
#})

# ======== BIOMOD running options
# este paso es rapido
#system.time({ # 72 seconds
#  myBiomodModelOut <- BIOMOD_Modeling(myBiomodData,
#                                      models            = c('GLM'),
#                                      models.options    = myBiomodOption,
#                                      NbRunEval         = 1,
#                                      DataSplit         = 70,
#                                      Prevalence        = NULL,
#                                      VarImport         = 0,
#                                      models.eval.meth  = c('TSS'),
#                                      SaveObj           = TRUE,
#                                      rescal.all.models = TRUE,
#                                      do.full.models    = FALSE,
#                                      modeling.id       = paste(myRespName,"_autocor",sep="")
#                                      )
#})
#rm(list=c("myBiomodData","myBiomodOption"))

# ======== Extract information from the models
# (Cargar si se ha salido del espacio de trabajo)
#load("/mnt/workspace_cluster_6/haachicanoy/Helianthus.tuberosus/Helianthus.tuberosus.Helianthus_tuberosus_PAsre.out")
#myBiomodModelOut <- Helianthus.tuberosus.Helianthus_tuberosus_PAsre.out
#myBiomodModelOut

#load("/home/haachicanoy/Helianthus.tuberosus/Helianthus.tuberosus.Helianthus_tuberosus_final2.models.out")
#myBiomodModelOut <- Helianthus.tuberosus.Helianthus_tuberosus_final2.models.out
#myBiomodModelOut

# ======== read training data
load(myBiomodModelOut@calib.lines@link) # ID puntos de entrenamiento
load(myBiomodModelOut@formated.input.data@link) # Datos e información

#nData <- data.frame(data@coord)
#names(nData) <- c("lon","lat")
#nData['occ'] <- 0
#nData[(1:sum(!is.na(data@data.species))),3] <- 1
#nData[150:170,]

#calib.lines <- as.data.frame(calib.lines)
#for(i in 1:ncol(calib.lines)){
#  nData[paste('IDtrain_PA',i,sep="")] <- as.numeric(calib.lines[,i])
#}

#rm(list=c("i","calib.lines"))

# ======== cantidad de datos
# (presencias + pseudo-ausencias utilizadas para el ajuste del modelo)
# debe ser de 10000 + número de presencias
# por cada conjunto de pseudo ausencias

sum_col <- numeric(10)
for(i in 1:10){
  sum_col[i] <- sum(!is.na(PA_data[paste('IDtrain_PA',i,sep="")]))
}
sum_col

rm(list=c("i","sum_col"))

# ======== Spatial sorting bias function

ssbFun <- function(data, count){
  
  require(biomod2); require(dismo); require(BIOMOD)
  require(usdm); require(stringr)
  
  seedList <- c(1234,1235,1236,1237,1238,1239,1240,1241,1242,1243)
  dataSets <- list()
  summarySel <- list() # lista para guardar: SSB, semillas para escoger el conjunto de PA
  ssbPA <- numeric(10)
  ssbPA_c <- numeric(10)
  for(i in 1:10){
    ### Partición de datos
    # Presence training
    p_train <- data.frame(data$lon[data$occ==1 & data[paste('IDtrain_PA',i,sep="")]==1],
                          data$lat[data$occ==1 & data[paste('IDtrain_PA',i,sep="")]==1])
    names(p_train) <- c("lon","lat")
    p_train <- p_train[!is.na(p_train$lon),]
    p_train <- p_train[!duplicated(p_train),]
    
    # Presence testing
    p_test <- data.frame(data$lon[data$occ==1 & data[paste('IDtrain_PA',i,sep="")]==0],
                         data$lat[data$occ==1 & data[paste('IDtrain_PA',i,sep="")]==0])
    names(p_test) <- c("lon","lat")
    p_test <- p_test[!is.na(p_test$lon),]
    p_test <- p_test[!duplicated(p_test),]
    
    # Pseudo-absences training
    a_train <- data.frame(data$lon[data$occ==0 & data[paste('IDtrain_PA',i,sep="")]==1],
                          data$lat[data$occ==0 & data[paste('IDtrain_PA',i,sep="")]==1])
    names(a_train) <- c("lon","lat")
    a_train <- a_train[!is.na(a_train$lon),]
    a_train <- a_train[!duplicated(a_train),]
    
    # Pseudo-absences testing
    a_test <- data.frame(data$lon[data$occ==0 & data[paste('IDtrain_PA',i,sep="")]==0],
                         data$lat[data$occ==0 & data[paste('IDtrain_PA',i,sep="")]==0])
    names(a_test) <- c("lon","lat")
    a_test <- a_test[!is.na(a_test$lon),]
    a_test <- a_test[!duplicated(a_test),]
    
    ### Sesgo de selección espacial original
    SSB <- ssb(p_test, a_test, p_train)
    ssbPA[i] <- SSB[,1]/SSB[,2]
    
    ### Selección de conjuntos de pseudo-ausencias reducidos PAIRWISE DISTANCE (utilizar semillas)
    PA_total <- rbind(a_train, a_test)
    set.seed(seedList[i])
    pointSel <- pwdSample(fixed=p_test, sample=PA_total, reference=p_train,
                          n=((7*count)/(0.30*count))) # 7 veces las presencias
    coordSel <- as.vector(na.exclude(as.vector(pointSel)))
    dataSel <- PA_total[coordSel,]
    rownames(dataSel) <- 1:nrow(dataSel)
    dataSets[[i]] <- dataSel
    
    ### Sesgo de selección espacial corregido
    SSB_c <- ssb(p_test, dataSel[,1:2], p_train)
    ssbPA_c[i] <- SSB_c[,1]/SSB_c[,2]
    
  }
  
  summarySel <- list(SSB=ssbPA,PA_dataSets=dataSets,cSSB=ssbPA_c)
  return(summarySel)
}

# running function
count <- nrow(myRespXY)
system.time({ # 40 seconds
  dt_sel <- ssbFun(PA_data,count)
})
rm(list=c("count"))

### =============== AJUSTAR LOS MODELOS

### Toda la información

PA_data # Info de pseudo ausencias
rm(list=c("myBiomodData","PA_data","ssbFun"))

nData <- cbind(nData,data@data.env.var)
nData <- as.data.frame(nData)

# la idea sería realizar varias repeticiones 5 ?????

# por set de PA seleccionado [esto debe ir en un ciclo y una función]

maxdir <- "/mnt/workspace_cluster_6/haachicanoy/lib"
maxdir <- "/home/haachicanoy"

for(i in 1:length(dsSel[[2]])){

  id_PA <- as.numeric(row.names(dsSel[[2]][[i]]))
  presences <- nData[1:sum(!is.na(data@data.species)),c("occ","lon","lat","bio_1","bio_2",
                                                        "bio_3","bio_4","bio_5","bio_6",
                                                        "bio_7","bio_8","bio_9","bio_10",
                                                        "bio_11","bio_12","bio_13","bio_14",
                                                        "bio_15","bio_16","bio_17","bio_18",
                                                        "bio_19")]
  
  pabsences <- nData[id_PA,c("occ","lon","lat","bio_1","bio_2",
                              "bio_3","bio_4","bio_5","bio_6",
                              "bio_7","bio_8","bio_9","bio_10",
                              "bio_11","bio_12","bio_13","bio_14",
                              "bio_15","bio_16","bio_17","bio_18",
                              "bio_19")]
  
  rm(id_PA)
  
  presences <- presences[complete.cases(presences),]; rownames(presences) <- 1:nrow(presences)
  pabsences <- pabsences[complete.cases(pabsences),]; rownames(pabsences) <- 1:nrow(pabsences)
  pabsences$occ[pabsences$occ==0] = NA
  
  # Training and testing sets 70/30 respectively (fijar la semilla)
  set.seed(1234); sp_sel <- sample(1:nrow(presences),size=round((nrow(presences)*.30),0))
  set.seed(1234); pa_sel <- sample(1:nrow(pabsences),size=round((nrow(pabsences)*.30),0))
  
  spp_tr <- presences[-sp_sel,]; rownames(spp_tr) <- 1:nrow(spp_tr)
  spp_te <- presences[sp_sel,]; rownames(spp_te) <- 1:nrow(spp_te)
  
  pab_tr <- pabsences[-pa_sel,]; rownames(pab_tr) <- 1:nrow(pab_tr)
  pab_te <- pabsences[pa_sel,]; rownames(pab_te) <- 1:nrow(pab_te)
  
  rm(list=c("sp_sel","pa_sel"))
  tr_data <- rbind(spp_tr,pab_tr); rm(list=c("spp_tr","pab_tr","presences"))
  te_data <- rbind(spp_te,pab_te); rm(list=c("spp_te","pab_te","pabsences"))
  
  names(tr_data)[1] <- myRespName
  names(te_data)[1] <- myRespName
  rownames(tr_data) <- 1:nrow(tr_data)
  rownames(te_data) <- 1:nrow(te_data)
  
  # Toca hacer este paso porque al formatear los datos en biomod2 surge un
  # error al incluir en la variable de respuesta valores "NA", que indican
  # las pseudo-ausencias utilizadas para la evaluación del modelo
  tr_data[is.na(tr_data),myRespName] <- 0
  te_data[is.na(te_data),myRespName] <- 0
  #te_data[which(te_data[,1]==0),myRespName] <- NA
  
  system.time({ # 1 second
  myBiomodData <- BIOMOD_FormatingData(resp.var  = tr_data[myRespName],
                                       resp.xy   = tr_data[,2:3],
                                       expl.var  = tr_data[,2:22],
                                       eval.resp.var = te_data[myRespName],
                                       eval.resp.xy = te_data[,2:3],
                                       eval.expl.var = te_data[,2:22],
                                       resp.name = myRespName,
                                       PA.nb.rep = 0,
                                       PA.strategy = "none",
                                       na.rm = FALSE)
  })
  
  system.time({ # 1 second
  myBiomodOption <- BIOMOD_ModelingOptions(GLM = list(myFormula=Helianthus.tuberosus~lon+lat+bio_1+bio_2+bio_3+bio_4+bio_5+bio_6+bio_7+bio_8+bio_9+bio_10+bio_11+bio_12+bio_13+bio_14+bio_15+bio_16+bio_17+bio_18+bio_19, family="binomial"),
                                           ANN = list(NbCV=5, maxit=500),
                                           RF  = list(do.classif=TRUE,ntree=100),
                                           MAXENT = list(path_to_maxent.jar = maxdir,
                                                         maximumiterations = 500,
                                                         visible = TRUE,
                                                         linear = TRUE,
                                                         quadratic = TRUE,
                                                         product = TRUE,
                                                         threshold = TRUE,
                                                         hinge = TRUE,
                                                         lq2lqptthreshold = 80,
                                                         l2lqthreshold = 10,
                                                         hingethreshold = 15,
                                                         beta_threshold = -1,
                                                         beta_categorical = -1,
                                                         beta_lqp = -1,
                                                         beta_hinge = -1,
                                                         defaultprevalence = 0.5))
  })
  
  system.time({ # 10 seconds
  myBiomodModelOut <- BIOMOD_Modeling(myBiomodData,
                                      models            = c('GLM','ANN','RF','MAXENT'), # ,
                                      models.options    = myBiomodOption,
                                      NbRunEval         = 5,
                                      DataSplit         = 100,
                                      Prevalence        = NULL,
                                      models.eval.meth  = c('TSS','KAPPA','ROC'),
                                      SaveObj           = TRUE,
                                      rescal.all.models = TRUE,
                                      do.full.models    = FALSE,
                                      modeling.id       = paste(myRespName,"_final2",sep=""))
  })

}

fix(BIOMOD_Modeling); getAnywhere(.Biomod.Models.loop); getAnywhere(.Biomod.Models)

if (em.by %in% c("PA_dataset", "PA_dataset+algo", "PA_dataset+repet")) {
  obs <- get_formal_data(modeling.output, "resp.var")
  expl <- get_formal_data(modeling.output, "expl.var")
  if (head(unlist(strsplit(assemb, "_")), 1) != "AllData") {
    kept_cells <- get_formal_data(modeling.output)@PA[, 
                                                      paste(head(unlist(strsplit(assemb, "_")), 1))]
    obs <- obs[kept_cells]
    expl <- expl[kept_cells, , drop = F]
  }
}

# Las repeticiones las realiza teniendo en cuenta diferentes selecciones de
# datos de entrenamiento y validación, al especificar el 100% de datos usados
# para entrenar el modelo se emplea la totalidad de la información. Por tanto,
# el conjunto de evaluación previamente especificado solo servira de base para
# ejecutar la evaluación del modelo.

# Cuento con problemas al correr los modelos maxent

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #
# Modelos

path.models <- "/mnt/workspace_cluster_6/haachicanoy/Helianthus.tuberosus/models/Helianthus_tuberosus_final2"
load(paste(path.models,"/Helianthus.tuberosus_AllData_Full_GLM",sep=""))
glm_adj <- get_formal_model(Helianthus.tuberosus_AllData_Full_GLM)

path.models <- "/mnt/workspace_cluster_6/haachicanoy/Helianthus.tuberosus/models/Helianthus_tuberosus_final2"
load(paste(path.models,"/Helianthus.tuberosus_AllData_Full_ANN",sep=""))
ann_adj <- get_formal_model(Helianthus.tuberosus_AllData_Full_ANN)

source('https://gist.github.com/fawda123/6860630/raw/e50fc6ef30b8269660b4e65aeec7ce02beb9b551/lek_fun.r')
lek.fun(ann_adj)

path.models <- "/mnt/workspace_cluster_6/haachicanoy/Helianthus.tuberosus/models/Helianthus_tuberosus_final2"
load(paste(path.models,"/Helianthus.tuberosus_AllData_Full_RF",sep=""))
rf_adj <- get_formal_model(Helianthus.tuberosus_AllData_Full_RF)

plot(rf_adj$votes)


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

# Spatial autocorrelation for residuals training model
tr_dist <- as.matrix(dist(x=cbind(tr_data$lon,tr_data$lat),method='euclidean'))
tr_dist_inv <- 1/tr_dist
diag(tr_dist_inv) <- 0

moranI <- gearymoran(bilis=tr_dist_inv,X=glm_adj$residuals,nrepet=1000,alter="two-sided")
Moran.I(glm_adj$residuals,tr_dist_inv) # Moran.I(as.numeric(glm_adj$fitted.values),tr_dist_inv)
as.numeric(predict(tmodel,spp_tr[,3:ncol(spp_tr)]))
stats::predict(glm_adj,te_data[,c(1,4:ncol(te_data))],type="link")

# Datos de evaluación del modelo
proba <- 1/(1+exp(-as.numeric(stats::predict(glm_adj,te_data[,c(1,4:ncol(te_data))],type="link"))))

te_dist <- as.matrix(dist(x=cbind(te_data$lon,te_data$lat),method='euclidean'))
te_dist_inv <- 1/te_dist
diag(te_dist_inv) <- 0

Moran.I(proba,te_dist_inv)

####
library(rgl)
data(volcano)
z <- proba # Exaggerate the relief
x <- 10 * (1:nrow(z)) # 10 meter spacing (S to N)
y <- 10 * (1:ncol(z)) # 10 meter spacing (E to W)
zlim <- range(z)
zlen <- zlim[2] - zlim[1] + 1
colorlut <- terrain.colors(zlen,alpha=0) # height color lookup table
col <- colorlut[ z-zlim[1]+1 ] # assign colors to heights for each point
open3d()
rgl.surface(x=te_data$lon, y=te_data$lat, z=proba, alpha=0.75, back="lines")
####

# Cargar un modelo en particular
path.models <- "/mnt/workspace_cluster_6/haachicanoy/Helianthus.tuberosus.1/models/Helianthus_tuberosus_final1"
load(paste(path.models,"/Helianthus.tuberosus.1_AllData_RUN1_ANN",sep=""))
load(paste(path.models,"/Helianthus.tuberosus.1_AllData_RUN2_ANN",sep=""))
load(paste(path.models,"/Helianthus.tuberosus.1_AllData_RUN3_ANN",sep=""))
load(paste(path.models,"/Helianthus.tuberosus.1_AllData_RUN4_ANN",sep=""))
load(paste(path.models,"/Helianthus.tuberosus.1_AllData_RUN5_ANN",sep=""))

ann1 <- get_formal_model(Helianthus.tuberosus.1_AllData_RUN1_ANN)
ann2 <- get_formal_model(Helianthus.tuberosus.1_AllData_RUN2_ANN)
ann3 <- get_formal_model(Helianthus.tuberosus.1_AllData_RUN3_ANN)
ann4 <- get_formal_model(Helianthus.tuberosus.1_AllData_RUN4_ANN)
ann5 <- get_formal_model(Helianthus.tuberosus.1_AllData_RUN5_ANN)

glm1 <- get_formal_model(Helianthus.tuberosus.1_AllData_RUN1_GLM)
glm2 <- get_formal_model(Helianthus.tuberosus.1_AllData_RUN2_GLM)
glm3 <- get_formal_model(Helianthus.tuberosus.1_AllData_RUN3_GLM)
glm4 <- get_formal_model(Helianthus.tuberosus.1_AllData_RUN4_GLM)
glm5 <- get_formal_model(Helianthus.tuberosus.1_AllData_RUN5_GLM)

residuals(glm1)
residuals(glm2)
residuals(glm3)
residuals(glm4)
residuals(glm5)

### 

### proyecciones

# Al parecer este código no genera mapas en formato tiff
myBiomodProj <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                  new.env = myExpl,
                                  proj.name = "current",
                                  selected.models = 'all',
                                  binary.meth = "TSS",
                                  compress = "gzip",
                                  clamping.mask = TRUE,
                                  output.format = ".tif")

### Con GLM

# Leer el archivo que contiene el modelo final
fit_file <- "/mnt/workspace_cluster_6/haachicanoy/Helianthus.tuberosus/Helianthus.tuberosus.Helianthus_tuberosus_final2.models.out"
sp_mOut <- get(load(fit_file))

# Cargar un modelo especifico
tmodel <- get(BIOMOD_LoadModels(sp_mOut,models="GLM"))

# determinar en que celdas existen valores
prj_data <- as.data.frame(xyFromCell(myResp.ras,which(!is.na(myResp.ras[]))))
prj_data$CELL <- which(!is.na(myResp.ras[]))

# extraer información ambiental
bio_stk <- myExpl
prj_data <- cbind(prj_data,extract(bio_stk,prj_data[,c("x","y")]))
rm(bio_stk)

# omitir celdas que por alguna razón contienen valores faltantes
prj_data$NAs <- apply(prj_data,1,FUN=function(x) {nac <- length(which(is.na(x))); return(nac)})
prj_data <- prj_data[which(prj_data$NAs == 0),]
prj_data$NAs <- NULL

# proyectar
names(prj_data)[1:2] <- c("lon","lat")
prjVect <- predict(tmodel,prj_data)
prj_rs <- raster(myResp.ras)
prj_rs[prj_data$CELL] <- prjVect

plot(prj_rs)

### Con ANN

fit_file <- "/mnt/workspace_cluster_6/haachicanoy/Helianthus.tuberosus/Helianthus.tuberosus.Helianthus_tuberosus_final2.models.out"
sp_mOut <- get(load(fit_file))

tmodel <- get(BIOMOD_LoadModels(sp_mOut,models="ANN")) # Cargar un modelo especifico

# proyectar
prjVect <- predict(tmodel,prj_data)
prj_rs <- raster(myResp.ras)
prj_rs[prj_data$CELL] <- prjVect

plot(prj_rs)

### Con Random Forest

fit_file <- "/mnt/workspace_cluster_6/haachicanoy/Helianthus.tuberosus/Helianthus.tuberosus.Helianthus_tuberosus_final2.models.out"
sp_mOut <- get(load(fit_file))

tmodel <- get(BIOMOD_LoadModels(sp_mOut,models="RF")) # Cargar un modelo especifico

# proyectar
prjVect <- predict(tmodel,prj_data)
prj_rs <- raster(myResp.ras)
prj_rs[prj_data$CELL] <- prjVect

plot(prj_rs)

prj_rs <- writeRaster(prj_rs,"/home/haachicanoy/sdm_proj.tiff",datatype="FLT4S", format="GTiff", overwrite=TRUE, options="COMPRESS=NONE")

# Ensemble modelling


myBiomodEM <- BIOMOD_EnsembleModeling(
  modeling.output = myBiomodModelOut,
  chosen.models = 'all',
  em.by='all',
  eval.metric = c('ROC'),
  eval.metric.quality.threshold = c(0.7),
  prob.mean = T,
  prob.cv = T,
  prob.ci = T,
  prob.ci.alpha = 0.05,
  prob.median = T,
  committee.averaging = T,
  prob.mean.weight = T,
  prob.mean.weight.decay = 'proportional')

myBiomodEM
get_evaluations(myBiomodEM)


# proyeccion con la mediana

EMmodel <- get(BIOMOD_LoadModels(myBiomodEM,models="EMmedian"))

# determinar en que celdas existen valores
prj_data <- as.data.frame(xyFromCell(myResp.ras,which(!is.na(myResp.ras[]))))
prj_data$CELL <- which(!is.na(myResp.ras[]))

# extraer información ambiental
bio_stk <- myExpl
prj_data <- cbind(prj_data,extract(bio_stk,prj_data[,c("x","y")]))
rm(bio_stk)

# omitir celdas que por alguna razón contienen valores faltantes
prj_data$NAs <- apply(prj_data,1,FUN=function(x) {nac <- length(which(is.na(x))); return(nac)})
prj_data <- prj_data[which(prj_data$NAs == 0),]
prj_data$NAs <- NULL

# proyectar
prjVect <- predict(EMmodel,prj_data)
prj_rs <- raster(myResp.ras)
prj_rs[prj_data$CELL] <- prjVect

png("/mnt/workspace_cluster_6/haachicanoy/EMmedian.png",width=2048,height=2048,res=300,pointsize=13)
plot(prj_rs)
dev.off()

#prj_rs <- writeRaster(prj_rs,"/home/haachicanoy/sdm_proj.tiff",datatype="FLT4S", format="GTiff", overwrite=TRUE, options="COMPRESS=NONE")

# proyeccion con la media

EMmodel <- get(BIOMOD_LoadModels(myBiomodEM,models="EMmean"))

# determinar en que celdas existen valores
prj_data <- as.data.frame(xyFromCell(myResp.ras,which(!is.na(myResp.ras[]))))
prj_data$CELL <- which(!is.na(myResp.ras[]))

# extraer información ambiental
bio_stk <- myExpl
prj_data <- cbind(prj_data,extract(bio_stk,prj_data[,c("x","y")]))
rm(bio_stk)

# omitir celdas que por alguna razón contienen valores faltantes
prj_data$NAs <- apply(prj_data,1,FUN=function(x) {nac <- length(which(is.na(x))); return(nac)})
prj_data <- prj_data[which(prj_data$NAs == 0),]
prj_data$NAs <- NULL

# proyectar
prjVect <- predict(EMmodel,prj_data)
prj_rs <- raster(myResp.ras)
prj_rs[prj_data$CELL] <- prjVect

png("/mnt/workspace_cluster_6/haachicanoy/EMmean.png",width=2048,height=2048,res=300,pointsize=13)
plot(prj_rs)
dev.off()

#prj_rs <- writeRaster(prj_rs,"/home/haachicanoy/sdm_proj.tiff",datatype="FLT4S", format="GTiff", overwrite=TRUE, options="COMPRESS=NONE")

# proyeccion con la committe averaging
# EMca, EMwmean

EMmodel <- get(BIOMOD_LoadModels(myBiomodEM,models="EMca"))

# proyectar
prjVect <- predict(EMmodel,prj_data)
prj_rs <- raster(myResp.ras)
prj_rs[prj_data$CELL] <- prjVect

png("/mnt/workspace_cluster_6/haachicanoy/EMca.png",width=2048,height=2048,res=300,pointsize=13)
plot(prj_rs)
dev.off()

# proyeccion con la weigth mean

EMmodel <- get(BIOMOD_LoadModels(myBiomodEM,models="EMwmean"))

# proyectar
prjVect <- predict(EMmodel,prj_data)
prj_rs <- raster(myResp.ras)
prj_rs[prj_data$CELL] <- prjVect

png("/mnt/workspace_cluster_6/haachicanoy/EMwmean.png",width=2048,height=2048,res=300,pointsize=13)
plot(prj_rs)
dev.off()





# Todos los datos incluyendo presencias y pseudo-ausencias
allData <- cbind(data@coord,data@data.env.var)

# Filas del primer conjunto de PA
rows_1 <- as.numeric(row.names(res[[2]][[1]]))
# Total de PA_1 con información ambiental
data_PA <- allData[rows_1,]
data_PA <- data_PA[complete.cases(data_PA),]
# Total de presencias con información ambiental
data_sp <- allData[1:p_count,]
#data_1 <- rbind(data_env_pres, data_1)
#data_1["occ"] <- NA
#data_1[1:p_count,"occ"] <- 1

# Datos para evaluación del modelo
set.seed(1234); id_pa <- sample(1:nrow(data_PA),size=round((nrow(data_PA)*.30),0))
set.seed(1234); id_sp <- sample(1:nrow(data_sp),size=round((nrow(data_sp)*.30),0))

testn_pa <- data_PA[id_pa,]; row.names(testn_pa) <- 1:nrow(testn_pa)
testn_sp <- data_PA[id_sp,]; row.names(testn_sp) <- 1:nrow(testn_sp)

testn_pa["y"] <- 0
testn_sp["y"] <- 1

train_pa <- data_PA[-id_pa,]; row.names(train_pa) <- 1:nrow(train_pa)
train_sp <- data_sp[-id_sp,]; row.names(train_sp) <- 1:nrow(train_sp)

train_pa["y"] <- 0
train_sp["y"] <- 1

trainAll <- rbind(train_sp,train_pa)
testnAll <- rbind(testn_sp,testn_pa)

library(nnet)
ANNmodel <- nnet(y~bio_1+bio_2+bio_4+bio_5+bio_6+
                   bio_8+bio_9+bio_10+bio_11+bio_12+
                   bio_13+bio_14+bio_15+bio_16+bio_17+
                   bio_18+bio_19,size=10,data=trainAll)

predict(ANNmodel, testnAll, type = "class")

...

data_10 <- algo

save.image()

# ================================================= #
### Medir autocorrelación
# ================================================= #

# Packages
require(ape); require(spdep); require(geoR); require(gstat); require(MASS); require(ade4)

# Bajo aleatoriedad: uniformes
q1 <- matrix(runif(36), nrow = 6)
q2 <- matrix(runif(36), nrow = 6)

mantel.test(q1, q2, graph = TRUE,
            main = "Mantel test: a random example with 6 X 6 matrices",
            xlab = "Mantel-statistic", ylab = "Density",
            sub = "The vertical line shows the observed M-statistic")

# Bajo aleatoriedad: normales
q1 <- matrix(rnorm(36, 10, 1), nrow = 6)
q2 <- matrix(rnorm(36, 30, 3), nrow = 6)
mantel.test(q1, q2, graph = TRUE,
            main = "Mantel test: a random example with 6 X 6 matrices",
            xlab = "Mantel-statistic", ylab = "Density",
            sub = "The vertical line shows the observed M-statistic")

# Bajo no aleatoriedad: normales
grid<-expand.grid(seq(0,9),seq(0,9))
plot(grid)
distancia<-dist(grid,diag=T,upper=T)
distancia<-as.matrix(distancia)
covarianza<-4*exp(-3*distancia/6)
covarianza<-as.matrix(covarianza)
medias.cte<-rep(20,100)
medias.ncte<-2*grid[,1]+3*grid[,2]
normal.cte<-mvrnorm(1,medias.cte,covarianza)
dist_cte=as.matrix(dist(normal.cte,diag=T,upper=T))
normal.ncte<-mvrnorm(1,medias.ncte,covarianza)
dist_ncte=as.matrix(dist(normal.ncte,diag=T,upper=T))
mantel.test(distancia, dist_cte, graph = TRUE,
            main = "Mantel test: a simulatied data set with spatial correlation",
            xlab = "M-statistic", ylab = "Density",
            sub = "The vertical line shows the observed M-statistic")

mantel.test(distancia, dist_ncte, graph = TRUE,
            main = "Mantel test: a simulated data set with spatial correlation",
            xlab = "M-statistic", ylab = "Density",
            sub = "The vertical line shows the observed M-statistic")

#### Paquete APE

ozone<-read.table("http://www.ats.ucla.edu/stat/r/faq/ozone.csv", sep=",", header=T)
head(ozone, 10)

ozone.dists <- as.matrix(dist(cbind(ozone$Lon, ozone$Lat))) # Matriz de distancias
ozone.dists.inv <- 1/ozone.dists # Matriz de distancias invertidas
diag(ozone.dists.inv) <- 0
ozone.dists.inv[1:5, 1:5]

Moran.I(ozone$Av8top, ozone.dists.inv) # Variable de interés + matriz de distancias investidas

### Paquete ADE4: Mide la autocorrelación espacial en todas las variables
# del modelo

# Moran asume homogeneidad sobre la zona de estudio
# (lo cual es dificil de suponer)

# Con los datos de entrenamiento
trainDist <- as.matrix(dist(x=cbind(trainAll$lon,trainAll$lat),method='euclidean'))
trainDist_inv <- 1/trainDist
diag(trainDist_inv) <- 0

ade4_autocorr <- gearymoran(bilis=trainDist_inv,X=trainAll[,3:19],nrepet=1000,alter="two-sided")

# Con los datos de prueba
testnAllFilter <- testnAll[!duplicated(testnAll[,1:2]),]
testDist <- as.matrix(dist(x=cbind(testnAllFilter$lon,testnAllFilter$lat),method='euclidean'))
testDist_inv <- 1/testDist
diag(testDist_inv) <- 0

ade4_autocorrTest <- gearymoran(bilis=testDist_inv,X=testnAllFilter[,3:19],nrepet=1000,alter="two-sided")

# ================================================= #

# Variograms for to explore spatial autocorrelation in explanatory variables
library(geoR)

### -----------------------------------------------
### Training data

trainAll
dists <- dist(trainAll[,1:2])
summary(as.vector(dists))
breaks = seq(0, 40, l = 11)

varList <- names(trainAll)[3:17]

smv <- list()

for(var in varList){
  v1 <- variog(coords = trainAll[,1:2], data = trainAll[,var], breaks = breaks)
  v1.summary <- cbind(c(1:10), v1$v, v1$n)
  colnames(v1.summary) <- c("lag", "semi-variance", "pairs")
  smv[[var]] <- v1.summary
  Sys.sleep(1) # time in seconds
  plot(v1, type = "b", main = paste(var,sep=""))
  rm(v1); rm(v1.summary)
}

### -----------------------------------------------
### Testing data

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
  Sys.sleep(1) # time in seconds
  plot(v1, type = "b", main = paste(var,sep=""))
  rm(v1); rm(v1.summary)
}

# ================================================= #

require(ape)
mantel_res <- mantel.test(testAll[,1:12],nperm=1000,graph=TRUE,alternative="two.sided")


# ================================================= #

# Presences training
#pt1 <- data.frame(mycoord$lon[mycoord$occ==1 & mycoord$id.train1==1],
#                  mycoord$lat[mycoord$occ==1 & mycoord$id.train1==1])
#names(pt1) <- c('lon','lat'); dim(pt1)

# Presences testing
#pp1 <- data.frame(mycoord$lon[mycoord$occ==1 & mycoord$id.train1==0],
#                  mycoord$lat[mycoord$occ==1 & mycoord$id.train1==0])
#names(pp1) <- c('lon','lat'); dim(pp1)

# Absences training
#at1 <- data.frame(mycoord$lon[mycoord$occ==0 & mycoord$id.train1==1],
#                  mycoord$lat[mycoord$occ==0 & mycoord$id.train1==1])
#names(at1) <- c('lon','lat')

#at1 <- at1[!is.na(at1$lon),]; dim(at1)

# Absences testing
#ap1 <- data.frame(mycoord$lon[mycoord$occ==0 & mycoord$id.train1==0],
#                  mycoord$lat[mycoord$occ==0 & mycoord$id.train1==0])
#names(ap1) <- c('lon','lat')

#ap1 <- ap1[!is.na(ap1$lon),]; dim(ap1)

# ======== Spatial sorting bias [Hijmans, 2012]
# Presences testing, Absences testing, Presences training

#sb1 <- ssb(pp1, ap1, pt1)
#sb1[,1]/sb1[,2]

# ======== Pairwise distance for correct SSB

# Cantidad de presencias utilizadas en la validación (25% de la información total)
# 0.25*sum(!is.na(myResp))

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

# sample.p.a1 <- pwdSample(fixed=pp1, sample=ap1, reference=pt1)

# Esta muestra se escoge sobre las pseudo-ausencias del conjunto de validación
# por lo cual es problable que sea un número reducido

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

#total.pa <- rbind(at1,ap1)
#sample.p.a2 <- pwdSample(fixed=pp1, sample=total.pa, reference=pt1,
#                         n=(2*sum(!is.na(myResp)))/(0.25*sum(!is.na(myResp))))

# Esta muestra se escoge sobre el total de pseudo-ausencias (conjuntos de 
# entrenamiento y validación) seleccionado a partir del método ambiental 
# de generación de pseudo ausencias SRE en el paso previo, aquí se tiene 
# en cuenta el número de puntos seleccionados por el método final este 
# puede ser igual al número de presencias, el doble, triple, ...
# y se determina por n=(2*sum(!is.na(myResp)))/(0.25*sum(!is.na(myResp)))
# Aquí se calcula el doble de conjuntos de pseudo-ausencias en relación
# al número de presencias.
# Esto indica que por cada punto de presencia (validación) teniendo en cuenta
# la distancia a los puntos de presencia (entrenamiento) se escogen un determinado
# número de pseudo ausencias para completar el número de pseudo ausencias deseado.

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

# Extraer los puntos seleccionados como pseudo-ausencias a partir del proceso
# de dos pasos
#mysample <- as.vector(na.exclude(as.vector(sample.p.a2)))

# Extraer las coordenadas de los puntos seleccionados por el procedimiento de
# Pairwise distance
#p.a.select2 <- total.pa[mysample,]

# ======== Sesgo de selección espacial para las pseudo-ausencias generadas
#sb <- ssb(pp1, p.a.select2[,1:2], pt1)
#sb[,1]/sb[,2]

# ======== readjust models
oc.dir   <- "/curie_data2/ncastaneda/cluster_variables/pruebas_aposteriori_20052013/occurrences"
myRespXY <- read.csv(paste(oc.dir,"/",myRespName,".csv",sep=""),header=T,sep=",")
myRespXY <- myRespXY[,1:3]
myRespXY <- subset(myRespXY[,2:3],myRespXY$Taxon=="Helianthus_tuberosus")
dups     <- duplicated(myRespXY[])
myRespXY <- myRespXY[!dups,]
myRespXY <- myRespXY[,1:2] # 178 registros

myRespXY_n <- rbind(myRespXY,p.a.select2[,1:2])

myResp <- extract(x=myResp.ras, y=myRespXY)
myResp_n <- c(myResp,rep(0,length(p.a.select2[,1])))

myRespXY_n <- cbind(myRespXY_n,myResp_n)
myRespXY_n <- myRespXY_n[complete.cases(myRespXY_n),] # 157 registros
colnames(myRespXY_n) <- c("lon","lat",myRespName)

myRespXY_n <- myRespXY_n[,1:2]

myResp_n <- c(rep(1,sum(!is.na(myResp))),rep(NA,dim(p.a.select2)[1]))

# ======== load explanatory variables
myExpl <- c("/curie_data2/ncastaneda/cluster_variables/pruebas_aposteriori_20052013/clim/Helianthus_tuberosus-clim/bio_1.asc",
            "/curie_data2/ncastaneda/cluster_variables/pruebas_aposteriori_20052013/clim/Helianthus_tuberosus-clim/bio_2.asc",
            "/curie_data2/ncastaneda/cluster_variables/pruebas_aposteriori_20052013/clim/Helianthus_tuberosus-clim/bio_4.asc",
            "/curie_data2/ncastaneda/cluster_variables/pruebas_aposteriori_20052013/clim/Helianthus_tuberosus-clim/bio_5.asc",
            "/curie_data2/ncastaneda/cluster_variables/pruebas_aposteriori_20052013/clim/Helianthus_tuberosus-clim/bio_6.asc",
            "/curie_data2/ncastaneda/cluster_variables/pruebas_aposteriori_20052013/clim/Helianthus_tuberosus-clim/bio_8.asc",
            "/curie_data2/ncastaneda/cluster_variables/pruebas_aposteriori_20052013/clim/Helianthus_tuberosus-clim/bio_9.asc",
            "/curie_data2/ncastaneda/cluster_variables/pruebas_aposteriori_20052013/clim/Helianthus_tuberosus-clim/bio_10.asc",
            "/curie_data2/ncastaneda/cluster_variables/pruebas_aposteriori_20052013/clim/Helianthus_tuberosus-clim/bio_11.asc",
            "/curie_data2/ncastaneda/cluster_variables/pruebas_aposteriori_20052013/clim/Helianthus_tuberosus-clim/bio_12.asc",
            "/curie_data2/ncastaneda/cluster_variables/pruebas_aposteriori_20052013/clim/Helianthus_tuberosus-clim/bio_13.asc",
            "/curie_data2/ncastaneda/cluster_variables/pruebas_aposteriori_20052013/clim/Helianthus_tuberosus-clim/bio_14.asc",
            "/curie_data2/ncastaneda/cluster_variables/pruebas_aposteriori_20052013/clim/Helianthus_tuberosus-clim/bio_15.asc",
            "/curie_data2/ncastaneda/cluster_variables/pruebas_aposteriori_20052013/clim/Helianthus_tuberosus-clim/bio_16.asc",
            "/curie_data2/ncastaneda/cluster_variables/pruebas_aposteriori_20052013/clim/Helianthus_tuberosus-clim/bio_17.asc",
            "/curie_data2/ncastaneda/cluster_variables/pruebas_aposteriori_20052013/clim/Helianthus_tuberosus-clim/bio_18.asc",
            "/curie_data2/ncastaneda/cluster_variables/pruebas_aposteriori_20052013/clim/Helianthus_tuberosus-clim/bio_19.asc"
)

myExpl <- lapply(myExpl,raster)
myExpl_n <- stack(myExpl)

# ======== BIOMOD2 formating data
myBiomodData_n <- BIOMOD_FormatingData(resp.var = myResp_n,
                                       resp.xy = myRespXY_n,
                                       expl.var = myExpl_n,
                                       resp.name = myRespName,
                                       PA.nb.rep = 0,
                                       na.rm = FALSE)

# ======== BIOMOD2 modeling options
myBiomodOption_n <- BIOMOD_ModelingOptions(GLM = list(myFormula=Helianthus.tuberosus~bio_1+bio_2+bio_4+bio_5+bio_6+bio_8+bio_9+bio_10+bio_11+bio_12+bio_13+bio_14+bio_15+bio_16+bio_17+bio_18, family='binomial'),
                                           ANN = list(NbCV=10, maxit=500),
                                           RF = list(do.classif=TRUE,ntree=100),
                                           MAXENT = list( path_to_maxent.jar = getwd(),
                                                          maximumiterations = 500,
                                                          visible = FALSE,
                                                          linear = TRUE,
                                                          quadratic = TRUE,
                                                          product = TRUE,
                                                          threshold = TRUE,
                                                          hinge = TRUE,
                                                          lq2lqptthreshold = 80,
                                                          l2lqthreshold = 10,
                                                          hingethreshold = 15,
                                                          beta_threshold = -1,
                                                          beta_categorical = -1,
                                                          beta_lqp = -1,
                                                          beta_hinge = -1,
                                                          defaultprevalence = 0.5)
                                                          )

# ======== ajustar los modelos
myBiomodModelOut_n <- BIOMOD_Modeling(myBiomodData_n,
                                      models = c('GLM','ANN','RF','MAXENT'),
                                      models.options = myBiomodOption_n,
                                      NbRunEval=1,
                                      DataSplit=75,
                                      Prevalence=0.5,
                                      VarImport=0,
                                      models.eval.meth = c('TSS','KAPPA','ROC'),
                                      SaveObj = TRUE,
                                      rescal.all.models = TRUE,
                                      do.full.models = FALSE,
                                      modeling.id = paste(myRespName,"_final",sep="")
                                      )

# ======== leer modelos ajustados
load("/home/haachicanoy/Helianthus.tuberosus/Helianthus.tuberosus.Helianthus_tuberosus_final.models.out")
myBiomodModelOut_n <- Helianthus.tuberosus.Helianthus_tuberosus_final.models.out
myBiomodModelOut_n

# ======== construir las proyecciones
myBiomodProj <- BIOMOD_Projection(modeling.output = myBiomodModelOut_n,
                                    new.env = myExpl_n,
                                    proj.name = "current",
                                    selected.models = 'all',
                                    binary.meth = "TSS",
                                    compress = "xz",
                                    clamping.mask = F,
                                    output.format = "GTiff")

#png(paste(myRespName,'_mxt.png',sep=''),width=4096,height=2048,res=300,pointsize=13)
#plot(myBiomodProj, str.grep = 'MAXENT')
#dev.off()

# ======== load projections
load("/home/haachicanoy/Helianthus.tuberosus/proj_current/Helianthus.tuberosus.current.projection.out")
myBiomodProj <- Helianthus.tuberosus.current.projection.out
myBiomodProj

myCurrentProj <- get_predictions(myBiomodProj)

########################################################
# STEP 2
########################################################

# --------------------------------- Test de hipótesis sobre las variables
# Explicativas para realizar la selección de las mismas a ingresar en cada
# modelo

# Realizar test de hipótesis sobre los puntos de presencia vs las pseudo-ausencias
# seleccionadas idenficando las variables en donde se observen las diferencias
# más marcadas

expl_var <- extract(x=myExpl_n, y=myRespXY_n)
expl_var <- as.data.frame(expl_var)

hist(expl_var$bio_8[-(1:157)],prob=F,main="Bio 1",add=F,border=NA,col=2)
hist(expl_var$bio_8[1:157],prob=F,main="Bio 1",add=T,border=NA,col=1)

imp.bios <- 0
for(i in 1:ncol(expl_var)){
  
  x <- expl_var[1:sum(!is.na(myResp_n)),i]
  y <- expl_var[-(1:sum(!is.na(myResp_n))),i]
  
  imp.bios[i] <- wilcox.test(x,y,alternative="two.sided",mu=0,paired=F)$p.value
  
}

names(imp.bios) <- names(expl_var)

#png("Wilcox.png",width=2048,height=2048,res=300,pointsize=13)
par(mar=c(4.5,6,4,2))
barplot(sort(imp.bios),,horiz=T,las=1,border=F,
        xlab='p-value (Wilcoxon Test)',main="")
abline(v=0.05,col=2,lty=3)
#dev.off()

res_imp <- data.frame(names(expl_df),imp.bios)
names(res_imp) <- c("Variables","p_value")

########################################################
# STEP 3
########################################################

# --------------------------------- Selección de variables explicativas
# omitiendo las variables que presenten una alta colinealidad con el resto
# de variables a incluir en el modelo

v_select <- vifstep(expl_var,th=10)

# Selección final

res_vif <- v_select@results
sel_var <- merge(res_imp,res_vif,by.x="Variables")

# retornar vector con los nombres de las variables candidatas

as.vector(sel_var$Variables[which(sel_var$p_value < 0.05)])

########################################################
# STEP 4
########################################################

load()

?nnet()

for(n in nUnits){
  
  fit <- nnet(v_resp ~ bio_i+bio_j+bio_k+..., size=n)
  error[n] <- sum(fit$residuals^2)
  
}

n_sel <- which.max(error)
n_sel # number of hidden unities to neural network

?glm()

?randomForest()

?MAXENT()

# ================== Ensemble modeling ================== #

myBiomodEM <- BIOMOD_EnsembleModeling(modeling.output = myBiomodModelOut_n,
                                      chosen.models = 'all',
                                      em.by='algo',
                                      eval.metric = c('TSS'),
                                      #eval.metric.quality.threshold = c(0.7),
                                      prob.mean = T,
                                      prob.cv = F,
                                      prob.ci = F,
                                      prob.ci.alpha = 0.05,
                                      prob.median = T,
                                      committee.averaging = T,
                                      prob.mean.weight = T,
                                      prob.mean.weight.decay = 'proportional'
                                      )

# ================== Get evaluations ================== #

get_evaluations(myBiomodEM)

# ================== Ensemble forecasting ================== #

# tuve problemas con esta parte

myBiomodEF <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEM,
                                         projection.output = myBiomomodProj)

plot(myBiomodEF)


