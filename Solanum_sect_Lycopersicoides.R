# H. Achicanoy
# CIAT, 2013

# ----------------------------------------------------------------------------- #
# ----------------------------------------------------------------------------- #
# ----------------------------------------------------------------------------- #
#                         Solanum (sect lycopersicoides) models
# ----------------------------------------------------------------------------- #
# ----------------------------------------------------------------------------- #
# ----------------------------------------------------------------------------- #

# ======== load packages
library(biomod2); library(dismo); library(BIOMOD); library(usdm); library(stringr)

# ======== call code
crop <- "tomato"
src.dir <- paste("/curie_data2/ncastaneda/gap-analysis/gap_",crop,"/_scripts",sep="") # !!! change accordingly !!!
source(paste(src.dir,"/000.zipRead.R",sep=""))

# ======== set working directory
wd <- "/curie_data2/ncastaneda/gap-analysis/gap_tomatoNHM"
dir <- "/curie_data2/ncastaneda/gap-analysis/gap_tomato"
spp <- "259"

# setwd("/home/haachicanoy")
setwd(wd)

# ======== load native area
# na.dir <- "/curie_data2/ncastaneda/cluster_variables/pruebas_aposteriori_20052013/nareas"
na.dir <- paste(dir,"/biomod_modeling/native-areas/asciigrids/",spp,sep="")
# myResp.ras2 <- list.files(na.dir, pattern=".asc", full.names=T)
myResp.ras <- zipRead(na.dir,"narea.asc.gz")

# ======== read response variable
# myRespName <- "Helianthus_tuberosus"
myRespName <- spp

oc.dir   <- paste(dir,"/occurrence_files/",sep="")
myRespXY <- read.csv(paste(oc.dir,"/",myRespName,".csv",sep=""),header=T,sep=",")
myRespXY <- myRespXY[,1:3]
myRespXY <- subset(myRespXY[,2:3],myRespXY$Taxon==myRespName)
dups     <- duplicated(myRespXY[])
myRespXY <- myRespXY[!dups,]
myRespXY <- myRespXY[,1:2] # 178 registros

myResp   <- extract(x=myResp.ras, y=myRespXY)
myRespXY <- cbind(myRespXY, myResp)

# coordenadas de puntos dentro del área nativa
myRespXY <- myRespXY[complete.cases(myRespXY),] # 157 registros
colnames(myRespXY) <- c("lon","lat",myRespName)

myResp   <- myRespXY[,3]
myRespXY <- SpatialPoints(myRespXY[,1:2])

# ======== read explanatory variables
# myExpl <- c("/curie_data2/ncastaneda/cluster_variables/pruebas_aposteriori_20052013/clim/Helianthus_tuberosus-clim/bio_1.asc",
#             "/curie_data2/ncastaneda/cluster_variables/pruebas_aposteriori_20052013/clim/Helianthus_tuberosus-clim/bio_2.asc",
#             "/curie_data2/ncastaneda/cluster_variables/pruebas_aposteriori_20052013/clim/Helianthus_tuberosus-clim/bio_4.asc",
#             "/curie_data2/ncastaneda/cluster_variables/pruebas_aposteriori_20052013/clim/Helianthus_tuberosus-clim/bio_5.asc",
#             "/curie_data2/ncastaneda/cluster_variables/pruebas_aposteriori_20052013/clim/Helianthus_tuberosus-clim/bio_6.asc",
#             "/curie_data2/ncastaneda/cluster_variables/pruebas_aposteriori_20052013/clim/Helianthus_tuberosus-clim/bio_8.asc",
#             "/curie_data2/ncastaneda/cluster_variables/pruebas_aposteriori_20052013/clim/Helianthus_tuberosus-clim/bio_9.asc",
#             "/curie_data2/ncastaneda/cluster_variables/pruebas_aposteriori_20052013/clim/Helianthus_tuberosus-clim/bio_10.asc",
#             "/curie_data2/ncastaneda/cluster_variables/pruebas_aposteriori_20052013/clim/Helianthus_tuberosus-clim/bio_11.asc",
#             "/curie_data2/ncastaneda/cluster_variables/pruebas_aposteriori_20052013/clim/Helianthus_tuberosus-clim/bio_12.asc",
#             "/curie_data2/ncastaneda/cluster_variables/pruebas_aposteriori_20052013/clim/Helianthus_tuberosus-clim/bio_13.asc",
#             "/curie_data2/ncastaneda/cluster_variables/pruebas_aposteriori_20052013/clim/Helianthus_tuberosus-clim/bio_14.asc",
#             "/curie_data2/ncastaneda/cluster_variables/pruebas_aposteriori_20052013/clim/Helianthus_tuberosus-clim/bio_15.asc",
#             "/curie_data2/ncastaneda/cluster_variables/pruebas_aposteriori_20052013/clim/Helianthus_tuberosus-clim/bio_16.asc",
#             "/curie_data2/ncastaneda/cluster_variables/pruebas_aposteriori_20052013/clim/Helianthus_tuberosus-clim/bio_17.asc",
#             "/curie_data2/ncastaneda/cluster_variables/pruebas_aposteriori_20052013/clim/Helianthus_tuberosus-clim/bio_18.asc",
#             "/curie_data2/ncastaneda/cluster_variables/pruebas_aposteriori_20052013/clim/Helianthus_tuberosus-clim/bio_19.asc")

myExpl <- c("/curie_data2/ncastaneda/gap-analysis/gap_tomato/biomod_modeling/current-clim/bio_1.asc",
            "/curie_data2/ncastaneda/gap-analysis/gap_tomato/biomod_modeling/current-clim/bio_2.asc",
            "/curie_data2/ncastaneda/gap-analysis/gap_tomato/biomod_modeling/current-clim/bio_4.asc",
            "/curie_data2/ncastaneda/gap-analysis/gap_tomato/biomod_modeling/current-clim/bio_5.asc",
            "/curie_data2/ncastaneda/gap-analysis/gap_tomato/biomod_modeling/current-clim/bio_6.asc",
            "/curie_data2/ncastaneda/gap-analysis/gap_tomato/biomod_modeling/current-clim/bio_8.asc",
            "/curie_data2/ncastaneda/gap-analysis/gap_tomato/biomod_modeling/current-clim/bio_9.asc",
            "/curie_data2/ncastaneda/gap-analysis/gap_tomato/biomod_modeling/current-clim/bio_10.asc",
            "/curie_data2/ncastaneda/gap-analysis/gap_tomato/biomod_modeling/current-clim/bio_11.asc",
            "/curie_data2/ncastaneda/gap-analysis/gap_tomato/biomod_modeling/current-clim/bio_12.asc",
            "/curie_data2/ncastaneda/gap-analysis/gap_tomato/biomod_modeling/current-clim/bio_13.asc",
            "/curie_data2/ncastaneda/gap-analysis/gap_tomato/biomod_modeling/current-clim/bio_14.asc",
            "/curie_data2/ncastaneda/gap-analysis/gap_tomato/biomod_modeling/current-clim/bio_15.asc",
            "/curie_data2/ncastaneda/gap-analysis/gap_tomato/biomod_modeling/current-clim/bio_16.asc",
            "/curie_data2/ncastaneda/gap-analysis/gap_tomato/biomod_modeling/current-clim/bio_17.asc",
            "/curie_data2/ncastaneda/gap-analysis/gap_tomato/biomod_modeling/current-clim/bio_18.asc",
            "/curie_data2/ncastaneda/gap-analysis/gap_tomato/biomod_modeling/current-clim/bio_19.asc")

myExpl <- lapply(myExpl,raster)
myExpl <- raster::stack(myExpl)

# ======== BIOMOD formating data
# Generación de pseudo-ausencias mediante el método "sre"
# producir un set cuantioso de pseudo-ausencias

myBiomodData <- BIOMOD_FormatingData(resp.var = myRespXY,
                                     expl.var = myExpl,
                                     resp.name = myRespName,
                                     PA.strategy = 'sre',
                                     PA.nb.rep = 3,
                                     PA.nb.absences = 10000,
                                     PA.sre.quant = 0.10
)


myBiomodData <- BIOMOD_FormatingData(resp.var       = myRespXY[,3],
                                     expl.var       = myExpl,
                                     resp.xy        = myRespXY[,1:2],
                                     resp.name      = myRespName,
                                     PA.strategy    = 'random',
                                     PA.nb.rep      = 10,
                                     PA.nb.absences = 10000,
                                     PA.sre.quant   = 0.10
)

# ======== BIOMOD modeling options
# myBiomodOption <- BIOMOD_ModelingOptions(GLM = list(myFormula=Helianthus_tuberosus~bio_5+bio_10, family='binomial'))
myBiomodOption <- BIOMOD_ModelingOptions(GLM = list(myFormula=spp~bio_5+bio_10, family='binomial'))

# ======== BIOMOD running options
myBiomodModelOut <- BIOMOD_Modeling(myBiomodData,
                                    models = c('GLM'),
                                    models.options = myBiomodOption,
                                    NbRunEval = 1,
                                    DataSplit = 75,
                                    Prevalence = 0.5,
                                    VarImport = 0,
                                    models.eval.meth = c('KAPPA','TSS','ROC'),
                                    SaveObj = TRUE,
                                    rescal.all.models = TRUE,
                                    do.full.models = FALSE,
                                    modeling.id = paste(myRespName,"_sre",sep=""))

# ======== Extract information from the models
load("/home/haachicanoy/Helianthus.tuberosus/Helianthus.tuberosus.Helianthus_tuberosus_sre.models.out")

myBiomodModelOut <- Helianthus.tuberosus.Helianthus_tuberosus_sre.models.out
myBiomodModelOut

# ======== read training data
load(myBiomodModelOut@calib.lines@link)
load(myBiomodModelOut@formated.input.data@link)

mycoord <- data.frame(data@coord)
names(mycoord) <- c("lon","lat")
mycoord['occ'] <- 0
mycoord[(1:sum(myResp,na.rm=T)),3] <- 1

calib.lines <- as.data.frame(calib.lines)

for(i in 1:3){ # npa_rep es el número de conjuntos de pseudo-ausencias generados
  mycoord[paste('id.train',i,sep="")] <- as.numeric(calib.lines[,i])
}

# ======== data training 1

# Presences training
pt1 <- data.frame(mycoord$lon[mycoord$occ==1 & mycoord$id.train1==1],
                  mycoord$lat[mycoord$occ==1 & mycoord$id.train1==1])
names(pt1) <- c('lon','lat'); dim(pt1)

# Presences testing
pp1 <- data.frame(mycoord$lon[mycoord$occ==1 & mycoord$id.train1==0],
                  mycoord$lat[mycoord$occ==1 & mycoord$id.train1==0])
names(pp1) <- c('lon','lat'); dim(pp1)

# Absences training
at1 <- data.frame(mycoord$lon[mycoord$occ==0 & mycoord$id.train1==1],
                  mycoord$lat[mycoord$occ==0 & mycoord$id.train1==1])
names(at1) <- c('lon','lat')

at1 <- at1[!is.na(at1$lon),]; dim(at1)

# Absences testing
ap1 <- data.frame(mycoord$lon[mycoord$occ==0 & mycoord$id.train1==0],
                  mycoord$lat[mycoord$occ==0 & mycoord$id.train1==0])
names(ap1) <- c('lon','lat')

ap1 <- ap1[!is.na(ap1$lon),]; dim(ap1)

# ======== Spatial sorting bias [Hijmans, 2012]
# Presences testing, Absences testing, Presences training

sb1 <- ssb(pp1, ap1, pt1)
sb1[,1]/sb1[,2]

# ======== Pairwise distance for correct SSB

# Cantidad de presencias utilizadas en la validación (25% de la información total)
0.25*sum(!is.na(myResp))

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

# sample.p.a1 <- pwdSample(fixed=pp1, sample=ap1, reference=pt1)

# Esta muestra se escoge sobre las pseudo-ausencias del conjunto de validación
# por lo cual es problable que sea un número reducido

# :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: #

total.pa <- rbind(at1,ap1)
sample.p.a2 <- pwdSample(fixed=pp1, sample=total.pa, reference=pt1,
                         n=(2*sum(!is.na(myResp)))/(0.25*sum(!is.na(myResp))))

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
mysample <- as.vector(na.exclude(as.vector(sample.p.a2)))

# Extraer las coordenadas de los puntos seleccionados por el procedimiento de
# Pairwise distance
p.a.select2 <- total.pa[mysample,]

# ======== Sesgo de selección espacial para las pseudo-ausencias generadas
sb <- ssb(pp1, p.a.select2[,1:2], pt1)
sb[,1]/sb[,2]

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
                                           ANN = list(NbCV=5, maxit=500),
                                           RF = list(do.classif=TRUE,ntree=100),
                                           MAXENT = list( path_to_maxent.jar = getwd(),
                                                          maximumiterations = 200,
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
                                  output.format = ".grd")

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


# Muestra completa presencia + pseudo ausencias
mycoord.s

# Seleccionar las coordenadas de los sitios escogidos a partir de la distancia
# pairwise
p.select <- mycoord.s[sample.p.a[complete.cases(sample.p.a)],]
pp1[30,]

png('pa_select.png',width=2048,height=2048,res=300,pointsize=13)
plot(myResp.ras,xlab='lon',ylab='lat',xlim=c(-130,-50),ylim=c(20,60)) 
points(p.select[,1:2],pch=20,cex=0.6)
points(pp1[30,],pch=20,cex=0.6,col=2)
dev.off()

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






















