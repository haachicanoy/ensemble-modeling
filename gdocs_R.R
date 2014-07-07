# Expert analysis
# H. Achicanoy
# CIAT, 2014

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
# Paths
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
# Load packages & source code
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

require(XML)
require(RCurl)
require(RGoogleDocs)

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
# Read data online
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

connection   = getGoogleDocsConnection(getGoogleAuth("harold22010@gmail.com","welcometohell1",service="wise"))
spreadsheets = getDocs(connection, ssl.verifypeer=FALSE)
sheets       = names(spreadsheets)

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
# Classify data
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

# Priority crops
pList <- c("Avena","Cajanus","Cicer","Daucus","Eggplant",
           "Eleusine","Helianthus","Hordeum","Ipomoea",
           "Lathyrus","Lens","Malus","Medicago","Musa",
           "Pennisetum","Phaseolus","Pisum","Potato",
           "Rice","Secale","Sorghum","Vicia","Vigna",
           "Wheat")

# List of priority taxa (Utilizar CWR_GA_Data)
taxaPrior = getWorksheets(spreadsheets[["PtaxaofPcrops_all"]],connection)
comp      = taxaPrior$"names_codes_new455"; # names_codes_newest(445)_2014_1
taxaPrior = sheetAsMatrix(comp,header=T,as.data.frame=T,trim=T)
rm(comp)

# Non priority crops
npList <- c("allium","asparagus","beet","brassica","breadfruit",
            "cacao","capsicum","cassava","citrus","cocoyam",
            "cotton","cucumis","dioscorea","grape","groundnut",
            "lettuce","maize","mango","millet_panicum",
            "millet_setaria","papaya","pear","pineapple",
            "prunus","quinoa","safflower","soy","spinach",
            "squash","strawberry","sugar_cane","tomato","vigna",
            "watermelon")

# List of non priority taxa (Utilizar NonP_PrioritiesTables)
taxaNonPrior = getWorksheets(spreadsheets[["CWR_non-priority_crops_processing"]],connection)
comp         = taxaNonPrior$"non_priority_taxa-2013-09-13"
taxaNonPrior = sheetAsMatrix(comp,header=T,as.data.frame=T,trim=T)
taxaNonPrior = taxaNonPrior[,1:12]
rm(comp)

# Publications crops
publs <- c("Helianthus(Paper)")

### =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= ###
### Scores experts evaluation
### =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= ###

g_shts = sheets[grep(pattern="CWR_Expert_Evaluation_", sheets)]
g_shts = unique(g_shts); g_shts = sort(g_shts)
# g_shts   = unlist(strsplit(g_shts, split=" "))
# g_shts   = g_shts[seq(1,length(g_shts),by=2)]

# =-=-= Priority crops
gp_shts = g_shts[unlist(sapply(pList,function(x){grep(pattern=paste("_",x," ",sep=""),g_shts,fixed=T)}))] #"\\", "$"

# =-=-= Non priority
gnp_shts = g_shts[unlist(sapply(npList,function(x){grep(pattern=paste("_",x," ",sep=""),g_shts,fixed=T)}))] #"\\"

# =-=-= Publications crops
pub_shts = g_shts[unlist(sapply(publs,function(x){grep(pattern=paste("_",x," ",sep=""),g_shts,fixed=T)}))] #"\\", "$"

# Organizar la información
g_crop = unlist(strsplit(g_shts, split=" "))
g_crop = g_crop[seq(1,length(g_crop),by=2)]
g_crop = unlist(strsplit(g_crop, split="CWR_Expert_Evaluation_",fixed=T))
g_crop = g_crop[seq(2,length(g_crop),by=2)]

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
# =-=-=-=-= General experts evaluation =-=-=-=-=- #
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

require(gdata)

count = 1
for(i in 1:length(g_shts)){
  
  cat("\n")
  
  cat("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= \n")
  cat("Procesing information for",g_shts[i],"\n")
  cat("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= \n")
  
  cat("\n")
  
  # Crear conexión con el archivo en gdocs
  con1 = getWorksheets(spreadsheets[[g_shts[i]]],connection)
  comp = con1$"Form Responses"
  con2 = try(sheetAsMatrix(comp,header=T,as.data.frame=T,trim=T),silent=T)
  
  # Revisar si hay datos (nrow > 0)
  if(is.data.frame(con2)){ # & nrow(con2) > 0
    
    # Revisar si los datos digitados por cada experto son validos
    aux = 1
    for(nr in 1:nrow(con2)){
      
      info = as.vector(t(con2[nr,])) # Extraer la información de cada experto
      responses = c("High Priority: 0","No Need for Further Collection: 10","Strongly Disagree","Disagree","Neutral","Agree","Strongly Agree",paste(1:9,sep=""))
      
      if(!any(info %in% responses)){
        
        cat("The information typed in the row",nr,"isn't valid \n")
        cat("Storing information \n")
        
        if(aux == 1){
          rows = nr
        } else {
          rows = c(rows,nr)
          rows = as.numeric(rows)
        }
        
        aux = aux + 1
        
      } else {
        cat("The information typed in the row",nr,"is valid \n")
      }
      
    }
    
    cat("\n")
    cat("\n")
    rm(aux); rm(nr)
    rm(info); rm(responses)
    
    # Si existe información invalida, eliminar los registros
    if(exists("rows")){
      cat("Delete information \n")
      cat("\n")
      cat("\n")
      con2 = con2[-rows,]
      rm(rows)
    }
    
    rownames(con2) = 1:nrow(con2)
    
    # Etiquetar la información de cada experto
    for(j in 1:nrow(con2)){
      con2$Expert[j] = paste("Expert_",j,sep="")
    }; rm(j)
    
    ### =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= ###
    ### TEXTUAL DATA
    ### =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= ###
    
    ### Data by crop
    # Expert_id, Expert_nm, Position,Institution, e-mail, Notes, 
    # Additional_taxa, Comments_gap_analysis, Priority_crop, Crop_code
    
    data_crop = data.frame(Expert_id       = con2$Expert,
                           Expert_nm       = con2$Name,
                           Position        = con2$"Position title",
                           Institution     = con2$Institution,
                           e_mail          = con2$Email,
                           Notes           = con2$Notes,
                           Additional_taxa = if(is.null(con2$Additional_Taxa)){NA} else{con2$Additional_Taxa},
                           Eval_gap_scores = if(is.null(con2$"Evaluation of Gap Analysis Scores")){NA} else{con2$"Evaluation of Gap Analysis Scores"},
                           Comments_gap_analysis = if(is.null(con2$"Comments on Gap Analysis Scores")){NA} else{con2$"Comments on Gap Analysis Scores"})
    
    #data_crop = data.frame(Expert_id       = if(is.null(con2$Expert)){NA} else{con2$Expert},
    #                       Expert_nm       = if(is.null(con2$Name)){NA} else{con2$Name},
    #                       Position        = if(is.null(con2$"Position title")){NA} else{con2$"Position title"},
    #                       Institution     = if(is.null(con2$Institution)){NA} else{con2$Institution},
    #                       e_mail          = if(is.null(con2$Email)){NA} else{con2$Email},
    #                       Notes           = if(is.null(con2$Notes)){NA} else{con2$Notes},
    #                       Additional_taxa = if(is.null(con2$Additional_Taxa)){NA} else {con2$Additional_Taxa},
    #                       Comments_gap_analysis = if(is.null(con2$"Comments on Gap Analysis Scores")){NA} else{con2$"Comments on Gap Analysis Scores"})
    
    data_crop$Priority_crop = if(g_crop[i] %in% pList){1}else{0}
    data_crop$Crop_code = g_crop[i]
    
    if(count == 1){
      text_data = data_crop
    } else {
      text_data = rbind(text_data,data_crop)
    }
    
    ### =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= ###
    ### QUANTITATIVE DATA
    ### =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= ###
    
    ### Data by taxon
    # Priority_crop, Crop_code, Taxon, Expert_id, Expert_nm,
    # Position, Comparable, Contextual, Evaluation
    
    # List of taxa by crop
    if(g_crop[i] %in% pList){ # PRIORITY CROPS
      if(g_crop[i] == "Phaseolus"){ # Phaseolus case
        taxList = sort(c(as.character(taxaPrior$Taxon_final[taxaPrior$Crop_code==tolower("Bean")]),
                         as.character(taxaPrior$Taxon_final[taxaPrior$Crop_code==tolower("Lima_bean")])))
      } else {
        if(g_crop[i] == "Vicia"){   # Vicia case
          taxList = sort(c(as.character(taxaPrior$Taxon_final[taxaPrior$Crop_code==tolower("Faba_bean")]),
                           as.character(taxaPrior$Taxon_final[taxaPrior$Crop_code==tolower("Vetch")])))
        } else {
          if(g_crop[i] == "Wheat"){ # Wheat case
            taxList = sort(as.character(taxaPrior$Taxon_final[taxaPrior$Crop_code==tolower("Triticum")]))
          } else {
            if(g_crop[i] == "Vigna"){
              taxList = sort(c(as.character(taxaPrior$Taxon_final[taxaPrior$Crop_code==tolower("Bambara")]),
                               as.character(taxaPrior$Taxon_final[taxaPrior$Crop_code==tolower("Cowpea")])))
            } else {
              taxList = sort(as.character(taxaPrior$Taxon_final[taxaPrior$Crop_code==tolower(g_crop[i])]))
            }
          }
        }
      }
    } else {
      if(g_crop[i] %in% npList){ # NON PRIORITY CROPS
        if(g_crop[i] == "dioscorea"){ # Dioscorea case
          taxList = sort(c(as.character(taxaNonPrior$taxon[taxaNonPrior$crop_code==tolower("yam_lagos")]),
                           as.character(taxaNonPrior$taxon[taxaNonPrior$crop_code==tolower("yam_water")]),
                           as.character(taxaNonPrior$taxon[taxaNonPrior$crop_code==tolower("yam_whiteguinea")]),
                           "Dioscorea_cayennensis_subsp._rotundata"))
          taxList = unique(taxList)
        } else {
          if(g_crop[i] == "soy"){
            taxList = sort(as.character(taxaNonPrior$taxon[taxaNonPrior$crop_code==tolower("soybean")]))
          } else {
            taxList = sort(as.character(taxaNonPrior$taxon[taxaNonPrior$crop_code==tolower(g_crop[i])]))
          }
        }
      } else {
        if(g_crop[i] %in% publs){ # PUBLICATION CROPS
          taxList = sort(c(as.character(taxaPrior$Taxon_final[taxaPrior$Crop_code==tolower("Helianthus")]),"Helianthus_pauciflorus_subsp._pauciflorus","Helianthus_pauciflorus_subsp._subrhomboideus","Helianthus_winteri"))
        }
      }
    }
    
    taxList = as.vector(taxList)
    cols    = as.vector(colnames(con2)) # Names of columns in data set
    
    # Intersect names of taxa
    mtch = unlist(sapply(taxList, function(x){grep(pattern=paste("[",x,"]",sep=""), cols, fixed=TRUE)}))
    
    if(is.null(dim(mtch))){
      mtch = unlist(sapply(names(mtch),FUN=function(x){substr(x, 1, nchar(x)-1)}))
      mtch = sort(unique(mtch))
    } else {
      mtch = colnames(mtch)
    }
    
    #--- Comparable data
    comp = matchcols(con2, with=paste("Comparable Expert Priority Score [",mtch,"]",sep=""), method="or", fixed=T)
    names(comp) = NULL; comparable = con2[,comp]; rm(comp)
    
    comparable = t(comparable); colnames(comparable) = con2$Expert
    comparable[which(comparable=="High Priority: 0")] <- "0"
    comparable[which(comparable=="No Need for Further Collection: 10")] <- "10"
    options(warn=-1)
    class(comparable) = "numeric"
    
    comparable = as.data.frame(comparable)
    comparable = stack(comparable, select=colnames(comparable))
    colnames(comparable) = c("Comparable","Expert_id")
    
    #--- Contextual
    cont = matchcols(con2, with=paste("Contextual Expert Priority Score [",mtch,"]",sep=""), method="or", fixed=T)
    names(cont) = NULL; contextual = con2[,cont]; rm(cont)
    
    contextual = t(contextual); colnames(contextual) = con2$Expert
    contextual[which(contextual=="High Priority: 0")] <- "0"
    contextual[which(contextual=="No Need for Further Collection: 10")] <- "10"
    options(warn=-1)
    class(contextual) = "numeric"
    
    contextual = as.data.frame(contextual)
    contextual = stack(contextual, select=colnames(contextual))
    colnames(contextual) = c("Contextual","Expert_id")
    
    #--- Evaluation
    eval = matchcols(con2, with=paste("Evaluation of Gap Analysis Scores [",mtch,"]",sep=""), method="or", fixed=T)
    names(eval) = NULL; evaluation = con2[,eval]; rm(eval)
    
    evaluation$Expert = con2$Expert
    evaluation$time = 1
    
    evaluation = reshape(evaluation, direction="long", v.names="Evaluation", idvar="Expert", varying=1:length(mtch))
    evaluation = evaluation[order(evaluation$Expert),]
    rownames(evaluation) = 1:nrow(evaluation)
    evaluation$time = NULL
    
    #--- Position
    ne = 1
    for(j in 1:nrow(con2)){
      if(ne == 1){
        position = rep(con2$"Position title"[j],length(mtch))
        position = as.character(position)
      } else {
        position = c(position, as.character(rep(con2$"Position title"[j],length(mtch))))
      }
      ne = ne + 1
    }; rm(j); rm(ne)
    
    #--- Name
    ne = 1
    for(j in 1:nrow(con2)){
      if(ne == 1){
        name = rep(con2$Name[j],length(mtch))
        name = as.character(name)
      } else {
        name = c(name, as.character(rep(con2$Name[j],length(mtch))))
      }
      ne = ne + 1
    }; rm(j); rm(ne)
    
    data_taxon = data.frame(Expert_id     = comparable$Expert_id,
                            Expert_nm     = name,
                            Position      = position,
                            Comparable    = comparable$Comparable,
                            Contextual    = contextual$Contextual,
                            Evaluation    = evaluation$Evaluation,
                            Priority_crop = if(g_crop[i] %in% pList){1}else{0},
                            Crop_code     = g_crop[i],
                            Taxon         = rep(mtch,nrow(con2)))
    
    rm(list=c("name","position","comparable","contextual","evaluation"))
    
    if(count == 1){
      quan_data = data_taxon
    } else {
      quan_data = rbind(quan_data,data_taxon)
    }
    
    rm(cols); rm(taxList)
    
    ### =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= ###
    
    count = count + 1
    
    ### =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= ###
    
  } else {
    cat("Awaiting for experts evaluation \n")
    cat("\n")
    cat("\n")
    rm(con2)
  }
  
}

rm(con1); rm(con2); rm(count)
rm(data_crop); rm(data_taxon)
rm(i); rm(mtch)

sumData <- data.frame(Crop_code = tolower(unique(quan_data$Crop_code)),
                      nTaxa     = unlist(sapply(as.character(unique(quan_data$Crop_code)),FUN=function(x){length(unique(quan_data$Taxon[quan_data$Crop_code==x]))})),
                      nExperts  = unlist(sapply(as.character(unique(quan_data$Crop_code)),FUN=function(x){length(unique(quan_data$Expert_id[quan_data$Crop_code==x]))})),
                      Priority  = unlist(sapply(as.character(unique(quan_data$Crop_code)),FUN=function(x){unique(quan_data$Priority_crop[quan_data$Crop_code==x])})))
rownames(sumData) <- 1:nrow(sumData)

sumData

### =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= ###
### Maps experts evaluation
### =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= ###

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#
#=-=-=-=-=-=- Maps experts evaluation =-=-=-=-=-#
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=#

m_shts   = sheets[grep(pattern="CWR_Map_Evaluation_", sheets)]
# m_shts   = unlist(strsplit(m_shts, split=" "))
# m_shts   = m_shts[seq(1,length(m_shts),by=2)]

#=-=-= Priority crops
mp_shts  = m_shts[unlist(sapply(pList,function(x){grep(pattern=paste("_",x," ",sep=""),m_shts)}))] #"\\"
#=-=-= Non priority
mnp_shts = m_shts[unlist(sapply(npList,function(x){grep(pattern=paste("_",x," ",sep=""),m_shts)}))] #"\\"




test = getWorksheets(spreadsheets[["CWR_Expert_Evaluation_Daucus (Responses)"]],connection)
names(test)
sup = test$"Form Responses"
x = sheetAsMatrix(sup,header=T,as.data.frame=T,trim=T)






# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
# Arrange data
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
