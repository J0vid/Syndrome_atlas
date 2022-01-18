library(plumber)
library(future)
library(Jovid)
library(dplyr)
library(promises)
library(Morpho)
library(Rvcg)
future::plan("multicore")

setwd("~/shiny/shinyapps/Syndrome_model/")
# setwd("/data/Syndrome_model_data/")
# save(atlas, d.meta.combined, front.face, PC.eigenvectors, synd.lm.coefs, synd.mshape, PC.scores, synd.mat, file = "data.Rdata")
load("data.Rdata")
load("modules_400PCs.Rdata")
# load("modules_PCA.Rdata")
eye.index <- as.numeric(read.csv("eye_small.csv", header = F)) +1 # eye.index <- as.numeric(read.csv("~/Desktop/eye_lms.csv", header = F)) +1
load("texture_300PCs.Rdata")
texture.coefs <- lm(texture.pca$x[,1:300] ~ d.meta.combined$Sex + d.meta.combined$Age + d.meta.combined$Age^2 + d.meta.combined$Age^3 + d.meta.combined$Syndrome + d.meta.combined$Age:d.meta.combined$Syndrome)$coef
texture.pcs <-  texture.pca$rotation[,1:300]
texture.mean <- texture.pca$center

tmp.mesh <- atlas
# synd.mshape <- d.registered$mshape

d.meta.combined$Sex <- as.numeric(d.meta.combined$Sex == "F")
d.meta.combined$Syndrome <- factor(d.meta.combined$Syndrome, levels = unique(d.meta.combined$Syndrome))

num_pcs <- 200
PC.eigenvectors <- PC.eigenvectors[,1:num_pcs]
PC.scores <- PC.scores[,1:num_pcs]
meta.lm <- lm(PC.scores[,1:num_pcs] ~ d.meta.combined$Sex + d.meta.combined$Age + d.meta.combined$Age^2 + d.meta.combined$Age^3 + d.meta.combined$Syndrome + d.meta.combined$Age:d.meta.combined$Syndrome)
synd.lm.coefs <- meta.lm$coefficients


predshape.lm <- function(fit, datamod, PC, mshape){
  dims <- dim(mshape)
  mat <- model.matrix(datamod)
  pred <- mat %*% fit

  predPC <- (PC %*% t(pred))
  out <- mshape + matrix(predPC, dims[1], dims[2], byrow = F)

  return(out * 1e10)
}

predPC.lm <- function(fit, datamod){
  mat <- model.matrix(datamod)
  pred <- mat %*% fit
  
  return(pred * 1e10)
}

predtexture.lm <- function(fit, datamod, PC, mshape, gestalt_combo = NULL){
  dims <- dim(mshape)
  mat <- model.matrix(datamod)
  pred <- mat %*% fit
  names <- as.matrix(model.frame(datamod))
  names <- apply(names, 1, paste, collapse = "_")
  names <- gsub(" ", "", names)
  predPC <- t(PC %*% t(pred))
  out <- mshape + predPC
  
  predicted.texture3d <- row2array3d(out)[,,1]
  
  if(is.null(gestalt_combo) == F){
    final.texture <-  2*(predicted.texture3d) + t(col2rgb(atlas$material$color))
  } else {final.texture <- predicted.texture3d}
  
  #scale values
  maxs <- apply(final.texture, 2, max)
  mins <- apply(final.texture, 2, min)
  additive.texture <- scale(final.texture, center = mins, scale = maxs - mins)
  hex.mean <- rgb(additive.texture, maxColorValue = 1)
 
  return(hex.mean)
}

#* @apiTitle Syndrome model API

#* generate atlas prediction
#* @param selected.sex predicted sex effect
#* @param selected.age predicted age effect
#* @param selected.synd predicted syndrome effect
#* @get /predshape
function(selected.sex = "Female", selected.age = 12, selected.synd = "Achondroplasia") {
  selected.synd <- factor(selected.synd, levels = levels(d.meta.combined$Syndrome))
  if(selected.sex == "Female"){selected.sex <-1
  } else if(selected.sex == "Male"){selected.sex <- 0} 
  selected.age <- as.numeric(selected.age)
  
  datamod <- ~ selected.sex + selected.age + selected.age^2 + selected.age^3 + selected.synd + selected.age:selected.synd
 future_promise({
  predicted.shape <- predshape.lm(synd.lm.coefs, datamod, PC.eigenvectors, synd.mshape)
 })
  
}

#* generate atlas prediction
#* @param selected.sex predicted sex effect
#* @param selected.age predicted age effect
#* @param selected.synd predicted syndrome effect
#* @get /predPC
function(selected.sex = "Female", selected.age = 12, selected.synd = "Achondroplasia") {
  # selected.synd <- factor(selected.synd, levels = levels(d.meta.combined$Syndrome))
  if(selected.sex == "Female"){selected.sex <-1
  } else if(selected.sex == "Male"){selected.sex <- 0} 
  selected.age <- as.numeric(selected.age)
  
  predicted.shape <- matrix(NA, nrow = length(unique(d.meta.combined$Syndrome)), ncol = 2)
  future_promise({
    for(i in 1:nrow(predicted.shape)){
      selected.synd <- factor(levels(d.meta.combined$Syndrome)[i], levels = levels(d.meta.combined$Syndrome))
      datamod <- ~ selected.sex + selected.age + selected.age^2 + selected.age^3 + selected.synd + selected.age:selected.synd
    predicted.shape[i,] <- predPC.lm(synd.lm.coefs, datamod)[1:2]
    }
    predicted.shape
  })
  
}


#* generate atlas prediction
#* @param selected.sex predicted sex effect
#* @param selected.age predicted age effect
#* @param selected.synd predicted syndrome effect
#* @get /predshape_mesh
function(selected.sex = "Female", selected.age = 12, selected.synd = "Achondroplasia", res) {
  selected.synd <- factor(selected.synd, levels = levels(d.meta.combined$Syndrome))
  if(selected.sex == "Female"){selected.sex <-1
  } else if(selected.sex == "Male"){selected.sex <- 0} 
  selected.age <- as.numeric(selected.age)
  
  datamod <- ~ selected.sex + selected.age + selected.age^2 + selected.age^3 + selected.synd + selected.age:selected.synd
  print(selected.synd)
  print(selected.sex)
  
  future_promise({
  predicted.shape <- predshape.lm(synd.lm.coefs, datamod, PC.eigenvectors, synd.mshape)
  
  tmp.mesh <- atlas
  tmp.mesh$vb[-4,] <- t(predicted.shape)
  
  Rvcg::vcgSmooth(tmp.mesh)

  })
  
}

#* comparison plot
#* @param comp_age age for syndrome comparison
#* @param comp_sex sex for syndrome comparison
#* @param reference reference syndrome
#* @param synd_comp compared syndrome
#* @param facial_subset what part of the face to analyze
#* @get /synd_comp
#pick two synds, construct mesh estimates, get subsets, register mesh subsets, project and score
function(comp_age = 12, comp_sex = "Female", reference = "Unaffected Unrelated", synd_comp = "Costello Syndrome", facial_subset = ""){

  #synd
  selected.synd <- factor(synd_comp, levels = levels(d.meta.combined$Syndrome))
  if(comp_sex == "Female"){selected.sex <-1
  } else if(comp_sex == "Male"){selected.sex <- 0} 
  selected.age <- as.numeric(comp_age)
  
  future({
  datamod <- ~ selected.sex + selected.age + selected.age^2 + selected.age^3 + selected.synd + selected.age:selected.synd
  synd.pred <- predshape.lm(synd.lm.coefs, datamod, PC.eigenvectors, synd.mshape)
  
  #ref
  selected.synd <- factor(reference, levels = levels(d.meta.combined$Syndrome))
  datamod <- ~ selected.sex + selected.age + selected.age^2 + selected.age^3 + selected.synd + selected.age:selected.synd
  ref.pred <- predshape.lm(synd.lm.coefs, datamod, PC.eigenvectors, synd.mshape) 
  
  #calculate syndrome severity scores for selected syndrome
  #calculate score for the whole face
  S <- synd.lm.coefs[grepl(pattern = synd_comp, rownames(synd.lm.coefs)),][1,]
  Snorm <- S/sqrt(sum(S^2))
  syndscores.main <- PC.scores %*% Snorm
  
  syndscores.df <- data.frame(Syndrome = d.meta.combined$Syndrome, face.score = syndscores.main)
  syndscores.wholeface <- syndscores.df%>%
    group_by(Syndrome) %>%
    summarise(face_score = mean(face.score))
  
  # calculate score for the selected subregion
  if(is.null(facial_subset)) selected.node <- 1 else if(facial_subset == ""){
    selected.node <- 1} else{
      selected.node <- as.numeric(facial_subset)}
  
  if(selected.node > 1){
    subregion.coefs <- manova(get(paste0(tolower(colnames(modules)[selected.node]), ".pca"))$x ~ d.meta.combined$Sex + d.meta.combined$Age + d.meta.combined$Age^2 + d.meta.combined$Age^3 + d.meta.combined$Syndrome + d.meta.combined$Age:d.meta.combined$Syndrome)$coef
    S <- subregion.coefs[grepl(pattern = synd_comp, rownames(subregion.coefs)),][1,]
    print(rownames(subregion.coefs)[grepl(pattern = synd_comp, rownames(subregion.coefs))])
    Snorm <- S/sqrt(sum(S^2))
    syndscores.main <- get(paste0(tolower(colnames(modules)[selected.node]), ".pca"))$x %*% Snorm
    
    syndscores.df <- data.frame(Syndrome = d.meta.combined$Syndrome, module.score = syndscores.main)
    syndscores.module <- syndscores.df%>%
      group_by(Syndrome) %>%
      summarise(face_score = mean(module.score))
  } else{syndscores.module <- syndscores.wholeface}
  
  list(reference_pred = ref.pred, syndrome_pred = synd.pred, wholeface_scores = syndscores.wholeface, subregion_scores = syndscores.module)
  
  }) #end future
  
  
}

#* generate atlas prediction for texture (vertex colors)
#* @param selected.sex predicted sex effect
#* @param selected.age predicted age effect
#* @param selected.synd predicted syndrome effect
#* @get /predtexture
function(selected.sex = "Female", selected.age = 12, selected.synd = "Achondroplasia", gestalt_combo = NULL) {
  selected.synd <- factor(selected.synd, levels = levels(d.meta.combined$Syndrome))
  if(selected.sex == "Female"){selected.sex <-1
  } else if(selected.sex == "Male"){selected.sex <- 0} 
  selected.age <- as.numeric(selected.age)
  
  datamod <- ~ selected.sex + selected.age + selected.age^2 + selected.age^3 + selected.synd + selected.age:selected.synd
  
  
  future_promise({
  predicted.texture <- predtexture.lm(texture.coefs[,1:300], datamod, texture.pcs, texture.mean, gestalt_combo = gestalt_combo)
  })
  
}


#* generate atlas morphtarget
#* @param selected.sex predicted sex effect
#* @param selected.synd predicted syndrome effect
#* @param selected.severity one of Mild, Typical, or Severe
#* @param min_age morphtarget min
#* @param max_age morphtarget max
#* @param severity_sd standard deviation of the severity scores
#* @get /gestalt_morphtarget
function(selected.sex = "Female", selected.synd = "Unaffected Unrelated", selected.severity = "Typical", min_age = 1, max_age = 20, severity_sd = .02, selected.color = "Generic") {
  selected.synd <- factor(selected.synd, levels = levels(d.meta.combined$Syndrome))
  if(selected.sex == "Female"){selected.sex <- 1.5
  } else if(selected.sex == "Male"){selected.sex <- -.5}
  
  #severity math####
  S <- matrix(synd.lm.coefs[grepl(pattern = selected.synd, rownames(synd.lm.coefs)),], nrow = 1, ncol = ncol(PC.eigenvectors))
  Snorm <- S/sqrt(sum(S^2))
  
  minage <- min_age
  nframes <- max_age
  
  if(selected.severity == "Mild"){selected.severity <- -1.5 * severity_sd} else if(selected.severity == "Severe"){selected.severity <- 1.5 * severity_sd} else if(selected.severity == "Typical"){selected.severity <- 0}
future_promise({
  values.shape <- matrix(NA, ncol = 166131 * 3, nrow = 2)
  values_col <- matrix(NA, ncol = 166131 * 3, nrow = 2)
  for(i in 1:nrow(values.shape)){
    selected.age <- as.numeric(c(min_age, max_age)[i])
    
    main.res <- 1e10 * matrix(t(PC.eigenvectors %*% t(selected.severity * Snorm)), dim(synd.mshape)[1], dim(synd.mshape)[2])
    
    datamod <- ~ selected.sex + selected.age + selected.age^2 + selected.age^3 + selected.synd + selected.age:selected.synd
    predicted.shape <- predshape.lm(synd.lm.coefs, datamod, PC.eigenvectors, synd.mshape)
    
    tmp.mesh$vb[-4,] <- t(predicted.shape + main.res)
    final.shape <- vcgSmooth(tmp.mesh)
    
    shape.tmp <- array(t(final.shape$vb[-4,])[atlas$it, ], dim = c(166131, 3, 1))
    values.shape[i,] <- geomorph::two.d.array(shape.tmp)
    
    if(selected.color == "Generic") col.tmp <- array(t(col2rgb(atlas$material$color[atlas$it])), dim = c(166131, 3, 1))/255
    if(selected.color == "Lightgrey") col.tmp <- array(211/255, dim = c(166131, 3, 1))
    if(selected.color == "Generic + Gestalt") col.tmp <- array(predtexture.lm(texture.coefs, datamod, texture.pcs, texture.mean, gestalt_combo = T)[atlas$it], dim = c(166131, 3, 1))
    values_col[i,] <- geomorph::two.d.array(col.tmp)
  }
  
  combined_values <- rbind(as.numeric(t(cbind(values.shape[1,], values_col[1,]))), as.numeric(t(cbind(values.shape[2,], values_col[2,]))))
  return(combined_values)
})
  
}

#* comparison plot
#* @param comp_age age for syndrome comparison
#* @param comp_sex sex for syndrome comparison
#* @param reference reference syndrome
#* @param synd_comp compared syndrome
#* @param facial_subset what part of the face to analyze
#* @get /comparison_morphtarget
function(selected.sex = "Female", selected.synd = "Unaffected Unrelated", synd_comp = "Achondroplasia", selected.severity = "Typical", min_age = 1, max_age = 20, severity_sd = .02, facial_subregion = 1) {
  selected.synd <- factor(selected.synd, levels = levels(d.meta.combined$Syndrome))
  synd_comp <- factor(synd_comp, levels = levels(d.meta.combined$Syndrome))
  if(selected.sex == "Female"){selected.sex <- 1.5
  } else if(selected.sex == "Male"){selected.sex <- -.5}
  
  facial_subregion <- as.numeric(facial_subregion)
  #severity math####
  S <- matrix(synd.lm.coefs[grepl(pattern = selected.synd, rownames(synd.lm.coefs)),], nrow = 1, ncol = ncol(PC.eigenvectors))
  Snorm <- S/sqrt(sum(S^2))
  
  minage <- min_age
  nframes <- max_age
  
  if(selected.severity == "Mild"){selected.severity <- -1.5 * severity_sd} else if(selected.severity == "Severe"){selected.severity <- 1.5 * severity_sd} else if(selected.severity == "Typical"){selected.severity <- 0}
  
  future_promise({
  values <- matrix(NA, ncol = 166131 * 3, nrow = 3)
  values_col <- matrix(NA, ncol = 166131 * 3, nrow = 3)
  
  for(i in 1:nrow(values)){
    selected.age <- as.numeric(c(min_age, mean(c(as.numeric(min_age), as.numeric(max_age))), max_age)[i])
    
    main.res <- 1e10 * matrix(t(PC.eigenvectors %*% t(selected.severity * Snorm)), dim(synd.mshape)[1], dim(synd.mshape)[2])
    
    datamod <- ~ selected.sex + selected.age + selected.age^2 + selected.age^3 + selected.synd + selected.age:selected.synd
    predicted.shape <- predshape.lm(synd.lm.coefs, datamod, PC.eigenvectors, synd.mshape)
    
    tmp.mesh$vb[-4,] <- t(predicted.shape + main.res)
    final.shape <- vcgSmooth(tmp.mesh)
    
    shape.tmp <- array(t(final.shape$vb[-4,])[atlas$it, ], dim = c(166131, 3, 1))
    values[i,] <- geomorph::two.d.array(shape.tmp)
    
    datamod_comp <- ~ selected.sex + selected.age + selected.age^2 + selected.age^3 + synd_comp + selected.age:synd_comp
    predicted.shape <- predshape.lm(synd.lm.coefs, datamod_comp, PC.eigenvectors, synd.mshape)
    
    tmp.mesh$vb[-4,] <- t(predicted.shape + main.res)
    final.shape2 <- vcgSmooth(tmp.mesh)
    
    #if a facial subregion is selected, let's only color the region and leave everything else grey
    if(facial_subregion == 1){
    col.tmp <- array(t(col2rgb(meshDist(final.shape, final.shape2, plot = F)$cols[atlas$it])), dim = c(166131, 3, 1))/255
    } else{
      node.code <- c("posterior_mandible" = 2, "nose" = 3,"anterior_mandible" = 4, "brow" = 5, "zygomatic" = 6, "premaxilla" = 7)
      tmp.meshdist <- meshDist(final.shape, final.shape2, plot = F)$cols
      tmp.meshdist[modules[,names(node.code)[node.code == facial_subregion]] == F] <- "#d3d3d3"
      
      col.tmp <- array(t(col2rgb(tmp.meshdist[atlas$it])), dim = c(166131, 3, 1))/255
    }
    values_col[i,] <- geomorph::two.d.array(col.tmp)
  }
  
  combined_values <- rbind(as.numeric(t(cbind(values[1,], values_col[1,]))), as.numeric(t(cbind(values[2,], values_col[2,]))), as.numeric(t(cbind(values[3,], values_col[3,]))))
  return(combined_values)
  })
  
  
}


#* get similarity scores for whole face and selected subregion
#* @param reference reference syndrome
#* @param synd_comp compared syndrome
#* @param facial_subset what part of the face to analyze
#* @get /similarity_scores
function(reference = "Unaffected Unrelated", synd_comp = "Costello Syndrome", facial_subregion = 1){
  
  selected.synd <- factor(synd_comp, levels = levels(d.meta.combined$Syndrome))

  # future_promise({

    #calculate syndrome severity scores for selected syndrome
    #calculate score for the whole face
    S <- synd.lm.coefs[grepl(pattern = synd_comp, rownames(synd.lm.coefs)),][1,]
    Snorm <- S/sqrt(sum(S^2))
    syndscores.main <- PC.scores %*% Snorm
    
    syndscores.df <- data.frame(Syndrome = d.meta.combined$Syndrome, face.score = syndscores.main)
    syndscores.wholeface <- syndscores.df%>%
      group_by(Syndrome) %>%
      summarise(face_score = mean(face.score))
    
    # calculate score for the selected subregion
    if(is.null(facial_subregion)) selected.node <- 1 else if(facial_subregion == 1){
      selected.node <- 1} else{
        selected.node <- as.numeric(facial_subregion)}
    
    if(selected.node > 1){
      node.code <- c("posterior_mandible" = 2, "nose" = 3,"anterior_mandible" = 4, "brow" = 5, "zygomatic" = 6, "premaxilla" = 7)
      subregion.coefs <- manova(get(paste0(tolower(names(node.code)[node.code == selected.node]), ".pca"))$x ~ d.meta.combined$Sex + d.meta.combined$Age + d.meta.combined$Age^2 + d.meta.combined$Age^3 + d.meta.combined$Syndrome + d.meta.combined$Age:d.meta.combined$Syndrome)$coef
      S <- subregion.coefs[grepl(pattern = synd_comp, rownames(subregion.coefs)),][1,]
      
      Snorm <- S/sqrt(sum(S^2))
      syndscores.main <- get(paste0(tolower(names(node.code)[node.code == selected.node]), ".pca"))$x %*% Snorm
      
      syndscores.df <- data.frame(Syndrome = d.meta.combined$Syndrome, module.score = syndscores.main)
      syndscores.module <- syndscores.df%>%
        group_by(Syndrome) %>%
        summarise(face_score = mean(module.score))
    } else{syndscores.module <- syndscores.wholeface}
    
    list(wholeface_scores = syndscores.wholeface, subregion_scores = syndscores.module)
    
  # }) #end future
  
}