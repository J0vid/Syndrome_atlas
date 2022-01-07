library(plumber)
library(future)
library(Jovid)
library(dplyr)
library(Rvcg)
library(promises)
future::plan("multisession")

setwd("~/shiny/shinyapps/Syndrome_model/")
# setwd("/srv/shiny-server/testing_ground/")
# save(atlas, d.meta.combined, front.face, PC.eigenvectors, synd.lm.coefs, synd.mshape, PC.scores, synd.mat, file = "data.Rdata")
load("2k_topo_objects.Rdata")
# load("modules_PCA.Rdata")
# eye.index <- as.numeric(read.csv("~/Desktop/eye_lms.csv", header = F)) +1
# load("~/shiny/shinyapps/Syndrome_model/FB2_texture_PCA.Rdata")
# texture.coefs <- lm(texture.pca$x[,1:300] ~ d.meta.combined$Sex + d.meta.combined$Age + d.meta.combined$Age^2 + d.meta.combined$Age^3 + d.meta.combined$Syndrome + d.meta.combined$Age:d.meta.combined$Syndrome)$coef
# texture.pcs <-  texture.pca$rotation[,1:300]
# texture.mean <- texture.pca$center

atlas <- new_atlas
tmp.mesh <- atlas
synd.mshape <- d.registered$mshape
PC.eigenvectors <- d.registered$PCs[,1:200]
d.meta.combined$Sex <- as.numeric(d.meta.combined$Sex == "F")
d.meta.combined$Syndrome <- factor(d.meta.combined$Syndrome, levels = unique(d.meta.combined$Syndrome))

num_pcs <- 200
meta.lm <- lm(d.registered$PCscores[,1:num_pcs] ~ d.meta.combined$Sex + d.meta.combined$Age + d.meta.combined$Age^2 + d.meta.combined$Age^3 + d.meta.combined$Syndrome + d.meta.combined$Age:d.meta.combined$Syndrome)
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

# predtexture.lm <- function(fit, datamod, PC, mshape, gestalt_combo = NULL){
#   dims <- dim(mshape)
#   mat <- model.matrix(datamod)
#   pred <- mat %*% fit
#   names <- as.matrix(model.frame(datamod))
#   names <- apply(names, 1, paste, collapse = "_")
#   names <- gsub(" ", "", names)
#   predPC <- t(PC %*% t(pred))
#   out <- mshape + predPC
#   
#   predicted.texture3d <- row2array3d(out)[,,1]
#   
#   if(is.null(gestalt_combo) == F){
#     final.texture <-  2*(predicted.texture3d) + t(col2rgb(atlas$material$color))
#   } else {final.texture <- predicted.texture3d}
#   
#   
#   #scale values
#   maxs <- apply(final.texture, 2, max)
#   mins <- apply(final.texture, 2, min)
#   additive.texture <- scale(final.texture, center = mins, scale = maxs - mins)
#   hex.mean <- rgb(additive.texture, maxColorValue = 1)
#   
#   # dimnames(out)[[3]] <- names
#   # rownames(pred) <- names
#   return(hex.mean)
# }


#* @apiTitle Syndrome model API

#* generate atlas prediction
#* @param selected.sex predicted sex effect
#* @param selected.age predicted age effect
#* @param selected.synd predicted syndrome effect
#* @get /predshape
function(selected.sex = "Female", selected.synd = "Achondroplasia", selected.severity = "typical", min_age = 1, max_age = 20, severity_sd = .02) {
  selected.synd <- factor(selected.synd, levels = levels(d.meta.combined$Syndrome))
  if(selected.sex == "Female"){selected.sex <- 2
  } else if(selected.sex == "Male"){selected.sex <- -1}
  
  #severity math####
  S <- matrix(synd.lm.coefs[grepl(pattern = selected.synd, rownames(synd.lm.coefs)),], nrow = 1, ncol = num_pcs)
  Snorm <- S/sqrt(sum(S^2))
  
  minage <- min_age
  nframes <- max_age
  
  if(selected.severity == "mild"){selected.severity <- -2 * severity_sd} else if(selected.severity == "severe"){selected.severity <- 2 * severity_sd} else if(selected.severity == "typical"){selected.severity <- 0}
  
 future_promise({
   values <- matrix(NA, ncol = 11628 * 3, nrow = length(minage:nframes))
   
   for(i in 1:nrow(values)){
     selected.synd <- factor(selected.synd, levels = levels(d.meta.combined$Syndrome))
     selected.age <- i
     
     main.res <- 1e10 * matrix(t(PC.eigenvectors %*% t(selected.severity * Snorm)), dim(synd.mshape)[1], dim(synd.mshape)[2])
     
     datamod <- ~ selected.sex + selected.age + selected.age^2 + selected.age^3 + selected.synd + selected.age:selected.synd
     predicted.shape <- predshape.lm(synd.lm.coefs, datamod, PC.eigenvectors, synd.mshape)
     
     tmp.mesh$vb[-4,] <- t(predicted.shape + main.res)
     final.shape <- t(vcgSmooth(tmp.mesh)$vb[-4,])
     
     shape.tmp <- array(final.shape[as.numeric(atlas$it), ], dim = c(11628, 3, 1))
     values[i,] <- geomorph::two.d.array(shape.tmp)
   }
   values
 })
  
}

#* generate atlas prediction
#* @param selected.sex predicted sex effect
#* @param selected.age predicted age effect
#* @param selected.synd predicted syndrome effect
#* @get /predPC
function(selected.sex = "Female", selected.age = 12, selected.synd = "Achondroplasia") {
  selected.synd <- factor(selected.synd, levels = levels(d.meta.combined$Syndrome))
  if(selected.sex == "Female"){selected.sex <-1
  } else if(selected.sex == "Male"){selected.sex <- 0} 
  selected.age <- as.numeric(selected.age)
  
  predicted.shape <- matrix(NA, nrow = 50, ncol = 200)
  future_promise({
    for(i in 1:50){
      selected.age <- i
      datamod <- ~ selected.sex + selected.age + selected.age^2 + selected.age^3 + selected.synd + selected.age:selected.synd
    predicted.shape[i,] <- predPC.lm(synd.lm.coefs, datamod)
    }
    predicted.shape
  })
  
}


#* generate atlas prediction
#* @param selected.sex predicted sex effect
#* @param selected.age predicted age effect
#* @param selected.synd predicted syndrome effect
#* @serializer htmlwidget
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
  
  plot3d(Rvcg::vcgSmooth(tmp.mesh), aspect = "iso", col = "lightgrey", specular = 1)
  rglwidget()

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







