
#duncan's ex

open3d()
shade3d(cube3d(), col="red")
objid <- ids3d()$id
xyz <- rgl.attrib(objid, "vertices")

# Move the first 2 vertices

values <- matrix(NA, 100, 6)
theta <- seq(0, 2*pi, len = 100)
for (i in 1:100) {
  xyz1 <- rotate3d(xyz[1:2,], x = 0, y = 0, z = 1, angle = theta[i])
  values[i, 1:3] <- xyz1[1,]
  values[i, 4:6] <- xyz1[2,]
}



widget <- rglwidget(width = 500, height = 500) %>%
  playwidget(vertexControl(values = values, 
                           vertices = rep(1:2, each=3),
                           attributes = rep(c("x", "y", "z"), 2),
                           objid = objid,
                           param = 1:100))

widget


#do it with mesh####
library(Morpho)
library(rgl)
library(Rvcg)

setwd("~/shiny/shinyapps/Syndrome_model/")
load("data.Rdata")

predshape.lm <- function(fit, datamod, PC, mshape){
  dims <- dim(mshape)
  mat <- model.matrix(datamod)
  pred <- mat %*% fit
  
  predPC <- (PC %*% t(pred))
  out <- mshape + matrix(predPC, dims[1], dims[2], byrow = F)
  
  return(out * 1e10)
}



close3d()

selected.synd <- factor("Achondroplasia", levels = levels(d.meta.combined$Syndrome))
selected.sex <-1
selected.age <- 1
datamod <- ~ selected.sex + selected.age + selected.age^2 + selected.age^3 + selected.synd + selected.age:selected.synd
predicted.shape <- predshape.lm(synd.lm.coefs, datamod, PC.eigenvectors, synd.mshape)/1e10

atlas$vb[-4,] <- t(predicted.shape)

open3d()
shade3d(vcgSmooth(atlas), aspect = "iso", col = "lightgrey", specular = 1)
objid <- ids3d()$id
xyz <- rgl.attrib(objid[1], "vertices")
dim(xyz)
#calc age vector values

selected.synd <- "Achondroplasia"
selected.sex <- "Female"
selected.age <- 1

raw_api_res <- httr::GET(url = paste0("http://localhost:6352", "/predshape"),
                         query = list(selected.sex = selected.sex, selected.age = selected.age, selected.synd = selected.synd),
                         encode = "json")

parsed_resp  <- jsonlite::fromJSON(httr::content(raw_api_res, "text"))/1e10

selected.age <- 58

raw_api_res <- httr::GET(url = paste0("http://localhost:6352", "/predshape"),
                         query = list(selected.sex = selected.sex, selected.age = selected.age, selected.synd = selected.synd),
                         encode = "json")

parsed_resp2  <- jsonlite::fromJSON(httr::content(raw_api_res, "text"))/1e10

values <- two.d.array(abind::abind(parsed_resp[as.numeric(atlas$it),], parsed_resp2[as.numeric(atlas$it),], along = 3))

rglwidget(width = 500, height = 500) %>%
  playwidget(vertexControl(values = values, 
                           vertices = rep(1:nrow(xyz), each = 3),
                           attributes = rep(c("x", "y", "z"), nrow(xyz)),
                           objid = objid), step = 1/58)



#try out new topology####
library(Morpho)
library(rgl)
library(geomorph)
library(Rvcg)

predshape.lm <- function(fit, datamod, PC, mshape){
  dims <- dim(mshape)
  mat <- model.matrix(datamod)
  pred <- mat %*% fit
  
  predPC <- (PC %*% t(pred))
  out <- mshape + matrix(predPC, dims[1], dims[2], byrow = F)
  
  return(out)
}

load("~/Downloads/2k_topo_objects.Rdata")

atlas <- new_atlas
#only processed through first 425 inds so far
# d.lm.sparse <- d.lm.sparse[,,1:425]
# d.meta.combined <- d.meta.combined[1:425,]
d.meta.combined$Sex <- as.numeric(d.meta.combined$Sex == "F")
d.meta.combined$Syndrome <- factor(d.meta.combined$Syndrome, levels = unique(d.meta.combined$Syndrome))

# d.registered <- procSym(d.lm.sparse)

num_pcs <- 200
meta.lm <- lm(d.registered$PCscores[,1:num_pcs] ~ d.meta.combined$Sex + d.meta.combined$Age + d.meta.combined$Age^2 + d.meta.combined$Age^3 + d.meta.combined$Syndrome + d.meta.combined$Age:d.meta.combined$Syndrome)
synd.lm.coefs <- meta.lm$coefficients

selected.synd <- factor("Treacher Collins Syndrome", levels = levels(d.meta.combined$Syndrome))
selected.sex <-1
selected.age <- 1
datamod <- ~ selected.sex + selected.age + selected.age^2 + selected.age^3 + selected.synd + selected.age:selected.synd
predicted.shape <- predshape.lm(synd.lm.coefs, datamod, d.registered$PCs[,1:num_pcs], d.registered$mshape)

atlas$vb[-4,] <- t(predicted.shape)

open3d()
shade3d(vcgSmooth(atlas), aspect = "iso", col = "lightgrey", specular = 1)
objid <- ids3d()$id
xyz <- rgl.attrib(objid[1], "vertices")
dim(xyz)

age_vec <- c(1,50)

values <- matrix(NA, ncol = nrow(xyz) * 3, nrow = 2)
for(i in 1:2){
  selected.age <- age_vec[i]
  selected.sex <- -2
  datamod <- ~ selected.sex + selected.age + selected.age^2 + selected.age^3 + selected.synd + selected.age:selected.synd
  predicted.shape <- predshape.lm(synd.lm.coefs, datamod, d.registered$PCs[,1:num_pcs], d.registered$mshape)
  
  shape.tmp <- array(predicted.shape[as.numeric(atlas$it), ], dim = c(nrow(xyz), 3, 1))
  values[i,] <- geomorph::two.d.array(shape.tmp)
}

#severity morph target####
S <- matrix(synd.lm.coefs[grepl(pattern = "Treacher Collins Syndrome", rownames(synd.lm.coefs)),], nrow = 1, ncol = num_pcs)
Snorm <- S/sqrt(sum(S^2))

syndscores.main <- d.registered$PCscores[,1:num_pcs] %*% t(Snorm)
severity_range <- seq(round(-1*sd(syndscores.main), 2), round(1*sd(syndscores.main), 2), length.out = 10)

selected.synd <- factor("Treacher Collins Syndrome", levels = levels(d.meta.combined$Syndrome))

syndonly_lm <- lm(d.registered$PCscores[,1:num_pcs] ~ d.meta.combined$Syndrome)$coef

datamod <- ~ selected.synd 
predicted.shape <- predshape.lm(syndonly_lm, datamod, d.registered$PCs[,1:num_pcs], d.registered$mshape)

values2 <- matrix(NA, ncol = nrow(xyz) * 3, nrow = 10)

for(i in 1:nrow(values2)){
  
  selected.severity <- severity_range[i] #.01
  
  main.res <- matrix(t(d.registered$PCs[,1:num_pcs] %*% t(selected.severity * Snorm)), dim(d.registered$mshape)[1], dim(d.registered$mshape)[2])
  
  #add severity residuals
  severity_shape <-  predicted.shape + main.res
  
  shape.tmp <- array(severity_shape[as.numeric(atlas$it), ], dim = c(nrow(xyz), 3, 1))
  values2[i,] <- geomorph::two.d.array(shape.tmp)
  
}

control <- vertexControl(values = values[,], 
              vertices = rep(1:nrow(xyz), each = 3),
              attributes = rep(c("x", "y", "z"), nrow(xyz)),
              objid = objid)

control2 <- vertexControl(values = values2, 
                          vertices = rep(1:nrow(xyz), each = 3),
                          attributes = rep(c("x", "y", "z"), nrow(xyz)),
                          objid = objid,
                          param = 1:nrow(values2))

rglwidget(width = 500, height = 500) %>%
  playwidget(control, step = .001)


#combined controller hack####


open3d()
shade3d(vcgSmooth(atlas), aspect = "iso", col = "lightgrey", specular = 1)
objid <- ids3d()$id
xyz <- rgl.attrib(objid[1], "vertices")
dim(xyz)

nframes <- 50

values3 <- matrix(NA, ncol = nrow(xyz) * 3, nrow = 10 * nframes)
for(i in 0:nframes){
  selected.age <- i
  selected.sex <- -2
  datamod <- ~ selected.sex + selected.age + selected.age^2 + selected.age^3 + selected.synd + selected.age:selected.synd
  predicted.shape <- predshape.lm(synd.lm.coefs, datamod, d.registered$PCs[,1:num_pcs], d.registered$mshape)
  
for(j in 1:10){
  
  selected.severity <- severity_range[j] #.01
  
  main.res <- matrix(t(d.registered$PCs[,1:num_pcs] %*% t(selected.severity * Snorm)), dim(d.registered$mshape)[1], dim(d.registered$mshape)[2])
  
  #add severity residuals
  severity_shape <-  predicted.shape + main.res
  
  shape.tmp <- array(severity_shape[as.numeric(atlas$it), ], dim = c(nrow(xyz), 3, 1))
  values3[i*(10) + j,] <- geomorph::two.d.array(shape.tmp)
  
}
}
control3 <- vertexControl(values = na.omit(values3),
                          vertices = rep(1:nrow(xyz), each = 3),
                          attributes = rep(c("x", "y", "z"), nrow(xyz)),
                          objid = objid)

rglwidget(width = 500, height = 500) %>%
  playwidget(control3)


#try out new topology for heatmaps####
library(Morpho)
library(rgl)
library(geomorph)
library(Rvcg)

predshape.lm <- function(fit, datamod, PC, mshape){
  dims <- dim(mshape)
  mat <- model.matrix(datamod)
  pred <- mat %*% fit
  
  predPC <- (PC %*% t(pred))
  out <- mshape + matrix(predPC, dims[1], dims[2], byrow = F)
  
  return(out)
}

load("~/Downloads/2k_topo_objects.Rdata")

atlas <- new_atlas
#only processed through first 425 inds so far
# d.lm.sparse <- d.lm.sparse[,,1:425]
# d.meta.combined <- d.meta.combined[1:425,]
d.meta.combined$Sex <- as.numeric(d.meta.combined$Sex == "F")
d.meta.combined$Syndrome <- factor(d.meta.combined$Syndrome, levels = unique(d.meta.combined$Syndrome))

# d.registered <- procSym(d.lm.sparse)

num_pcs <- 200
meta.lm <- lm(d.registered$PCscores[,1:num_pcs] ~ d.meta.combined$Sex + d.meta.combined$Age + d.meta.combined$Age^2 + d.meta.combined$Age^3 + d.meta.combined$Syndrome + d.meta.combined$Age:d.meta.combined$Syndrome)
synd.lm.coefs <- meta.lm$coefficients

selected.synd <- factor("Achondroplasia", levels = levels(d.meta.combined$Syndrome))
selected.sex <-1
selected.age <- 1
datamod <- ~ selected.sex + selected.age + selected.age^2 + selected.age^3 + selected.synd + selected.age:selected.synd
predicted.shape <- predshape.lm(synd.lm.coefs, datamod, PC.eigenvectors, synd.mshape)

atlas$vb[-4,] <- t(predicted.shape)

open3d()
shade3d(vcgSmooth(atlas), aspect = "iso", col = "lightgrey", specular = 1)
objid <- ids3d()$id
xyz <- rgl.attrib(objid[1], "vertices")
dim(xyz)


#severity morph target####
S <- matrix(synd.lm.coefs[grepl(pattern = "Treacher Collins Syndrome", rownames(synd.lm.coefs)),], nrow = 1, ncol = num_pcs)
Snorm <- S/sqrt(sum(S^2))

syndscores.main <- d.registered$PCscores[,1:num_pcs] %*% t(Snorm)
severity_range <- seq(round(-1*sd(syndscores.main), 2), round(1*sd(syndscores.main), 2), length.out = 10)

selected.synd <- factor("Treacher Collins Syndrome", levels = levels(d.meta.combined$Syndrome))

syndonly_lm <- lm(d.registered$PCscores[,1:num_pcs] ~ d.meta.combined$Syndrome)$coef

selected.synd <- factor(selected.synd, levels = levels(d.meta.combined$Syndrome))
if(selected.sex == "Female"){selected.sex <- 2
} else if(selected.sex == "Male"){selected.sex <- -1}


#example syndrome comparison morphtarget API call####
selected.synd <- factor("Crouzon Syndrome", levels = levels(d.meta.combined$Syndrome))
selected.sex <-1
selected.age <- 1
datamod <- ~ selected.sex + selected.age + selected.age^2 + selected.age^3 + selected.synd + selected.age:selected.synd
predicted.shape <- predshape.lm(synd.lm.coefs, datamod, PC.eigenvectors, synd.mshape)

atlas$vb[-4,] <- t(predicted.shape)

open3d()
shade3d(vcgSmooth(atlas), aspect = "iso", col = md1$cols, specular = 1)
objid <- ids3d()$id
xyz <- rgl.attrib(objid[1], "vertices")
dim(xyz)

api_test <- function(selected.sex = "Female", selected.synd = "Unaffected Unrelated", synd_comp = "Achondroplasia", selected.severity = "typical", min_age = 1, max_age = 20, severity_sd = .02) {
  selected.synd <- factor(selected.synd, levels = levels(d.meta.combined$Syndrome))
  synd_comp <- factor(synd_comp, levels = levels(d.meta.combined$Syndrome))
  if(selected.sex == "Female"){selected.sex <- 1.5
  } else if(selected.sex == "Male"){selected.sex <- -.5}
  
  #severity math####
  S <- matrix(synd.lm.coefs[grepl(pattern = selected.synd, rownames(synd.lm.coefs)),], nrow = 1, ncol = ncol(PC.eigenvectors))
  Snorm <- S/sqrt(sum(S^2))
  
  minage <- min_age
  nframes <- max_age
  
  if(selected.severity == "mild"){selected.severity <- -1.5 * severity_sd} else if(selected.severity == "severe"){selected.severity <- 1.5 * severity_sd} else if(selected.severity == "typical"){selected.severity <- 0}
  

    values <- matrix(NA, ncol = nrow(xyz) * 3, nrow = 2)
    values_col <- matrix(NA, ncol = nrow(xyz) * 3, nrow = 2)
    
    for(i in 1:nrow(values)){
      selected.synd <- factor(selected.synd, levels = levels(d.meta.combined$Syndrome))
      selected.age <- c(min_age, max_age)[i]
      
      main.res <- 1e10 * matrix(t(PC.eigenvectors %*% t(selected.severity * Snorm)), dim(synd.mshape)[1], dim(synd.mshape)[2])
      
      datamod <- ~ selected.sex + selected.age + selected.age^2 + selected.age^3 + selected.synd + selected.age:selected.synd
      predicted.shape <- predshape.lm(synd.lm.coefs, datamod, PC.eigenvectors, synd.mshape)
      
      tmp.mesh$vb[-4,] <- t(predicted.shape + main.res)
      final.shape <- vcgSmooth(tmp.mesh)
      
      shape.tmp <- array(t(final.shape$vb[-4,])[atlas$it, ], dim = c(nrow(xyz), 3, 1))
      values[i,] <- geomorph::two.d.array(shape.tmp)
      
      datamod_comp <- ~ selected.sex + selected.age + selected.age^2 + selected.age^3 + synd_comp + selected.age:synd_comp
      predicted.shape <- predshape.lm(synd.lm.coefs, datamod_comp, PC.eigenvectors, synd.mshape)
      
      tmp.mesh$vb[-4,] <- t(predicted.shape + main.res)
      final.shape2 <- vcgSmooth(tmp.mesh)
      
      col.tmp <- array(t(col2rgb(meshDist(final.shape, final.shape2, plot = F)$cols[atlas$it])), dim = c(nrow(xyz), 3, 1))/255
      values_col[i,] <- geomorph::two.d.array(col.tmp)
    }
    
    combined_values <- rbind(as.numeric(t(cbind(values[1,], values_col[1,]))), as.numeric(t(cbind(values[2,], values_col[2,]))))
    
    return(combined_values)

}

test <- api_test(selected.synd = "Crouzon Syndrome", max_age = 22, selected.severity = "mild")


control_combined <- vertexControl(values = test, 
                         vertices = c(rep(1:nrow(xyz), each = 6)),
                         attributes = c(rep(c("x", "red", "y", "green", "z", "blue"), nrow(xyz))),
                         objid = objid)

rglwidget(width = 500, height = 500) %>%
  playwidget(control_combined, step = .001)


#example gestalt morphtarget API call####
selected.synd <- factor("Crouzon Syndrome", levels = levels(d.meta.combined$Syndrome))
selected.sex <-1
selected.age <- 1
datamod <- ~ selected.sex + selected.age + selected.age^2 + selected.age^3 + selected.synd + selected.age:selected.synd
predicted.shape <- predshape.lm(synd.lm.coefs, datamod, PC.eigenvectors, synd.mshape)

atlas$vb[-4,] <- t(predicted.shape)

open3d()
shade3d(vcgSmooth(atlas), aspect = "iso", specular = 1)
objid <- ids3d()$id
xyz <- rgl.attrib(objid[1], "vertices")
dim(xyz)

gestalt_test <- function(selected.sex = "Female", selected.synd = "Unaffected Unrelated", selected.severity = "typical", min_age = 1, max_age = 20, severity_sd = .02) {
  selected.synd <- factor(selected.synd, levels = levels(d.meta.combined$Syndrome))
  synd_comp <- factor(synd_comp, levels = levels(d.meta.combined$Syndrome))
  if(selected.sex == "Female"){selected.sex <- 1.5
  } else if(selected.sex == "Male"){selected.sex <- -.5}
  
  #severity math####
  S <- matrix(synd.lm.coefs[grepl(pattern = selected.synd, rownames(synd.lm.coefs)),], nrow = 1, ncol = ncol(PC.eigenvectors))
  Snorm <- S/sqrt(sum(S^2))
  
  minage <- min_age
  nframes <- max_age
  
  if(selected.severity == "mild"){selected.severity <- -1.5 * severity_sd} else if(selected.severity == "severe"){selected.severity <- 1.5 * severity_sd} else if(selected.severity == "typical"){selected.severity <- 0}
  
  values <- matrix(NA, ncol = nrow(xyz) * 3, nrow = 2)
  for(i in 1:nrow(values)){
    selected.synd <- factor(selected.synd, levels = levels(d.meta.combined$Syndrome))
    selected.age <- c(min_age, max_age)[i]
    
    main.res <- 1e10 * matrix(t(PC.eigenvectors %*% t(selected.severity * Snorm)), dim(synd.mshape)[1], dim(synd.mshape)[2])
    
    datamod <- ~ selected.sex + selected.age + selected.age^2 + selected.age^3 + selected.synd + selected.age:selected.synd
    predicted.shape <- predshape.lm(synd.lm.coefs, datamod, PC.eigenvectors, synd.mshape)
    
    tmp.mesh$vb[-4,] <- t(predicted.shape + main.res)
    final.shape <- vcgSmooth(tmp.mesh)
    
    shape.tmp <- array(t(final.shape$vb[-4,])[atlas$it, ], dim = c(nrow(xyz), 3, 1))
    values[i,] <- geomorph::two.d.array(shape.tmp)
    
  }
  return(values)
}

test2 <- gestalt_test(selected.synd = "Crouzon Syndrome", max_age = 22, selected.severity = "severe")

control <- vertexControl(values = test2,
                         vertices = rep(1:nrow(xyz), each = 3),
                         attributes = rep(c("x", "y", "z"), nrow(xyz)),
                         objid = objid)

rglwidget(width = 500, height = 500) %>%
  playwidget(control, step = .001)



#


#do it with mesh####
library(Morpho)
library(rgl)
library(Rvcg)

setwd("~/shiny/shinyapps/Syndrome_model/")
load("data.Rdata")

predshape.lm <- function(fit, datamod, PC, mshape){
  dims <- dim(mshape)
  mat <- model.matrix(datamod)
  pred <- mat %*% fit
  
  predPC <- (PC %*% t(pred))
  out <- mshape + matrix(predPC, dims[1], dims[2], byrow = F)
  
  return(out * 1e10)
}



close3d()

selected.synd <- factor("Achondroplasia", levels = levels(d.meta.combined$Syndrome))
selected.sex <-1
selected.age <- 1
datamod <- ~ selected.sex + selected.age + selected.age^2 + selected.age^3 + selected.synd + selected.age:selected.synd
predicted.shape <- predshape.lm(synd.lm.coefs, datamod, PC.eigenvectors, synd.mshape)/1e10

atlas$vb[-4,] <- t(predicted.shape)

open3d()
shade3d(vcgSmooth(atlas), aspect = "iso", col = "lightgrey", specular = 1)
objid <- ids3d()$id
xyz <- rgl.attrib(objid[1], "vertices")
dim(xyz)
#calc age vector values

age.vec <- c(1, 35)
values3d <- array(NA, dim = c(nrow(xyz), 3, 2))
for(i in 1:2){
  selected.age <- age.vec[i]
  
  datamod <- ~ selected.sex + selected.age + selected.age^2 + selected.age^3 + selected.synd + selected.age:selected.synd
  values3d[,,i] <- (predshape.lm(synd.lm.coefs, datamod, PC.eigenvectors, synd.mshape)/1e10)[as.numeric(atlas$it),]
}

num <- which(duplicated(as.numeric(atlas$it)) == F)
values <- two.d.array(values3d[num,,])

rglwidget(width = 500, height = 500) %>%
  playwidget(vertexControl(values = values, 
                           vertices = rep(num, each = 3),
                           attributes = rep(c("x", "y", "z"), length(num)),
                           objid = objid,
                           param = 1:2), step = 1/58)





