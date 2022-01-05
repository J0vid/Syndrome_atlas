
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
atlas <- file2mesh("~/whoami3_2k.ply")

selected.synd <- factor("Achondroplasia", levels = levels(d.meta.combined$Syndrome))
selected.sex <-1
selected.age <- 1
datamod <- ~ selected.sex + selected.age + selected.age^2 + selected.age^3 + selected.synd + selected.age:selected.synd
predicted.shape <- predshape.lm(synd.lm.coefs, datamod, PC.eigenvectors, synd.mshape)/1e10

sample2k <- sample(1:27903, 1982)

atlas$vb[-4,] <- t(predicted.shape[sample2k,])

open3d()
shade3d(vcgSmooth(atlas), aspect = "iso", col = "lightgrey", specular = 1)
objid <- ids3d()$id
xyz <- rgl.attrib(objid[1], "vertices")
dim(xyz)
#calc age vector values

selected.synd <- "Achondroplasia"
selected.sex <- "Female"
selected.age <- i

raw_api_res <- httr::GET(url = paste0("http://localhost:6352", "/predPC"),
                         query = list(selected.sex = selected.sex, selected.age = selected.age, selected.synd = selected.synd),
                         encode = "json")

parsed_resp  <- jsonlite::fromJSON(httr::content(raw_api_res, "text"))/1e10

nframes <- 50

values <- matrix(NA, ncol = nrow(xyz) * 3, nrow = nframes)
for(i in 1:nframes){
  
  predicted.shape <- synd.mshape + matrix(PC.eigenvectors %*% (parsed_resp[i,]), nrow = 27903, ncol = 3, byrow = F) 
  
  # tmp.mesh$vb[-4,] <- t(predicted.shape)
  # vertices <- t(asEuclidean2(tmp.mesh$vb))
  
  shape.tmp <- array(predicted.shape[sample2k,][as.numeric(atlas$it), ], dim = c(nrow(xyz), 3, 1))
  values[i,] <- geomorph::two.d.array(shape.tmp)
  print(i)
}

widget <- rglwidget(width = 500, height = 500) %>%
  playwidget(vertexControl(values = values[,], 
                           vertices = rep(1:nrow(xyz), each = 3),
                           attributes = rep(c("x", "y", "z"), nrow(xyz)),
                           objid = objid), start = 1, stop = 20, rate = .5)

widget

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

d.registered <- procSym(d.lm.sparse)

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

nframes <- 50

values <- matrix(NA, ncol = nrow(xyz) * 3, nrow = nframes)
for(i in 1:nframes){
  selected.age <- i
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

control <- vertexControl(values = values, 
              vertices = rep(1:nrow(xyz), each = 3),
              attributes = rep(c("x", "y", "z"), nrow(xyz)),
              objid = objid)

control2 <- vertexControl(values = values2, 
                          vertices = rep(1:nrow(xyz), each = 3),
                          attributes = rep(c("x", "y", "z"), nrow(xyz)),
                          objid = objid,
                          param = 1:nrow(values2))

rglwidget(width = 500, height = 500) %>%
  playwidget(control2)


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





