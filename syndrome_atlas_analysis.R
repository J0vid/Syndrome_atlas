library(Morpho)

setwd("/mnt/Hallgrimsson/Users/Jovid/Andlit/")

load("Dense_combined_starter_PC_outliers_duplicates_removed.Rdata")

r.registered <- procSym(d.lms.combined)

#what do I need to keep for the synd app?

PC.eigenvectors <- r.registered$PCs
PC.scores <- r.registered$PCscores
synd.mshape <- r.registered$mshape
synd.names <- as.factor(d.meta.combined$Syndrome)

# save(PC.eigenvectors, PC.scores, synd.mshape, synd.names, file = "syndromic_shiny.Rdata")

#test out scoring
model.data <- r.registered$rotated[,,]
model.mshape <- r.registered$mshape[,]
residuals <- sweep(model.data, c(1,2), model.mshape, "-")
residuals2d <- two.d.array(residuals)

mf.mod <- manova(PC.scores ~ d.meta.combined$Syndrome == "Achondroplasia")
# mf.mod <- manova(residuals2d ~ ns.cov$Sex[ns.sample], ns.cov$Age + poly(ns.cov$Age, 2))

#score individuals
S <- mf.mod$coefficients[2,]
Snorm <- S/sqrt(sum(S^2))
syndscores <- PC.scores %*% Snorm

scale10 <- function(x, max.range = 10) ((x-min(x))/(max(x)-min(x))) * max.range

boxplot((syndscores)  ~ d.meta.combined$Syndrome == "Achondroplasia")

#figures setup####

library(dsai)
library(Morpho)
library(dplyr)
library(rgl)
# setwd("/mnt/Hallgrimsson/Users/Jovid/Andlit/")
setwd("~/shiny/shinyapps/Syndrome_model/")
# save(atlas, d.meta.combined, front.face, PC.eigenvectors, PC.scores, synd.mshape, file = "syndromic_shiny.Rdata")
load("syndromic_shiny.Rdata")
# atlas <- file2mesh("hailey_reg.ply")
atlas <- file2mesh("Full.ply")
#supp table 1, demographics####
#for each syndrome, what's the age range, sex dist and severity distribution####

#tidyverse count sex
#fix the two incorrect Turner Synd M
d.meta.combined$Sex[d.meta.combined$Syndrome == "Turner Syndrome"] <- "F"
# make the reference class non-syndromic
# d.meta.combined$Syndrome <- relevel(d.meta.combined$Syndrome, ref = "Unaffected Unrelated")
d.meta.combined %>% 
  count(Sex, Syndrome, sort = TRUE) %>%
  View()

#get ranges for age by syndrome
d.meta.combined %>%
  group_by(Syndrome) %>%
  summarize(min_age = min(Age), max_age = max(Age))

#get scores, min and max per syndrome
d.meta.combined2 <- cbind(d.meta.combined, PC.scores)
d.meta.combined2$Syndrome <- relevel(d.meta.combined2$Syndrome, ref = "Unaffected Unrelated")
sex.scores <- rep(NA, nrow(d.meta.combined))
man.coefs <- manova(as.matrix(d.meta.combined2[,grep(colnames(d.meta.combined2), pattern = "PC")]) ~ d.meta.combined2$Sex + d.meta.combined2$Age + d.meta.combined2$Age^2 + d.meta.combined2$Age^3 + d.meta.combined$Syndrome + d.meta.combined2$Age:d.meta.combined2$Syndrome)$coef
rownames(man.coefs)[1] <- "Unaffected Unrelated"
for(i in unique(d.meta.combined$Syndrome)[]){
  #calculate score for the main effect
  S <- (as.numeric(man.coefs[grepl(pattern = i, rownames(man.coefs)),]))[1:500]
  Snorm <- S/sqrt(sum(S^2))
  sex.scores[d.meta.combined$Syndrome == i] <- (PC.scores %*% Snorm)[d.meta.combined$Syndrome == i]
}

d.meta.combined2 <- data.frame(d.meta.combined, sex.scores)
d.meta.combined2 %>%
  group_by(Syndrome) %>%
  summarise(min_score = min(sex.scores), max_score = max(sex.scores)) %>%
  View()


#fig1 Visualizations of syndromic morphology####
#fig 1A, range of syndromic morphologies for an example syndrome: Achondroplasia####
selected.sex <- .5
selected.age <- mean(d.meta.combined$Age[d.meta.combined$Syndrome == "Achondroplasia"])
selected.synd <- factor("Achondroplasia", levels = c("Unaffected Unrelated", "Achondroplasia"))
Achondroplasians.index <- which(d.meta.combined$Syndrome == "Achondroplasia" | d.meta.combined$Syndrome == "Unaffected Unrelated")
Achondroplasians.scores <- PC.scores[Achondroplasians.index,]
Achondroplasians.meta <- d.meta.combined[Achondroplasians.index,]
Achondroplasians.meta$Syndrome <- factor(Achondroplasians.meta$Syndrome, levels = c("Unaffected Unrelated", "Achondroplasia"))

Achondroplasians.coefs <- manova(Achondroplasians.scores ~ Achondroplasians.meta$Sex + Achondroplasians.meta$Age + Achondroplasians.meta$Age^2 + Achondroplasians.meta$Age^3 + Achondroplasians.meta$Syndrome + Achondroplasians.meta$Age:Achondroplasians.meta$Syndrome)$coef
Achondroplasians.datamod <- ~ selected.sex + selected.age + selected.age^2 + selected.age^3 + selected.synd + selected.age:selected.synd

#calculate score for the main effect
S <- Achondroplasians.coefs[4,] #Achondroplasians.coefs[grepl(pattern = selected.synd, rownames(Achondroplasians.coefs)),][1,]
Snorm <- S/sqrt(sum(S^2))
syndscores <- Achondroplasians.scores %*% Snorm

#plot scores of syndrome versus non-syndrome
boxplot((syndscores)  ~ Achondroplasians.meta$Syndrome)

#math check
predicted.shape <- predshape.lm(Achondroplasians.coefs[,1:500], Achondroplasians.datamod, PC.eigenvectors[,1:500], synd.mshape)
main.res <- matrix(t(PC.eigenvectors %*% (-.012 * Snorm)), dim(synd.mshape)[1], dim(synd.mshape)[2])
interaction.res <- matrix(t(PC.eigenvectors %*% (min(syndscores2[Achondroplasians.meta$Syndrome == "Achondroplasia"]) * Snorm2)), dim(synd.mshape)[1], dim(synd.mshape)[2])

#3d viz
# atlas <- file2mesh("baking/jordan_registered.ply")
Achondroplasia.mesh <- atlas
Achondroplasia.mesh$vb[-4,] <-  t(predicted.shape + main.res)# + interaction.res)
# samplereg <- sample(1:27000, 100)
Achondroplasia.max <- rotmesh.onto(Achondroplasia.mesh, t(Achondroplasia.mesh$vb[-4,samplereg]),  t(atlas$vb[-4,samplereg]), scale = T)$mesh
mesh2ply(Achondroplasia.max, filename = "Achondroplasia.max")
Achondroplasia.max <- file2mesh("Achondroplasia.max.ply")
plot3d(Achondroplasia.max, aspect = "iso", specular = "black", color = "lightgrey")


#Achondroplasia pheno####
#take the original prediction, mean Achondroplasia and add Shape from estimated PC scores
predicted.shape <- predshape.lm(Achondroplasians.coefs[,1:500], Achondroplasians.datamod, PC.eigenvectors[,1:500], synd.mshape)

#standard pred
#3d viz
# atlas <- file2mesh("baking/jordan_registered.ply")
Achondroplasia.mean <- atlas
Achondroplasia.mean$vb[-4,] <-  t(predicted.shape[,])
samplereg <- sample(1:27000, 100)
Achondroplasia.mean <- rotmesh.onto(Achondroplasia.mean, t(Achondroplasia.mean$vb[-4,samplereg]),  t(atlas$vb[-4,samplereg]), scale = T)$mesh
mesh2ply(Achondroplasia.mean, filename = "Achondroplasia_mean")
Achondroplasia.mean <- file2mesh("Achondroplasia_mean.ply")

# plot3d(Achondroplasia.mean, aspect = "iso", specular = 1, color = "lightgrey")
par3d(windowRect = c(0,0, 950,1000), userMatrix = lat.face, zoom = .75)
plot3d(Achondroplasia.mean, aspect = "iso", specular = 1, color = "lightgrey", add =F, axes = F, box = F, xlab = "", ylab = "", zlab = "")

#screenshot
system(paste0("screencapture -R 20,60,900,950 ~/shiny/shinyapps/Syndrome_model/Achond_mean_lat.png"))

# light3d(viewpoint.rel = T)
# rglwidget()

# mesh2ply(Achondroplasia.mean, "ns_mean")

#most main/int effect Achondroplasia ind
main.res <- matrix(t(PC.eigenvectors %*% (.05 * Snorm)), dim(synd.mshape)[1], dim(synd.mshape)[2])

#3d viz
# atlas <- file2mesh("baking/jordan_registered.ply")
Achondroplasia.mesh <- atlas
Achondroplasia.mesh$vb[-4,] <-  t(predicted.shape + main.res)
# samplereg <- sample(1:27000, 100)
Achondroplasia.max <- rotmesh.onto(Achondroplasia.mesh, t(Achondroplasia.mesh$vb[-4,samplereg]),  t(atlas$vb[-4,samplereg]), scale = T)$mesh
mesh2ply(Achondroplasia.max, filename = "Achondroplasia.max")
Achondroplasia.max <- file2mesh("Achondroplasia.max.ply")

# plot3d(Achondroplasia.mean, aspect = "iso", specular = 1, color = "lightgrey")
par3d(windowRect = c(0,0, 950,1000), userMatrix = diag.face, zoom = .75)
plot3d(Achondroplasia.max, aspect = "iso", specular = 1, color = "lightgrey", add =F, axes = F, box = F, xlab = "", ylab = "", zlab = "")

#screenshot
system(paste0("screencapture -R 20,60,900,950 ~/shiny/shinyapps/Syndrome_model/Achond_max_diag.png"))

#diagnostic overlay
plot3d(Achondroplasia.max, aspect = "iso", specular = "black", color = "lightgrey")
plot3d(Achondroplasia.mean, add = T, col = 2)

#heatmap
meshDist(Achondroplasia.mean, Achondroplasia.max)
system(paste0("screencapture -R 20,60,900,950 ~/shiny/shinyapps/Syndrome_model/Achond_severity_heatmap_front.png"))

#least main/int effect Achondroplasia ind
main.res <- matrix(t(PC.eigenvectors %*% (-.011 * Snorm)), dim(synd.mshape)[1], dim(synd.mshape)[2])

#3d viz
# atlas <- file2mesh("baking/jordan_registered.ply")
Achondroplasia.mesh <- atlas
Achondroplasia.mesh$vb[-4,] <-  t(predicted.shape + main.res)
# samplereg <- sample(1:27000, 100)
Achondroplasia.min <- rotmesh.onto(Achondroplasia.mesh, t(Achondroplasia.mesh$vb[-4,samplereg]),  t(atlas$vb[-4,samplereg]), scale = T)$mesh
mesh2ply(Achondroplasia.min, filename = "Achondroplasia_min")
Achondroplasia.min <- file2mesh("Achondroplasia_min.ply")
par3d(windowRect = c(0,0, 950,1000), userMatrix = front.face, zoom = .75)
plot3d(Achondroplasia.min, aspect = "iso", specular = 1, color = "lightgrey", add =F, axes = F, box = F, xlab = "", ylab = "", zlab = "")

#screenshot
system(paste0("screencapture -R 20,60,900,950 ~/shiny/shinyapps/Syndrome_model/Achond_min_lat.png"))


plot3d(Achondroplasia.max, aspect = "iso", specular = "black", color = "red", alpha = .5)
plot3d(Achondroplasia.min, aspect = "iso", specular = "black", color = "lightgrey", add = T, alpha = .5)

#heatmap
meshDist(Achondroplasia.mean, Achondroplasia.min)
system(paste0("screencapture -R 20,60,900,950 ~/shiny/shinyapps/Syndrome_model/Achond_leastseverity_heatmap_front.png"))

#fig1B, image warping to generate a 2D gestalt####
#FB2 nager screenshots

library(Morpho)
library(rgl)

#load metadata
setwd("~/shiny/shinyapps/Syndrome_model/")
load("data.Rdata")

#plot 1st face, define front pose
setwd("~/Downloads/FB2_baked/")
faces <- list.files(".", pattern = "*.ply")

first.mesh <- file2mesh(faces[1], readcol = T)
par3d(windowRect = c(0,0, 950,1000), userMatrix = front.face, zoom = .65)
plot3d(first.mesh, aspect = "iso", lit = F, xlab = "", ylab = "", zlab = "", specular = 1)

#see if all faces are already registered
tmp.mesh <- file2mesh(faces[200], readcol = T)
par3d(windowRect = c(0,0, 950,1000), userMatrix = front.face, zoom = .65)
plot3d(tmp.mesh, aspect = "iso", lit = F, xlab = "", ylab = "", zlab = "", specular = 1)

#not registered, so let's do that
samplereg <- sample(1:27000, 100)

#affine register to first mesh####
pose <- "lat.face2"
achond.scans <- d.meta.combined$scan[d.meta.combined$Syndrome == "Achondroplasia"]
for(i in achond.scans){
  tmp.mesh <- file2mesh(i, readcol = T)
  tmp.mesh <- rotmesh.onto(tmp.mesh, t(tmp.mesh$vb[-4,samplereg]), t(first.mesh$vb[-4,samplereg]), scale = T)$mesh
  
  
  par3d(windowRect = c(0,0, 950,1000), userMatrix = get(pose), zoom = .75)
  plot3d(tmp.mesh, aspect = "iso", lit = F, xlab = "", ylab = "", zlab = "", specular = 1, axes = F, box = F)
  
  #screenshot
  system(paste0("screencapture -R 20,60,900,950 ~/Downloads/FB2_baked_screenshots/", i, "_", pose,".png"))
  
}

#non-affine registration and screenshots####
#register to estimated average face generated from the atlas model
dense.achond <- file2mesh("~/shiny/shinyapps/Syndrome_model/FB2_baked_screenshots/Achondroplasia/achond_mean.ply")

pose <- "lat.face2"
achond.scans <- d.meta.combined$scan[d.meta.combined$Syndrome == "Achondroplasia"]
for(i in achond.scans){
  tmp.mesh <- file2mesh(i, readcol = T)
  tmp.mesh <- rotmesh.onto(tmp.mesh, t(tmp.mesh$vb[-4,samplereg]), t(first.mesh$vb[-4,samplereg]), scale = T)$mesh
  tmp.mesh <- tps3d(tmp.mesh, t(tmp.mesh$vb[-4,samplereg]), t(dense.achond$vb[-4,samplereg]))
  
  
  par3d(windowRect = c(0,0, 950,1000), userMatrix = get(pose), zoom = .75)
  plot3d(tmp.mesh, aspect = "iso", lit = F, xlab = "", ylab = "", zlab = "", specular = 1, axes = F, box = F)
  
  #screenshot
  system(paste0("screencapture -R 20,60,900,950 ~/shiny/shinyapps/Syndrome_model/FB2_baked_screenshots/Achondroplasia/tps/", i, "_", pose,".png"))
  
}



#load in screenshots and average them
library(EBImage)
library(Colormesh)
library(imager)

face.images <- list.files("~/shiny/shinyapps/Syndrome_model/FB2_baked_screenshots/Achondroplasia/tps/", pattern = "*lat.face2.png")
test.image <- image_reader("~/shiny/shinyapps/Syndrome_model/FB2_baked_screenshots/Achondroplasia/tps/", image.names = "1502241408312nd.ply_front.face.png")
tmp.array <- array(NA, dim = c(nrow(test.image), ncol(test.image), 3, length(face.images)))

for(i in 1:length(face.images)){
  tmp.image <- image_reader("~/shiny/shinyapps/Syndrome_model/FB2_baked_screenshots/Achondroplasia/tps/", image.names = face.images[i])
  tmp.array[,,,i] <- tmp.image[,,1:3]
  print(i)
}

gestalt <- as.cimg(rowMeans(tmp.array, dims = 3))

plot(gestalt)

save.image(gestalt, "~/shiny/shinyapps/Syndrome_model/FB2_baked_screenshots/Achondroplasia/nonlinear_lat_achond.png")

#fig1C, wireframe comparisons####
#get 65 landmark set
# sparse.atlas <- file2mesh("~/Downloads/atlas_sparse.ply")
# atlas.lms <- read.table("~/Downloads/landmarks_65.txt")

sparse.atlas <- file2mesh("~/shiny/shinyapps/Syndrome_model/Full.ply")
atlas.lms <- read.mpp("~/shiny/shinyapps/Syndrome_model/Full_picked_points2.pp")

#load achond mean, register atlas mesh, transform lms, create achond wireframe
Achondroplasia.mesh <- file2mesh("FB2_baked_screenshots/Achondroplasia/achond_mean.ply")
samplereg <- sample(1:27000, 100)
sparse.atlas <- rotmesh.onto(sparse.atlas, t(sparse.atlas$vb[-4,samplereg]), t(Achondroplasia.mesh$vb[-4,samplereg]), scale = T)

atlas.lms2 <- applyTransform(atlas.lms, sparse.atlas$trafo)#atlas.lms %*% sparse.atlas$trafo[-4,]

achond.lms <- tps3d(atlas.lms2, refmat = sparse.atlas$yrot, tarmat = t(Achondroplasia.mesh$vb[-4,samplereg]))

par3d(windowRect = c(0,0, 950,1000), userMatrix =  diag.face2, zoom = .75)

plot3d(sparse.atlas$mesh, aspect = "iso", specular = "black", color = "lightgrey", alpha = .5, box = F, axes = F, xlab = "", ylab = "", zlab = "")
# plot3d(Achondroplasia.mesh, add = T, col = 2)
# spheres3d(atlas.lms, col = 2)
# text3d(atlas.lms2, texts = 1:65, col = 2)

lines3d(atlas.lms2[c(1,4, 16,15,17,39,16,4,42,41,43,65, 42,4,5,34,6,60,5,NA,53,52,51,54,64,13,12,11,36,24,8,50,62,11,8,7,6, NA,24,9,50, NA,27,26,25,28,38,13, NA,34,21,6,47,60, NA,21,22,6,5,6,48,47, NA, 14, 1, 40, 2,14, NA, 30,19,29, 30, NA,56,45,55,56),], lwd = 3)
lines3d(achond.lms[c(1,4, 16,15,17,39,16,4,42,41,43,65, 42,4,5,34,6,60,5,NA,53,52,51,54,64,13,12,11,36,24,8,50,62,11,8,7,6, NA,24,9,50, NA,27,26,25,28,38,13, NA,34,21,6,47,60, NA,21,22,6,5,6,48,47, NA, 14, 1, 40, 2,14, NA, 30,19,29, 30, NA,56,45,55,56),], col = 2, lwd = 3)
# rglwidget(width = 1000, height = 1000)

#screenshot
system(paste0("screencapture -R 20,60,900,950 ~/Downloads/achond_wireframe_diag.png"))






#fig 2, gestalts, heatmaps, morphospace?####
#fig 2A, 4 y/o female nager####
selected.sex <- 1
selected.age <- 4
selected.synd <- factor("Nager Syndrome", levels = c("Unaffected Unrelated", "Nager Syndrome"))
nagerns.index <- which(d.meta.combined$Syndrome == "Nager Syndrome" | d.meta.combined$Syndrome == "Unaffected Unrelated")
nagerns.scores <- PC.scores[nagerns.index,]
nagerns.meta <- d.meta.combined[nagerns.index,]
nagerns.meta$Syndrome <- factor(nagerns.meta$Syndrome, levels = c("Unaffected Unrelated", "Nager Syndrome"))

nagerns.coefs <- manova(nagerns.scores ~ nagerns.meta$Sex + nagerns.meta$Age + nagerns.meta$Age^2 + nagerns.meta$Age^3 + nagerns.meta$Syndrome + nagerns.meta$Age:nagerns.meta$Syndrome)$coef
nagerns.datamod <- ~ selected.sex + selected.age + selected.age^2 + selected.age^3 + selected.synd + selected.age:selected.synd

#calculate score for the main effect
S <- nagerns.coefs[4,]
Snorm <- S/sqrt(sum(S^2))
syndscores <- nagerns.scores %*% Snorm

#plot scores of syndrome versus non-syndrome
boxplot((syndscores)  ~ nagerns.meta$Syndrome)


#nager pheno####
#take the original prediction, 5 y/o M Nager and add Shape from estimated PC scores
predicted.shape <- predshape.lm(nagerns.coefs[,1:500], nagerns.datamod, PC.eigenvectors[,1:500], synd.mshape)

#most main/int effect nager ind
main.res <- matrix(t(PC.eigenvectors %*% (min(syndscores[nagerns.meta$Syndrome == "Nager Syndrome"]) * Snorm)), dim(synd.mshape)[1], dim(synd.mshape)[2])

#3d viz
nager.mesh <- atlas
nager.mesh$vb[-4,] <-  t(predicted.shape + main.res) 
nager.mesh <- rotmesh.onto(nager.mesh, t(nager.mesh$vb[1:3,samplereg]), t(atlas$vb[1:3,samplereg]), scale = T)$mesh
plot3d(nager.mesh, aspect = "iso", specular = "black", lit = F)


#generic shiny code####
pred.mesh <- atlas 
selected.sex <- 1
selected.age <- 20
selected.synd <- factor("Treacher Collins Syndrome", levels = levels(d.meta.combined$Syndrome))
datamod <- ~ selected.sex + selected.age + selected.age^2 + selected.age^3 + selected.synd + selected.age:selected.synd

synd.lm.coefs <- lm(PC.scores ~ d.meta.combined$Sex + d.meta.combined$Age + d.meta.combined$Age^2 + d.meta.combined$Age^3 + d.meta.combined$Syndrome + d.meta.combined$Age:d.meta.combined$Syndrome)$coef

predicted.shape <- predshape.lm(synd.lm.coefs, datamod, PC.eigenvectors, synd.mshape)


#synd
pred.mesh$vb[-4,] <-  t(predicted.shape)
rigid.points <- sample(1:27000, 100)
test <- rotmesh.onto(pred.mesh, t(pred.mesh$vb[1:3,rigid.points]), t(atlas$vb[1:3,rigid.points]), scale = T)

plot3d(test$mesh, aspect = "iso", axes = F, box = F, xlab = "", ylab = "", zlab = "", lit = F, specular = 1)
plot3d(test$mesh, aspect = "iso", axes = F, box = F, xlab = "", ylab = "", zlab = "", col = "lightgrey", specular = 1)






#trying out claes style spectral clustering of lms####

setwd("/mnt/Hallgrimsson/Users/Jovid/Andlit/")
load("Dense_combined_starter_PC_outliers_duplicates_removed.Rdata")

#get point-wise correlations


atlas.seg <- kmeans(PC.eigenvectors[,1:500], 2)
seg.colors <- atlas.seg$cluster[seq(1,length(atlas.seg$cluster), 3)]

plot3d(atlas, col = seg.colors, aspect = "iso", specular = 1)




#rolling our own segments####
eyes <- as.numeric(read.csv("../Perception/eyeslms", header = F)) + 1
nose <- as.numeric(read.csv("../Perception/noselms.txt", header = F)) + 1
jaw <- as.numeric(read.csv("../Perception/jawlms", header = F)) + 1
chin <- as.numeric(read.csv("../Perception/chinlms", header = F)) + 1
mouth <- as.numeric(read.csv("../Perception/mouthlms", header = F)) + 1
cheeks <- as.numeric(read.csv("../Perception/cheeklms", header = F)) + 1
forehead <- as.numeric(read.csv("../Perception/foreheadlms", header = F)) + 1
brow <- c(as.numeric(read.csv("../Perception/brow_left", header = F)) + 1, as.numeric(read.csv("../Perception/brow_right", header = F)) + 1)
face <- 1:ncol(atlas$vb)

for(i in 1:8){
  mod.names <- c("face", "eyes", "nose", "jaw", "chin", "mouth", "cheeks", "forehead")
  colorvec<- rep("lightgrey", ncol(atlas$vb))
  colorvec[get(mod.names[i])] <- "#a6192e"
  
  mesh2ply(atlas, "tmpmesh")
  atlas <- file2mesh("tmpmesh.ply", readcol = T)
  par3d(windowRect = c(0,0, 950,1000), zoom = .75)
  plot3d(atlas, col = colorvec, specular = 1, aspect = "iso", box = F, axes = F, xlab = "", ylab = "", zlab = "")
  
  #screenshot
  system(paste0("screencapture -R 20,60,900,950 ~/Desktop/mod_images/mod", i,".png"))
}

#create index matrix
modules <- matrix(0, nrow = nrow(d.lms.combined), ncol = 8)
colnames(modules) <- c("Whole Face", "Nose", "Eyes", "Jaw", "Chin", "Mouth", "Cheeks", "Forehead")
modules[,1] <- 1
modules[nose,2] <- 1
modules[eyes,3] <- 1
modules[jaw,4] <- 1
modules[chin,5] <- 1
modules[mouth,6] <- 1
modules[cheeks,7] <- 1
modules[forehead,8] <- 1


synd.mat <- array(NA, dim = c(nrow(d.lms.combined), 3, length(unique(d.meta.combined$Syndrome))))

for(i in 1:dim(synd.mat)[3]){
  
  selected.sex <- 1
  selected.age <- 15
  selected.synd <- factor(unique(d.meta.combined$Syndrome)[i], levels = levels(d.meta.combined$Syndrome))
  datamod <- ~ selected.sex + selected.age + selected.age^2 + selected.age^3 + selected.synd + selected.age:selected.synd
  
  synd.mat[,,i] <- predshape.lm(synd.lm.coefs, datamod, PC.eigenvectors, synd.mshape)
  
}

#morphospace tests####
selected.node <- 2
module.names <- colnames(modules)

full.morpho <- procdist.array(synd.mat, synd.mat[,,unique(d.meta.combined$Syndrome) == "Nager Syndrome"])

sub.morpho <- procdist.array(synd.mat[modules[,selected.node] == 1,,], synd.mat[modules[,selected.node] == 1,, unique(d.meta.combined$Syndrome) == "Nager Syndrome"])



#calculate syndrome severity scores for selected syndrome
#calculate score for the whole face
S <- synd.lm.coefs[grepl(pattern = "Nager Syndrome", rownames(synd.lm.coefs)),][1,]
Snorm <- S/sqrt(sum(S^2))
syndscores.main <- PC.scores %*% Snorm

syndscores.df <- data.frame(Syndrome = d.meta.combined$Syndrome, face.score = syndscores.main)
syndscores.wholeface <- syndscores.df%>%
  group_by(Syndrome) %>%
  summarise(face_score = mean(face.score))

#calculate score for the selected subregion
#for each subregion, we'll need a pca just for those lms
# c("Whole Face", "Nose", "Eyes", "Jaw", "Chin", "Mouth", "Cheeks", "Forehead")
# nose.pca <- prcomp(array2row3d(d.lms.combined[nose,,]))
# eyes.pca <- prcomp(array2row3d(d.lms.combined[eyes,,]))
# jaw.pca <-  prcomp(array2row3d(d.lms.combined[jaw,,]))
# chin.pca <-  prcomp(array2row3d(d.lms.combined[chin,,]))
# mouth.pca <-  prcomp(array2row3d(d.lms.combined[mouth,,]))
# cheeks.pca <-  prcomp(array2row3d(d.lms.combined[cheeks,,]))
# forehead.pca <-  prcomp(array2row3d(d.lms.combined[forehead,,]))

load("modules_PCA.Rdata")

nose.coefs <- manova(nose.pca$x ~ d.meta.combined$Sex + d.meta.combined$Age + d.meta.combined$Age^2 + d.meta.combined$Age^3 + d.meta.combined$Syndrome + d.meta.combined$Age:d.meta.combined$Syndrome)$coef
S <- nose.coefs[grepl(pattern = "Nager Syndrome", rownames(nose.coefs)),][1,]
Snorm <- S/sqrt(sum(S^2))
syndscores.main <- nose.pca$x %*% Snorm

boxplot(syndscores.main ~ d.meta.combined$Syndrome == "Nager Syndrome")

# S <- synd.lm.coefs[grepl(pattern = "Nager Syndrome", rownames(synd.lm.coefs)),][1, modules[,selected.node] == 1]
# Snorm <- S/sqrt(sum(S^2))
# syndscores.main <- PC.scores %*% Snorm

syndscores.df <- data.frame(Syndrome = d.meta.combined$Syndrome, module.score = syndscores.main)
syndscores.module <- syndscores.df%>%
  group_by(Syndrome) %>%
  summarise(face_score = mean(module.score))


morpho.df <- data.frame(full = syndscores.wholeface$face_score, sub = syndscores.module$face_score, Syndrome = levels(d.meta.combined$Syndrome))

highlight_df <- morpho.df[sort(morpho.df$sub, index.return = T)$ix,]
#to do
#highlight 5 closest syndromes with names
#maybe all points should be names, with 5 highlighted
p <- ggplot(morpho.df) +
  geom_point(aes(y = full, x = sub), color = "white") +
  geom_point(data = highlight_df[1:5,], aes(y = full, x = sub), color = "white", fill = "#fca311", shape = 21, size = 3) +
  geom_label_repel(data = highlight_df[c(1:5, 85:89),], aes(y = full, x = sub, label = Syndrome),
                   force = 2,
                   size = 4,
                   color = c(rep("#fca311", 5), rep("white", 5)),
                   fill = "#1a1a1a") +
  ylab(paste0("Full face Nager Syndrome score")) +
  xlab(paste0(module.names[selected.node], " Nager Syndrome score")) +
  # ylim(-.25,3.5) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "#1a1a1a"),
        plot.background = element_rect(fill = "#1a1a1a"),
        axis.text = element_text(color = "white"),
        axis.title = element_text(color = "white"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "#404040"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "#404040")
  )


g <- ggplotGrob(p)
bg <- g$grobs[[1]]
round_bg <- roundrectGrob(x=bg$x, y=bg$y, width=bg$width, height=bg$height,
                          r=unit(0.1, "snpc"),
                          just=bg$just, name=bg$name, gp=bg$gp, vp=bg$vp)
g$grobs[[1]] <- round_bg

plot(g)






#does the model tend to overshoot and create people more syndromic than their counterparts?####

predmoresevere <- rep(NA, nrow(d.meta.combined))

for(i in 1:nrow(d.meta.combined)){
selected.age <- d.meta.combined$Age[i]
selected.sex <- if(d.meta.combined$Sex[i] == "M"){ 1} else {2} 
selected.synd <- factor(d.meta.combined$Syndrome[i], levels = levels(d.meta.combined$Syndrome))

datamod <- ~ selected.sex + selected.age + selected.age^2 + selected.age^3 + selected.synd + selected.age:selected.synd

predicted.shape <- predshape.lm(synd.lm.coefs, datamod, PC.eigenvectors, synd.mshape)

#project prediction into PC space
S <- as.numeric(synd.lm.coefs[grepl(pattern = d.meta.combined$Syndrome[i], rownames(synd.lm.coefs)),])[1:500] #Achondroplasians.coefs[grepl(pattern = selected.synd, rownames(Achondroplasians.coefs)),][1,]
Snorm <- S/sqrt(sum(S^2))
syndscores <- PC.scores %*% Snorm

predpcscore <- getPCscores(predicted.shape, PC.eigenvectors, synd.mshape)[,1:dim(PC.scores)[2]] %*% Snorm

#who is more distant from the mean syndscore?
meanss <- mean(syndscores[d.meta.combined$Syndrome == d.meta.combined$Syndrome[i]])

distreal <- sqrt((syndscores[i] - meanss)^2)
distsynth <- sqrt((predpcscore - meanss)^2)

predmoresevere[i] <- if(distsynth > distreal){1} else{0} 
}

sum(predmoresevere)/2748

which(predmoresevere)

#plot scores of syndrome versus non-syndrome
boxplot((syndscores)  ~ d.meta.combined$Syndrome == d.meta.combined$Syndrome[i])










