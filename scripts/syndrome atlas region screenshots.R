#new subregion screenshots####
library(Morpho)
library(rgl)
library(Jovid)
library(Rvcg)

#load data####
setwd("~/shiny/shinyapps/Syndrome_model/")
# save(atlas, d.meta.combined, front.face, PC.eigenvectors, synd.lm.coefs, synd.mshape, PC.scores, synd.mat, file = "data.Rdata")
load("data.Rdata")

#load new modules####
painted_atlas <- file2mesh("syndrome_atlas_regions.ply", readcol = T)
plot3d(painted_atlas, aspect = "iso", xlab = "", ylab = "", zlab = "", box = F, axes = F, specular = 1)

unique(painted_atlas$material$color)

posterior_mandible <- painted_atlas$material$color == "#BA5EA9"
brow <- painted_atlas$material$color == "#BD0004" | painted_atlas$material$color == "#C10004"
eye <- painted_atlas$material$color == "#000000"
zygomatic<- painted_atlas$material$color == "#EC813F" | painted_atlas$material$color == "#EB814A"
anterior_mandible <- painted_atlas$material$color == "#4E2CBB"
premaxilla <- painted_atlas$material$color == "#24707C"
nose <- painted_atlas$material$color == "#F0F51C"
face <- 1:27903

#generate plots####
newmean <- atlas 
newmean$vb[-4,] <- t(synd.mshape)
newmean <- vcgSmooth(newmean)
plot3d(newmean, aspect = "iso", xlab = "", ylab = "", zlab = "", box = F, axes = F, col = "lightgrey", specular = 1, alpha = .6)

points3d(t(newmean$vb[-4, nose]), col = 2)

for(i in 1:7){
  mod.names <- c("face", "posterior_mandible", "nose", "anterior_mandible", "brow", "zygomatic", "premaxilla")
  colorvec<- rep("lightgrey", ncol(newmean$vb))
  colorvec[get(mod.names[i])] <- "#a6192e"
  
  # mesh2ply(newmean, "tmpmesh")
  # newmean <- file2mesh("tmpmesh.ply", readcol = T)
  par3d(windowRect = c(0,0, 950,1000), zoom = .75, userMatrix = front.face2)
  plot3d(newmean, col = colorvec, specular = 1, aspect = "iso", box = F, axes = F, xlab = "", ylab = "", zlab = "")
  
  #screenshot
  system(paste0("screencapture -R 20,60,900,950 ~/Desktop/mod_images2/mod", i,".png"))
}


#modules.Rdata####
#PCA for all regions, data.frame for lm indices in the same order as the images
load("/Users/jovid/Downloads/Dense_combined_starter_PC_outliers_duplicates_removed.Rdata")

posterior_mandible.pca <- prcomp(two.d.array(d.lms.combined[posterior_mandible,,]))
nose.pca <- prcomp(two.d.array(d.lms.combined[nose,,]))
brow.pca <- prcomp(two.d.array(d.lms.combined[brow,,]))
zygomatic.pca <- prcomp(two.d.array(d.lms.combined[zygomatic,,]))
premaxilla.pca <- prcomp(two.d.array(d.lms.combined[premaxilla,,]))
anterior_mandible.pca <- prcomp(two.d.array(d.lms.combined[anterior_mandible,,]))

modules <- data.frame(posterior_mandible, nose, anterior_mandible, brow, zygomatic, premaxilla)











