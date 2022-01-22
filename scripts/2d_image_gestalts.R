#FB2 achondroplasia screenshots

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
plot3d(tmp.mesh, aspect = "iso", lit = F, xlab = "", ylab = "", zlab = "", specular = 1)

#see if all faces are already registered
tmp.mesh <- file2mesh(faces[200], readcol = T)
par3d(windowRect = c(0,0, 950,1000), userMatrix = front.face, zoom = .65)
plot3d(tmp.mesh, aspect = "iso", lit = F, xlab = "", ylab = "", zlab = "", specular = 1)

#not registered, so let's do that
samplereg <- sample(1:27000, 100)

#affine register to first mesh####
pose <- "diag.face"
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
dense.achond <- file2mesh("~/Downloads/FB2_baked_screenshots/achond_mean.ply")

pose <- "front.face"
achond.scans <- d.meta.combined$scan[d.meta.combined$Syndrome == "Achondroplasia"]
for(i in achond.scans){
  tmp.mesh <- file2mesh(i, readcol = T)
  tmp.mesh <- rotmesh.onto(tmp.mesh, t(tmp.mesh$vb[-4,samplereg]), t(first.mesh$vb[-4,samplereg]), scale = T)$mesh
  tmp.mesh <- tps3d(tmp.mesh, t(tmp.mesh$vb[-4,samplereg]), t(dense.achond$vb[-4,samplereg]))
  
  
  par3d(windowRect = c(0,0, 950,1000), userMatrix = get(pose), zoom = .75)
  plot3d(tmp.mesh, aspect = "iso", lit = F, xlab = "", ylab = "", zlab = "", specular = 1, axes = F, box = F)
  
  #screenshot
  system(paste0("screencapture -R 20,60,900,950 ~/Downloads/FB2_baked_screenshots/tps/", i, "_", pose,".png"))
  
}



#load in screenshots and average them
library(EBImage)
library(Colormesh)
library(imager)

face.images <- list.files("~/Downloads/FB2_baked_screenshots/tps/", pattern = "*front.face.png")
test.image <- image_reader("~/Downloads/FB2_baked_screenshots/tps/", image.names = "1502241408312nd.ply_front.face.png")
tmp.array <- array(NA, dim = c(nrow(test.image), ncol(test.image), 3, length(face.images)))

for(i in 1:length(face.images)){
  tmp.image <- image_reader("~/Downloads/FB2_baked_screenshots/tps/", image.names = face.images[i])
  tmp.array[,,,i] <- tmp.image[,,1:3]
  print(i)
}

gestalt <- as.cimg(rowMeans(tmp.array, dims = 3))

plot(gestalt)

save.image(gestalt, "~/Downloads/FB2_baked_screenshots/nonlinear_front_achond.png")




