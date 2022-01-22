library(opencv)
library(Morpho)
library(Rvcg)
library(rgl)

#load mesh####
test.mesh <- file2mesh("~/shiny/shinyapps/Perception/to_bake/jovid_baked.ply", readcol = T)


#generate screenshot####
par3d(userMatrix = front.face, zoom = .75, windowRect = c(0,0, 950,1000))
plot3d(test.mesh, aspect = "iso", specular = 1, axes = F, box = F, xlab = "", ylab = "", zlab = "")

img.save <- rglwidget(snapshot=TRUE)
system(paste0("mv ", img.save[1], " ", "jovid.png"))

#get landmarks from the screenshot####

unconf <- ocv_read('jovid.png')
facemask <- ocv_facemask(unconf)
attr(facemask, 'faces')

keypoints <- ocv_keypoints(unconf, method = "FAST", type = "TYPE_9_16", threshold = 40)

plot(imager::load.image("jovid.png"))
points(keypoints$x, keypoints$y, col = 2)

#try raycasting those landmarks at the mesh


