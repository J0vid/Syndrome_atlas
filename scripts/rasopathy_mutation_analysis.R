#mutation analysis
library(Jovid)
library(Morpho)
library(geomorph)


#load syndromic data
setwd("Downloads/hammond_lms/")
claes_data <- read.csv("MetadataAndLandmarks24-Jul-2020.txt")
claes_data <- claes_data[is.na(Age) == F,]

#register
claes_lms <- claes_data[,21:215]

claes_reg <- procSym(row2array3d(claes_lms))

#standardize for age + sex
Age <- claes_data$age
Age2 <- claes_data$age^2
Age3 <- claes_data$age^3

# gdf for regression components
pr_gm_s <- geomorph.data.frame(shape= claes_reg$rotated[,,is.na(Age) == F], Age= Age[is.na(Age) == F], Age_2=Age2[is.na(Age) == F],Age_3=Age3[is.na(Age) == F], Sex= as.factor(claes_data$sexNumeric)[is.na(Age) == F])

# regression
age_sex_lm_s <- procD.lm(shape ~ Age + Age_2 + Age_3 + Sex , data = pr_gm_s, iter = 1)

# regression residuals  
resids <- age_sex_lm_s$residuals

# add back the grandmean and convert to 2d array
# mean_sym <- apply(sym_dat, 2, mean)
# resid_shapes <- sweep(resids, MARGIN = 2, STATS = mean_sym, FUN = "+")

# add back the average 20 year old face and convert to 2d array
mu_age <- 10
mean_sym <- mu_age * age_sex_lm_s$coefficients[2,] + mu_age^2 * age_sex_lm_s$coefficients[3,] + mu_age^3 * age_sex_lm_s$coefficients[4,] + age_sex_lm_s$coefficients[1,]
resid_shapes <- sweep(resids, MARGIN = 2, STATS = mean_sym, FUN = "+")


#filter registered data
claes_ras <- claes_data[claes_data$superLabel == "Rasopathy",]

#load mutation data
library(readxl)
Metadata_ras <- read_excel("Metadata_export rasopathy 6.9.2020.xlsx", sheet = "data avail - needed", skip = 1)
View(Metadata_ras)

#cross reference confirmed mutations with landmark data
BRAF <- as.character(na.omit(Metadata_ras$Title[Metadata_ras$`Short Molecular Diag` == "BRAF"]))
PTPN <- as.character(na.omit(Metadata_ras$Title[Metadata_ras$`Short Molecular Diag` == "PTPN11"]))
HRAS <- as.character(na.omit(Metadata_ras$Title[Metadata_ras$`Short Molecular Diag` == "HRAS"]))
NF1 <- as.character(na.omit(Metadata_ras$Title[Metadata_ras$`Short Molecular Diag` == "NF1"]))

combined.BRAF <- unique(c(which(claes_data$CalgarySubjectID %in% BRAF), which(claes_data$genotype == "BRAF")))
combined.PTPN <- unique(c(which(claes_data$CalgarySubjectID %in% PTPN), which(claes_data$genotype == "PTPN11")))
combined.HRAS <- unique(c(which(claes_data$CalgarySubjectID %in% HRAS), which(claes_data$genotype == "HRAS")))
combined.NF1 <- unique(c(which(claes_data$CalgarySubjectID %in% NF1), which(claes_data$genotype == "NF1")))



BRAF.mean <- arrayspecs(t(as.matrix(colMeans(resid_shapes[combined.BRAF,]))), 65, 3)[,,1]
PTPN.mean <- arrayspecs(t(as.matrix(colMeans(resid_shapes[combined.PTPN,]))), 65, 3)[,,1]
HRAS.mean <- arrayspecs(t(as.matrix(colMeans(resid_shapes[combined.HRAS,]))), 65, 3)[,,1]
NF1.mean <- arrayspecs(t(as.matrix(colMeans(resid_shapes[combined.NF1,]))), 65, 3)[,,1]

#get the atlas in this space
atlas <- file2mesh("atlas.ply")
atlas.lms <- as.matrix(read.table("landmarks_65.txt"))

plot3d(atlas, col = "lightgrey", alpha = .3, aspect = "iso", axes = F, specular = 1, xlab = "", ylab = "", zlab = "")
points3d(atlas.lms, col = 2)

BRAF.mesh <- tps3d(atlas, atlas.lms, BRAF.mean, scale = T)
PTPN.mesh <- tps3d(atlas, atlas.lms, PTPN.mean, scale = T)
HRAS.mesh <- tps3d(atlas, atlas.lms, HRAS.mean, scale = T)
NF1.mesh <- tps3d(atlas, atlas.lms, NF1.mean, scale = T)

plot3d(BRAF.mesh, col = "lightgrey", alpha = .3, aspect = "iso", axes = F, specular = 1, xlab = "", ylab = "", zlab = "")
points3d(BRAF.mean, col = 2)

#overall means for comparison
Cardiofaciocutaneous.mean <- arrayspecs(t(as.matrix(colMeans(resid_shapes[claes_data$diagnosis == "Cardiofaciocutaneous",]))), 65, 3)[,,1]
Noonan.mean <- arrayspecs(t(as.matrix(colMeans(resid_shapes[claes_data$diagnosis == "Noonan",]))), 65, 3)[,,1]
Costello.mean <- arrayspecs(t(as.matrix(colMeans(resid_shapes[claes_data$diagnosis == "Costello",]))), 65, 3)[,,1]
Neurofibromatosis.mean <- arrayspecs(t(as.matrix(colMeans(resid_shapes[claes_data$diagnosis == "Neurofibromatosis",]))), 65, 3)[,,1]

Cardiofaciocutaneous.mesh <- tps3d(atlas, atlas.lms, Cardiofaciocutaneous.mean, scale = T)
Noonan.mesh <- tps3d(atlas, atlas.lms, Noonan.mean, scale = T)
Costello.mesh <- tps3d(atlas, atlas.lms, Costello.mean, scale = T)
Neurofibromatosis.mesh <- tps3d(atlas, atlas.lms, Neurofibromatosis.mean, scale = T)

plot3d(Cardiofaciocutaneous.mesh, col = "lightgrey", alpha = .3, aspect = "iso", axes = F, specular = 1, xlab = "", ylab = "", zlab = "")
points3d(Cardiofaciocutaneous.mean, col = 2)
shade3d(BRAF.mesh, col = 2, alpha = .3)

#difference between overall mean and confirmed diagnosis mean
meshDist(Cardiofaciocutaneous.mesh, BRAF.mesh)

meshDist(Noonan.mesh, PTPN.mesh)

meshDist(Costello.mesh, HRAS.mesh)

meshDist(Neurofibromatosis.mesh, NF1.mesh)




#pathway position####
ns.mean <- arrayspecs(t(as.matrix(colMeans(resid_shapes[claes_data$diagnosis == "Control",]))), 65, 3)[,,1]
ns.mesh <- tps3d(atlas, atlas.lms, ns.mean, scale = T)

#NRAS
nras.mesh <- tps3d(atlas, atlas.lms, arrayspecs(t(as.matrix(resid_shapes[claes_data$genotype == "NRAS mutation",])), 65, 3)[,,1], scale = T)
mek.mesh <- tps3d(atlas, atlas.lms, arrayspecs(t(as.matrix(colMeans(resid_shapes[claes_data$genotype == "MEK1" | claes_data$genotype == "MEK2",]))), 65, 3)[,,1], scale = T)
raf1.mesh <- tps3d(atlas, atlas.lms, arrayspecs(t(as.matrix(resid_shapes[claes_data$genotype == "RAF1 ",])), 65, 3)[,,1], scale = T)

meshDist(ns.mesh, BRAF.mesh)
meshDist(ns.mesh, PTPN.mesh)
meshDist(ns.mesh, HRAS.mesh)
meshDist(ns.mesh, NF1.mesh)

meshDist(ns.mesh, nras.mesh)
meshDist(ns.mesh, mek.mesh)
meshDist(ns.mesh, raf1.mesh)






