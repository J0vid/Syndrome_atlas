#syndrome atlas figures

#fig 1####
#syndrome heterogeneity####
#subset data to just selected synd and non-syndromic
selected.synd <- factor("Nager Syndrome", levels = c("Unaffected Unrelated", "Nager Syndrome"))
nagerns.index <- which(d.meta.combined$Syndrome == "Nager Syndrome" | d.meta.combined$Syndrome == "Unaffected Unrelated")
nagerns.scores <- PC.scores[nagerns.index,]
nagerns.meta <- d.meta.combined[nagerns.index,]
nagerns.meta$Syndrome <- factor(nagerns.meta$Syndrome, levels = c("Unaffected Unrelated", "Nager Syndrome"))

nagerns.coefs <- manova(nagerns.scores ~ nagerns.meta$Sex + nagerns.meta$Age + nagerns.meta$Age^2 + nagerns.meta$Age^3 + nagerns.meta$Syndrome + nagerns.meta$Age:nagerns.meta$Syndrome)$coef
nagerns.datamod <- ~ selected.sex + selected.age + selected.age^2 + selected.age^3 + selected.synd + selected.age:selected.synd

#calculate score for the main effect
S <- nagerns.coefs[grepl(pattern = selected.synd, rownames(nagerns.coefs)),][1,]
Snorm <- S/sqrt(sum(S^2))
sexscores <- nagerns.scores %*% Snorm

#calculate scores for the interaction effect
S2 <- nagerns.coefs[grepl(pattern = selected.synd, rownames(nagerns.coefs)),][2,]
Snorm2 <- S2/sqrt(sum(S2^2))
sexscores2 <- nagerns.scores %*% Snorm2

#plot scores of syndrome versus non-syndrome
boxplot((sexscores)  ~ nagerns.meta$Syndrome)
boxplot((sexscores2)  ~ nagerns.meta$Syndrome)

#nager pheno####
#take the original prediction, 5 y/o M Nager and add Shape from estimated PC scores
predicted.shape <- predshape.lm(nagerns.coefs, nagerns.datamod, PC.eigenvectors, synd.mshape)

#most main/int effect nager ind
main.res <- matrix(t(PC.eigenvectors %*% (max(sexscores) * S)), dim(synd.mshape)[1], dim(synd.mshape)[2])
interaction.res <- matrix(t(PC.eigenvectors %*% (min(sexscores2) * S2)), dim(synd.mshape)[1], dim(synd.mshape)[2])
#3d viz
nager.mesh <- pred.mesh
nager.mesh$vb[-4,] <-  t(predicted.shape + main.res + interaction.res) 
plot3d(nager.mesh, aspect = "iso", specular = "black", lit = F)

#least main/int effect nager ind
main.res <- matrix(t(PC.eigenvectors %*% (min(sexscores) * S)), dim(synd.mshape)[1], dim(synd.mshape)[2])
interaction.res <- matrix(t(PC.eigenvectors %*% (max(sexscores2) * S2)), dim(synd.mshape)[1], dim(synd.mshape)[2])
#3d viz
nager.mesh <- pred.mesh
nager.mesh$vb[-4,] <-  t(predicted.shape + main.res + interaction.res) 
plot3d(nager.mesh, aspect = "iso", specular = "black", lit = F)

#ns pheno####
#take the original prediction, 5 y/o M Nager and add Shape from estimated PC scores
selected.synd <- factor("Unaffected Unrelated", levels = c("Unaffected Unrelated", "Nager Syndrome"))
predicted.shape <- predshape.lm(nagerns.coefs, nagerns.datamod, PC.eigenvectors, synd.mshape)

#most main/int effect ns ind
main.res <- matrix(t(PC.eigenvectors %*% (max(sexscores) * S)), dim(synd.mshape)[1], dim(synd.mshape)[2])
interaction.res <- matrix(t(PC.eigenvectors %*% (min(sexscores2) * S2)), dim(synd.mshape)[1], dim(synd.mshape)[2])
#3d viz
nager.mesh <- pred.mesh
nager.mesh$vb[-4,] <-  t(predicted.shape + main.res + interaction.res) 
plot3d(nager.mesh, aspect = "iso", specular = "black", lit = F)

#least main/int effect nager ind
main.res <- matrix(t(PC.eigenvectors %*% (min(sexscores) * S)), dim(synd.mshape)[1], dim(synd.mshape)[2])
interaction.res <- matrix(t(PC.eigenvectors %*% (max(sexscores2) * S2)), dim(synd.mshape)[1], dim(synd.mshape)[2])
#3d viz
nager.mesh <- pred.mesh
nager.mesh$vb[-4,] <-  t(predicted.shape + main.res + interaction.res) 
plot3d(nager.mesh, aspect = "iso", specular = "black", lit = F)









#fit model
mf.mod <- manova(PC.scores ~ d.meta.combined$Sex + d.meta.combined$Age + d.meta.combined$Age^2 + d.meta.combined$Age^3 + d.meta.combined$Syndrome + d.meta.combined$Age:d.meta.combined$Syndrome)
# mf.mod <- manova(residuals2d ~ ns.cov$Sex[ns.sample], ns.cov$Age + poly(ns.cov$Age, 2))

#score individuals
#find columns for the selected syndrome
selected.synd <- factor("Nager Syndrome", levels = levels(d.meta.combined$Syndrome))
synd.lm.coefs[grepl(pattern = selected.synd, rownames(synd.lm.coefs)),1:2]# <- synd.lm.coefs.severity[grepl(pattern = selected.synd, rownames(synd.lm.coefs)),] + (input$severity/100) * synd.lm.coefs.severity[grepl(pattern = selected.synd, rownames(synd.lm.coefs)),]

#calculate score for the main effect
S <- synd.lm.coefs[grepl(pattern = selected.synd, rownames(synd.lm.coefs)),][1,]
Snorm <- S/sqrt(sum(S^2))
sexscores <- PC.scores[,1:100] %*% Snorm

#calculate scores for the interaction effect
S2 <- synd.lm.coefs[grepl(pattern = selected.synd, rownames(synd.lm.coefs)),][2,]
Snorm2 <- S2/sqrt(sum(S2^2))
sexscores2 <- PC.scores[,1:100] %*% Snorm2

#plot scores of syndrome versus non-syndrome
boxplot((sexscores)  ~ d.meta.combined$Syndrome == selected.synd)
boxplot((sexscores2)  ~ d.meta.combined$Syndrome == selected.synd)

#plot scores of nager vs treacher collins
ntc.index <- which(d.meta.combined$Syndrome == selected.synd | d.meta.combined$Syndrome == "Treacher Collins Syndrome")
View(d.meta.combined[ntc.index,])

boxplot((sexscores[ntc.index])  ~ d.meta.combined$Syndrome[ntc.index] == selected.synd)
boxplot((sexscores2[ntc.index])  ~ d.meta.combined$Syndrome[ntc.index] == selected.synd)

#plot scores of nager vs non-syndromic
ntc.index <- which(d.meta.combined$Syndrome == selected.synd | d.meta.combined$Syndrome == "Unaffected Unrelated")
View(d.meta.combined[ntc.index,])

boxplot((sexscores[ntc.index])  ~ d.meta.combined$Syndrome[ntc.index] == selected.synd)
boxplot((sexscores2[ntc.index])  ~ d.meta.combined$Syndrome[ntc.index] == selected.synd)

#what's a standard deviation for each score type
sd.main <- sd(sexscores[d.meta.combined$Syndrome == selected.synd])
sd.interaction <- sd(sexscores2[d.meta.combined$Syndrome == selected.synd])

#take the original prediction, 5 y/o M Nager and add Shape from estimated PC scores
#don't add the mean again, just the residual shape from the estimated PC scores
predicted.shape <- predshape.lm(synd.lm.coefs, datamod, PC.eigenvectors[,1:100], synd.mshape)

#most main/int effect nager ind
main.res <- matrix(t(PC.eigenvectors[,1:100] %*% (max(sexscores) * S)), dim(synd.mshape)[1], dim(synd.mshape)[2])
interaction.res <- matrix(t(PC.eigenvectors[,1:100] %*% (min(sexscores2) * S2)), dim(synd.mshape)[1], dim(synd.mshape)[2])
#synd
pred.mesh2$vb[-4,] <-  t(predicted.shape + main.res + interaction.res) 
plot3d(pred.mesh, aspect = "iso", specular = "black", lit = F)

#least main/int effect nager ind
main.res <- matrix(t(PC.eigenvectors[,1:100] %*% (min(sexscores) * S)), dim(synd.mshape)[1], dim(synd.mshape)[2])
interaction.res <- matrix(t(PC.eigenvectors[,1:100] %*% (max(sexscores2) * S2)), dim(synd.mshape)[1], dim(synd.mshape)[2])
#synd
pred.mesh$vb[-4,] <-  t(synd.mshape + main.res + interaction.res) 
plot3d(pred.mesh, aspect = "iso", specular = "black", lit = F)


meshDist(pred.mesh, pred.mesh2, displace = T, alpha = .5)
