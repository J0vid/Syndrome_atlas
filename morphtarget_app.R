library(shiny)
library(rgl)
library(Morpho)

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

ui <- (fluidPage(
  selectInput("synd", label = "Syndrome", choices = sort(levels(d.meta.combined$Syndrome)), selected = "Achondroplasia"),
  selectInput("sex", label = "Sex", choices = c("Female", "Male")),
  sliderInput("age", label = "Age", min = 1, max = 70, value = 12, step = 1, animate = animationOptions(loop = T, interval = 40)),
  selectInput("severity", label = "Severity", choices = c("mild", "typical", "severe"), selected = "typical"),
  actionButton("update", label = "Update"),
  playwidgetOutput("control"),
  rglwidgetOutput("wdg")
))

server <- function(input, output, session) {
  options(rgl.useNULL = TRUE)
  save <- options(rgl.inShiny = TRUE)
  on.exit(options(save))
  
  close3d()
  
  #synd reactive####
  outVar <- reactive({
    #get selected syndrome age range
    age.range <- d.meta.combined$Age[d.meta.combined$Syndrome == input$synd]
    
    #calculate syndrome severity scores for selected syndrome
    #calculate score for the main effect
    S <- synd.lm.coefs[grepl(pattern = input$synd, rownames(synd.lm.coefs)),]#[1,]
    Snorm <- S/sqrt(sum(S^2))
    syndscores.main <- d.registered$PCscores[,1:num_pcs] %*% t(Snorm)
    # syndscores.main <- PC.scores %*% Snorm

    return(list(age.range, syndscores.main[d.meta.combined$Syndrome == input$synd]))
  })

  doutVar <- outVar #debounce(outVar, 2000)

  observe({
    slider_min <-  round(abs(min(doutVar()[[1]])))
    slider_max <- round(max(doutVar()[[1]]))
    updateSliderInput(session, "age", label = "Age", min = slider_min, max = slider_max, value = mean(doutVar()[[1]]))

  })
  
  selected.synd <- factor("Achondroplasia", levels = levels(d.meta.combined$Syndrome))
  selected.sex <-1
  selected.age <- 19 
  
  datamod <- ~ selected.sex + selected.age + selected.age^2 + selected.age^3 + selected.synd + selected.age:selected.synd
  predicted.shape <- predshape.lm(synd.lm.coefs, datamod, d.registered$PCs[,1:num_pcs], d.registered$mshape)

  atlas$vb[-4,] <- t(predicted.shape)
  tmp.mesh <- atlas
  
  open3d()
  shade3d(vcgSmooth(atlas), aspect = "iso", col = "lightgrey", specular = 1)
  objid <- ids3d()$id
  xyz <- rgl.attrib(objid[1], "vertices")
  dim(xyz)
  
  
  morph_target <- eventReactive(input$update, {
  #age morph target####
  minage <- round(abs(min(doutVar()[[1]])))
  nframes <- round(max(doutVar()[[1]]))
  if(input$sex == "Female"){selected.sex <- 2
  } else if(input$sex == "Male"){selected.sex <- -1}
  
  #severity math####
  S <- matrix(synd.lm.coefs[grepl(pattern = input$synd, rownames(synd.lm.coefs)),], nrow = 1, ncol = num_pcs)
  Snorm <- S/sqrt(sum(S^2))
  
  syndscores.main <- doutVar()[[2]]
  
  
  if(input$severity == "mild"){selected.severity <- -2 * sd(syndscores.main)} else if(input$severity == "severe"){selected.severity <- 2 * sd(syndscores.main)} else if(input$severity == "typical"){selected.severity <- 0}

  values <- matrix(NA, ncol = nrow(xyz) * 3, nrow = length(minage:nframes))
values.sev <- matrix(NA, ncol = nrow(xyz) * 3, nrow = 10)
  
  for(i in 1:nrow(values)){
    selected.synd <- factor(input$synd, levels = levels(d.meta.combined$Syndrome))
    selected.age <- i

    main.res <- matrix(t(d.registered$PCs[,1:num_pcs] %*% t(selected.severity * Snorm)), dim(d.registered$mshape)[1], dim(d.registered$mshape)[2])
    # selected.severity <- "mild"
    main.res <- matrix(t(d.registered$PCs[,1:num_pcs] %*% t(selected.severity * Snorm)), dim(d.registered$mshape)[1], dim(d.registered$mshape)[2])
    
    datamod <- ~ selected.sex + selected.age + selected.age^2 + selected.age^3 + selected.synd + selected.age:selected.synd
    predicted.shape <- predshape.lm(synd.lm.coefs, datamod, d.registered$PCs[,1:num_pcs], d.registered$mshape)
    
    tmp.mesh$vb[-4,] <- t(predicted.shape + main.res)
    final.shape <- t(vcgSmooth(tmp.mesh)$vb[-4,])
    
    shape.tmp <- array(final.shape[as.numeric(atlas$it), ], dim = c(nrow(xyz), 3, 1))
    values[i,] <- geomorph::two.d.array(shape.tmp)

  }

  control <- vertexControl(values = values,
                           vertices = rep(1:nrow(xyz), each = 3),
                           attributes = rep(c("x", "y", "z"), nrow(xyz)),
                           objid = objid)

 
  scene <- scene3d()
  # rgl.close()
  # return(list(scene, control, control2))
  return(list(scene, control))

})


  output$wdg <- renderRglwidget({
    rglwidget(morph_target()[[1]], controllers = c("control"))
  })


  output$control <- renderPlaywidget({
    playwidget("wdg", morph_target()[[2]],
               respondTo = "age")
  })

#    output$control2 <- renderPlaywidget({
#   playwidget("wdg", morph_target()[[3]],
#              respondTo = "severity")
# })

}



shinyApp(ui = ui, server = server)