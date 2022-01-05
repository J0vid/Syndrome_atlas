library(shiny)
library(rgl)
library(Morpho)


ui <- (fluidPage(
  selectInput("synd", label = "Syndrome", choices = sort(levels(d.meta.combined$Syndrome)), selected = "Achondroplasia"),
  selectInput("sex", label = "Sex", choices = c("Female", "Male")),
  sliderInput("age", label = "Age", min = 1, max = 70, value = 12, step = 1),
  sliderInput("severity", label = "Severity", min = 0, max = 4, value = 2, step = 1),
  # sliderInput("secret_slider", label = "", min = 0, max = 10, value = 0, step = 1),
  playwidgetOutput("control"),
  playwidgetOutput("control2"),
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
    updateSliderInput(session, "age", label = "Age", min = slider_min, max = slider_max)
    
    syndscores.main <- doutVar()[[2]] * 100
    #updateSliderInput(session, "severity", label = "Syndrome severity", min = round(-1*sd(syndscores.main)), max = round(1*sd(syndscores.main)), step = round(diff(range(syndscores.main))/10), value = 0)
    updateSliderInput(session, "severity", label = "Syndrome severity", min = 0, max = 4, step = 1, value = 2)
    # updateSliderInput(session, "secret_slider", min = 0, max = (slider_max - slider_min) * 5, step = 1, value = input$severity + (input$age*5))
  })
  
  selected.synd <- factor("Achondroplasia", levels = levels(d.meta.combined$Syndrome))
  selected.sex <-1
  selected.age <- 1 #round(abs(min(doutVar()[[1]])))
  
  datamod <- ~ selected.sex + selected.age + selected.age^2 + selected.age^3 + selected.synd + selected.age:selected.synd
  predicted.shape <- predshape.lm(synd.lm.coefs, datamod, d.registered$PCs[,1:num_pcs], d.registered$mshape)

  atlas$vb[-4,] <- t(predicted.shape)
  
  open3d()
  shade3d(vcgSmooth(atlas), aspect = "iso", col = "lightgrey", specular = 1)
  objid <- ids3d()$id
  xyz <- rgl.attrib(objid[1], "vertices")
  dim(xyz)
  print(currentSubscene3d())
  
  morph_target <- reactive({
  #age morph target####
  minage <- round(abs(min(doutVar()[[1]])))
  nframes <- round(max(doutVar()[[1]]))
  if(input$sex == "Female"){selected.sex <- 2
  } else if(input$sex == "Male"){selected.sex <- -1}

  values <- matrix(NA, ncol = nrow(xyz) * 3, nrow = length(minage:nframes))

  for(i in 1:nrow(values)){
    selected.synd <- factor(input$synd, levels = levels(d.meta.combined$Syndrome))
    selected.age <- i

    datamod <- ~ selected.sex + selected.age + selected.age^2 + selected.age^3 + selected.synd + selected.age:selected.synd
    predicted.shape <- predshape.lm(synd.lm.coefs, datamod, d.registered$PCs[,1:num_pcs], d.registered$mshape)

    shape.tmp <- array(predicted.shape[as.numeric(atlas$it), ], dim = c(nrow(xyz), 3, 1))
    values[i,] <- geomorph::two.d.array(shape.tmp)

  }

  control <- vertexControl(values = values,
                           vertices = rep(1:nrow(xyz), each = 3),
                           attributes = rep(c("x", "y", "z"), nrow(xyz)),
                           objid = objid)

  #severity morph target####
  S <- matrix(synd.lm.coefs[grepl(pattern = input$synd, rownames(synd.lm.coefs)),], nrow = 1, ncol = num_pcs)
  Snorm <- S/sqrt(sum(S^2))

  syndscores.main <- doutVar()[[2]]
  severity_range <- seq(round(-1*sd(syndscores.main), 2), round(1*sd(syndscores.main), 2), length.out = 10)

  selected.synd <- factor(input$synd, levels = levels(d.meta.combined$Syndrome))
  selected.age <- 0#input$age
  # syndonly_lm <- lm(d.registered$PCscores[,1:num_pcs] ~ d.meta.combined$Syndrome)$coef

  datamod <- ~ selected.sex + selected.age + selected.age^2 + selected.age^3 + selected.synd + selected.age:selected.synd
  predicted.shape <- predshape.lm(synd.lm.coefs, datamod, d.registered$PCs[,1:num_pcs], d.registered$mshape)

  values2 <- matrix(NA, ncol = nrow(xyz) * 3, nrow = 10)

  for(i in 1:nrow(values2)){

    selected.severity <- severity_range[i] #.01

    main.res <- matrix(t(d.registered$PCs[,1:num_pcs] %*% t(selected.severity * Snorm)), dim(d.registered$mshape)[1], dim(d.registered$mshape)[2])

    #add severity residuals
    severity_shape <-  predicted.shape + main.res

    shape.tmp <- array(severity_shape[as.numeric(atlas$it), ], dim = c(nrow(xyz), 3, 1))
    values2[i,] <- geomorph::two.d.array(shape.tmp)

  }

  # control2 <- vertexControl(values = values2,
  #                          vertices = rep(1:nrow(xyz), each = 3),
  #                          attributes = rep(c("x", "y", "z"), nrow(xyz)),
  #                          objid = objid
  #                          )


  scene <- scene3d()
  # rgl.close()
  # return(list(scene, control, control2))
  return(list(scene, control))

})


  output$wdg <- renderRglwidget({
    rglwidget(elementId="plot3drgl2", morph_target()[[1]], controllers = c("control"))
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