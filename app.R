#syndrome gestalts w/age flexibility and face submission/classification
library(shiny)
library(Morpho)
library(geomorph)
library(rgl)
library(shinycssloaders)
library(Jovid)
library(Rvcg)
library(visNetwork)
library(plotly)
library(ggrepel)
library(shinyBS)
library(grid)
library(mesheR)
library(RvtkStatismo)
options(shiny.maxRequestSize=300*1024^2)

# meta.lm, d.registered
 setwd("~/shiny/Syndrome_atlas/")
 test.atlas <- file2mesh("atlas.ply")
 atlas.points <- read.table("atlas_lm_5.txt")
 # atlas <- test.atlas
 a.man.lm <- as.matrix(atlas.points)
 atlas.lms <- as.matrix(atlas.points)
 # setwd("/srv/shiny-server/testing_ground/")
# save(atlas, d.meta.combined, front.face, PC.eigenvectors, synd.lm.coefs, synd.mshape, PC.scores, synd.mat, file = "data.Rdata")

load("data.Rdata")
load("modules_PCA.Rdata")
eye.index <- as.numeric(read.csv("eye_lms.csv", header = F)) + 1
load("FB2_texture_PCA.Rdata")
 
 predshape.lm <- function(fit, datamod, PC, mshape){
   dims <- dim(mshape)
   mat <- model.matrix(datamod)
   pred <- mat %*% fit
   names <- as.matrix(model.frame(datamod))
   names <- apply(names, 1, paste, collapse = "_")
   names <- gsub(" ", "", names)
   predPC <- t(PC %*% t(pred))
   out <- mshape + matrix(predPC, dims[1], dims[2])
   # dimnames(out)[[3]] <- names
   # rownames(pred) <- names
   return(out)
 }
 
 predtexture.lm <- function(fit, datamod, PC, mshape){
   dims <- dim(mshape)
   mat <- model.matrix(datamod)
   pred <- mat %*% fit
   names <- as.matrix(model.frame(datamod))
   names <- apply(names, 1, paste, collapse = "_")
   names <- gsub(" ", "", names)
   predPC <- t(PC %*% t(pred))
   out <- mshape + predPC
   
   # dimnames(out)[[3]] <- names
   # rownames(pred) <- names
   return(out)
 }

ui <- fluidPage(
  
  titlePanel(""),
  sidebarLayout(
    sidebarPanel(
      conditionalPanel(condition = "input.Atlas_tabs=='Gestalt'",
                       selectInput("synd", label = "Syndrome", choices = sort(levels(d.meta.combined$Syndrome)), selected = "Costello Syndrome"),
                       selectInput("sex", label = "Sex", choices = c("Female", "Male")),
      # sliderInput("age", label = "Age", min = 1, max = 50, step = 1, value = 12),
                      uiOutput("age_slider"),
                      uiOutput("severity_slider"),
                      bsButton("texture_help", label = "", icon = icon("question"), style = "default", size = "extra-small"),
                      selectInput("texture", label = "Texture type", choices = c("Generic", "Lightgrey", "Gestalt", "Generic + Gestalt")),
                      bsPopover(id = "texture_help", title = "Plot types",
                                content = paste0("Texture (mesh color) can be really helpful for understanding syndromic faces. By default, we provide a generic synthetic texture. If it does not fit the face features well, you can remove it by selecting lightgrey. You can also view syndromic texture or add syndromic aspects of texture to the generic mesh, but these feature are experimental and may not work well."
                                ),
                                placement = "right", 
                                trigger = "hover", 
                                options = list(container = "body")
                      ),
                      actionButton("Update", "Update Gestalt", icon("refresh"))),
      conditionalPanel(condition = "input.Atlas_tabs=='Comparisons'",
                      selectInput("reference", label = "Reference", choices = sort(levels(d.meta.combined$Syndrome)), selected = "Unaffected Unrelated"),
                      selectInput("synd_comp", label = "Syndrome", choices = sort(levels(d.meta.combined$Syndrome)), selected = "Costello Syndrome"),
                      selectInput("sex", label = "Sex", choices = c("Female", "Male")),
                      # uiOutput("age_slider"),
                      sliderInput("transparency", label = "Mesh transparency", min = 0, max = 1, step = .1, value = 1),
                      checkboxInput("displace", label = "Plot vectors?", value = F),
                      hr(),
                      bsButton("q1", label = "", icon = icon("question"), style = "default", size = "extra-small"),
                      selectInput("score_plot", label = "Plot type", choices = c("Similarity", "Raw score")),
                      bsPopover(id = "q1", title = "Plot types",
                                content = paste0("When syndromes are very severe, they may outscore the selected syndrome. In these cases, it may be more informative to look at syndromes that are similar in the magnitude of their effects."
                                ),
                                placement = "right", 
                                trigger = "hover", 
                                options = list(container = "body")
                      )),
      conditionalPanel(condition = "input.Atlas_tabs=='Submitted face'",
                       fileInput("file1", "",
                                 multiple = T,
                                 accept = c(".obj", ".ply", ".png", ".pp", ".csv")),
                       splitLayout(cellWidths = c("50%", "50%"),
                                   textInput("current_age", label = "Age (years)", value = "", placeholder = "ie. 9.5"),
                                   selectInput("current_sex", label = "Sex", selected = "Male", choices = c("Male", "Female"))
                       ),
                       checkboxGroupInput("submitted_heatmap", label = "Morphology", choices = "Compare to syndrome?"),
                       conditionalPanel(condition = "input.submitted_heatmap == 'Compare to syndrome?'",
                       selectInput("submitted_comp", label = "Syndrome", choices = sort(levels(d.meta.combined$Syndrome)), selected = "Costello Syndrome")),
                       ),
      width = 3, 
      tags$head(
        tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
      )),
    
    mainPanel(
      tabsetPanel(id = "Atlas_tabs",#img(src='uc_logo.jpg', align = "right", height = 70 * 1.15, width = 90 * 1.25),
        tabPanel("Gestalt", br(),
                 withSpinner(rglwidgetOutput("gestalt", width = "65vw", height="80vh"), type = 6, color = "#fca311")),
        tabPanel("Comparisons", br(),
                 withSpinner(rglwidgetOutput("comparison", width = "65vw", height="80vh"), type = 6, color = "#fca311"),
                 br(),
                 HTML("<center><font style=\"color: red; font-size: xx-large;\"> Bigger </font><font style=\"color: #fcfcfc; font-size: xx-large;\"> | </font><font style=\"color: lightgreen; font-size: xx-large;\"> Similar </font><font style=\"color: #fcfcfc; font-size: xx-large;\"> | </font><font style=\"color: blue; font-size: xx-large;\"> Smaller </font></center>"),
                 br(), plotOutput("morphospace")),
        tabPanel("Segments", br(), HTML("<p style=\"color:black;\">Select a facial partition using the dropdown menu or by clicking the bubbles directly.</p>"), 
                 visNetworkOutput("network", height="80vh")),
        tabPanel("Submitted face", br(),
                 withSpinner(rglwidgetOutput("submitted_face", width = "65vw", height="80vh"), type = 6, color = "#fca311"),
                 br(),
                 plotlyOutput("posterior_scree")),
        # tabPanel("Morphospace", br(), plotOutput("morphospace")),
        tabPanel("About", br(), HTML("<p style=\"color:black;\">This app aims to help clinical geneticists better understand the characteristic craniofacial features of various genetic syndromes. There are 3 sections to this app and here is my description of how they work. Here are the people that made this app possible.</p>"))#, 
        #includeHTML("~/shiny/shinyapps/Syndrome_gestalts/about_test.html"), 
        #includeCSS("~/shiny/shinyapps/Syndrome_gestalts/test.css"))
        
      )
    )
  )
)



server <- function(input, output) {
  # atlas <- file2mesh("whoami.ply", readcol = T)
  
  # atlas <- file2mesh("whoami2.ply", readcol = T)
  pred.mesh <- atlas 
  pred.mesh2 <- atlas
  texture.coefs <- lm(texture.pca$x[,1:100] ~ d.meta.combined$Sex + d.meta.combined$Age + d.meta.combined$Age^2 + d.meta.combined$Age^3 + d.meta.combined$Syndrome + d.meta.combined$Age:d.meta.combined$Syndrome)$coef
  # synd.lm.coefs <- lm(PC.scores ~ d.meta.combined$Sex + d.meta.combined$Age + d.meta.combined$Age^2 + d.meta.combined$Age^3 + d.meta.combined$Syndrome + d.meta.combined$Age:d.meta.combined$Syndrome)$coef
  
  #process reactive####
  outVar <- reactive({
    #get selected syndrome age range
    age.range <- d.meta.combined$Age[d.meta.combined$Syndrome == input$synd]
    
    #calculate syndrome severity scores for selected syndrome
    #calculate score for the main effect
    S <- synd.lm.coefs[grepl(pattern = input$synd, rownames(synd.lm.coefs)),][1,]
    Snorm <- S/sqrt(sum(S^2))
    syndscores.main <- PC.scores %*% Snorm
    
    return(list(age.range, syndscores.main[d.meta.combined$Syndrome == input$synd]))
  }) 
  
  doutVar <- outVar#debounce(outVar, 2000)
  
  output$age_slider <- renderUI({
    sliderInput('dynamic_age', 'Age', min = if(min(na.omit(doutVar()[[1]])) < 1){1} else{round(min(na.omit(doutVar()[[1]])))}, max = round(max(na.omit(doutVar()[[1]]))), value = 12)
    #Diagnostic slider: sliderInput('dynamic_age', 'Age', min = 1, max = 25, value = 12)
  })
  
  output$severity_slider<- renderUI({
  syndscores.main <- doutVar()[[2]]  
  #what's the current synd age range
  # age.range <- doutVar()[[1]]
  #filter to scores for only nearby age ranges
  # syndscores.main <- syndscores.main[age.range > (input$dynamic_age - 5) & age.range < (input$dynamic_age + 5)]
  sliderInput("severity", label = "Syndrome severity", min = round(-1*sd(syndscores.main), 2), max = round(1*sd(syndscores.main), 2), step = round(diff(range(syndscores.main))/10, 2), value = 0) 
  })
  
  debounced.inputs <- reactive({
    return(list(input$dynamic_age, input$severity)) 
  }) %>% debounce(1)
  
  
  pred.synd <- eventReactive(input$Update, {
    selected.synd <- factor(input$synd, levels = levels(d.meta.combined$Syndrome))
    if(input$sex == "Female"){selected.sex <-1
    } else if(input$sex == "Male"){selected.sex <- 0} 
    # selected.age <- input$age
    
    if(is.null(debounced.inputs()[[1]])){ selected.age <- 12 } else{
      selected.age <- debounced.inputs()[[1]]
    }
    
    # if(is.null(debounced.inputs()[[2]])){ selected.severity <- 0} else{
    selected.severity <- debounced.inputs()[[2]]
    # }
    # synd.lm.coefs.severity <- synd.lm.coefs
    #old severity code: synd.lm.coefs.severity[grepl(pattern = input$synd, rownames(synd.lm.coefs)),] <- synd.lm.coefs.severity[grepl(pattern = input$synd, rownames(synd.lm.coefs)),] + (input$severity/100) * synd.lm.coefs.severity[grepl(pattern = input$synd, rownames(synd.lm.coefs)),]
    datamod <- ~ selected.sex + selected.age + selected.age^2 + selected.age^3 + selected.synd + selected.age:selected.synd
    predicted.shape <- predshape.lm(synd.lm.coefs, datamod, PC.eigenvectors, synd.mshape)
    
    # S <- synd.lm.coefs[grepl(pattern = selected.synd, rownames(synd.lm.coefs)),][1,]
    # Snorm <- S/sqrt(sum(S^2))
    # syndscores.main <- PC.scores %*% Snorm
    # 
    # #solve for max in min syndromic severities
    # 
    # #if the mean is higher for the selected syndrome, then max is max
    # meanscore.ns <- mean(syndscores.main[d.meta.combined$Syndrome == "Unaffected Unrelated"])
    # meanscore.synd <- mean(syndscores.main[d.meta.combined$Syndrome == selected.synd])
    # 
    # 
    
    return(list(predicted.shape, selected.synd, selected.age, selected.severity))
    
  })
  
  #module selection####
  output$network <- renderVisNetwork({
    # 63 modules
    # nodes <- data.frame(id = 1:63, label = 1:63, shape = "image", image = "https://genopheno.ucalgary.ca/presentations/Child_health_data_science_Aponte/images/atlas.png", value = c(rep(120,21), rep(5,21), rep(2,21)))
    # edges <- data.frame(from = rep(1:31, each = 2), to = 2:63)#c(2,3,4,5,6,7,8,9,10,11,12,13,14,15))
    # 
    # # 31 modules
    # nodes <- data.frame(id = 1:31, label = 1:31, shape = "image", image = "https://genopheno.ucalgary.ca/presentations/Child_health_data_science_Aponte/images/atlas.png", value = c(rep(10,21), rep(5,10)))
    # edges <- data.frame(from = rep(1:15, each = 2), to = 2:31)#c(2,3,4,5,6,7,8,9,10,11,12,13,14,15))
    # 15 modules
    path_to_images <- "https://raw.githubusercontent.com/J0vid/Syndrome_atlas/master/mod_images/" #"https://genopheno.ucalgary.ca/presentations/Child_health_data_science_Aponte/mod_images/"
    module.names <- c("Whole face", "Nose", "Mandible/Sphenoid/Frontal", "Upper lip", "Nasal/Maxilla subset", "Cheek/Mandible", "Sphenoid/Frontal", "Nasolabial", "Philtrum", "Lateral Nasal", "Nose", "Cheek", "Chin/Mandible", "Frontal", "Orbital/Temporal/Sphenoid")
    module.names[1:8] <- c("Whole Face", "Nose", "Eyes", "Jaw", "Chin", "Mouth", "Cheeks", "Forehead")
    nodes <- data.frame(id = 1:8, label = module.names[1:8], title = module.names[1:8], shape = "circularImage", color = "#fca311", color.highlight = "#fca311", image = paste0(path_to_images, "mod", 1:8, ".png"))#, value = c(rep(130,7), rep(130,8)))
    edges <- data.frame(from = 1, to = 2:8)
    
    visNetwork(nodes, edges)  %>%
      visInteraction(hover = F,
                     dragNodes = T,
                     dragView = T) %>%
      visOptions(nodesIdSelection = list(enabled = TRUE, 
                                         selected = "1"),
                 highlightNearest = F
      ) %>% 
      visNodes(size = 50,
               shapeProperties = list(useBorderWithImage = TRUE),
               borderWidthSelected = 7,
               borderWidth = 5,
               font = list(color = "#fca311",
                           size = 17)) %>%
      visEdges(shadow = F,
               arrows =list(to = list(enabled = TRUE, scaleFactor = .5)),
               color = list(color = "#f0c173", highlight = "#fca311")) %>%
      visEvents(type = "once", startStabilizing = "function() {
            this.moveTo({scale:1.35})}") %>%
      visPhysics(stabilization = FALSE) %>%
      visLayout(randomSeed = 12)# %>%
      #visHierarchicalLayout()
    
  })
  
  
  output$shiny_return <- renderText({
    paste("Current node selection : ", input$network_selected)
  })
  
  
  syndrome_comps <- reactive({
    
    selected.sex <- 1.5
    selected.age <- 10
    
    selected.synd <- factor(input$synd_comp, levels = levels(d.meta.combined$Syndrome))
    
    # datamod <- ~ selected.sex + selected.age + selected.age^2 + selected.age^3 + selected.synd + selected.age:selected.synd
    
    # predicted.shape <- predshape.lm(synd.lm.coefs, datamod, PC.eigenvectors, synd.mshape)
    
    #synd
    pred.mesh$vb[-4,] <-  t(synd.mat[,,unique(d.meta.combined$Syndrome) == input$synd_comp])
    rigid.points <- sample(1:27000, 100)
    reduce.mesh <-  sample(1:27000, 1000)
    test <- rotmesh.onto(pred.mesh, t(pred.mesh$vb[1:3,rigid.points]), t(atlas$vb[1:3,rigid.points]))
    # test <- rmVertex(test$mesh, index = reduce.mesh, keep = T)
    #ref
    pred.mesh$vb[-4,] <-  t(synd.mat[,,unique(d.meta.combined$Syndrome) == input$reference])
    test2 <- rotmesh.onto(pred.mesh, t(pred.mesh$vb[1:3,rigid.points]), t(atlas$vb[1:3,rigid.points]))
    # test2 <- rmVertex(test2$mesh, index = reduce.mesh, keep = T)
    
    #calculate syndrome severity scores for selected syndrome
    #calculate score for the whole face
    S <- synd.lm.coefs[grepl(pattern = input$synd_comp, rownames(synd.lm.coefs)),][1,]
    Snorm <- S/sqrt(sum(S^2))
    syndscores.main <- PC.scores %*% Snorm
    
    syndscores.df <- data.frame(Syndrome = d.meta.combined$Syndrome, face.score = syndscores.main)
    syndscores.wholeface <- syndscores.df%>%
      group_by(Syndrome) %>%
      summarise(face_score = mean(face.score))
    
    # calculate score for the selected subregion
    if(is.null(input$network_selected)) selected.node <- 1 else if(input$network_selected == ""){
      selected.node <- 1} else{
        selected.node <- as.numeric(input$network_selected)}

    if(selected.node > 1){
    subregion.coefs <- manova(get(paste0(tolower(colnames(modules)[selected.node]), ".pca"))$x ~ d.meta.combined$Sex + d.meta.combined$Age + d.meta.combined$Age^2 + d.meta.combined$Age^3 + d.meta.combined$Syndrome + d.meta.combined$Age:d.meta.combined$Syndrome)$coef
    S <- subregion.coefs[grepl(pattern = input$synd_comp, rownames(subregion.coefs)),][1,]
    print(rownames(subregion.coefs)[grepl(pattern = input$synd_comp, rownames(subregion.coefs))])
    Snorm <- S/sqrt(sum(S^2))
    syndscores.main <- get(paste0(tolower(colnames(modules)[selected.node]), ".pca"))$x %*% Snorm

    syndscores.df <- data.frame(Syndrome = d.meta.combined$Syndrome, module.score = syndscores.main)
    syndscores.module <- syndscores.df%>%
      group_by(Syndrome) %>%
      summarise(face_score = mean(module.score))
    } else{syndscores.module <- syndscores.wholeface}
    
    
    return(list(test2$mesh, test$mesh, syndscores.wholeface, syndscores.module))
  })

  
  output$gestalt <- renderRglwidget({
    pdf(NULL)
    dev.off()
    
    clear3d()
    #visualize full model estimates####
     if(input$sex == "Female"){selected.sex <-1
     } else if(input$sex == "Male"){selected.sex <- 0} 
    # selected.age <- input$age
    
    # if(is.null(debounced.inputs()[[1]])){ selected.age <- 12 } else{
    selected.age <- pred.synd()[[3]]
    # }
    # selected.sex <- 1
    # selected.age <- 6
    
    
    # if(is.null(input$synd)){ selected.synd <- "Costello Syndrome" } else{
      selected.synd <- pred.synd()[[2]]
    # }
    
    datamod <- ~ selected.sex + selected.age + selected.age^2 + selected.age^3 + selected.synd + selected.age:selected.synd
    #calculate score for the main effect
    S <- synd.lm.coefs[grepl(pattern = selected.synd, rownames(synd.lm.coefs)),][1,]
    Snorm <- S/sqrt(sum(S^2))
    # syndscores.main <- PC.scores %*% Snorm
    # 
    # #solve for max in min syndromic severities
    # 
    # #if the mean is higher for the selected syndrome, then max is max
    # meanscore.ns <- mean(syndscores.main[d.meta.combined$Syndrome == "Unaffected Unrelated"])
    # meanscore.synd <- mean(syndscores.main[d.meta.combined$Syndrome == selected.synd])
    
    # if(is.null(pred.synd()[[4]])){ selected.severity <- meanscore.synd } else{
      selected.severity <- pred.synd()[[4]]
    # }
    print(selected.severity)
      
    # if(meanscore.synd > meanscore.ns){main.res <- matrix(t(PC.eigenvectors %*% (selected.severity * Snorm)), dim(synd.mshape)[1], dim(synd.mshape)[2])} else{main.res <- matrix(t(PC.eigenvectors %*% (selected.severity * Snorm)), dim(synd.mshape)[1], dim(synd.mshape)[2])} 
      main.res <- matrix(t(PC.eigenvectors %*% (selected.severity * Snorm)), dim(synd.mshape)[1], dim(synd.mshape)[2])
  
    #add severity residuals
    pred.mesh$vb[-4,] <-  t(pred.synd()[[1]] + main.res) #+ interaction.res)
    rigid.points <- sample(1:27000, 100)
    test <- rotmesh.onto(pred.mesh, t(pred.mesh$vb[1:3,rigid.points]), t(atlas$vb[1:3,rigid.points]))
    
    # nonsynd.mean <- modules[[selected.node]]
    # nonsynd.mean$vb[1:3,] <- synd.mat[, segmentation[,selected.node] == T,d.meta.combined$Syndrome == input$reference]
    # 
    # 
    # synd.mean <- modules[[selected.node]]
    # synd.mean$vb[1:3,] <- synd.mat[, segmentation[,selected.node] == T, d.meta.combined$Syndrome == input$synd]
    
    # par3d(userMatrix = matrix(c(.998,-.005,.0613,0,.0021,.999,.045,0,-.061,-.045,.997,0,0,0,0,1),ncol = 4,nrow = 4))
    # par3d(userMatrix = front.face)

    # 
    # 
    #   par3d(userMatrix = front.face)
    #   
    #   synd.ages <- data.frame(age = d.meta.combined$Age[which(d.meta.combined$Syndrome == input$synd)], index = which(d.meta.combined$Syndrome == input$synd))
    #   
    #   synd.textures <- synd.ages$index[sort(abs(synd.ages$age - if(is.null(input$dynamic_age)){12}else{input$dynamic_age}), index.return = T)$ix][1:2]
    #   
    #   print(d.meta.combined[synd.textures,])
      
      # gen.scan <- generate.face(texture.pca, num.pcs = 2000, mesh = test$mesh, draws = synd.textures)
      # plot3d(gen.scan, aspect = "iso", axes = F, box = F, xlab = "", ylab = "", zlab = "", lit = F, specular = 1)
    
    if(input$texture == "Generic"){ 
      mesh.color <- test$mesh$material$color
      lit <- F}
    if(input$texture == "Lightgrey"){
      mesh.color <- "lightgrey"
      lit <- T}
    if(input$texture == "Gestalt"){ 
      #texture simulation####
      datamod <- ~ selected.sex + selected.age + selected.age^2 + selected.age^3 + selected.synd + selected.age:selected.synd
      predicted.texture <- predtexture.lm(texture.coefs, datamod, texture.pca$rotation[,1:100], texture.pca$center)
      
      predicted.texture3d <- row2array3d(predicted.texture)[,,1]
      
      final.texture <-  (predicted.texture3d)
      
      #scale values
      maxs <- apply(final.texture, 2, max)
      mins <- apply(final.texture, 2, min)
      additive.texture <- scale(final.texture, center = mins, scale = maxs - mins)
      hex.mean <- rgb(additive.texture, maxColorValue = 1)
      test$mesh$material$color[-eye.index] <- hex.mean[-eye.index]
      mesh.color <- test$mesh$material$color
      lit <- F}
    if(input$texture == "Generic + Gestalt"){ 

      #texture simulation####
      datamod <- ~ selected.sex + selected.age + selected.age^2 + selected.age^3 + selected.synd + selected.age:selected.synd
      predicted.texture <- predtexture.lm(texture.coefs, datamod, texture.pca$rotation[,1:100], texture.pca$center)
      
      predicted.texture3d <- row2array3d(predicted.texture)[,,1]
      
      starting.texture <- matrix(NA, nrow = length(atlas$material$color), ncol = 3)
      for(j in 1:length(atlas$material$color)) starting.texture[j,] <- col2rgb(atlas$material$color[j])
      
      final.texture <- starting.texture
      # final.texture[brow.index,] <- starting.texture[brow.index,] + predicted.texture3d[brow.index,]
      # final.texture[-brow.index,1] <- final.texture[-brow.index,1]  + colMeans(predicted.texture3d[brow.index,])[1]
      # final.texture[-brow.index,2] <- final.texture[-brow.index,2]  + colMeans(predicted.texture3d[brow.index,])[2]
      # final.texture[-brow.index,3] <- final.texture[-brow.index,3]  + colMeans(predicted.texture3d[brow.index,])[3]
      final.texture <-  2*(predicted.texture3d) + starting.texture
      
      
      #scale values
      maxs <- apply(final.texture, 2, max)
      mins <- apply(final.texture, 2, min)
      additive.texture <- scale(final.texture, center = mins, scale = maxs - mins)
      hex.mean <- rgb(additive.texture, maxColorValue = 1)
      test$mesh$material$color[-eye.index] <- hex.mean[-eye.index]
      mesh.color <- test$mesh$material$color
      lit <- F}
    
      plot3d(test$mesh, aspect = "iso", axes = F, box = F, xlab = "", ylab = "", zlab = "", lit = lit, specular = 1, col = mesh.color)
      # plot3d(atlas, aspect = "iso", axes = F, box = F, xlab = "", ylab = "", zlab = "", lit = F)
      par3d(zoom = .83)
      # bg3d("#e5e5e5")
      bg3d("#fcfcfc")
      
    par3d(userMatrix = front.face)
    rglwidget()
    
  })
  
  output$comparison <- renderRglwidget({
    pdf(NULL)
    dev.off()
    
    clear3d()
    
    if(is.null(input$network_selected)) selected.node <- 1 else if(input$network_selected == ""){
      selected.node <- 1} else{
        selected.node <- as.numeric(input$network_selected)} 
    
    
    if(input$displace & input$synd_comp != input$reference & selected.node == 1){
      mD.synd <- meshDist(syndrome_comps()[[1]], syndrome_comps()[[2]], plot = F, scaleramp = F, displace = input$displace, alpha = input$transparency)
      a <- render(mD.synd, displace = T, alpha = input$transparency)
    } else if(input$displace == F & input$synd_comp != input$reference & selected.node == 1){
      
      mD.synd <- meshDist(syndrome_comps()[[1]], syndrome_comps()[[2]], plot = F, scaleramp = F, displace = input$displace, alpha = input$transparency)
      a <- render(mD.synd, displace = F, alpha = input$transparency)
      
    }
    
    if(selected.node > 1){
     
      mD.synd <- meshDist(syndrome_comps()[[1]], syndrome_comps()[[2]], plot = F, scaleramp = F, displace = input$displace, alpha = input$transparency)
      mD.synd <- rmVertex(mD.synd$colMesh, index = which(modules[,selected.node] == 1), keep = T)
      # a <- render(mD.synd, displace = F, alpha = input$transparency)
      plot3d(mD.synd, aspect = "iso", axes = F, box = F, xlab = "", ylab = "", zlab = "", lit = T, specular = 1)
      shade3d(syndrome_comps()[[1]], col = "lightgrey", alpha = .2, specular = 1)
    }
    
    par3d(userMatrix = front.face, zoom = .83)
    aspect3d("iso")
    rglwidget()
    
  })
  
  output$morphospace <- renderPlot({

    if(is.null(input$network_selected)) selected.node <- 1 else if(input$network_selected == ""){
      selected.node <- 1} else{
        selected.node <- as.numeric(input$network_selected)}

    module.names <- colnames(modules)
#
#     full.morpho <- procdist.array(synd.mat, synd.mat[,,d.meta.combined$Syndrome == input$synd_comp])
#
#     sub.morpho <- procdist.array(synd.mat[modules[,selected.node] == 1,,], synd.mat[modules[,selected.node] == 1,, d.meta.combined$Syndrome == input$synd_comp])
    if(input$score_plot == "Raw score"){
        morpho.df <- data.frame(full = syndrome_comps()[[3]]$face_score, sub = syndrome_comps()[[4]]$face_score, Syndrome = levels(d.meta.combined$Syndrome))
        #sort to get which synds to highlight
        highlight_df <- morpho.df[sort(morpho.df$sub, index.return = T)$ix,]
        highlighted_points <- geom_point(data = highlight_df[85:89,], aes(y = full, x = sub), color = "black", fill = "#a6192e", shape = 21, size = 3) 
        tmp_xlab <- xlab(paste0(module.names[selected.node], " score for ", input$synd_comp))
        tmp_ylab <- ylab(paste0("Whole face score for ", input$synd_comp))
        tmp_labels <- geom_label_repel(data = highlight_df[c(1:5, 85:89),], aes(y = full, x = sub, label = Syndrome),
                                       force = 2,
                                       size = 4,
                                       color = c(rep("black", 5), rep("#a6192e", 5)),
                                       fill = "white")
    }
    
    if(input$score_plot == "Similarity"){
      full_similarity <- sqrt((syndrome_comps()[[3]]$face_score - syndrome_comps()[[3]]$face_score[levels(d.meta.combined$Syndrome) == input$synd_comp])^2)
      module_similarity <- sqrt((syndrome_comps()[[4]]$face_score - syndrome_comps()[[4]]$face_score[levels(d.meta.combined$Syndrome) == input$synd_comp])^2)
      
      morpho.df <- data.frame(full = full_similarity, sub = module_similarity, Syndrome = levels(d.meta.combined$Syndrome))
      #sort to get which synds to highlight
      highlight_df <- morpho.df[sort(morpho.df$sub, index.return = T)$ix,]
      highlighted_points <- geom_point(data = highlight_df[1:5,], aes(y = full, x = sub), color = "black", fill = "#a6192e", shape = 21, size = 3) 
      tmp_xlab <- xlab(paste0(module.names[selected.node], " similarity to ", input$synd_comp))
      tmp_ylab <- ylab(paste0("Whole face similarity to ", input$synd_comp))
      tmp_labels <- geom_label_repel(data = highlight_df[c(1:5, 85:89),], aes(y = full, x = sub, label = Syndrome),
                       force = 2,
                       size = 4,
                       color = c(rep("#a6192e", 5), rep("black", 5)),
                       fill = "white")
    }
    #to do
    #highlight 5 closest syndromes with names
    #maybe all points should be names, with 5 highlighted
    p <- ggplot(morpho.df) +
      geom_point(aes(y = full, x = sub), color = "black") +
      highlighted_points +
      tmp_labels +
      tmp_ylab +
      tmp_xlab +
      # ylim(-.25,3.5) +
      theme_bw() +
      theme(plot.background = element_rect(fill = "#fcfcfc"),
            axis.text = element_text(color = "black"),
            axis.title = element_text(color = "black"),
            panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                            colour = "grey65"),
            panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                            colour = "grey65"), 
            panel.background = element_blank(),
            plot.margin = margin(20, 20, 20, 20)
      )
    print(selected.node)
    g <- ggplotGrob(p)
    bg <- g$grobs[[1]]
    round_bg <- roundrectGrob(x=bg$x, y=bg$y, width=bg$width, height=bg$height,
                              r=unit(0.1, "snpc"),
                              just=bg$just, name=bg$name, gp=bg$gp, vp=bg$vp)
    g$grobs[[1]] <- round_bg
    
    plot(g)
    

  })
  
  
  mesh.et.lms <- eventReactive(input$file1, {
    print(input$file1)
    # file.mesh <- file2mesh("/mnt/Hallgrimsson/Users/Jovid/Andlit/baking/DavidAponte.ply")
    # file.mesh <- file2mesh("~/shiny/Syndrome_atlas/GCS_1.ply")
    file.mesh <- file2mesh(input$file1$datapath[grepl("*.ply", input$file1$datapath)])
    # file.name <- substr(input$file1$name, start = 1, stop = nchar(input$file1$name) - 4)
    file.name <- "test"
    # file.lms <- read.mpp("/mnt/Hallgrimsson/Users/Jovid/Andlit/baking/DavidAponte_picked_points.pp")
    file.lms <- read.mpp(input$file1$datapath[grepl("*.pp", input$file1$datapath)])
    # file.lms <- read.mpp("/mnt/Hallgrimsson/Users/Jovid/FB2_ML/sylvia_picked_points.pp")
    # file.lms <- read.mpp("~/shiny/Syndrome_atlas/GCS_1_picked_points.pp")
    
    tmp.fb <- rotmesh.onto(file.mesh, refmat = file.lms, tarmat = a.man.lm, scale = T)
    gp.fb <- tmp.fb$yrot
    
    tmp.fb <- tmp.fb$mesh
    
    Kernels <- IsoKernel(0.1,atlas)
    mymod <- statismoModelFromRepresenter(test.atlas,kernel=Kernels,ncomp=100)
    postDef <- posteriorDeform(mymod, tmp.fb, modlm = atlas.lms, samplenum = 1000)
    
    withProgress(message = 'Fitting mesh', value = 0, {
      
      for (i in 1:7) {
        # Increment the progress bar, and update the detail text.
        incProgress(1/7, detail = c(rep("Shaping, thinking, doing...", 4), "Coffee break?", "Non-rigid deformation", "Last few measurements")[i])
        
        if(i < 4) postDef <- posteriorDeform(mymod, tmp.fb, modlm = a.man.lm, tarlm = gp.fb, samplenum = 1000, reference = postDef)
        
        if(i > 4){
          postDefFinal <- postDef
          postDefFinal <- posteriorDeform(mymod, tmp.fb, modlm=atlas.lms, samplenum = 3000, reference = postDefFinal, deform = T, distance = 3)
        }
      }
      
    })
    
    # atlas.map <- c(6532,6632,4355,4188, 5368, 5629, 5373, 5452,5156, 5404, 4845, 6740, 6715, 6537, 3660, 2527, 2210, 4358, 2303, 4625, 3950, 4516, 5328, 4144, 5709, 821, 617, 1765, 6306, 1693, 1195, 2677, 1778, 3845, 4006, 5332, 2537, 6500, 3089, 2488, 3441, 2685, 2233, 4283, 2382, 4881, 4115, 4603, 5258, 3974, 1086, 807, 649, 1807, 1804, 6068, 5768, 2531, 1615, 4030, 4157, 5166, 2426, 3208, 6557)
    
    return(list(postDefFinal, file.name, file.mesh))
  })
  
  
  output$submitted_face <- renderRglwidget({
    pdf(NULL)
    dev.off()
    
    par3d(userMatrix = diag(4), zoom = .75)
    bg3d(color = "#e5e5e5")
    plot3d(mesh.et.lms()[[1]], col = "lightgrey", axes = F, specular = 1, xlab = "", ylab = "", zlab = "", aspect = "iso")  
    rglwidget()
    
    #if heatmap, register face, make comp same age?, compare
    
  })
  
  output$posterior_scree <- renderPlotly({
    sample1k <- sample(1:27903, 1000)
    #register landmarks to the space
    registered.mesh <- t(rotmesh.onto(mesh.et.lms()[[1]], t(mesh.et.lms()[[1]]$vb[-4, sample1k]), FB2.mean[sample1k,], scale = T)$mesh$vb[-4,])
    #testing: 
    # test.mesh <- file2mesh("/mnt/Hallgrimsson/Users/Jovid/Andlit/baking/davida_registered.ply")
    # registered.mesh <- t(rotmesh.onto(test.mesh, t(test.mesh$vb[-4, sample1k]), FB2.mean[sample1k,], scale = T)$mesh$vb[-4,])
    #project mesh lms into FB2 PC space
    projected.mesh <- getPCscores(registered.mesh, FB2.vectors, FB2.mean)
    
    #remove the fitted values, add residuals to the coefficients
    # age <- 31
    
    if(input$current_age == ""){ age <- 20
    } else{age <- as.numeric(input$current_age)}
    sex <- as.numeric(input$current_sex == "Male")
  
    fitted.ind <- age * FB2.coefs[3,] + age^2 * FB2.coefs[4,] + age^3 * FB2.coefs[5,] + FB2.coefs[1,] + sex * FB2.coefs[2,]

    projected.residuals <- projected.mesh - fitted.ind
    
    #classify individual's scores using the model
    colnames(projected.residuals) <- colnames(FB2.scores)
    
    posterior.distribution <- predict(HDRDA.mod, newdata = rbind(projected.residuals, projected.residuals), type = "prob")$post[1,]
    
    posterior.distribution <- sort(posterior.distribution, decreasing = T)
    
    #used to be part of plot.df: ID = as.factor(1:10), 
    plot.df <- data.frame(Probs = round(as.numeric(posterior.distribution[1:10]), digits = 4), Syndrome = as.factor(names(posterior.distribution[1:10])))
    plot.df$Syndrome <- as.character(plot.df$Syndrome)
    plot.df$Syndrome[plot.df$Syndrome == "Unrelated Unaffected"] <- "Non-syndromic"
    
    plot_ly(data = plot.df, x = ~Syndrome, y = ~Probs, type = "bar", color = I("grey"), hoverinfo = paste0("Syndrome: ", "x", "<br>", "Probability: ", "y")) %>%
      layout(xaxis = list(tickvals = gsub("_", " ", plot.df$Syndrome), tickangle = 45, ticktext = c(Syndrome = plot.df$Syndrome, Probability = plot.df$Probs), title = "<b>Syndrome</b>"),
             yaxis = list(title = "<b>Class probability</b>"),
             paper_bgcolor='#e5e5e5',
             margin = list(b = 125, l = 50, r = 100)
      )
  })
  
}

shinyApp(ui = ui, server = server)











