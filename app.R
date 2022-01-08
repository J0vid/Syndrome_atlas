#syndrome gestalts w/age flexibility
library(shiny)
library(Morpho)
# library(geomorph)
library(rgl)
library(shinycssloaders)
library(Jovid)
library(Rvcg)
library(visNetwork)
library(plotly)
library(ggrepel)
library(shinyBS)
library(grid)

# meta.lm, d.registered
 setwd("~/shiny/shinyapps/Syndrome_model/")
 # setwd("/srv/shiny-server/testing_ground/")
# save(atlas, d.meta.combined, front.face, PC.eigenvectors, synd.lm.coefs, synd.mshape, PC.scores, synd.mat, file = "data.Rdata")
 load("data.Rdata")
 load("modules_PCA.Rdata")

 eye.index <- as.numeric(read.csv("~/shiny/shinyapps/Syndrome_model/lm_indices/eye_small.csv", header = F)) +1#as.numeric(read.csv("~/Desktop/eye_lms.csv", header = F)) +1

# eye.index <- as.numeric(read.csv("~/Desktop/eye_lms.csv", header = F)) +1


 # load("~/shiny/shinyapps/Syndrome_model/FB2_texture_PCA.Rdata")
 
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
 
 render <- function(x,...) UseMethod("render")
 
 render.meshDist <- function(x,from=NULL,to=NULL,steps=NULL,ceiling=NULL,uprange=NULL,tol=NULL,tolcol=NULL,rampcolors=NULL,NAcol=NULL,displace=FALSE,shade=TRUE,sign=NULL,add=FALSE,scaleramp=NULL,titleplot="Distance in mm",...) {
   clost <- x$clost
   dists <- x$dists
   distsOrig <- dists
   colorall <- x$cols
   colramp <- x$colramp
   params <- x$params
   distqual <- x$distqual
   if (!is.null(tolcol))
     tolcol <- colorRampPalette(tolcol)(1)
   if (!add) {
     if (rgl.cur() !=0)
       rgl.clear()
   }
   if (!is.null(from) || !is.null(to) || !is.null(uprange) ||  !is.null(tol)  ||  !is.null(sign) || !is.null(steps) || !is.null(rampcolors) || !is.null(NAcol) || !is.null(tolcol) || !is.null(scaleramp)) {
     neg=FALSE
     colMesh <- x$colMesh
     if(is.null(steps))
       steps <- x$params$steps
     if (is.null(rampcolors))
       rampcolors <- x$params$rampcolors
     if (is.null(NAcol))
       NAcol <- x$params$NAcol
     if (is.null(tolcol))
       tolcol <- x$params$tolcol
     if (is.null(tol))
       tol <- x$params$tol
     if(is.null(sign))
       sign <- x$params$sign
     if (!sign) {
       distsOrig <- dists
       dists <- abs(dists)
     }
     if(is.null(ceiling))
       ceiling <- x$params$ceiling
     if(is.null(uprange))
       uprange <- x$params$uprange
     
     if (is.null(from)) {
       mindist <- min(dists)
       if (sign && mindist < 0 ) {
         from <- quantile(dists,probs=(1-uprange)) 
         neg <- TRUE            
       } else {
         from <- 0
       }             
     }
     if (is.null(scaleramp))
       scaleramp <- x$params$scaleramp
     
     if (from < 0)
       neg <- TRUE
     if (is.null(to))
       to <- quantile(dists,probs=uprange)    
     if(ceiling)
       to <- ceiling(to)
     
     to <- to+1e-10
     #ramp <- blue2green2red(maxseq*2)
     ramp <- colorRampPalette(rampcolors)(steps-1)
     colseq <- seq(from=from,to=to,length.out=steps)
     coldif <- colseq[2]-colseq[1]
     if (neg && sign) {
       
       negseq <- length(which(colseq<0))
       poseq <- steps-negseq
       maxseq <- max(c(negseq,poseq))
       if (scaleramp) {
         ramp <- colorRampPalette(rampcolors)(maxseq*2)
         ramp <- ramp[c(maxseq-negseq+1):(maxseq+poseq)]
         
       }
       else
         ramp <- colorRampPalette(rampcolors)(steps-1)
       distqual <- ceiling(((dists+abs(from))/coldif)+1e-14)
       #distqual[which(distqual < 1)] <- steps+10
     } else if (from > 0) {
       distqual <- ceiling(((dists-from)/coldif)+1e-14)
     } else {
       distqual <- ceiling((dists/coldif)+1e-14)
     }
     distqual[which(distqual < 1)] <- steps+10
     colorall <- ramp[distqual]
     if (!is.null(tol)) {
       if ( length(tol) < 2 ) {
         if (sign) {
           tol <- c(-tol,tol)
         } else {
           tol <- c(0,tol)
         }
       }
       good <- which(abs(dists) < tol[2])
       colorall[good] <- tolcol
     }
     colfun <- function(x){x <- colorall[x];return(x)}
     colMesh$material$color <- colorall
     colMesh$material$color[is.na(colMesh$material$color)] <- NAcol
     #colMesh$material$color <- matrix(colfun(colMesh$it),dim(colMesh$it))
     colramp <- list(1,colseq, matrix(data=colseq, ncol=length(colseq),nrow=1),col=ramp,useRaster=T,ylab=titleplot,xlab="",xaxt="n")
   } else {
     if (is.null(tol))
       tol <- x$params$tol
     colramp <- x$colramp
     colMesh <- x$colMesh
   }
   if (is.null(tolcol))
     tolcol <- x$params$tolcol
   
   if (shade)
     shade3d(vcgUpdateNormals(colMesh),specular="black",...)
   if (displace) {
     dismesh <- colMesh
     vl <- dim(colMesh$vb)[2]
     dismesh$vb <- cbind(colMesh$vb,rbind(clost,1))
     dismesh$it <- rbind(1:vl,1:vl,(1:vl)+vl)
     dismesh$material$color <- colorall
     dismesh$normals <- cbind(dismesh$normals, dismesh$normals)
     wire3d(dismesh,lit=FALSE)
   }
   diffo <- ((colramp[[2]][2]-colramp[[2]][1])/2)
   image(colramp[[1]],colramp[[2]][-1]-diffo,t(colramp[[3]][1,-1])-diffo,col=colramp[[4]],useRaster=TRUE,ylab=titleplot,xlab="",xaxt="n")
   if (!is.null(tol)) {
     if (sum(abs(tol)) != 0)
       image(colramp[[1]],c(tol[1],tol[2]),matrix(c(tol[1],tol[2]),1,1),col=tolcol,useRaster=TRUE,add=TRUE)
   }
   params <- list(steps=steps,from=from,to=to,uprange=uprange,ceiling=ceiling,sign=sign,tol=tol,rampcolors=rampcolors,NAcol=NAcol,tolcol=tolcol)
   out <- list(colMesh=colMesh,dists=distsOrig,cols=colorall,colramp=colramp,params=params,distqual=distqual,clost=clost)
   
   class(out) <- "meshDist"
   invisible(out)
 }

ui <- fluidPage(
  
  titlePanel(""),
  sidebarLayout(
    sidebarPanel(
      conditionalPanel(condition = "input.Atlas_tabs=='Gestalt'",
                       selectInput("synd", label = "Syndrome", choices = sort(levels(d.meta.combined$Syndrome)), selected = "Costello Syndrome"),
                       selectInput("sex", label = "Sex", choices = c("Female", "Male")),
                      uiOutput("age_slider"),
                      uiOutput("severity_slider"),
                      bsButton("texture_help", label = "", icon = icon("question"), style = "default", size = "extra-small"),
                      selectInput("texture", label = "Texture type", choices = c("Generic", "Lightgrey", "Gestalt", "Generic + Gestalt"), selected = "Lightgrey"),
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
                      selectInput("comp_sex", label = "Sex", choices = c("Female", "Male")),
                      sliderInput("comp_age", label = "Age", min = 1, max = 40, step = 1, value = 12),
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
                      ),
                      actionButton("update_comp", "Update Comparison", icon("refresh"))),
      width = 3, 
      tags$head(
        tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
      )
    ),
    
    mainPanel(
      tabsetPanel(id = "Atlas_tabs",#img(src='uc_logo.jpg', align = "right", height = 70 * 1.15, width = 90 * 1.25),
        tabPanel("Gestalt", br(),
                 withSpinner(rglwidgetOutput("gestalt", width = "65vw", height="80vh"), type = 6, color = "#fca311")),
        tabPanel("Comparisons", br(),
                 withSpinner(rglwidgetOutput("comparison", width = "65vw", height="80vh"), type = 6, color = "#fca311"),
                 br(),
                 HTML("<center><font style=\"color: red; font-size: xx-large;\"> Bigger </font><font style=\"color: #fcfcfc; font-size: xx-large;\"> | </font><font style=\"color: lightgreen; font-size: xx-large;\"> Similar </font><font style=\"color: #fcfcfc; font-size: xx-large;\"> | </font><font style=\"color: blue; font-size: xx-large;\"> Smaller </font></center>"),
                 br(), plotOutput("morphospace"),
                 br(), 
                 HTML("<h3 style=\"color:black; text-align:center\">Select a facial partition by clicking the bubbles directly and then update the comparison.</h2>"),
                 visNetworkOutput("network", height="80vh")),
        # tabPanel("Segments", br(), HTML("<p style=\"color:black;\">Select a facial partition using the dropdown menu or by clicking the bubbles directly.</p>"), 
        #          visNetworkOutput("network", height="80vh")),
        # tabPanel("Morphospace", br(), plotOutput("morphospace")),
        tabPanel("Submitted face", br(),
                 withSpinner(rglwidgetOutput("submitted_face", width = "65vw", height="80vh"), type = 6, color = "#fca311"),
                 br(),
                 plotlyOutput("posterior_scree")),
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
  # texture.coefs <- lm(texture.pca$x[,1:100] ~ d.meta.combined$Sex + d.meta.combined$Age + d.meta.combined$Age^2 + d.meta.combined$Age^3 + d.meta.combined$Syndrome + d.meta.combined$Age:d.meta.combined$Syndrome)$coef
  # synd.lm.coefs <- lm(PC.scores ~ d.meta.combined$Sex + d.meta.combined$Age + d.meta.combined$Age^2 + d.meta.combined$Age^3 + d.meta.combined$Syndrome + d.meta.combined$Age:d.meta.combined$Syndrome)$coef
  
  #synd reactive####
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
    selected.synd <- input$synd#factor(input$synd, levels = levels(d.meta.combined$Syndrome))
    selected.sex <- input$sex
    # selected.age <- input$age
    
    if(is.null(debounced.inputs()[[1]])){ selected.age <- 12 } else{
      selected.age <- debounced.inputs()[[1]]
    }
    
    selected.severity <- debounced.inputs()[[2]]
 
    raw_api_res <- httr::GET(url = paste0("http://localhost:6352", "/predshape"),
                             query = list(selected.sex = selected.sex, selected.age = selected.age, selected.synd = selected.synd),
                             encode = "json")
    
    predicted.shape  <- jsonlite::fromJSON(httr::content(raw_api_res, "text"))/1e6
    
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
    path_to_images <- "https://raw.githubusercontent.com/J0vid/Syndrome_atlas/main/mod_images/" #"https://genopheno.ucalgary.ca/presentations/Child_health_data_science_Aponte/mod_images/"
    module.names <- c("Whole face", "Nose", "Mandible/Sphenoid/Frontal", "Upper lip", "Nasal/Maxilla subset", "Cheek/Mandible", "Sphenoid/Frontal", "Nasolabial", "Philtrum", "Lateral Nasal", "Nose", "Cheek", "Chin/Mandible", "Frontal", "Orbital/Temporal/Sphenoid")
    module.names[1:7] <- c("face", "posterior_mandible", "nose", "anterior_mandible", "brow", "zygomatic", "premaxilla") #c("Whole Face", "Nose", "Eyes", "Jaw", "Chin", "Mouth", "Cheeks")
    nodes <- data.frame(id = 1:7, label = module.names[1:7], title = module.names[1:7], shape = "circularImage", color = "#fca311", color.highlight = "#fca311", image = paste0(path_to_images, "mod", 1:7, ".png"))#, value = c(rep(130,7), rep(130,8)))
    edges <- data.frame(from = 1, to = 2:7)
    
    visNetwork(nodes, edges)  %>%
      visInteraction(hover = F,
                     dragNodes = T,
                     dragView = T) %>%
      visOptions(nodesIdSelection = list(enabled = T, 
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
  
  
  # output$shiny_return <- renderText({
  #   paste("Current node selection : ", input$network_selected)
  # })
  
  #syndrome comparison reactive####
  syndrome_comps <- eventReactive(input$update_comp, {
    comp_age <- input$comp_age
    comp_sex <- input$comp_sex
    reference <- input$reference
    synd_comp <- input$synd_comp
    
    if(is.null(input$network_selected)) selected.node <- 1 else if(input$network_selected == ""){
      selected.node <- 1} else{
        selected.node <- as.numeric(input$network_selected)}
    
    print(paste0("current node reactive: ", input$network_selected))
    
    raw_api_res <- httr::GET(url = paste0("http://localhost:6352", "/synd_comp"),
                             query = list(comp_age = comp_age, comp_sex = comp_sex, reference = reference, synd_comp = synd_comp, facial_subset = selected.node),
                             encode = "json")
    
    json2list <- jsonlite::fromJSON(httr::content(raw_api_res, "text"))
    
    synd_mesh <- atlas
    ref_mesh <- atlas
    
    synd_mesh$vb[-4,] <- t(json2list$syndrome_pred)
    ref_mesh$vb[-4,] <- t(json2list$reference_pred)
    
    return(list(ref_mesh, synd_mesh, json2list$wholeface_scores, json2list$subregion_scores))
    
    
  })

  
  output$gestalt <- renderRglwidget({
    pdf(NULL)
    dev.off()
    
    clear3d()
    #visualize full model estimates####
     # if(input$sex == "Female"){selected.sex <-1
     # } else if(input$sex == "Male"){selected.sex <- 0} 
    selected.sex <- input$sex
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
      print(range(pred.synd()[[1]]))
      print(range(main.res * 1e4))
    pred.mesh$vb[-4,] <-  t(pred.synd()[[1]] + main.res*1e4) #+ interaction.res)
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
      raw_api_res <- httr::GET(url = paste0("http://localhost:6352", "/predtexture"),
                               query = list(selected.sex = selected.sex, selected.age = selected.age, selected.synd = selected.synd),
                               encode = "json")
      
      hex.mean  <- jsonlite::fromJSON(httr::content(raw_api_res, "text"))
      test$mesh$material$color[-eye.index] <- hex.mean[-eye.index]
      mesh.color <- test$mesh$material$color
      lit <- F
      }
    if(input$texture == "Generic + Gestalt"){ 

      #texture simulation####
     
      raw_api_res <- httr::GET(url = paste0("http://localhost:6352", "/predtexture"),
                               query = list(selected.sex = selected.sex, selected.age = selected.age, selected.synd = selected.synd, gestalt_combo = T),
                               encode = "json")
      
      hex.mean  <- jsonlite::fromJSON(httr::content(raw_api_res, "text"))
      test$mesh$material$color[-eye.index] <- hex.mean[-eye.index]
      mesh.color <- test$mesh$material$color
      lit <- F
      }
    
      plot3d(vcgSmooth(test$mesh), aspect = "iso", axes = F, box = F, xlab = "", ylab = "", zlab = "", lit = lit, specular = 1, col = mesh.color)
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
    
    print(paste0("current node: ", selected.node))
    
    if(input$displace & input$synd_comp != input$reference & selected.node == 1){
      mD.synd <- meshDist(syndrome_comps()[[1]], syndrome_comps()[[2]], plot = F, scaleramp = F, displace = input$displace, alpha = input$transparency)
      a <- render(mD.synd, displace = T, alpha = input$transparency)
    } else if(input$displace == F & input$synd_comp != input$reference & selected.node == 1){
      # print(syndrome_comps())
      mD.synd <- meshDist(syndrome_comps()[[1]], syndrome_comps()[[2]], plot = F, scaleramp = F)
      a <- render(mD.synd, alpha = input$transparency)
      
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
    # print(selected.node)
    g <- ggplotGrob(p)
    bg <- g$grobs[[1]]
    round_bg <- roundrectGrob(x=bg$x, y=bg$y, width=bg$width, height=bg$height,
                              r=unit(0.1, "snpc"),
                              just=bg$just, name=bg$name, gp=bg$gp, vp=bg$vp)
    g$grobs[[1]] <- round_bg
    
    plot(g)
    

  })
  
  
  
}

shinyApp(ui = ui, server = server)











