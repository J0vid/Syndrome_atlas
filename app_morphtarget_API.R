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

 load("2k_topo_objects.Rdata")
 
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
                       sliderInput("age", label = "Age", min = 1, max = 70, value = 12, step = 1, animate = animationOptions(loop = T, interval = 40)),
                       selectInput("severity", label = "Severity", choices = c("mild", "typical", "severe"), selected = "typical"),
                       playwidgetOutput("control"),
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
                 withSpinner(rglwidgetOutput("gestalt", width="75vw", height="80vh"), type = 6, color = "#fca311")),
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
  predicted.shape <- predshape.lm(synd.lm.coefs, datamod, d.registered$PCs[,1:200], d.registered$mshape)
  
  atlas$vb[-4,] <- t(predicted.shape)
  tmp.mesh <- atlas
  
  open3d()
  shade3d(vcgSmooth(atlas), aspect = "iso", col = "lightgrey", specular = 1)
  objid <- ids3d()$id
  xyz <- rgl.attrib(objid[1], "vertices")
  dim(xyz)
  
  
  morph_target <- eventReactive(input$Update, {
    #age morph target####
    min_age <- round(abs(min(doutVar()[[1]])))
    max_age <- round(max(doutVar()[[1]]))
   
    #api call
    selected.synd <- input$synd
    selected.sex <- input$sex
    selected.severity <- input$severity
    
    raw_api_res <- httr::GET(url = paste0("http://localhost:6352", "/predshape"),
                             query = list(selected.sex = selected.sex, selected.synd = selected.synd, selected.severity = selected.severity, min_age = min_age, max_age = max_age),
                             encode = "json")
    
    values  <- jsonlite::fromJSON(httr::content(raw_api_res, "text"))
    
    control <- vertexControl(values = values,
                             vertices = rep(1:nrow(xyz), each = 3),
                             attributes = rep(c("x", "y", "z"), nrow(xyz)),
                             objid = objid)
    
    scene <- scene3d()
    # rgl.close()
    # return(list(scene, control, control2))
    return(list(scene, control))
    
  })
  
  
  output$gestalt <- renderRglwidget({
    rglwidget(morph_target()[[1]], controllers = c("control"))
  })
  
  
  output$control <- renderPlaywidget({
    playwidget("gestalt", morph_target()[[2]],
               respondTo = "age")
  })
  
# 
#   #module selection####
#   output$network <- renderVisNetwork({
#     # 63 modules
#     # nodes <- data.frame(id = 1:63, label = 1:63, shape = "image", image = "https://genopheno.ucalgary.ca/presentations/Child_health_data_science_Aponte/images/atlas.png", value = c(rep(120,21), rep(5,21), rep(2,21)))
#     # edges <- data.frame(from = rep(1:31, each = 2), to = 2:63)#c(2,3,4,5,6,7,8,9,10,11,12,13,14,15))
#     # 
#     # # 31 modules
#     # nodes <- data.frame(id = 1:31, label = 1:31, shape = "image", image = "https://genopheno.ucalgary.ca/presentations/Child_health_data_science_Aponte/images/atlas.png", value = c(rep(10,21), rep(5,10)))
#     # edges <- data.frame(from = rep(1:15, each = 2), to = 2:31)#c(2,3,4,5,6,7,8,9,10,11,12,13,14,15))
#     # 15 modules
#     path_to_images <- "https://raw.githubusercontent.com/J0vid/Syndrome_atlas/main/mod_images/" #"https://genopheno.ucalgary.ca/presentations/Child_health_data_science_Aponte/mod_images/"
#     module.names <- c("Whole face", "Nose", "Mandible/Sphenoid/Frontal", "Upper lip", "Nasal/Maxilla subset", "Cheek/Mandible", "Sphenoid/Frontal", "Nasolabial", "Philtrum", "Lateral Nasal", "Nose", "Cheek", "Chin/Mandible", "Frontal", "Orbital/Temporal/Sphenoid")
#     module.names[1:7] <- c("face", "posterior_mandible", "nose", "anterior_mandible", "brow", "zygomatic", "premaxilla") #c("Whole Face", "Nose", "Eyes", "Jaw", "Chin", "Mouth", "Cheeks")
#     nodes <- data.frame(id = 1:7, label = module.names[1:7], title = module.names[1:7], shape = "circularImage", color = "#fca311", color.highlight = "#fca311", image = paste0(path_to_images, "mod", 1:7, ".png"))#, value = c(rep(130,7), rep(130,8)))
#     edges <- data.frame(from = 1, to = 2:7)
#     
#     visNetwork(nodes, edges)  %>%
#       visInteraction(hover = F,
#                      dragNodes = T,
#                      dragView = T) %>%
#       visOptions(nodesIdSelection = list(enabled = T, 
#                                          selected = "1"),
#                  highlightNearest = F
#       ) %>% 
#       visNodes(size = 50,
#                shapeProperties = list(useBorderWithImage = TRUE),
#                borderWidthSelected = 7,
#                borderWidth = 5,
#                font = list(color = "#fca311",
#                            size = 17)) %>%
#       visEdges(shadow = F,
#                arrows =list(to = list(enabled = TRUE, scaleFactor = .5)),
#                color = list(color = "#f0c173", highlight = "#fca311")) %>%
#       visEvents(type = "once", startStabilizing = "function() {
#             this.moveTo({scale:1.35})}") %>%
#       visPhysics(stabilization = FALSE) %>%
#       visLayout(randomSeed = 12)# %>%
#       #visHierarchicalLayout()
#     
#   })
#   
#   
#   # output$shiny_return <- renderText({
#   #   paste("Current node selection : ", input$network_selected)
#   # })
#   
#   #syndrome comparison reactive####
#   syndrome_comps <- eventReactive(input$update_comp, {
#     comp_age <- input$comp_age
#     comp_sex <- input$comp_sex
#     reference <- input$reference
#     synd_comp <- input$synd_comp
#     
#     if(is.null(input$network_selected)) selected.node <- 1 else if(input$network_selected == ""){
#       selected.node <- 1} else{
#         selected.node <- as.numeric(input$network_selected)}
#     
#     print(paste0("current node reactive: ", input$network_selected))
#     
#     raw_api_res <- httr::GET(url = paste0("http://localhost:6352", "/synd_comp"),
#                              query = list(comp_age = comp_age, comp_sex = comp_sex, reference = reference, synd_comp = synd_comp, facial_subset = selected.node),
#                              encode = "json")
#     
#     json2list <- jsonlite::fromJSON(httr::content(raw_api_res, "text"))
#     
#     synd_mesh <- atlas
#     ref_mesh <- atlas
#     
#     synd_mesh$vb[-4,] <- t(json2list$syndrome_pred)
#     ref_mesh$vb[-4,] <- t(json2list$reference_pred)
#     
#     return(list(ref_mesh, synd_mesh, json2list$wholeface_scores, json2list$subregion_scores))
#     
#     
#   })
# 
#   output$comparison <- renderRglwidget({
#     pdf(NULL)
#     dev.off()
#     
#     clear3d()
#     
#     if(is.null(input$network_selected)) selected.node <- 1 else if(input$network_selected == ""){
#       selected.node <- 1} else{
#         selected.node <- as.numeric(input$network_selected)} 
#     
#     print(paste0("current node: ", selected.node))
#     
#     if(input$displace & input$synd_comp != input$reference & selected.node == 1){
#       mD.synd <- meshDist(syndrome_comps()[[1]], syndrome_comps()[[2]], plot = F, scaleramp = F, displace = input$displace, alpha = input$transparency)
#       a <- render(mD.synd, displace = T, alpha = input$transparency)
#     } else if(input$displace == F & input$synd_comp != input$reference & selected.node == 1){
#       # print(syndrome_comps())
#       mD.synd <- meshDist(syndrome_comps()[[1]], syndrome_comps()[[2]], plot = F, scaleramp = F)
#       a <- render(mD.synd, alpha = input$transparency)
#       
#     }
#     
#     if(selected.node > 1){
#      
#       mD.synd <- meshDist(syndrome_comps()[[1]], syndrome_comps()[[2]], plot = F, scaleramp = F, displace = input$displace, alpha = input$transparency)
#       mD.synd <- rmVertex(mD.synd$colMesh, index = which(modules[,selected.node] == 1), keep = T)
#       # a <- render(mD.synd, displace = F, alpha = input$transparency)
#       plot3d(mD.synd, aspect = "iso", axes = F, box = F, xlab = "", ylab = "", zlab = "", lit = T, specular = 1)
#       shade3d(syndrome_comps()[[1]], col = "lightgrey", alpha = .2, specular = 1)
#     }
#     
#     par3d(userMatrix = front.face, zoom = .83)
#     aspect3d("iso")
#     rglwidget()
#     
#   })
#   
#   output$morphospace <- renderPlot({
# 
#     if(is.null(input$network_selected)) selected.node <- 1 else if(input$network_selected == ""){
#       selected.node <- 1} else{
#         selected.node <- as.numeric(input$network_selected)}
# 
#     module.names <- colnames(modules)
# #
# #     full.morpho <- procdist.array(synd.mat, synd.mat[,,d.meta.combined$Syndrome == input$synd_comp])
# #
# #     sub.morpho <- procdist.array(synd.mat[modules[,selected.node] == 1,,], synd.mat[modules[,selected.node] == 1,, d.meta.combined$Syndrome == input$synd_comp])
#     if(input$score_plot == "Raw score"){
#         morpho.df <- data.frame(full = syndrome_comps()[[3]]$face_score, sub = syndrome_comps()[[4]]$face_score, Syndrome = levels(d.meta.combined$Syndrome))
#         #sort to get which synds to highlight
#         highlight_df <- morpho.df[sort(morpho.df$sub, index.return = T)$ix,]
#         highlighted_points <- geom_point(data = highlight_df[85:89,], aes(y = full, x = sub), color = "black", fill = "#a6192e", shape = 21, size = 3) 
#         tmp_xlab <- xlab(paste0(module.names[selected.node], " score for ", input$synd_comp))
#         tmp_ylab <- ylab(paste0("Whole face score for ", input$synd_comp))
#         tmp_labels <- geom_label_repel(data = highlight_df[c(1:5, 85:89),], aes(y = full, x = sub, label = Syndrome),
#                                        force = 2,
#                                        size = 4,
#                                        color = c(rep("black", 5), rep("#a6192e", 5)),
#                                        fill = "white")
#     }
#     
#     if(input$score_plot == "Similarity"){
#       full_similarity <- sqrt((syndrome_comps()[[3]]$face_score - syndrome_comps()[[3]]$face_score[levels(d.meta.combined$Syndrome) == input$synd_comp])^2)
#       module_similarity <- sqrt((syndrome_comps()[[4]]$face_score - syndrome_comps()[[4]]$face_score[levels(d.meta.combined$Syndrome) == input$synd_comp])^2)
#       
#       morpho.df <- data.frame(full = full_similarity, sub = module_similarity, Syndrome = levels(d.meta.combined$Syndrome))
#       #sort to get which synds to highlight
#       highlight_df <- morpho.df[sort(morpho.df$sub, index.return = T)$ix,]
#       highlighted_points <- geom_point(data = highlight_df[1:5,], aes(y = full, x = sub), color = "black", fill = "#a6192e", shape = 21, size = 3) 
#       tmp_xlab <- xlab(paste0(module.names[selected.node], " similarity to ", input$synd_comp))
#       tmp_ylab <- ylab(paste0("Whole face similarity to ", input$synd_comp))
#       tmp_labels <- geom_label_repel(data = highlight_df[c(1:5, 85:89),], aes(y = full, x = sub, label = Syndrome),
#                        force = 2,
#                        size = 4,
#                        color = c(rep("#a6192e", 5), rep("black", 5)),
#                        fill = "white")
#     }
#     #to do
#     #highlight 5 closest syndromes with names
#     #maybe all points should be names, with 5 highlighted
#     p <- ggplot(morpho.df) +
#       geom_point(aes(y = full, x = sub), color = "black") +
#       highlighted_points +
#       tmp_labels +
#       tmp_ylab +
#       tmp_xlab +
#       # ylim(-.25,3.5) +
#       theme_bw() +
#       theme(plot.background = element_rect(fill = "#fcfcfc"),
#             axis.text = element_text(color = "black"),
#             axis.title = element_text(color = "black"),
#             panel.grid.major = element_line(size = 0.5, linetype = 'solid',
#                                             colour = "grey65"),
#             panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
#                                             colour = "grey65"), 
#             panel.background = element_blank(),
#             plot.margin = margin(20, 20, 20, 20)
#       )
#     # print(selected.node)
#     g <- ggplotGrob(p)
#     bg <- g$grobs[[1]]
#     round_bg <- roundrectGrob(x=bg$x, y=bg$y, width=bg$width, height=bg$height,
#                               r=unit(0.1, "snpc"),
#                               just=bg$just, name=bg$name, gp=bg$gp, vp=bg$vp)
#     g$grobs[[1]] <- round_bg
#     
#     plot(g)
#     
# 
#   })
#   
  
  
}

shinyApp(ui = ui, server = server)










