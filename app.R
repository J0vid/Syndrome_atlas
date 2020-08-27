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


# load("/srv/shiny-server/data_segmented.Rdata")
load("~/shiny/shinyapps/Syndrome_gestalts/data_segmented.Rdata")

#define new morpho plotting method until he releases it on cran

#' render <- function(x,...) UseMethod("render")
#' 
#' #' @rdname render
#' #' @method render meshDist
#' #' @export
#' render.meshDist <- function(x,from=NULL,to=NULL,steps=NULL,ceiling=NULL,uprange=NULL,tol=NULL,tolcol=NULL,rampcolors=NULL,NAcol=NULL,displace=FALSE,shade=TRUE,sign=NULL,add=FALSE,scaleramp=NULL,...) {
#'   clost <- x$clost
#'   dists <- x$dists
#'   distsOrig <- dists
#'   colorall <- x$cols
#'   colramp <- x$colramp
#'   params <- x$params
#'   distqual <- x$distqual
#'   if (!is.null(tolcol))
#'     tolcol <- colorRampPalette(tolcol)(1)
#'   if (!add) {
#'     if (rgl.cur() !=0)
#'       rgl.clear()
#'   }
#'   if (!is.null(from) || !is.null(to) || !is.null(uprange) ||  !is.null(tol)  ||  !is.null(sign) || !is.null(steps) || !is.null(rampcolors) || !is.null(NAcol) || !is.null(tolcol) || !is.null(scaleramp)) {
#'     neg=FALSE
#'     colMesh <- x$colMesh
#'     if(is.null(steps))
#'       steps <- x$params$steps
#'     if (is.null(rampcolors))
#'       rampcolors <- x$params$rampcolors
#'     if (is.null(NAcol))
#'       NAcol <- x$params$NAcol
#'     if (is.null(tolcol))
#'       tolcol <- x$params$tolcol
#'     if (is.null(tol))
#'       tol <- x$params$tol
#'     if(is.null(sign))
#'       sign <- x$params$sign
#'     if (!sign) {
#'       distsOrig <- dists
#'       dists <- abs(dists)
#'     }
#'     if(is.null(ceiling))
#'       ceiling <- x$params$ceiling
#'     if(is.null(uprange))
#'       uprange <- x$params$uprange
#' 
#'     if (is.null(from)) {
#'       mindist <- min(dists)
#'       if (sign && mindist < 0 ) {
#'         from <- quantile(dists,probs=(1-uprange))
#'         neg <- TRUE
#'       } else {
#'         from <- 0
#'       }
#'     }
#'     if (is.null(scaleramp))
#'       scaleramp <- x$params$scaleramp
#' 
#'     if (from < 0)
#'       neg <- TRUE
#'     if (is.null(to))
#'       to <- quantile(dists,probs=uprange)
#'     if(ceiling)
#'       to <- ceiling(to)
#' 
#'     to <- to+1e-10
#'     #ramp <- blue2green2red(maxseq*2)
#'     ramp <- colorRampPalette(rampcolors)(steps-1)
#'     colseq <- seq(from=from,to=to,length.out=steps)
#'     coldif <- colseq[2]-colseq[1]
#'     if (neg && sign) {
#' 
#'       negseq <- length(which(colseq<0))
#'       poseq <- steps-negseq
#'       maxseq <- max(c(negseq,poseq))
#'       if (scaleramp) {
#'         ramp <- colorRampPalette(rampcolors)(maxseq*2)
#'         ramp <- ramp[c(maxseq-negseq+1):(maxseq+poseq)]
#' 
#'       }
#'       else
#'         ramp <- colorRampPalette(rampcolors)(steps-1)
#'       distqual <- ceiling(((dists+abs(from))/coldif)+1e-14)
#'       #distqual[which(distqual < 1)] <- steps+10
#'     } else if (from > 0) {
#'       distqual <- ceiling(((dists-from)/coldif)+1e-14)
#'     } else {
#'       distqual <- ceiling((dists/coldif)+1e-14)
#'     }
#'     distqual[which(distqual < 1)] <- steps+10
#'     colorall <- ramp[distqual]
#'     if (!is.null(tol)) {
#'       if ( length(tol) < 2 ) {
#'         if (sign) {
#'           tol <- c(-tol,tol)
#'         } else {
#'           tol <- c(0,tol)
#'         }
#'       }
#'       good <- which(abs(dists) < tol[2])
#'       colorall[good] <- tolcol
#'     }
#'     colfun <- function(x){x <- colorall[x];return(x)}
#'     colMesh$material$color <- matrix(colfun(colMesh$it),dim(colMesh$it))
#'     colMesh$material$color[is.na(colMesh$material$color)] <- NAcol
#'     #colMesh$material$color <- matrix(colfun(colMesh$it),dim(colMesh$it))
#'     colramp <- list(1,colseq, matrix(data=colseq, ncol=length(colseq),nrow=1),col=ramp,useRaster=T,ylab="Distance in mm",xlab="",xaxt="n")
#'   } else {
#'     if (is.null(tol))
#'       tol <- x$params$tol
#'     colramp <- x$colramp
#'     colMesh <- x$colMesh
#'   }
#'   if (is.null(tolcol))
#'     tolcol <- x$params$tolcol
#' 
#'   if (shade)
#'     shade3d(vcgUpdateNormals(colMesh),specular="black",meshColor="legacy",...)
#'   if (displace) {
#'     dismesh <- colMesh
#'     vl <- dim(colMesh$vb)[2]
#'     dismesh$vb <- cbind(colMesh$vb,rbind(clost,1))
#'     dismesh$it <- rbind(1:vl,1:vl,(1:vl)+vl)
#'     dismesh$material$color <- rbind(colorall,colorall,colorall)
#'     wire3d(dismesh,lit=FALSE,meshColor="legacy")
#'   }
#'   diffo <- ((colramp[[2]][2]-colramp[[2]][1])/2)
#'   image(colramp[[1]],colramp[[2]][-1]-diffo,t(colramp[[3]][1,-1])-diffo,col=colramp[[4]],useRaster=TRUE,ylab="Distance in mm",xlab="",xaxt="n")
#'   if (!is.null(tol)) {
#'     if (sum(abs(tol)) != 0)
#'       image(colramp[[1]],c(tol[1],tol[2]),matrix(c(tol[1],tol[2]),1,1),col=tolcol,useRaster=TRUE,add=TRUE)
#'   }
#'   params <- list(steps=steps,from=from,to=to,uprange=uprange,ceiling=ceiling,sign=sign,tol=tol,rampcolors=rampcolors,NAcol=NAcol,tolcol=tolcol)
#'   out <- list(colMesh=colMesh,dists=distsOrig,cols=colorall,colramp=colramp,params=params,distqual=distqual,clost=clost)
#' 
#'   class(out) <- "meshDist"
#'   invisible(out)
#' }

ui <- fluidPage(
  
  titlePanel(""),
  sidebarLayout(
    sidebarPanel(
      selectInput("reference", label = "Reference", choices = sort(synd.names), selected = "Control"),
      selectInput("synd", label = "Syndrome", choices = sort(synd.names), selected = sample(synd.names, 1)),
      sliderInput("transparency", label = "Mesh transparency", min = 0, max = 1, step = .1, value = .3),
      # sliderInput("module", label = "Module", min = 1, max = 63, step = 1, value = 1),
      checkboxInput("displace", label = "Plot vectors?", value = T),
      # checkboxInput("variance", label = "Show syndrome heterogeneity?", value = F),
      submitButton("Update", icon("refresh")),
      width = 3, 
      tags$head(
        tags$style("body {background-color: #1a1a1a; }
                   .well {background-color: #404040;}
                   .control-label {color: white}
                   .irs-min, .irs-max {color: white}
                   .checkbox {color: white}
                   .nav-tabs > li > a {
    background-color: transparent;
    color: black;
  }
  
  .nav-tabs > li > a:hover {
    background-color: #a6192e;
    border-color: transparent;
    color: white;
  }
       
  .tabbable > .nav > li[class=active] > a {
    background-color: #a6192e; 
    color:white;
    border-color: transparent;
  }
    
      .skin-blue .main-header .navbar {
    background-color: white;
    border-bottom: 1px solid #ddd;
      }
  
  div.vis-network{
  background-color: white; #1a1a1a
  } 
                   ")
      )
    ),
    
    mainPanel(
      tabsetPanel(#img(src='uc_logo.jpg', align = "right", height = 70 * 1.15, width = 90 * 1.25),
        tabPanel("Gestalt", HTML("<center><font style=\"color: red; font-size: xx-large;\"> Bigger </font><font style=\"color: white; font-size: xx-large;\"> | </font><font style=\"color: lightgreen; font-size: xx-large;\"> Similar </font><font style=\"color: white; font-size: xx-large;\"> | </font><font style=\"color: blue; font-size: xx-large;\"> Smaller </font></center>"),
                 withSpinner(rglwidgetOutput("gestalt", width = "75vw", height="95vh"), type = 6, color = "#fca311")),
        tabPanel("Segments", br(), HTML("<p style=\"color:white;\"><b>Select a facial partition using the dropdown menu or by clicking the bubbles directly</b></p>"), visNetworkOutput("network", height="95vh")),
        tabPanel("Morphospace", br(), plotOutput("morphospace")),
        tabPanel("About", br(), HTML("<p style=\"color:white;\"><b>This app aims to help clinical geneticists better understand the characteristic craniofacial features of various genetic syndromes. There are 3 sections to this app and here is my description of how they work. Here are the people that made this app possible.</b></p>"))#, 
                 #includeHTML("~/shiny/shinyapps/Syndrome_gestalts/about_test.html"), 
                 #includeCSS("~/shiny/shinyapps/Syndrome_gestalts/test.css"))
        
      )
    )
  )
)
   

server <- function(input, output) {

  # observe({
  #   if(input$variance > 0){
  #     updateCheckboxInput(session, "displace", label = "Plot vectors?", value = F)
  #   }
  # })
  
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
    path_to_images <- "https://genopheno.ucalgary.ca/presentations/Child_health_data_science_Aponte/mod_images/"
    module.names <- c("Whole face", "Maxilla/Nose", "Mandible/Sphenoid/Frontal", "Upper lip", "Nasal/Maxilla subset", "Cheek/Mandible", "Sphenoid/Frontal", "Nasolabial", "Philtrum", "Lateral Nasal", "Nose", "Cheek", "Chin/Mandible", "Frontal", "Orbital/Temporal/Sphenoid")
    nodes <- data.frame(id = 1:15, label = module.names, title = module.names, shape = "circularImage", color = "#fca311", color.highlight = "#fca311", image = paste0(path_to_images, "mod", 1:15, ".png"))#, value = c(rep(130,7), rep(130,8)))
    edges <- data.frame(from = rep(1:7, each = 2), to = 2:15)
    
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
            this.moveTo({scale:1.15})}") %>%
      visPhysics(stabilization = FALSE) %>%
      visLayout(randomSeed = 12)
    # visHierarchicalLayout()
  
    
  })
  
  
  output$shiny_return <- renderText({
    paste("Current node selection : ", input$network_selected)
  })
  
  output$gestalt <- renderRglwidget({
    pdf(NULL)
    dev.off()


    clear3d()
    #file loading
    # nonsynd.mean <- file2mesh(paste0("Downloads/Module1/", input$reference, ".obj"))
    # synd.mean <- file2mesh(paste0("Downloads/Module1/", input$synd, ".obj"))
    if(is.null(input$network_selected)) selected.node <- 1 else if(input$network_selected == ""){
      selected.node <- 1} else{
      selected.node <- as.numeric(input$network_selected)} 
    
    
    nonsynd.mean <- modules[[selected.node]]
    nonsynd.mean$vb[1:3,] <- synd.mat[, segmentation[,selected.node] == T,synd.names == input$reference]


    synd.mean <- modules[[selected.node]]
    synd.mean$vb[1:3,] <- synd.mat[, segmentation[,selected.node] == T, synd.names == input$synd]

     # par3d(userMatrix = matrix(c(.998,-.005,.0613,0,.0021,.999,.045,0,-.061,-.045,.997,0,0,0,0,1),ncol = 4,nrow = 4))
    par3d(userMatrix = front.face)
    par3d(zoom = .83)
    bg3d(color = "#1a1a1a")
    bg3d(color = "white")
    if(input$displace){
    mD.synd <- meshDist(nonsynd.mean, synd.mean, plot = F, scaleramp = F, displace = input$displace, alpha = input$transparency)
    a <- render(mD.synd, displace = T, alpha = input$transparency)
    } else{

      mD.synd <- meshDist(nonsynd.mean, synd.mean, plot = F, scaleramp = F, displace = input$displace, alpha = input$transparency)
      a <- render(mD.synd, displace = F, alpha = input$transparency)

    }

    if(selected.node > 1){
      module1.ref <- modules[[1]]
      module1.ref$vb[1:3,] <- synd.mat[,,synd.names == input$reference]
      shade3d(module1.ref, col = "lightgrey", alpha = .3, specular = 1)
    }

    if(input$synd == input$reference){
      plot3d(nonsynd.mean, col = "lightgrey", alpha = 1, specular = 1, box = F, axes = F, xlab = "", ylab = "", zlab = "")
      aspect3d("iso")
    }

    rglwidget()

   })
  
  output$morphospace <- renderPlot({
    
    if(is.null(input$network_selected)) selected.node <- 1 else if(input$network_selected == ""){
      selected.node <- 1} else{
        selected.node <- as.numeric(input$network_selected)} 
    
    module.names <- c("Whole face", "Maxilla/Nose", "Mandible/Sphenoid/Frontal", "Upper lip", "Nasal/Maxilla subset", "Cheek/Mandible", "Sphenoid/Frontal", "Nasolabial", "Philtrum", "Lateral Nasal", "Nose", "Cheek", "Chin/Mandible", "Frontal", "Orbital/Temporal/Sphenoid")
    
    full.morpho <- procdist.array(synd.mat, synd.mat[,,synd.names == input$synd])
    
    sub.morpho <- procdist.array(synd.mat[,segmentation[,selected.node] == T,], synd.mat[, segmentation[,selected.node] == T, synd.names == input$synd])

    # 
    # full.morpho <- procdist.array(synd.mat, synd.mat[,,synd.names == "Nager"])
    # 
    # sub.morpho <- procdist.array(synd.mat[,segmentation[,2] == T,], synd.mat[, segmentation[,2] == T, synd.names == "Nager"])
    
    morpho.df <- data.frame(full = full.morpho, sub = sub.morpho, synd.names = synd.names)[full.morpho != 0,]
    
    highlight_df <- morpho.df[sort(morpho.df$sub, index.return = T)$ix,]
    #to do
    #highlight 5 closest syndromes with names
    #maybe all points should be names, with 5 highlighted
    ggplot(morpho.df) +
      geom_point(aes(y = full, x = sub), color = "white") +
      geom_point(data = highlight_df[1:5,], aes(y = full, x = sub), color = "white", fill = "#fca311", shape = 21, size = 3) +
      geom_label_repel(data = highlight_df[c(1:5, 60:64),], aes(y = full, x = sub, label = synd.names),
                       force = 2,
                       size = 4,
                       color = c(rep("#fca311", 5), rep("white", 5)),
                       fill = "#1a1a1a") +
      ylab(paste0("Full face distance from ", input$synd)) +
      xlab(paste0(module.names[selected.node], " distance from ", input$synd)) +
      # ylim(-.25,3.5) +
      theme_bw() +
      theme(panel.background = element_rect(fill = "#1a1a1a"),
            plot.background = element_rect(fill = "#1a1a1a"),
            axis.text = element_text(color = "white"),
            axis.title = element_text(color = "white"),
            panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                            colour = "#404040"), 
            panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                            colour = "#404040")
      )
    
  })
  
}

shinyApp(ui = ui, server = server)