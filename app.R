library(shiny)
library(Morpho)
library(geomorph)
library(shinycssloaders)
library(Jovid)
library(Rvcg)

 #save(starting.mean, starting.lms, filtered.lms, file = "~/shiny/Syndrome_gestalts/data.Rdata")

# load("/srv/shiny-server/Syndrome_gestalts/data.Rdata")
#load("~/shiny/Syndrome_gestalts/data.Rdata")
load("/srv/shiny-server/Syndrome_gestalts/data_downsampled.Rdata")

#define new morpho plotting method until he releases it on cran

render <- function(x,...) UseMethod("render")

#' @rdname render
#' @method render meshDist
#' @export
render.meshDist <- function(x,from=NULL,to=NULL,steps=NULL,ceiling=NULL,uprange=NULL,tol=NULL,tolcol=NULL,rampcolors=NULL,NAcol=NULL,displace=FALSE,shade=TRUE,sign=NULL,add=FALSE,scaleramp=NULL,...) {
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
    colMesh$material$color <- matrix(colfun(colMesh$it),dim(colMesh$it))
    colMesh$material$color[is.na(colMesh$material$color)] <- NAcol
    #colMesh$material$color <- matrix(colfun(colMesh$it),dim(colMesh$it))
    colramp <- list(1,colseq, matrix(data=colseq, ncol=length(colseq),nrow=1),col=ramp,useRaster=T,ylab="Distance in mm",xlab="",xaxt="n")
  } else {
    if (is.null(tol))
      tol <- x$params$tol
    colramp <- x$colramp
    colMesh <- x$colMesh
  }
  if (is.null(tolcol))
    tolcol <- x$params$tolcol
  
  if (shade)
    shade3d(vcgUpdateNormals(colMesh),specular="black",meshColor="legacy",...)
  if (displace) {
    dismesh <- colMesh
    vl <- dim(colMesh$vb)[2]
    dismesh$vb <- cbind(colMesh$vb,rbind(clost,1))
    dismesh$it <- rbind(1:vl,1:vl,(1:vl)+vl)
    dismesh$material$color <- rbind(colorall,colorall,colorall)
    wire3d(dismesh,lit=FALSE,meshColor="legacy")
  }
  diffo <- ((colramp[[2]][2]-colramp[[2]][1])/2)
  image(colramp[[1]],colramp[[2]][-1]-diffo,t(colramp[[3]][1,-1])-diffo,col=colramp[[4]],useRaster=TRUE,ylab="Distance in mm",xlab="",xaxt="n")
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
      selectInput("reference", label = "Reference", choices = sort(unique(filtered.lms$Syndrome)), selected = "Non-syndromic"),
      selectInput("synd", label = "Syndrome", choices = sort(unique(filtered.lms$Syndrome)), selected = sample(unique(filtered.lms$Syndrome), 1)),
      sliderInput("transparency", label = "Mesh transparency", min = 0, max = 1, step = .1, value = .5),
      checkboxInput("displace", label = "Plot vectors?", value = T),
      checkboxInput("variance", label = "Show syndrome heterogeneity?", value = F),
      submitButton("Update", icon("refresh")),
      width = 3, 
      tags$head(
        tags$style("body {background-color: #1a1a1a; }
                   .well {background-color: #404040;}
                   .control-label {color: white}
                   .irs-min, .irs-max {color: white}
                   .checkbox {color: white}
                   ")
      )
    ),
    
    mainPanel(
      #img(src='uc_logo.jpg', align = "right", height = 70 * 1.15, width = 90 * 1.25),
      withSpinner(rglwidgetOutput("gestalt", width = "75vw", height="95vh"), type = 6, color = "#fca311")
    )
  )
)


server <- function(input, output, session) {

  observe({
    if(input$variance > 0){
      updateCheckboxInput(session, "displace", label = "Plot vectors?", value = F)
    }
  })
  
  output$gestalt <- renderRglwidget({
    pdf(NULL)
    dev.off()
    
    controls <- matrix(colMeans(filtered.lms[filtered.lms$Syndrome == input$reference,-1]), nrow = 65, byrow = T)
    
    # nonsynd.mean <- rotmesh.onto(starting.mean$mesh, refmat = as.matrix(starting.lms), tarmat = controls, scale = T, reflection = T)
    nonsynd.mean <- tps3d(starting.mean$mesh, starting.lms, controls)
    
    syndrome.mean <- matrix(colMeans(filtered.lms[filtered.lms$Syndrome == input$synd,-1]), nrow = 65, byrow = T)
    
    clear3d()
  
    # par3d(userMatrix = matrix(c(.998,-.005,.0613,0,.0021,.999,.045,0,-.061,-.045,.997,0,0,0,0,1),ncol =4,nrow = 4))
    par3d(userMatrix = matrix(c(-.017,-.999,-.022,0,-.999,-.016,-.03,0,.03,-.023,.999,0,0,0,0,1),ncol =4,nrow = 4))
    par3d(zoom = .7)
    bg3d(color = "#1a1a1a")
    synd.mesh <- tps3d(nonsynd.mean, controls, syndrome.mean)
    if(input$displace){
    mD.synd <- meshDist(nonsynd.mean, synd.mesh, plot = F, scaleramp = F, displace = input$displace, alpha = input$transparency)
    a <- render(mD.synd, displace = input$displace, alpha = input$transparency)
    } 
    
    if(input$variance){
    mD.synd <- meshDist(nonsynd.mean, synd.mesh, plot = F, scaleramp = F, displace = F)
    a <- render(mD.synd, displace = F, alpha = input$transparency)
    old.target <- row2array3d(na.omit(filtered.lms[filtered.lms$Syndrome == input$synd,-1]))
    sampler <- sample(1:dim(old.target)[3], size = 10)
    if(dim(old.target)[3] < 10) sampler <- 1:dim(old.target)[3]
    for(j in sampler){
      for(i in 1:65){
        arrow3d(controls[i,], old.target[i,,j], type = "lines", col = "red", barblen = 0)
      }
    }
    }
    
    rglwidget()

   })
  
}

shinyApp(ui = ui, server = server)