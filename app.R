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
library(sparsediscrim)

setwd("~/shiny/shinyapps/Syndrome_model/")
#setwd("/data/Syndrome_model_data/")
# save(atlas, d.meta.combined, front.face, PC.eigenvectors, synd.lm.coefs, synd.mshape, PC.scores, synd.mat, file = "data.Rdata")
 load("data.Rdata")
 load("module_indices.Rdata")
 load("pose.Rdata")
 #calculations at startup that should make it into the startup file
 hdrda.df <- data.frame(synd = d.meta.combined$Syndrome, PC.scores[,1:200])
 hdrda.mod <- hdrda(synd ~ ., data = hdrda.df)
 
 predshape.lm <- function(fit, datamod, PC, mshape){
   dims <- dim(mshape)
   mat <- model.matrix(datamod)
   pred <- mat %*% fit
   
   predPC <- (PC %*% t(pred))
   out <- mshape + matrix(predPC, dims[1], dims[2], byrow = F)
   
   return(out)
 }
 
 predtexture.lm <- function(fit, datamod, PC, mshape, gestalt_combo = NULL){
   dims <- dim(mshape)
   mat <- model.matrix(datamod)
   pred <- mat %*% fit
   names <- as.matrix(model.frame(datamod))
   names <- apply(names, 1, paste, collapse = "_")
   names <- gsub(" ", "", names)
   predPC <- t(PC %*% t(pred))
   out <- mshape + predPC
   
   predicted.texture3d <- row2array3d(out)[,,1]
   
   if(is.null(gestalt_combo) == F){
     final.texture <- 3 * (predicted.texture3d) + t(col2rgb(atlas$material$color))
   } else {final.texture <- predicted.texture3d}
   #scale values
   maxs <- apply(final.texture, 2, max)
   mins <- apply(final.texture, 2, min)
   additive.texture <- scale(final.texture, center = mins, scale = maxs - mins)
   # hex.mean <- rgb(additive.texture, maxColorValue = 1)
   return(additive.texture)
 }
 
ui <- fluidPage(
  tags$head(tags$style(HTML('.irs-from, .irs-to, .irs-min, .irs-max, .irs-single {
            visibility: hidden !important;
    }'))),
  titlePanel(""),
  sidebarLayout(
    sidebarPanel(
      conditionalPanel(condition = "input.Atlas_tabs=='Gestalt'",
                       selectInput("synd", label = "Syndrome", choices = sort(levels(d.meta.combined$Syndrome)), selected = "Costello Syndrome"),
                       sliderInput("age", label = "Age", min = 0, max = 1, value = 0, step = .035, ticks = F, animate = animationOptions(loop = T, interval = 1)),
                       tableOutput("age_data"),
                       selectInput("sex", label = "Sex", choices = c("Female", "Male")),
                       selectInput("severity", label = "Severity", choices = c("Mild", "Typical", "Severe"), selected = "Typical"),
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
                      actionButton("Update", "Update Gestalt", icon("sync"))),
      conditionalPanel(condition = "input.Atlas_tabs=='Comparisons'",
                      selectInput("reference", label = "Reference", choices = sort(levels(d.meta.combined$Syndrome)), selected = "Unaffected Unrelated"),
                      selectInput("synd_comp", label = "Syndrome", choices = sort(levels(d.meta.combined$Syndrome)), selected = "Costello Syndrome"),
                      sliderInput("comp_age", label = "Age", min = 0, max = 1, step = .035, value = .3, ticks = F, animate = animationOptions(loop = T, interval = 25)),
                      tableOutput("age_data_comp"),
                      selectInput("comp_sex", label = "Sex", choices = c("Female", "Male")),
                      selectInput("comp_severity", label = "Severity", choices = c("Mild", "Typical", "Severe"), selected = "Typical"),
                      playwidgetOutput("control_comp"),
                      hr(),
                      bsButton("q1", label = "", icon = icon("question"), style = "default", size = "extra-small"),
                      selectInput("score_plot", label = "Plot type", choices = c("Similarity", "Raw score")),
                      bsPopover(id = "q1", title = "Plot types",
                                content = paste0("When syndromes are very Severe, they may outscore the selected syndrome. In these cases, it may be more informative to look at syndromes that are similar in the magnitude of their effects."
                                ),
                                placement = "right", 
                                trigger = "hover", 
                                options = list(container = "body")
                      ),
                      actionButton("update_comp", "Update Comparison", icon("sync"))),
      conditionalPanel(condition = "input.Atlas_tabs=='Submitted face'",
                       fileInput("file1", "",
                                 multiple = T,
                                 accept = c(".obj", ".ply", ".png", ".pp", ".csv")
                                 ),
                       splitLayout(cellWidths = c("50%", "50%"),
                                   textInput("current_age", label = "Age (years)", value = "", placeholder = "ie. 9.5"),
                                   selectInput("current_sex", label = "Sex", selected = "Male", choices = c("Male", "Female"))
                                  ),
                       checkboxGroupInput("submitted_heatmap", label = "Morphology", choices = "Compare to syndrome?"),
                       conditionalPanel(condition = "input.submitted_heatmap == 'Compare to syndrome?'",
                                        selectInput("submitted_comp", label = "Syndrome", choices = sort(levels(d.meta.combined$Syndrome)), selected = "Costello Syndrome")
                                        )
                       ),
      width = 3, 
      tags$head(
        tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
      )
    ),
    
    mainPanel(
      tabsetPanel(id = "Atlas_tabs",#img(src='uc_logo.jpg', align = "right", height = 70 * 1.15, width = 90 * 1.25),
        tabPanel("Gestalt", br(),
                 withSpinner(rglwidgetOutput("gestalt", width = "65vw", height="55vh"), type = 6, color = "#AE80E6")),
        tabPanel("Comparisons", br(),
                 withSpinner(rglwidgetOutput("comparison", width = "65vw", height="55vh"), type = 6, color = "#AE80E6"),
                 br(),
                 HTML("<center><font style=\"color: red; font-size: xx-large;\"> Bigger </font><font style=\"color: #fcfcfc; font-size: xx-large;\"> | </font><font style=\"color: lightgreen; font-size: xx-large;\"> Similar </font><font style=\"color: #fcfcfc; font-size: xx-large;\"> | </font><font style=\"color: blue; font-size: xx-large;\"> Smaller </font></center>"),
                 br(), plotOutput("morphospace"),
                 br(), 
                 HTML("<h3 style=\"color:black; text-align:center\">Select a facial partition by clicking the bubbles directly and then update the comparison.</h2>"),
                 visNetworkOutput("network", height="80vh")),
        tabPanel("Submitted face", br(),
                 column(width = 12,
                 withSpinner(rglwidgetOutput("submitted_face", width = "65vw", height="55vh"), type = 6, color = "#AE80E6"),
                 ),
                 br(),
                 column(width = 12,
                  br(),
                 shinydashboard::box(plotlyOutput("posterior_scree"), width = 6),
                 br(),
                 shinydashboard::box(plotlyOutput("personal_morphospace"), width = 6)
                 )),
        tabPanel("About", br(), HTML("<p style=\"color:black;\">This app aims to help clinical geneticists better understand the characteristic craniofacial features of various genetic syndromes. There are 3 sections to this app and here is my description of how they work. Here are the people that made this app possible.</p>"))#, 
      )
    )
  )
)


server <- function(input, output, session) {
  
  options(rgl.useNULL = T)
  save <- options(rgl.inShiny = T)
  on.exit(options(save))
  
  
  # initial.scene <- reactive({
  # close3d()
  selected.synd <- factor("Unaffected Unrelated", levels = levels(d.meta.combined$Syndrome))
  selected.sex <-1
  selected.age <- 10 
  
  datamod <- ~ selected.sex + selected.age + selected.age^2 + selected.age^3 + selected.synd + selected.age:selected.synd
  predicted.shape <- predshape.lm(synd.lm.coefs, datamod, PC.eigenvectors, synd.mshape)
  
  atlas$vb[-4,] <- t(predicted.shape)
  tmp.mesh <- atlas
  
  open3d()
  # shinySetPar3d(userMatrix = front.face2, session = session)
  par3d(userMatrix = front.face2)
  shade3d(vcgSmooth(atlas), aspect = "iso", col = "lightgrey", specular = 1)
  objid <- ids3d()$id
  xyz <- rgl.attrib(objid[length(objid)], "vertices")
  # return(list(objid, xyz))
  # })
  
  #synd reactive####
  outVar <- reactive({
    #get selected syndrome age range
    age.range <- d.meta.combined$Age[d.meta.combined$Syndrome == input$synd]
    
    #calculate syndrome severity scores for selected syndrome
    #calculate score for the main effect
    S <- synd.lm.coefs[grepl(pattern = input$synd, rownames(synd.lm.coefs)),][1,]
    Snorm <- S/sqrt(sum(S^2))
    syndscores.main <- PC.scores %*% Snorm
    
    return(list(age.range, syndscores.main[d.meta.combined$Syndrome == input$synd], Snorm, S))
  }) 
  
  doutVar <- outVar #debounce(outVar, 2000)
  
  output$age_data <- renderTable({
    age.data <- range(doutVar()[[1]])
    
    data.frame("Min" = as.integer(abs(age.data[1])), "Current" = as.integer((diff(age.data) * input$age + 1)), "Max" = as.integer(age.data[2]), check.names = F)
  }, align = "c", width = "100%")
  
  observeEvent(input$synd, {
    # slider_min <-  round(abs(min(doutVar()[[1]])))
    # slider_max <- round(max(doutVar()[[1]]))
    # animation_length <- (slider_max - slider_min) / 2000
    # print(animation_length)
    # updateSliderInput(session, "age", label = "Age", min = slider_min, max = slider_max, value = mean(doutVar()[[1]]))
    # how many years/second do I want for the animation? 
    # updateSliderInput(session, "age", label = "Age", min = 0, max = 1, value = 0.29, step = animation_length)
    M.synds <- c("Klinefelter Syndrome", "XXYY")
    F.synds <- c("XXX", "Turner Syndrome", "Craniofrontonasal Dysplasia", "18p Deletion")
    if(input$synd %in% M.synds | input$synd %in% F.synds){
      if(input$synd %in% M.synds){ sex.choices <- c("Male")
      } else {sex.choices <- "Female"}
    updateSelectInput(session, "sex", "Sex", choices = sex.choices)
    }
  })
  
  morph_target <- eventReactive(input$Update, {
    #age morph target####
    min_age <- round(abs(min(doutVar()[[1]])))
    max_age <- round(max(doutVar()[[1]]))
 
    #api call
    selected.synd <- input$synd
    selected.sex <- input$sex
    selected.severity <- input$severity
    selected.color <- input$texture
    severity_sd <- sd(doutVar()[[2]])
    Snorm <-  doutVar()[[3]]
    S <- doutVar()[[4]]
    
    if(selected.severity == "Mild"){ selected.severity <- -1.5 * severity_sd} else if(selected.severity == "Severe"){selected.severity <- 1.5 * severity_sd} else if(selected.severity == "Typical"){selected.severity <- 0}
 
    severity.resids <- S * (selected.severity)
    print(severity.resids)
    # severity.resids <- matrix(0, dim(synd.mshape)[1], dim(synd.mshape)[2])
    
    raw_api_res <- httr::GET(url = paste0("http://localhost:3636", "/predPC"),
                             query = list(selected.sex = selected.sex, selected.synd = selected.synd, min_age = min_age, max_age = max_age, selected.color = "Lightgrey"),#input$texture),
                             encode = "json")
    
    parsed_res  <- jsonlite::fromJSON(httr::content(raw_api_res, "text", encoding = "UTF-8"))
    
    values_shape <- matrix(NA, ncol = xyz.num * 3, nrow = 2)
    values_col <- matrix(NA, ncol = xyz.num * 3, nrow = 2)
    
    for(i in 1:nrow(values_shape)){
    
      #transform PC scores
      predicted.shape <- showPC((severity.resids[1:200] + parsed_res[i,1:200])/1e10, PC.eigenvectors[,1:200], synd.mshape)
      
      tmp.mesh$vb[-4,] <- t(predicted.shape)
      final.shape <- vcgSmooth(tmp.mesh)
      
      shape.tmp <- array(t(final.shape$vb[-4,])[atlas$it, ], dim = c(xyz.num, 3, 1))
      values_shape[i,] <- geomorph::two.d.array(shape.tmp)
      
      if(selected.color == "Generic") col.tmp <- array(t(col2rgb(atlas$material$color[atlas$it])), dim = c(xyz.num, 3, 1))/255
      if(selected.color == "Lightgrey") col.tmp <- array(211/255, dim = c(xyz.num, 3, 1))
      if(selected.color == "Generic + Gestalt"){
        selected.synd <- factor(selected.synd, levels = levels(d.meta.combined$Syndrome))
        if(selected.sex == "Female"){selected.sex <-1.5
        } else if(selected.sex == "Male"){selected.sex <- -.5} 
        
        selected.age <- as.numeric(c(min_age, max_age))
        
        datamod <- ~ selected.sex + selected.age[i] + selected.age[i]^2 + selected.age[i]^3 + selected.synd + selected.age[i]:selected.synd
        
        col.tmp <- array(predtexture.lm(texture.coefs, datamod, texture.pcs, texture.mean, gestalt_combo = T)[atlas$it,], dim = c(xyz.num, 3, 1))
        # col.tmp <- array(t(col2rgb(atlas$material$color[atlas$it])), dim = c(xyz.num, 3, 1))/255
        }
      values_col[i,] <- geomorph::two.d.array(col.tmp)
    }
    print(values_col[1,1:5])
    combined_values <- rbind(as.numeric(t(cbind(values_shape[1,], values_col[1,]))), as.numeric(t(cbind(values_shape[2,], values_col[2,]))))
    
    print(dim(combined_values))
    
    control <- vertexControl(values = combined_values,
                             vertices = rep(1:nrow(xyz), each = 6),
                             attributes = c(rep(c("x", "red", "y", "green", "z", "blue"), nrow(xyz))),
                             objid = objid)
    
    scene <- scene3d()
    
    return(list(scene, control))
    
  })
  
  output$gestalt <- renderRglwidget({
    par3d(userMatrix = front.face2, zoom = .75)
    rglwidget(morph_target()[[1]], controllers = c("control"))
  })
  
  output$control <- renderPlaywidget({
    playwidget("gestalt", morph_target()[[2]],
               respondTo = "age", step = .01)
  })
  
  #facial subregion selection####
  output$network <- renderVisNetwork({

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
  
  #comparison reactive####
  outVar2 <- reactive({
    #get selected syndrome age range
    age.range <- d.meta.combined$Age[d.meta.combined$Syndrome == input$reference]
    
    #calculate syndrome severity scores for selected syndrome
    #calculate score for the main effect
    if(input$reference == "Unaffected Unrelated"){
      S <- synd.lm.coefs[1,]
    } else {S <- synd.lm.coefs[grepl(pattern = input$reference, rownames(synd.lm.coefs)),][1,]}
    
    Snorm <- S/sqrt(sum(S^2))
    syndscores.main <- PC.scores %*% Snorm
    
    return(list(age.range, syndscores.main[d.meta.combined$Syndrome == input$reference]))
  }) 
  
  output$age_data_comp <- renderTable({
    age.data <- range(outVar2()[[1]])
    
    data.frame("Min" = as.integer(abs(age.data[1])), "Current" = as.integer((diff(age.data) * input$comp_age + 1)), "Max" = as.integer(age.data[2]), check.names = F)
  }, align = "c", width = "100%")
  
  observeEvent(input$synd_comp, {
    # slider_min <-  round(abs(min(outVar2()[[1]])))
    # slider_max <- round(max(outVar2()[[1]]))
    # updateSliderInput(session, "age", label = "Age", min = slider_min, max = slider_max, value = mean(doutVar()[[1]]))
    # updateSliderInput(session, "comp_age", label = "Age", min = 0, max = 1, value = 0)
    M.synds <- c("Klinefelter Syndrome", "XXYY")
    F.synds <- c("XXX", "Turner Syndrome", "Craniofrontonasal Dysplasia", "18p Deletion")
    if(input$synd_comp %in% M.synds | input$synd_comp %in% F.synds){
      if(input$synd_comp %in% M.synds){ sex.choices <- c("Male")
      } else {sex.choices <- "Female"}
      updateSelectInput(session, "comp_sex", "Sex", choices = sex.choices)
    }
  })
  
  morph_target_comparison <- eventReactive(input$update_comp, {
    #age morph target####
    min_age <- round(abs(min(outVar2()[[1]])))
    max_age <- round(max(outVar2()[[1]]))
    
    #api call
    selected.synd <- input$reference
    synd_comp <- input$synd_comp
    selected.sex <- input$comp_sex
    selected.severity <- input$comp_severity
    
    if(is.null(input$network_selected)){
      selected.node <- 1 } else if(input$network_selected == ""){
           selected.node <- 1} else{
           selected.node <- as.numeric(input$network_selected)
           }
    
    raw_api_res <- httr::GET(url = paste0("http://localhost:3636", "/comparison_morphtarget"),
                             query = list(selected.sex = selected.sex, selected.synd = selected.synd, synd_comp = synd_comp, selected.severity = selected.severity, min_age = min_age, max_age = max_age, facial_subregion = selected.node),
                             encode = "json")
    
    values  <- jsonlite::fromJSON(httr::content(raw_api_res, "text", encoding = "UTF-8"))
    
    control_combined <- vertexControl(values = values, 
                                      vertices = c(rep(1:nrow(xyz), each = 6)),
                                      attributes = c(rep(c("x", "red", "y", "green", "z", "blue"), nrow(xyz))),
                                      objid = objid)
    
    scene <- scene3d()
    # rgl.close()
    # return(list(scene, control, control2))
    return(list(scene, control_combined))
  })
  
  output$comparison <- renderRglwidget({
    par3d(userMatrix = front.face2, zoom = .75)
    rglwidget(morph_target_comparison()[[1]], controllers = c("control_comp"))
  })
  
  output$control_comp <- renderPlaywidget({
    playwidget("comparison", morph_target_comparison()[[2]],
               respondTo = "comp_age", step = .01)
  })
  
  #syndrome comparison reactive####
  syndrome_comps <- eventReactive(input$update_comp, {
    reference <- input$reference
    synd_comp <- input$synd_comp
    
    if(is.null(input$network_selected)) selected.node <- 1 else if(input$network_selected == ""){
      selected.node <- 1} else{
        selected.node <- as.numeric(input$network_selected)}
    
    raw_api_res <- httr::GET(url = paste0("http://localhost:3636", "/similarity_scores"),
                             query = list(reference = reference, synd_comp = synd_comp, facial_subregion = selected.node),
                             encode = "json")
    
    json2list <- jsonlite::fromJSON(httr::content(raw_api_res, "text"))
    
    return(list(json2list$wholeface_scores, json2list$subregion_scores))
    
  })
  
  output$morphospace <- renderPlot({
    
    node.code <- c("posterior_mandible" = 2, "nose" = 3,"anterior_mandible" = 4, "brow" = 5, "zygomatic" = 6, "premaxilla" = 7)
    
    if(is.null(input$network_selected)) selected.node <- 1 else if(input$network_selected == ""){
      selected.node <- 1} else{
        selected.node <- as.numeric(input$network_selected)}

    module.names <- colnames(modules)

    if(input$score_plot == "Raw score"){
        morpho.df <- data.frame(full = syndrome_comps()[[1]]$face_score, sub = syndrome_comps()[[2]]$face_score, Syndrome = syndrome_comps()[[1]]$Syndrome)
        #sort to get which synds to highlight
        highlight_df <- morpho.df[sort(morpho.df$sub, index.return = T)$ix,]
        highlighted_points <- geom_point(data = highlight_df[85:89,], aes(y = full, x = sub), color = "black", fill = "#a6192e", shape = 21, size = 3)
        tmp_xlab <- xlab(paste0(names(node.code)[node.code == selected.node], " score for ", input$synd_comp))
        tmp_ylab <- ylab(paste0("Whole face score for ", input$synd_comp))
        tmp_labels <- geom_label_repel(data = highlight_df[c(1:5, 85:89),], aes(y = full, x = sub, label = Syndrome),
                                       force = 2,
                                       size = 4,
                                       color = c(rep("black", 5), rep("#a6192e", 5)),
                                       fill = "white")
    }

    if(input$score_plot == "Similarity"){
      full_similarity <- sqrt((syndrome_comps()[[1]]$face_score - syndrome_comps()[[1]]$face_score[levels(d.meta.combined$Syndrome) == input$synd_comp])^2)
      module_similarity <- sqrt((syndrome_comps()[[2]]$face_score - syndrome_comps()[[2]]$face_score[levels(d.meta.combined$Syndrome) == input$synd_comp])^2)

      morpho.df <- data.frame(full = full_similarity, sub = module_similarity, Syndrome = syndrome_comps()[[1]]$Syndrome)
      #sort to get which synds to highlight
      highlight_df <- morpho.df[sort(morpho.df$sub, index.return = T)$ix,]
      highlighted_points <- geom_point(data = highlight_df[1:5,], aes(y = full, x = sub), color = "black", fill = "#a6192e", shape = 21, size = 3)
      tmp_xlab <- xlab(paste0(names(node.code)[node.code == selected.node], " similarity to ", input$synd_comp))
      tmp_ylab <- ylab(paste0("Whole face similarity to ", input$synd_comp))
      tmp_labels <- geom_label_repel(data = highlight_df[c(1:5, 85:89),], aes(y = full, x = sub, label = Syndrome),
                       force = 2,
                       size = 4,
                       color = c(rep("#a6192e", 5), rep("black", 5)),
                       fill = "white")

    }
   
    p <- ggplot(morpho.df) +
      geom_point(aes(y = full, x = sub), color = "black") +
      highlighted_points +
      tmp_labels +
      tmp_ylab +
      tmp_xlab +
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
  
  #elements of the face_submission branch####
  mesh.et.lms <- reactive({
    ex_mesh <- vcgSmooth(file2mesh("~/shiny/shinyapps/Perception/to_bake/da_reg.ply"))
    sample1k <- sample(1:27903, 200)
    #register landmarks to the space
    registered.mesh <- rotmesh.onto(ex_mesh, t(ex_mesh$vb[-4, sample1k]), synd.mshape[sample1k,], scale = T)$mesh
    registered.mesh
    })
  
  output$submitted_face <- renderRglwidget({
    pdf(NULL)
    dev.off()
    
    ex_mesh <- mesh.et.lms()
    if(is.null(input$submitted_heatmap)){
      print(shinyGetPar3d("userMatrix", session))
    par3d(userMatrix = front.face2, zoom = .75)
    plot3d(ex_mesh, col = "lightgrey", axes = F, specular = 1, xlab = "", ylab = "", zlab = "", aspect = "iso")  
    rglwidget()
    } else if(input$submitted_heatmap == "Compare to syndrome?"){
      #if heatmap, register face, make comp same age?, compare
      selected.sex <- input$current_sex
      selected.age <- input$current_age
      selected.synd <- input$submitted_comp
      
      raw_api_res <- httr::GET(url = paste0("http://localhost:3636", "/predshape"),
                               query = list(selected.sex = selected.sex, selected.age = selected.age, selected.synd = selected.synd),
                               encode = "json")
      
      pred.comp <- ex_mesh
      pred.comp$vb[-4,] <- t(jsonlite::fromJSON(httr::content(raw_api_res, "text", encoding = "UTF-8")))/1e10
      
      par3d(userMatrix = front.face2, zoom = .75)
      meshDist(ex_mesh, pred.comp)
      rglwidget()
      
    }
    
  })
  
  output$posterior_scree <- renderPlotly({
    
    projected.mesh <- getPCscores(t(mesh.et.lms()$vb[-4,]), PC.eigenvectors, synd.mshape)[1:200]
    # projected.mesh <- getPCscores(t(registered.mesh$vb[-4,]), PC.eigenvectors, synd.mshape)[1:200]
    
    #classify individual's scores using the model
    # colnames(projected.mesh) <- colnames(PC.scores)
    
    posterior.distribution <- predict(hdrda.mod, newdata = projected.mesh, type = "prob")$post
     
    posterior.distribution <- sort(posterior.distribution, decreasing = T)

    #used to be part of plot.df: ID = as.factor(1:10),
    plot.df <- data.frame(Probs = round(as.numeric(posterior.distribution[1:10]), digits = 4), Syndrome = as.factor(names(posterior.distribution[1:10])))
    plot.df$Syndrome <- as.character(plot.df$Syndrome)
    plot.df$Syndrome[plot.df$Syndrome == "Unrelated Unaffected"] <- "Non-syndromic"

    plot_ly(data = plot.df, x = ~Syndrome, y = ~Probs, type = "bar", color = I("grey"), hoverinfo = paste0("Syndrome: ", "x", "<br>", "Probability: ", "y")) %>%
      layout(xaxis = list(tickvals = gsub("_", " ", plot.df$Syndrome), tickangle = 45, ticktext = c(Syndrome = plot.df$Syndrome, Probability = plot.df$Probs), title = "<b>Syndrome</b>"),
             yaxis = list(title = "<b>Class probability</b>"),
             paper_bgcolor='white',
             margin = list(b = 125, l = 50, r = 100)
      )
  })
  
  output$personal_morphospace <- renderPlotly({
    
    projected.mesh <- getPCscores(t(mesh.et.lms()$vb[-4,]), PC.eigenvectors, synd.mshape)[1:200]
    # projected.mesh <- getPCscores(t(registered.mesh$vb[-4,]), PC.eigenvectors, synd.mshape)[1:200]
    
    selected.sex <- input$current_sex
    selected.age <- input$current_age

    raw_api_res <- httr::GET(url = paste0("http://localhost:3636", "/predPC"),
                             query = list(selected.sex = selected.sex, selected.age = selected.age),
                             encode = "json")

    parsed_res <- jsonlite::fromJSON(httr::content(raw_api_res, "text"))/1e10
    
    plot.df <- data.frame(Syndrome = c("Submitted mesh", levels(d.meta.combined$Syndrome)), Scores = rbind(projected.mesh[1:2], parsed_res))
    plot.df$Syndrome <- as.character(plot.df$Syndrome)
    plot.df$Syndrome[plot.df$Syndrome == "Unrelated Unaffected"] <- "Non-syndromic"
    
    plot_ly(data = plot.df, x = ~round(Scores.1, digits = 3), y = ~round(Scores.2,digits = 3), text = ~Syndrome, type = "scatter", marker = list(size = 10), 
            mode = "markers+text",
            textposition = 'top center', color = I(c(2, rep("#AE80E6", length(unique(d.meta.combined$Syndrome))))), hoverinfo = paste0("PC1 score: ", "x", "<br>", "PC2 score: ", "y")) %>%
      layout(xaxis = list(title = "<b>PC1 score</b>"),
             yaxis = list(title = "<b>PC2 score</b>"),
             paper_bgcolor='white',
             margin = list(b = 25, l = 25, r = 25, t = 25)
      ) %>%
      style(borderRadius = '15px')
  })
  
}

shinyApp(ui = ui, server = server)


