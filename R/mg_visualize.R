## Andrea Hita
## JUN-2021

## Define Multi-Graph as igraph object
## -----------------------------------------------------------

mg_build  <- function(inM, ml, attr){
  
  ## Scale adjacency matrix
  D <- diag(inM)
  sM <- inM
  sM <- sM/diag(sM)
  sM[is.na(sM)] <- 0
  
  ## (If no multimappers, no graph is built)
  if(all(sM == 0)){return(0)}
  
  ## Generate multi-mappers graph
  g <- graph_from_adjacency_matrix(sM, weighted = TRUE, mode = "directed", diag = FALSE)
  E(g)$weight[E(g)$weight > 1] <- 1.0
  
  ## Define graph visualization parameters
  V(g)$weight <- log10(D)
  V(g)$size <- scales::rescale(log10(D), to = c(1,5))
  V(g)$frame.color <-"white"; V(g)$color <- "black"; 
  V(g)$feat <- as.character(ml[[attr]])
  V(g)$label <- ''
  V(g)$label.cex <- 0.55
  E(g)$arrow.mode <- 0;  E(g)$curved <- 0.2
  E(g)$width <- scales::rescale(E(g)$weight, to = c(0.1,2))
  E(g)$color <- 'grey50'
  
  ## Define multi-loci groups
  cl <- as.character(ml$community_id)
  cl[cl == ""] <- paste0('ucl-',c(1:sum(cl == "")))
  V(g)$clust <- cl
  V(g)$ml <- ml$community_flag == "True"
  
  return(g)}


## Generate Multi-Graph exploratory plots
## -----------------------------------------------------------

mg_plotset <- function(plotfile, g){
  
  gc <- delete_vertices(g, V(g)[!V(g)$ml])
  
  l <- layout_nicely(g)
  lc <- layout_nicely(gc)
  cl <- lapply(unique(V(g)$clust), function(x) V(g)[V(g)$clust == x])
  clc <- lapply(unique(V(gc)$clust), function(x) V(gc)[V(gc)$clust == x])
  
  ci <- factor(unlist(lapply(cl, function(x) length(x)==1)), levels = c(FALSE, TRUE))
  
  V(g)$color <- V(gc)$color <- 'black'
  
  png(paste0(plotfile,'.png'), res = 300, units = 'in', width = 6, height = 6)
  par(mar=c(0.25,0.25,0.25,0.25))
  plot(g, layout = l)
  dev.off()
  
  png(paste0(plotfile,'_cl.png'), res = 300, units = 'in', width = 6, height = 6)
  par(mar = c(0.25,0.25,0.25,0.25))
  plot(g, layout = l, mark.groups = cl, mark.border = c('white', rgb(1, 1, 1, alpha = 0))[ci],
       mark.col =  c(rgb(0.5, 0.5, 0.5, alpha = 0.15),  rgb(1, 1, 1, alpha = 0))[ci]) 
  dev.off()
  
  png(paste0(plotfile,'c.png'), res = 300, units = 'in', width = 6, height = 6)
  par(mar = c(0.25,0.25,0.25,0.25))
  plot(gc, layout = lc) 
  dev.off()
  
  png(paste0(plotfile,'c_cl.png'), res = 300, units = 'in', width = 6, height = 6)
  par(mar = c(0.25,0.25,0.25,0.25))
  plot(gc, layout = lc,  mark.groups = clc, mark.col = rgb(0.5, 0.5, 0.5, alpha = 0.2), mark.border = 'white') 
  dev.off()
  
  tmplabel <-  V(g)$label
  tmplabelc <- V(gc)$label
  V(g)$label <- V(g)$feat
  V(gc)$label <- V(gc)$feat
  V(g)$color <- V(gc)$color <- 'goldenrod2'
  
  png(paste0(plotfile,'_lab.png'), res = 300, units = 'in', width = 6, height = 6)
  par(mar=c(0.25,0.25,0.25,0.25))
  plot(g, layout = l)
  dev.off()
  
  png(paste0(plotfile,'_cl_lab.png'), res = 300, units = 'in', width = 6, height = 6)
  par(mar = c(0.25,0.25,0.25,0.25))
  ci <- factor(unlist(lapply(cl, function(x) length(x)==1)), levels = c(FALSE, TRUE))
  plot(g, layout = l, mark.groups = cl, mark.border = c('white', rgb(1, 1, 1, alpha = 0))[ci],
       mark.col =  c(rgb(0.5, 0.5, 0.5, alpha = 0.15),  rgb(1, 1, 1, alpha = 0))[ci]) 
  dev.off()
  
  png(paste0(plotfile,'c_lab.png'), res = 300, units = 'in', width = 6, height = 6)
  par(mar = c(0.25,0.25,0.25,0.25))
  plot(gc, layout = lc) 
  dev.off()
  
  png(paste0(plotfile,'c_cl_lab.png'), res = 300, units = 'in', width = 6, height = 6)
  par(mar = c(0.25,0.25,0.25,0.25))
  plot(gc, layout = lc,  mark.groups = clc, mark.col = rgb(0.5, 0.5, 0.5, alpha = 0.2), mark.border = 'white') 
  dev.off()
  
  V(g)$label <- tmplabel
  V(gc)$label <- tmplabelc
  
  if(!is.null(V(g)$customColor)){
    
    png(paste0(plotfile,'_combined_color.png'), res = 300, units = 'in', width = 20, height = 6)
    par(mfrow = c(1,4), mar = c(0.1, 0.1, 0.1, 0.1))
    plot(g, layout = l)
    V(g)$color <- V(g)$customColor
    plot(g, layout = l )
    plot(gc,  mark.groups = clc, mark.col = rgb(0.5, 0.5, 0.5, alpha = 0.2), mark.border = 'white', layout = lc)
    V(gc)$color <- V(gc)$customColor
    plot(gc, mark.groups = clc, mark.col = rgb(0.5, 0.5, 0.5, alpha = 0.2), mark.border = 'white', layout = lc)
    dev.off()
    
    par(mar=c(0.5,0.5,0.5,0.5))
    png(paste0(plotfile,'_color.png'), res = 300, units = 'in', width = 6, height = 6)
    par(mar = c(0.25,0.25,0.25,0.25))
    plot(g, layout = l)
    dev.off()
    
    png(paste0(plotfile,'_cl_color.png'), res = 300, units = 'in', width = 6, height = 6)
    par(mar = c(0.25,0.25,0.25,0.25))
    plot(g, layout = l, mark.groups = cl, mark.border = c('white', rgb(1, 1, 1, alpha = 0))[ci],
         mark.col =  c(rgb(0.5, 0.5, 0.5, alpha = 0.15),  rgb(1, 1, 1, alpha = 0))[ci]) 
    dev.off()
    
    png(paste0(plotfile,'c_color.png'), res = 300, units = 'in', width = 6, height = 6)
    par(mar = c(0.25,0.25,0.25,0.25))
    plot(gc, layout = lc) 
    dev.off()
    
    png(paste0(plotfile,'c_cl_color.png'), res = 300, units = 'in', width = 6, height = 6)
    par(mar = c(0.25,0.25,0.25,0.25))
    plot(gc, layout = lc,  mark.groups = clc, mark.col = rgb(0.5, 0.5, 0.5, alpha = 0.2), mark.border = 'white')
    dev.off()
    
    V(g)$label <- V(g)$feat
    V(gc)$label <- V(gc)$feat
    
    png(paste0(plotfile,'_combined_color_lab.png'), res = 300, units = 'in', width = 20, height = 6)
    par(mfrow = c(1,4), mar = c(0.1, 0.1, 0.1, 0.1))
    plot(g, layout = l)
    V(g)$color <- V(g)$customColor
    plot(g, layout = l )
    plot(gc,  mark.groups = clc, mark.col = rgb(0.5, 0.5, 0.5, alpha = 0.2), mark.border = 'white', layout = lc)
    V(gc)$color <- V(gc)$customColor
    plot(gc, mark.groups = clc, mark.col = rgb(0.5, 0.5, 0.5, alpha = 0.2), mark.border = 'white', layout = lc)
    dev.off()
    
    par(mar=c(0.5,0.5,0.5,0.5))
    png(paste0(plotfile,'_color_lab.png'), res = 300, units = 'in', width = 6, height = 6)
    par(mar = c(0.25,0.25,0.25,0.25))
    plot(g, layout = l)
    dev.off()
    
    png(paste0(plotfile,'_cl_color_lab.png'), res = 300, units = 'in', width = 6, height = 6)
    par(mar = c(0.25,0.25,0.25,0.25))
    plot(g, layout = l, mark.groups = cl, mark.border = c('white', rgb(1, 1, 1, alpha = 0))[ci],
         mark.col =  c(rgb(0.5, 0.5, 0.5, alpha = 0.15),  rgb(1, 1, 1, alpha = 0))[ci]) 
    dev.off()
    
    png(paste0(plotfile,'c_color_lab.png'), res = 300, units = 'in', width = 6, height = 6)
    par(mar = c(0.25,0.25,0.25,0.25))
    plot(gc, layout = lc) 
    dev.off()
    
    png(paste0(plotfile,'c_cl_color_lab.png'), res = 300, units = 'in', width = 6, height = 6)
    par(mar = c(0.25,0.25,0.25,0.25))
    plot(gc, layout = lc,  mark.groups = clc, mark.col = rgb(0.5, 0.5, 0.5, alpha = 0.2), mark.border = 'white')
    dev.off()
  }
}


mg_interactive <- function(mg, plotfile){
  
  vs <- V(mg) #Get node list
  vsize <- scales::rescale(vs$size, to = c(6,12))
  vlabel <- V(mg)$feat
  
  es <- as.data.frame(get.edgelist(mg))
  
  Nv <- length(vs) #number of nodes
  Ne <- length(es[1]$V1) #number of edges
  
  Ncl <- length(unique(V(mg)$clust))
  colP <-  c(rainbow(Ncl))[as.factor(V(mg)$clust)]
  
  #Coordinates for nodes
  L <- layout_nicely(mg)
  Xn <- L[,1]
  Yn <- L[,2]
  
  #Creates the nodes (plots the points)
  network <- plot_ly(x = ~Xn, y = ~Yn, #Node points
                     mode = 'markers', 
                     text = vlabel, 
                     hoverinfo = 'text',
                     color = colP,
                     size = vsize)
  
  #Create edges
  edge_shapes <- list()
  
  for(i in c(1:Ne)){
    
    v0 <- es[i,1]
    v1 <- es[i,2]
    
    edge_shape = list(
      type = 'line',
      line = list(color = 'gray', width = E(mg)$weight[i]),
      x0 = Xn[v0],
      y0 = Yn[v0],
      x1 = Xn[v1],
      y1 = Yn[v1]
    )
    
    edge_shapes[[i]] <- edge_shape
  }
  
  axis <- list(title = '', showgrid = FALSE, 
               showticklabels = FALSE, zeroline = FALSE)
  
  p <- layout(
    network,
    title = 'Networks & Plotly',
    shapes = edge_shapes,
    xaxis = axis,
    yaxis = axis,
    showlegend=FALSE)
  
  htmlwidgets::saveWidget(p, paste0(plotfile,'.html'))
  return(p)
}