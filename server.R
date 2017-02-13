library(shiny)
library(edgebundleR)

shinyServer(function(input, output) {
  
  build_plot <- function(in_file)
  {
    require(data.table)
    require(dplyr)
    require(igraph)
    
    # Data file name
    hbond_file <- paste('~/Desktop/Protein_analysis/hbonds/', in_file, '/all_hbond_red.txt', sep = '')
    map_file <- paste('~/Desktop/Protein_analysis/map/', in_file, '/glycan_map.txt', sep = '')
    total_frames <- 25000
    
    # Read in hbond data
    all_hbonds <- data.frame(fread(hbond_file, header = T))
    colnames(all_hbonds) <- c('Frame', 'bond')
    
    # Creat new cols with resnos
    all_hbonds$A <- as.numeric(gsub('.*_([0-9]+)@.*-.*_([0-9]+)@.*$', '\\1', all_hbonds$bond))
    all_hbonds$D <- as.numeric(gsub('.*_([0-9]+)@.*-.*_([0-9]+)@.*$', '\\2', all_hbonds$bond))
    
    # Remove residue-residue interactions
    all_hbonds <- subset(all_hbonds, A != D)
    
    # Read in map
    glycan_map <- read.csv(map_file, sep = "", stringsAsFactors = FALSE)
    
    # Match hbond data and map data
    all_hbonds$A_ref_n <- glycan_map$ref_n[match(all_hbonds$A, glycan_map$hetresno)]
    all_hbonds$D_ref_n <- glycan_map$ref_n[match(all_hbonds$D, glycan_map$hetresno)]
    
    # Remove intra-glycan interactions 
    all_hbonds <- subset(all_hbonds, A_ref_n != D_ref_n)
    
    # Order interaction  
    all_hbonds$large <- with(all_hbonds, ifelse(A_ref_n > D_ref_n, A, D))
    all_hbonds$small <- with(all_hbonds, ifelse(A_ref_n < D_ref_n, A, D))
    
    # Select only relavant cols
    all_hbonds <- subset(all_hbonds, select = c('Frame', 'small', 'large'))
    
    # Match hbond data and map data
    glycan_map_n <- names(glycan_map)
    names(glycan_map) <- paste('small', glycan_map_n, sep = '_')
    all_hbonds <- merge(all_hbonds, glycan_map, by.x = 'small', by.y = 'small_hetresno', all.x = T)
    all_hbonds$small_label <- paste(all_hbonds$small_mon, sprintf('%03s', all_hbonds$small_ref_s), sep = '_')
    
    names(glycan_map) <- paste('large', glycan_map_n, sep = '_')
    all_hbonds <- merge(all_hbonds, glycan_map, by.x = 'large', by.y = 'large_hetresno', all.x = T)
    all_hbonds$large_label <- paste(all_hbonds$large_mon, sprintf('%03s', all_hbonds$large_ref_s), sep = '_')
    
    # Remove duplicate inter glycan interactions  
    all_hbonds <- all_hbonds %>% 
      group_by(Frame, small_label, large_label) %>% 
      summarise(cross_mon = min(ifelse(small_mon != large_mon, 1, 0)))
    
    net_data <- all_hbonds %>%
      group_by(small_label, large_label) %>%
      summarise(frac = length(Frame)/total_frames, cross_mon = min(cross_mon))
    net_data$lty <- cut(net_data$frac, breaks = c(0,0.25,0.5,0.75,1), labels = c("1,1","3,5,1","3,3","1,0"))
    net_data$opacity <- 0.5
    
    vertices <- unique(c(net_data$small_label, net_data$large_label))
    vertice_names <- gsub('([ABC])_([0-9]+.*)', '\\2 \\1', vertices)
    d <- structure(list(ID = vertices,
                        Type = substr(vertice_names, regexpr("[^0]", vertice_names), nchar(vertice_names))),
                   .Names = c("ID", "vertice_names"), class = "data.frame",
                   row.names = c(NA, -length(vertices)))
    
    
    g <- graph.data.frame(net_data, directed = F, vertices = d)
    return(g)
  }  
  
  output$distPlot <- renderEdgebundle({edgebundle(build_plot(input$variable), padding = 100)})
})





