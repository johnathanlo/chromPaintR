GraphScaff<-function(chromdat, scaffname){
  require(ggplot2);require(ggpubr);require(RColorBrewer);require(dplyr)
  if(scaffname %in% names(chromdat$summary$`Database hits`)){
    print("Graphing...")
    
    dbhits <- chromdat$summary$`Database hits`[[scaffname]]
    qlocs <- chromdat$summary$`Locations on query scaffold`[[scaffname]]  
    dblocs <- chromdat$summary$`Locations on database`[[scaffname]]
    scaff.length <- as.numeric(chromdat$lengthdata[[2]][grep(scaffname, chromdat$lengthdata[[1]])])
    scaff <- data.frame(dbhits = dbhits, qlocs = qlocs, dblocs = dblocs)
    
    dbchroms <- unique(scaff$dbhits)
    colors <- brewer.pal(length(dbchroms), "Set3")
    names(colors) <- dbchroms
    colorScale <- scale_color_manual(name = "Chromosomes", values = colors)
    
    chromplot <- ggplot(data=scaff, aes(x = 1, y=qlocs, color = dbhits )) +
      geom_rect(data = scaff, aes(xmin = .9, xmax = 1.1, ymin = 0, ymax = scaff.length, fill = I("white")), color = NA) +
      geom_segment(aes(x= .9, y = qlocs, xend = 1.1, yend = qlocs), size = .5) +
      scale_x_continuous(name = "Scaffold", limits = c(0,2))+
      scale_y_continuous(name = "Position", labels = fancy_scientific)+
      colorScale +
      guides(colour = guide_legend(override.aes = list(size = 10)))
    chromplot
  }else{
    print("Scaffold not found.")
  }

}