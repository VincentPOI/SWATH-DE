volcanoPlot <- function(data , save = FALSE, fc=0.5, p=0.05 ){ 
  require(ggplot2)
  require(ggrepel)
  require(svglite)
  for (i in 1:length(data)){
  
  plotTitle <- substr(colnames(data[[i]])[2], 9 ,nchar(colnames(data[[i]])[2]))
  values <- as.data.frame(data[[i]])
  forplot <- data.frame(x=as.numeric(values[,4]), y=-log10(values[,3]), id=as.character(values[,1]))
  tmp <- forplot[as.numeric(forplot$y)>=-log10(fc) & abs(forplot$x)>fc,]
  p <- ggplot(forplot) + geom_point(aes(x, y, label= id , color = ifelse(y>=-log10(p) & abs(x)>=fc, "not signi", "FC")),show.legend = F) +
    scale_color_manual(values = c("blue", "red")) +
    geom_text_repel(data = subset(forplot, abs(forplot$x)>=fc & forplot$y>=-log10(p)),
                    aes(x,y,label = id),
                    size = 2) +
    geom_vline(xintercept = fc ) +
    geom_vline(xintercept = -fc) + 
    geom_hline(yintercept = -log10(p)) + 
    labs(title = plotTitle,x="log2(Fold-change)", y="-log10(P.Value)") + theme_bw() 
    print(p)
  
  
  if(save){
    ggsave(paste0("VolcanoPlot-",plotTitle,".svg"),plot = p)
  }
  }
}
