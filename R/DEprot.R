DEprot <- function(data,nbCond=NULL,nbRep=NULL,condName=NULL, normalization="SCALING",
                   pvalue.treshold = 0.05, fc.treshold = 1){

  ##applying different normalization depending on the normalization parameter
  normalization=toupper(normalization)
  if(identical(normalization, "NULL")){
  }
  else if(identical(normalization, "MEAN")){
    data<-scale(log2(data))
  }
  else if(identical(normalization, "QUANTILE")){
    library(preprocessCore)
    data <- log2(data)
    data <- normalize.quantiles(as.matrix(data))
  }
  else if(identical(normalization, "MEDIAN")){
    data<-log2(data)
    for (x in 1:length(data)){
      data[,x]<-(data[,x]-median(data[,x],na.rm = T))/(mad(data[,x],na.rm = T))
    }
  }
  else{
    data<-scale(log2(data))
  }

  col<-1
  ##cond is a list containing the columns number of each different conditions.
  cond<-list()
  ## the p value is stored in a 3d array, the two first dimmention beeing the comparated conditions and the third dimension the values for each protein.
  p.value<-array(dim = c(nbCond,nbCond,nrow(data)))
  ## the adjusted p value is stored in a 3d array, the two first dimmention beeing the comparated conditions and the third dimension the values for each protein.
  adjust.p.value<-array(dim = c(nbCond,nbCond,nrow(data)))
  ## the fold change is stored in a 3d array, the two first dimmention beeing the comparated conditions and the third dimension the values for each protein.
  fc<-array(dim = c(nbCond,nbCond,nrow(data)))

  ## resultsTable is a 2D matrix of data frames containing the differents values calculated for each comparison.
  resultsTable<-matrix(data.frame(), nrow = nbCond, ncol = nbCond)

  ## DEprotein is a 2D matrix of data frames containing the values for the significant proteins for each comparison.
  DEprotein<-matrix(data.frame(), nrow = nbCond, ncol = nbCond)

  ## if the parameter condName is null, or there is an error in the number of name in the vector, the name are set numbers.
  if(is.null(condName) | length(condName) != nbCond){
    condName<-c(1:nbCond)
  }

  ## Setting the cond list depending on the number of conditions and replicates in parameters
  for(i in 1:nbCond){
    cond[[i]]<-c(col:(col+nbRep-1))
    col<-col+nbRep
  }
  data<-as.data.frame(data[complete.cases(data),])

  ## Calcul of the p value and fold change for each condition, and setting the resultsTable for each comparisons
  for(i in 1:nbCond){
    for(j in 1:nbCond){
      ##student's test
      p.value[i,j,]<-apply(data,1,function(x){t.test(as.numeric(x[c(cond[[i]])]),as.numeric(x[c(cond[[j]])]), alternative = "t") $p.value})

      ## Benjamini hochberg correction of the p value
      adjust.p.value[i,j,]<-p.adjust(p.value[i,j,], method = "BH")
      ##fold change
      fc[i,j,]<-rowMeans(data[,cond[[i]]])-rowMeans(data[,cond[[j]]])

      resultsTable[[i,j]]<-data.frame(protein=c(row.names(data)))
      resultsTable[[i,j]][paste("p.value.",condName[i],".vs.",condName[j],sep="")]=c(p.value[i,j,])
      resultsTable[[i,j]][paste("adjust.p.value.",condName[i],".vs.",condName[j],sep="")]=c(adjust.p.value[i,j,])
      resultsTable[[i,j]][paste("fc.",condName[i],".vs.",condName[j],sep="")]=c(fc[i,j,])
    }
  }


  for(i in 1:nbCond){
    for(j in 1:nbCond){
      if(i!=j){

        ## Threshold
        N<-which(p.value[i,j,]<=pvalue.treshold & abs(fc[i,j,])>=fc.treshold)
        if(length(N)!=0){
          print(paste0("comparison ", condName[i], " against ", condName[j], " : ", nrow(data[N,]), " significantly expressed proteins."))
          #Volcano plot
          filename<-paste0("VolcanoPlot.",condName[i],".vs.",condName[j], '.svg')
          svg(filename)
          plot(fc[i,j,],
               -log10(adjust.p.value[i,j,]),
               main = 'Volcano plot',
               xlab =sprintf("fold change %s vs %s", condName[i],condName[j]),
               ylab ="-log10(adjusted p value)",
               abline(h = -log10(pvalue.treshold), v = c(-fc.treshold,fc.treshold)),
               pch = 20,
               col = ifelse(p.value[i,j,]<=pvalue.treshold & abs(fc[i,j,])>=fc.treshold, "red", "blue"))
          ##adding text on plot
          text(fc[i,j,][N], -log10(p.value[i,j,])[N], labels = row.names(data[N,]), cex= 0.7, pos = 2)
          dev.off()
          #saving significant protein into a dataframe (for each comparison)
          DEprotein[[i,j]]<-data.frame(protein=c(row.names(data[N,])))
          DEprotein[[i,j]][paste("p.value.",condName[i],".vs.",condName[j],sep="")]=c(p.value[i,j,N])
          DEprotein[[i,j]][paste("adjust.p.value.",condName[i],".vs.",condName[j],sep="")]=c(adjust.p.value[i,j,N])
          DEprotein[[i,j]][paste("fc.",condName[i],".vs.",condName[j],sep="")]=c(fc[i,j,N])
        }
      }
    }
  }

  results <-list( "normdata" = data,"dataprocessed"=resultsTable, "DEprotein"=DEprotein)
  return(results)
}
