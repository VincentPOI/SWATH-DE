

norm <- function(data, normalization="MEAN"){
  if(!is.data.frame(data) & !is.matrix(data)){stop("the data must be a data.frame or a matrix")}
  
  normalization=toupper(normalization)
  if(identical(normalization, "MEAN")){
    data<-scale(log2(data))
  }
  else if(identical(normalization, "QUANTILE")){
    data <- log2(data)
    nomCol <- colnames(data)
    data <- as.data.frame(normalize.quantiles(as.matrix(data)), row.names = row.names(data))
    colnames(data) <- nomCol
  }
  else if(identical(normalization, "MEDIAN")){
    data<-log2(data)
    for (x in 1:length(data)){
      data[,x]<-(data[,x]-median(data[,x],na.rm = T))/(mad(data[,x],na.rm = T))
    }
  }
  else{
    stop("normalization must be : 'MEAN', 'QUANTILE' or 'MEDIAN'")
  }
  return(data)
}