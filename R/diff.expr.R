diff.expr <- function( data , thresh_fc = 0.5, thresh_p = 0.05){
  if(!is.list(data)){stop("data must be a list of data frame generated with the function test.stat")}
  de <- data
  for(i in 1:length(de)){
    fc = as.data.frame(de[[i]])[,4]
    p = as.data.frame(de[[i]])[,3]
    dt <-as.data.frame(de[[i]])
    de[[i]] <- dt[which(p<=thresh_p & abs(fc)>=thresh_fc),]
  }
  return(de)
}