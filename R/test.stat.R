test.stat <- function(data, stat = "t.test" , design , contrast){
  if(!is.data.frame(data) & !is.matrix(data)){stop("the data must be a data.frame or a matrix")}
  if(!is.data.frame(design) & !is.matrix(design)){stop("the experience design must be a matrix generated from model.matrix")}
  if(!is.data.frame(contrast) & !is.matrix(contrast)){stop("the contrast matrix must be  a matrix generated from makeContrasts()")}


  nbComp<-ncol(contrast)
  fc <- list()
  p.value <- list()
  adjust.p.value<-list()
  listResults<-list()
  
  for (i in 1:nbComp){
    
    op<-row.names(contrast)[contrast[,i] == 1 ]
    ref<-row.names(contrast)[contrast[,i] == -1 ]
    
    samplesop <- row.names(design)[row(as.matrix(design[,op]))[design[,op]==1]]
    samplesref <- row.names(design)[row(as.matrix(design[,ref]))[design[,ref]==1]]
    
    colop <- which(colnames(data) %in% samplesop)
    colref <- which(colnames(data) %in% samplesref)
    
    p.value[[i]]<-apply(data,1,function(x){t.test(as.numeric(x[colref]),as.numeric(x[colop]), alternative = "t") $p.value})
    adjust.p.value[[i]]<-p.adjust(p.value[[i]], method = "BH")
    
    fc[[i]]<-rowMeans(data[,colop])-rowMeans(data[,colref])
    
    listResults[[i]]<-data.frame(protein=c(row.names(data)))
    listResults[[i]][paste("p.value.",colnames(contrast)[i],sep="")]=c(p.value[[i]])
    listResults[[i]][paste("adjust.p.value.",colnames(contrast)[i],sep="")]=c(adjust.p.value[[i]])
    listResults[[i]][paste("fc.",colnames(contrast)[i],sep="")]=c(fc[[i]])
    
  }
  return(listResults)
}

# truc <- SWATH.example.data
# 
# 
# truc <- norm(truc)
# 
# aaa <- test.stat(truc,exp.design = exp.design, contrast = contrast)
# 
# exp.design <- data.frame(samples = colnames(SWATH.example.data), condition = 1)
# exp.design$condition[1:10] = "control"
# exp.design$condition[11:20] = "lessone"
# exp.design$condition[21:30] = "onetofive"
# exp.design$condition[31:40] = "adult"
# 
# 
# 
# design <- model.matrix(~0 + exp.design$condition, data = exp.design)
# colnames(design) <- sort(unique(exp.design$condition))
# row.names(design) <- exp.design$samples
# 
# contrast <- makeContrasts(adult-control, lessone-control, onetofive-control, adult+lessone+onetofive-control, levels=design) 
# 
# 
# nbCond<-ncol(design)
# nbComp<-ncol(contrast)
# fc <- list()
# p.value <- list()
# adjust.p.value<-list()
# resultsTable<-list()
# 
# 
# for (i in 1:nbComp){
#       print(levels(contrast))
#       print(contrast[,i])
#       print("******")
#       print(row.names(contrast)[contrast[,i] == 1 ])
#       op<-row.names(contrast)[contrast[,i] == 1 ]
#       print(row.names(contrast)[contrast[,i] == -1 ])
#       ref<-row.names(contrast)[contrast[,i] == -1 ]
#       
#       print(row.names(design)[row(as.matrix(design[,op]))[design[,op]==1]])
#       samplesop <- row.names(design)[row(as.matrix(design[,op]))[design[,op]==1]]
#       print(row.names(design)[row(as.matrix(design[,ref]))[design[,ref]==1]])
#       samplesref <- row.names(design)[row(as.matrix(design[,ref]))[design[,ref]==1]]
#       
#       
#       print(which(colnames(truc) %in% samplesop))
#       colop <- which(colnames(truc) %in% samplesop)
#       print(which(colnames(truc) %in% samplesref))
#       colref <- which(colnames(truc) %in% samplesref)
#       print("")
#       print("")
# 
#       
#       p.value[[i]]<-apply(truc,1,function(x){t.test(as.numeric(x[colref]),as.numeric(x[colop]), alternative = "t") $p.value})
#       adjust.p.value[[i]]<-p.adjust(p.value[[i]], method = "BH")
# 
#       fc[[i]]<-rowMeans(truc[,colop])-rowMeans(truc[,colref])
#       
#       resultsTable[[i]]<-data.frame(protein=c(row.names(truc)))
#       resultsTable[[i]][paste("p.value.",colnames(contrast)[i],sep="")]=c(p.value[[i]])
#       resultsTable[[i]][paste("adjust.p.value.",colnames(contrast)[i],sep="")]=c(adjust.p.value[[i]])
#       resultsTable[[i]][paste("fc.",colnames(contrast)[i],sep="")]=c(fc[[i]])
# 
# }



# fc<-array(dim = c(nbCond,nbCond,nrow(truc)))
# p.value<-array(dim = c(nbCond,nbCond,nrow(truc))) 
# adjust.p.value<-array(dim = c(nbCond,nbCond,nrow(truc)))
# resultsTable<-matrix(data.frame(), nrow = nbCond, ncol = nbCond)
# for (i in unique(exp.design$condition)){
#   cptj<-1
#   for (j in unique(exp.design$condition)){
#     if(i != j){
#       
#       # print(i)
#        rowc1 <- which(exp.design$condition == i)
#        # print(rowc1)
#        rowc2 <- which(exp.design$condition == j)
#       # print(colnames(truc[,rowc1]))
#       # print(colnames(truc[,rowc2]))
#       # print(exp.design[rowc1,]$samples)
#       # print(exp.design[rowc2,]$samples)
#       print(rowc1)
#       print(which(colnames(truc) %in% as.character(exp.design[rowc1,]$samples)))
# 
#       
#       p.value[cpti,cptj,]<-apply(truc,1,function(x){t.test(as.numeric(x[c(which(exp.design$condition == i))]),as.numeric(x[c(which(exp.design$condition == j))]), alternative = "t") $p.value})
#       adjust.p.value[cpti,cptj,]<-p.adjust(p.value[cpti,cptj,], method = "BH")
#       
#       fc[cpti,cptj,]<-rowMeans(truc[,which(exp.design$condition == i)])-rowMeans(truc[,which(exp.design$condition == j)])
#     
#       resultsTable[[cpti,cptj]]<-data.frame(protein=c(row.names(truc)))
#       resultsTable[[cpti,cptj]][paste("p.value.",i,".vs.",j,sep="")]=c(p.value[cpti,cptj,])
#       resultsTable[[cpti,cptj]][paste("adjust.p.value.",i,".vs.",j,sep="")]=c(adjust.p.value[cpti,cptj,])
#       resultsTable[[cpti,cptj]][paste("fc.",i,".vs.",j,sep="")]=c(fc[cpti,cptj,])
#       
#       }
#     cptj<-cptj+1
#   }
#   cpti<-cpti+1
# }


# lima normalization :
# library(limma)
# library(edgeR)
# dge <- DGEList(truc)
# dge <- calcNormFactors(dge)
# y <- voom(dge, design)
# norm.expr <- y$E
# 
# boxplot(log2(truc), 
#         main="Distribution of unnormalised data",
#         xlab="",
#         ylab="log2 intensity",
#         las=2,cex.axis=0.8)  
# boxplot(norm.expr, 
#         main="Distribution of normalised data",
#         xlab="",
#         ylab="log2 normalised intensity",
#         las=2,cex.axis=0.8) 
