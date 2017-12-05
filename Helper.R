checklibs<-function(libs){
  toinstall <- libs[!(libs %in% installed.packages()[,"Package"])]
  if(length(toinstall) > 0){
    print(paste("Installing packages",toinstall))
    install.packages(toinstall, repos = 'http://cran.rstudio.com', dependencies=TRUE)
  }
  toinstall <- libs[!(libs %in% installed.packages()[,"Package"])]
  if(length(toinstall) > 0){
    print(paste("Trying to install packages",toinstall,"from bioconductor"))
    source("https://bioconductor.org/biocLite.R")
    biocLite(,suppressUpdates=TRUE)
    biocLite(toinstall,suppressUpdates=TRUE)
  }
  for(lib in libs){
    lapply(lib, library, character.only = TRUE)
  }
}

libs=c("TxDb.Hsapiens.UCSC.hg19.knownGene","GO.db","impute","preprocessCore",
       "org.Hs.eg.db","GenomeInfoDb","WGCNA","VariantAnnotation") #"stats","R.utils"
checklibs(libs)

#Define convenient helper functions
suppressWarnings(library(WGCNA,warn.conflicts = FALSE, quietly=TRUE,verbose=FALSE))
suppressWarnings(allowWGCNAThreads())
suppressWarnings(enableWGCNAThreads())
# Possibly use Hmisc::Rcorr

# Stupid as.array doesn't work as smoothly as one would hope, 
# so I created this convenince function instead
convert2Array=function(x){
  return(array(unlist(x), dim=dim(x), dimnames=dimnames(x)))
}

# convenience function for doing a single correlation
mycor=function(d,method="pearson"){
  x=d[,"prot"]
  y=d[,"rna"]
  ret=cor.test(x,y,method=method, exact=TRUE)
  ret=c(ret$estimate,ret$p.value)
  names(ret)=c("correlations","Pval")
  return(ret)
}

# Convenience function for doing all correlations 
doCorrAllPairs = function(x,y, name, xid, yid, method="pearson", adjustP="None"){
  # use apply to do the correlation over marg
  tmp="tmp"
  if(method=="pearson"){
    tmp <- corAndPvalue(t(x),t(y), method=method, use= "all.obs", use="pairwise.complete.obs",pearsonFallback="individual")
  }
  else if(method=="bicor"){
    tmp <- bicorAndPvalue(t(x),t(y), use="all.obs")
    names(tmp)=gsub("bicor","cor",names(tmp))
    }
  else{
    stop(paste("Correlation method",method,"not handled"))
  }
  ret=data.frame(expand.grid(rownames(tmp$p),colnames(tmp$p)))
  names(ret)=c(xid,yid)
  ret[,paste(name,"p",adjustP,sep=".")]=p.adjust(tmp$p,method=adjustP)
  ret[,paste(name,"cor",sep=".")]=c(tmp$cor)
  ret[,paste(name,"p",sep=".")]=c(tmp$p)
  return(ret)
}

# Convenience function for doing all correlations 
doCorrDiagPairs = function(data, marg=c(2), meth="pearson", meta=NULL, indexname=NULL){
  # use apply to do the correlation over marg
  ret=as.data.frame(t(apply(data,marg,mycor,method=meth)))
  #  ret=data.frame(correlations=apply(data,marg,mycor,method=meth))
  # add the meta data and make names a proper column
  if(!is.null(meta)){
    ret = merge(ret,meta, by="row.names")
    rownames(ret) = ret$Row.names
    ret$Row.names=NULL
    # if(!is.null(indexname)){
    #   names(ret)[1]=indexname
    # }
  }
  # else{
    if(!is.null(indexname)){
      order=c(indexname,names(ret))
      ret[,indexname]=row.names(ret)
      ret=ret[,order]
    }
  # }
  ret$Padjust=p.adjust(ret$Pval, method="BH")
  return(ret)
}

# get cv,sd,iqr of data
getStat = function(data, marg=c(2)){
  # use apply to do the correlation over dimension 'marg'
  stat=data.frame(mean=apply(data,marg,mean))
  stat$sd=apply(data,marg,sd)
  stat$cv=abs(stat$sd/stat$mean)
  stat$iqr=apply(data,marg,IQR)
  return(stat)
}

# Function for getting genomic positions from gene symbol or vice versa
# global variables 
Symbol2id <- as.list( org.Hs.egSYMBOL2EG )
id2pos <- as.list( org.Hs.egCHRLOC[mappedkeys(org.Hs.egCHRLOC)] )
txdb=TxDb.Hsapiens.UCSC.hg19.knownGene
id2Symbol <- rep( names(Symbol2id), sapply(Symbol2id, length) )
names(id2Symbol) <- unlist(Symbol2id)

# the funcitons
getGene=function(name,chr,pos){
  input = data.frame(id=name,
                     chr=paste("chr",chr,sep=""), 
                     pos=pos)
  target <- with(input, GRanges( seqnames = Rle(chr),
                                 ranges   = IRanges(pos, end=pos, names=id),
                                 strand   = Rle(strand("*")) ) )
  loc <- locateVariants(target, txdb, CodingVariants())
#  loc <- locateVariants(target, txdb, AllVariants())
  names(loc) <- NULL
  out <- as.data.frame(loc)
  out$chrpos <- names(target)[ out$QUERYID ]
  out <- out[ , c("chrpos","GENEID")]
  out <- unique(out[!is.na(out$GENEID),])
  
  # get table of gene id and gene symbol and index by gene id
  out$GENESYMBOL <- id2Symbol[ as.character(out$GENEID) ]
  
  return(out)
}

getGenePos=function(genesymb){
  ids = Symbol2id[genesymb]
  ret = unlist(id2pos[unlist(ids)])
  ret=ret[!duplicated(names(ret))]
  ret = rbind(ret,chr=gsub("^.*\\.","",names(ret)))
  t=id2Symbol[gsub("\\..*$","",colnames(ret))]
  names(t) = NULL
  colnames(ret)=t
  #  colnames(ret)=id2Symbol[gsub("\\..*$","",colnames(ret))]
  rownames(ret)[rownames(ret)=="ret"] ="pos"
  ret=t(ret)
  return(ret[ret[,"chr"]!="Y" & (!grepl("_",ret[,"chr"])),])
}

Logicalcombination=function(X,Y){
  ret=rep("4TrueFalse",length(X))
  for(i in seq(1,length(X))){
    x=X[i]
    y=Y[i]
    if(x==FALSE & y==FALSE){
      ret[i]="1FalseFalse"
    }
    if(x==TRUE & y==TRUE){
      ret[i]="2TrueTrue"
    }
    if(x==FALSE & y==TRUE){
      ret[i]="3FalseTrue"
    }
  }
  return(ret)
}
