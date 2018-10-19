 # Setup libraires
libs=c("readxl","matrixStats","data.table")
checklibs(libs)
toinstall <- libs[!(libs %in% installed.packages()[,"Package"])]
if(length(toinstall) > 0){
  print(paste("Installing packages",toinstall))
  install.packages(toinstall)
}
lapply(libs, library, character.only = TRUE)

indir="./Data/"
outdir="./R-files/"
dir.create(outdir, showWarnings=FALSE)

# Three convenience functions

# identifies rows mapping to same GeneName and replaces these by the row 
# containingtheir median values for each column
mymedian=function(x){
  ret=data.table(x)
  #setkey(ret, GeneName)
  ret=as.data.frame(ret[, lapply(.SD,median),by=GeneName, on=NULL])
  row.names(ret)=ret$GeneName
  ret$GeneName = NULL
  return(ret)
}

# Stupid as.array doesn't work as smoothly as one would hope, so created this function instead
convert2Array=function(x){
  if(is.null(dim(x))){  # for vectors
    return(array(unlist(x)))
  }
  else{ # for lists
    return(array(unlist(x), dim=dim(x), dimnames=sapply(dimnames(x),make.names,allow_=FALSE)))
  }
}

# Wrapper for reading excel files 
readArrayFromExcel = function(file, cols.to.remove=c(), row.name.col=c(), sheet=1){
  library(readxl)
  ret=as.data.frame(read_excel(file,sheet=sheet, col_names = TRUE))
  # Set column row.name.col as row names and convert to array
  if(!is.null(row.name.col)){
    row.names(ret)=ret[,row.name.col]
  }
  ret=ret[,!names(ret) %in% c(row.name.col,cols.to.remove)]
  ret=convert2Array(ret)
  return(ret)
}

readListFromExcel = function(file, sheet=1){
  library(readxl)
  raw = as.data.frame(read_excel(file,sheet=sheet, col_names = TRUE))
  ret = list()
  for(key in names(raw)){
    tmp=as.vector(raw[key])
    ret[[key]]=tmp[!is.na(tmp),]
  }
  return(ret)
}

#setwd("/Users/bengts/WABI/Lehtio/j_lehtio_1411")
getwd()

# PROTEIN DATA
##############
# Read the expression data extracted from the original preoteome data from 
# Henrik Johansson, Lehtio group
prot_data = readArrayFromExcel(paste(indir,"diffProtdata_extracted_from_9995_BC_45_proteome_Gene_symbol_centric_1FDR.xlsx", sep=""), row.name.col="GeneSymbol")
prot_data=apply(prot_data,c(1,2),log2)
head(prot_data, 5)
# Correct 
typos=c("OSL2U.0334","OSL2U.0407","OSL2U.0030","OSL2U.0484","OSL2U.0289","OSL2U.0364","OSL2U.0429")
newcolnames=colnames(prot_data)
for(typo in typos){
  newcolnames=gsub(typo,paste(typo,"T1",sep=""), newcolnames)
}
colnames(prot_data) = newcolnames

# MRNA ADATA
############
# First we need to udate the the original RNA data from Kristine Kleivi in 
# Oslo to GRCh37 gene symbols -- this will also update the CNV files, to 
# be used below. 
source("updating_gene_symbols.R")

# This puts the new data file in the outdir, so read the updated data from there
rna_data = read.table(paste(outdir,"main_qnorm_missimp_cent-hosp_log2_annot_updated_ids_to_Janne_updated.txt.gz",sep=""), header=TRUE, sep="\t", stringsAsFactors=FALSE)

# Correct erroneous value jan.27 with 1.3527 as per Miriam email
# Notice that this error caused the col to be coerced to factor and we need to
# convert back to double
rna_data$OSL2U.0383T1[rna_data$ProbeUID.SampleArray.gDetrendedSignal=="12629"] = 1.3527
rna_data$OSL2U.0383T1=as.double(rna_data$OSL2U.0383T1)

# remove duplicate according to Kristines email
rna_data$OSL2U.0219T1 = NULL
# remove superfluous columns
rna_data=rna_data[,!names(rna_data) %in% c("ProbeUID.SampleArray.gDetrendedSignal","ProbeName","SystematicName","chrom","seq_beg","seq_end","accessions")]
# replace the rows of individual probes mapping to the same genename by their median 
# row (for each column) and finally convert to array
head(rna_data, 5)
rna_data=convert2Array(mymedian(rna_data))

print(indir)
# META DATA
###########
# read in the tumour meta data, i.e., PAM50 classification of the tumours
meta_tumours=readArrayFromExcel(paste(indir,
                                      "151111_clin_info_connect_protoemics.xlsx",
                                      sep=""),
                                row.name.col= "OSL2_full_name")
# Correct 
typos=c("OSL2U.0334","OSL2U.0407","OSL2U.0030","OSL2U.0484","OSL2U.0289","OSL2U.0364","OSL2U.0429")
newrownames=rownames(meta_tumours)
for(typo in typos){
  newrownames=gsub(typo,paste(typo,"T1",sep=""), newrownames)
}
rownames(meta_tumours) = newrownames

# read in the gene meta data, i.e., a list of PAM50 genes
# pam50 = readListFromExcel(paste(indir,
#                                 "PAM50_and_Polyak_Cell_Snapshot-gene_lists.xlsx",
#                                 sep="/")) 
# head(pam50, 5)

meta_genes = list()
for(i in seq(1,3)){
  meta_genes = c(meta_genes,readListFromExcel(paste(indir,"Genelists-summary160606.xlsx",sep="/"), sheet=i))
}
# Add Henrik's last list
print(paste(indir,"KEGG_n_Hallmark_genes_for_mRNA-protein_corr.xlsx",sep="/"))
meta_genes = c(meta_genes,readListFromExcel(paste(indir,"KEGG_n_Hallmark_genes_for_mRNA-protein_corr.xlsx",sep="/"), sheet=1))

# Add Henrik's very last list
print(paste(indir,"COSMIC_n_BC_drivers.xlsx",sep="/"))
meta_genes = c(meta_genes,readListFromExcel(paste(indir,"COSMIC_n_BC_drivers.xlsx",sep="/"), sheet=1))
print(length(meta_genes))

# meta_genes["PAM50"] = pam50
newnames = names(meta_genes)
newnames=gsub(" ","_", newnames)
names(meta_genes) = newnames
names(meta_genes)

## Not included in this release
## # CNV data
## ##########
## # The CNV files were updated to GRCh37 together with the mRNA data above, so 
## # read in the CNV  data from outdir
## filename = paste(outdir,
##                  "Gene_probe_centric.LogR_tperc_ploidy_adj.OSL2_n331_updated.txt.gz",
##                  sep="/")
## # filename = paste(indir,
## #                  "Gene_probe_centric.LogR_tperc_ploidy_adj.OSL2_n331.txt.gz", 
## #                  sep="/")
## cnv=read.table(gzfile(filename), header=TRUE, sep="\t")
## # This data contains dupliate lines -- remove them!
## cnv=cnv[!duplicated(cnv,MARGIN=1),]
## #cnv$chr[cnv$chr==23] = "X"
## rownames(cnv) = cnv$probes #paste(cnv$chr,cnv$start,cnv$stop,sep=":")
## #This data contains "-Inf" values (=314) which are removed
## cnv = cnv[-which(apply(cnv=="-Inf",1,sum)>=1),]

# OUTPUT
########
# Bundle protein and rna data arrays together with some info in a list
exprdata=list(name="Tumour expressionData",
              readme=paste(
                "This list comprise proteome expression data and mRNA",
                "expression data for 45 breast cancer tumours as 2 separate",
                "arrays indexed by tumours (columns) and gene symbol (rows),",
                "as well as GW cnv data for the same tumours. It is created",
                "by the script formatData.R from the original data files",
                "diffProtdata_extracted_from_9995_BC_45_proteome_Gene_symbol_centric_1FDR.xlsx",
                "(protein data from Henrik Johansson, Lehtio group),",
                "main_qnorm_missimp_cent-hosp_log2_annot_updated_ids_to_Janne.txt",
                "(mRNA data, from Kristine Kleivi, Oslo group),",
                "PAM50_and_Polyak_Cell_Snapshot-gene_lists.xlsx and",
                "Genelists-summary.xlsx (Meta data from Henrik Johansson), and",
                "140219_OSLO2.n325.CNA.byGeneExpressionProbePositons.HKMV_updated.Sample.Ids_MRA_SJ_41samples.txt"
#                "(CNV data, from Miriam Ragle Aure)"
                ), 
              protexpr=prot_data, rnaexpr=rna_data, metatumour=meta_tumours,
              metagene=meta_genes) #, cnv=cnv) 
# Save the list for later access
save(exprdata,file=paste(outdir,"tumourExpressionData",sep="/"))

