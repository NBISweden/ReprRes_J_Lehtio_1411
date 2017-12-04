#Updating gene symbols in CNA/mRNA data files

# Loading files

indir = "DATA/"
outdir = "R-files/"
new_ids = read.table(paste(indir,"Agilent60k.Hg19.annotation.txt", sep =""), header = TRUE)

# mRNA data 

infilename = "main_qnorm_missimp_cent-hosp_log2_annot_updated_ids_to_Janne.txt"
outfilename = "main_qnorm_missimp_cent-hosp_log2_annot_updated_ids_to_Janne_updated.txt"

rna = read.table(paste(indir,infilename, sep="/"), header=TRUE, sep="\t")
idx = match(rna$ProbeName,new_ids$probes)
dt1 = new_ids[idx,]
rna$GeneName = as.character(rna$GeneName)
rna$GeneName[!is.na(dt1$gene)] = as.character(dt1$gene[!is.na(dt1$gene)])
#write.table(rna$GeneName, "rna_gene.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(rna, gzfile(paste(outdir,outfilename,".gz",sep="")), sep = "\t")

## # CNV data
## infilename ="Gene_probe_centric.LogR_tperc_ploidy_adj.OSL2_n331.txt.gz"
## outfilename = "Gene_probe_centric.LogR_tperc_ploidy_adj.OSL2_n331_updated.txt"
## cnv=read.table(paste(indir,infilename,sep="/"), header=TRUE, sep="\t") 

## # Matching

## idx = match(cnv$probes,new_ids$probes)
## dt = new_ids[idx,]
## cnv$gene = as.character(cnv$gene)

## #Excluding NA values

## cnv$gene[!is.na(dt$gene)] = as.character(dt$gene[!is.na(dt$gene)])
## #write.table(cnv$gene, "cnv_gene.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
## write.table(cnv, gzfile(paste(outdir, outfilename,".gz",sep="")), sep = "\t")

