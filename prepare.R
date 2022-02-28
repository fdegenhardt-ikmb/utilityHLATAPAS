#########################################################################
# ToCSV and toDOSAGE of output from Michigan Imputation Server HLA imputation
#########################################################################

# February 2022
# Author: Frauke Degenhardt
# Contact: f.degenhardt@ikmb.uni-kiel.de

###############################
# SETTINGS
###############################

library(data.table)
library(reshape)

options(stringsAsFactors=F)
###############################
# FUNCTIONS
###############################

###############################
# MAIN
###############################
# READ VCF FILE
data  = data.frame(fread("candidates.dose.vcf",sep="\t"))
samples = fread("samples.txt", h=F,sep="\t")
colnames(data) = samples

# STORE FIRST 9 COLUMNS OF data
info = data[,1:9]
info$ID=unlist(lapply(strsplit(info$ID, ":"), function(x){paste0(x[1],":",x[2])}))

# GET GENOTYPES
geno = data[,-(1:9)]
geno = apply(geno,1,function(x){gsub(":[0-9].*", "", x)})
geno = apply(cbind(info$ID, t(geno)),1, 
	         function(x){x=gsub("1",x[1],x); return(x[-1])}) # REPLACE 1 with ALLELE NAME

# CREATE PHASED OUTPUT
tmp = apply(geno, 1, function(x){x=x[x!="0|0"];  
			x = do.call(rbind,strsplit(x,"\\|")); 
			x= c(paste(x[,1], collapse=","), paste(x[,2], collapse=",")); 
			x = gsub(",0|0,", "", x); return(x)})
tmp = t(tmp)
tmp = melt(tmp)

# FORMAT
tmp$value = gsub("HLA_", "", gsub(",","-", tmp$value))
tmp = tmp[tmp$value!=0 & tmp$value!="",]
colnames(tmp)=c("ID","", "haplotype")

n = nchar(gsub("[0-9A-Z\\*:]", "", tmp$haplotype)) + 1

# MAKE DOSAGE
dosage_a = as.data.frame.matrix(table(tmp$haplotype, tmp$ID))

dosage_b = apply(data[,-(1:9)], 1, function(x){x=do.call(rbind,strsplit(as.matrix(x),":")); return(as.numeric(x[,2]))})
dosage_b = t(dosage_b)
rownames(dosage_b)=data$ID
colnames(dosage_b)=colnames(data)[-(1:9)]

dosage = rbind(dosage_a, dosage_b)

# MAKE CSV

id = rownames(dosage)
CSV = apply(dosage, 2, function(x){x=cbind(id[x!=0], x[x!=0])})
CSV = cbind(rep(names(CSV), times = unlist(lapply(CSV, nrow))), do.call(rbind,CSV))
colnames(CSV)=c("IID", "id","dosage")

CSV = data.frame(CSV)
# EXCLUDE ALLELES WITH LOW PROBABILITY
CSV= CSV[CSV$dosage > 0.8 & !(CSV$dosage > 1 & CSV$dosage < 1.8) ,]

write.table(CSV, "final_michigan_imputation_server.csv",quote=F, row.names=F, sep="\t")
save(dosage,file = "final_michigan_imputation_server.RData")
