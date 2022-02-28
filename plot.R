#########################################################################
# Visualization of output from Michigan Imputation Server HLA imputation
#########################################################################

# February 2022
# Author: Frauke Degenhardt
# Contact: f.degenhardt@ikmb.uni-kiel.de

###############################
# SETTINGS
###############################
library(ggplot2)
library(reshape)
###############################
# FUNCTIONS
###############################

# Calculate distribution of imputation R^2 values
calc = function(tmp){
out = do.call(rbind, 
			tapply(tmp$Rsq, paste(tmp$digit, tmp$locus), 
				function(x){x = cut(x, seq(0,1,0.01)); 
					    x = table(x); 
					    x = cumsum(x)}))
names_old= seq(0,1,0.01)
names_new = colnames(data.frame(out))
info = do.call(rbind, strsplit(rownames(out), " "))
colnames(info) = c("digit", "locus")
out = cbind(info, out)
out = melt(data.frame(out), id=c("digit", "locus"))
out$variable = names_old[match(out$variable, names_new)]
return(out)}


###############################
# MAIN
###############################
# Read in file with imputation R^2 values
tmp = read.table("info.txt", h=T)
print(head(tmp))
tmp$digit=nchar(gsub("[0-9A-Z_\\*]", "", tmp$SNP))+1
tmp$locus = gsub("\\*.*", "", tmp$SNP)


# Plot discribution of R^2values
pdf("qualityInfo.pdf", width=10)
out = calc(tmp)

# whole disgribution
p = ggplot(out, aes(x=variable, y=as.numeric(as.matrix(value)))) + 
    geom_line(aes(col=locus)) + 
    facet_wrap(~digit, scale="free") +  xlab("Rsq") + ylab("cumsum (%)") + ggtitle("all")
print(p)
# distribution only of alleles with minor allele frequency (MAF) > 1%
out = calc(tmp[tmp$ALT_Frq > 0.01,])
p = ggplot(out, aes(x=variable, y=as.numeric(as.matrix(value)))) + 
    geom_line(aes(col=locus)) + 
    facet_wrap(~digit, scale="free") + xlab("Rsq") + ylab("cumsum (%)") + ggtitle("MAF > 1%")
print(p)
dev.off()
