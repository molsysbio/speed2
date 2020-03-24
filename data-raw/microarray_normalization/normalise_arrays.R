# Todo:
# Calculate variance if only one replecate
# More than 1 sample

library(GEOquery);
library(annaffy);
library(preprocessCore);

geo_set<-read.table("set.dat");

##############################
# Download expression set
tryCatch(gds <- getGEO(as.character(geo_set[1,1]),GSEMatrix = TRUE), error= gds <- getGEO(as.character(geo_set[1,1]),GSEMatrix = FALSE))
if (length(grep('GDS.*',geo_set$V1))>0) {
  eset <- GDS2eSet(gds,do.log2 = geo_set[2,1]=="log");
} else {
  eset <- gds[[1]];
}
##############################
# Find out if and where Gene ID column is
geneIDcolumn <- (varLabels(featureData(eset))=="Gene.ID") | (varLabels(featureData(eset))=="ENTREZ_GENE_ID")

# Sanety Check: Are the Gene IDs available?

if (sum(geneIDcolumn)!=1) {
  print("######### THERE ARE NO GENE ID ANNOTATIONS!!!!");
}

# Extract Gene IDs from annotations
tmp <- pData(featureData(eset));
geneIDs <- tmp[,geneIDcolumn]

##############################
# Extract data matrix and do quantile normalisation
# Bolstad, B. M., Irizarry R. A., Astrand, M, and Speed, T. P. (2003)
# A Comparison of Normalization Methods for High Density Oligonucleotide Array Data Based on Bias and Variance.
# Bioinformatics 19(2) ,pp 185-193
expressiondata <- exprs(eset)
if (length(grep('GDS.*',geo_set$V1)) == 0) {
  if (geo_set[2,1]=="log") {
    expressiondata=log2(expressiondata);
  }
}
expressiondatanormalised <- normalize.quantiles(expressiondata,copy=TRUE)


# Sanety Check: Are the ids of annotation in the same order as in expression set?

if (sum(labels(tmp)[[1]]!=labels(expressiondata)[[1]])!=0) {
  print("######### Annotations seem to be in different order");
}

##############################
# Now we extract the indecies for the correct arrays

gsm_accessions<-read.table("gsm_accession.dat");

arrayindices=c();
idx=1:length(labels(expressiondata)[[2]]);

for(i in 1:dim(gsm_accessions)[1]) {
  arrayindices=c(arrayindices,idx[labels(expressiondata)[[2]]==as.character(gsm_accessions[i,1])]);
}

##############################
# Now average over same samples
# Calculate Variance
# Calculate and save z-calues

if (length(arrayindices[gsm_accessions$V2==1])==1) {
  Means=rowMeans(cbind(expressiondatanormalised[,arrayindices[gsm_accessions$V2==1]],expressiondatanormalised[,arrayindices[gsm_accessions$V2==2]]));
  Stddev=apply(cbind(expressiondatanormalised[,arrayindices[gsm_accessions$V2==1]],expressiondatanormalised[,arrayindices[gsm_accessions$V2==2]]),1,sd);
} else {
  Means=rowMeans(expressiondatanormalised[,arrayindices[gsm_accessions$V2==1]]);
  Stddev=apply(expressiondatanormalised[,arrayindices[gsm_accessions$V2==1]],1,sd);
} 
stddevmodel = loess(Stddev ~ Means, data.frame(Means,Stddev)[1:dim(as.matrix(Means))[1],])
stddevthroughmodel <- predict(stddevmodel,Means)
sorted=sort(Means,index.return=TRUE);


if (length(arrayindices[gsm_accessions$V2==1])==1) {
  foldchange=expressiondatanormalised[,arrayindices[gsm_accessions$V2==2]]-expressiondatanormalised[,arrayindices[gsm_accessions$V2==1]];
  Means2=Means;
} else  {
   foldchange=rowMeans(expressiondatanormalised[,arrayindices[gsm_accessions$V2==2]])-rowMeans(expressiondatanormalised[,arrayindices[gsm_accessions$V2==1]]);
  Means2=rowMeans(cbind(expressiondatanormalised[,arrayindices[gsm_accessions$V2==2]],expressiondatanormalised[,arrayindices[gsm_accessions$V2==1]]));
}
zvalues=foldchange/predict(stddevmodel,Means2);
zvaluessorted=sort(zvalues,index.return=TRUE,decreasing = TRUE)

zvaluestable=data.frame(geneIDs[zvaluessorted$ix],zvaluessorted$x);
colnames(zvaluestable)=c("GeneID","zvalue")
write.table(zvaluestable, file = "zvalues.dat",  quote = FALSE, sep = "\t",
                 eol = "\n", na = "NA", dec = ".", row.names = TRUE,
                 col.names = TRUE, qmethod = c("escape", "double"))
pdf("MAplot.pdf")
plot(Means2,Stddev,ylim=c(0,max(foldchange)));
lines(sorted$x,stddevthroughmodel[sorted$ix],col="red");
points(Means2[zvaluessorted$ix[abs(zvaluessorted$x)>2]],abs(foldchange[zvaluessorted$ix[abs(zvaluessorted$x)>2]]),col="blue");
points(Means2[zvaluessorted$ix[abs(zvaluessorted$x)>3]],abs(foldchange[zvaluessorted$ix[abs(zvaluessorted$x)>3]]),col="green");
dev.off()
