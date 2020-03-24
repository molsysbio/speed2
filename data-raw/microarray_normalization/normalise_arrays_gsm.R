
library(GEOquery);
library(annaffy);
library(preprocessCore);

geo_set<-read.table("set.dat");

gsm_accessions<-read.table("gsm_accession.dat");

gsmtmp={};
ids={};
values={};
for (i in 1:length(gsm_accessions$V1)) {
  gsmtmp=getGEO(as.character(gsm_accessions[i,1]));
  ids=cbind(ids,as.character(Table(gsmtmp)$ID_REF));
  values=cbind(values,as.numeric(Table(gsmtmp)$VALUE));
}


stopifnot(dim(ids)[2]==length(gsm_accessions$V1),
          dim(values)[2]==length(gsm_accessions$V1),
          dim(ids)[1]==dim(values)[1],
          dim(ids)[1]>0);

platform=Meta(gsmtmp)$platform_id;
gpl=Table(getGEO(platform));
if (sum(labels(gpl)[[2]]=="ORF_LIST")>0) {
  geneIDcolumn <- labels(gpl)[[2]]=="ORF_LIST";
} else {
  geneIDcolumn <- (labels(gpl)[[2]]=="GeneID") |(labels(gpl)[[2]]=="Gene.ID") | (labels(gpl)[[2]]=="ENTREZ_GENE_ID") |  (labels(gpl)[[2]]=="Entrez_Gene_ID");
}
# assert:
stopifnot(sum(geneIDcolumn)>0);

geneIDs <- gpl[geneIDcolumn]
probesets=as.character(gpl$ID)

values2={};
for (i in 1:length(gsm_accessions$V1)) {
  mymatch <- match(probesets,ids[,i]);
  values2=cbind(values2,values[mymatch,i]);
}

if (geo_set[2,1]=="log") {
  expressiondata=log2(values2);
} else {
  expressiondata=values2;
}

expressiondatanormalised <- normalize.quantiles(expressiondata,copy=TRUE)

if (sum(gsm_accessions$V2==1)==1) {
  Means=rowMeans(cbind(expressiondatanormalised[,gsm_accessions$V2==1],expressiondatanormalised[,gsm_accessions$V2==2]));
  Stddev=apply(cbind(expressiondatanormalised[,gsm_accessions$V2==1],expressiondatanormalised[,gsm_accessions$V2==2]),1,sd);
} else {
  Means=rowMeans(expressiondatanormalised[,gsm_accessions$V2==1]);
  Stddev=apply(expressiondatanormalised[,gsm_accessions$V2==1],1,sd);
}

# Keeps track of rows with problems...

stddevmodel = loess(Stddev ~ Means, data.frame(Means,Stddev)[!(is.na(Means)|is.na(Stddev)),])

if (sum(gsm_accessions$V2==1)==1) {
  foldchange=expressiondatanormalised[,gsm_accessions$V2==2]-expressiondatanormalised[,gsm_accessions$V2==1];
  Means2=Means;
} else  {
  if (sum(gsm_accessions$V2==2)==1) { 
    foldchange=expressiondatanormalised[,gsm_accessions$V2==2]-rowMeans(expressiondatanormalised[,gsm_accessions$V2==1]);
    Means2=rowMeans(cbind(expressiondatanormalised[,gsm_accessions$V2==2],expressiondatanormalised[,gsm_accessions$V2==1]));
  } else {
    foldchange=rowMeans(expressiondatanormalised[,gsm_accessions$V2==2])-rowMeans(expressiondatanormalised[,gsm_accessions$V2==1]);
    Means2=rowMeans(cbind(expressiondatanormalised[,gsm_accessions$V2==2],expressiondatanormalised[,gsm_accessions$V2==1]));
  }
}
notnas=!(is.na(Means2)|is.na(foldchange))

zvalues=foldchange[notnas]/predict(stddevmodel,Means2[notnas]);
zvaluessorted=sort(zvalues,index.return=TRUE,decreasing = TRUE);
goodgeneids=geneIDs[notnas,1];
goodmeans2=Means2[notnas];

zvaluestable=data.frame(goodgeneids[zvaluessorted$ix],zvaluessorted$x,goodmeans2[zvaluessorted$ix]);
colnames(zvaluestable)=c("GeneID","zvalue")
write.table(zvaluestable, file = "zvalues.dat",  quote = FALSE, sep = "\t",
                 eol = "\n", na = "NA", dec = ".", row.names = TRUE,
                 col.names = TRUE, qmethod = c("escape", "double"))
