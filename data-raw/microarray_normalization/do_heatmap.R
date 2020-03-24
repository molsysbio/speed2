x=read.table('output.RDATA')
xreduced=x[rowSums(is.na(x))!=dim(x)[2],colSums(is.na(x))!=dim(x)[1]];
xbinary=matrix(0,nrow=dim(xreduced)[1],ncol=dim(xreduced)[2]);
for (i in 1:dim(xreduced)[2]) {
  goodones=abs(xreduced[,i])>quantile(abs(xreduced[,i]),probs=c(0.90),na.rm=TRUE);
  xbinary[goodones&!is.na(goodones),i]=sign(xreduced[goodones&!is.na(xreduced[,i]),i]);
  
}

manhattan=matrix(0,nrow=dim(xreduced)[2],ncol=dim(xreduced)[2])
for (i in 1:dim(xreduced)[2]) {
  print(i);
  whichtocompare=!(is.na(xreduced)|is.na(xreduced)[,i]);
  for (j in i:dim(xreduced)[2]) {
    manhattan[i,j]=sum(abs(abs(xbinary[whichtocompare[,j],i])-abs(xbinary[whichtocompare[,j],j])))/sum(whichtocompare[,j]);
    manhattan[j,i]=manhattan[i,j];
  }
}
plot(hclust(as.dist(manhattan),method="ward"),labels=labels(xreduced)[[2]])
