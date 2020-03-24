minrowcol=matrix(pmin(rep(diag(A),length(diag(A))),rep(diag(A),each=length(diag(A)))),length(diag(A)),length(diag(A)));
plot(hclust(as.dist(1-A/minrowcol)),cex=3)
