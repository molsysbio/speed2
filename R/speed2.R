options(stringsAsFactors=FALSE)

Ua=-1
Ub=1
BAD_PATHWAYS = c("PI3K", "RAR", "Trail", "MAPK")  # pathways with insufficient data


####################################
######## INTERNAL FUNCTIONS ########
####################################

# fast conversion to Entrez id (also works for background list)
convertSymbolToEntrez = function(genes) {

  genes = as.character(genes, stringsAsFactors = FALSE)
  df_entrezId_joinEntrez = suppressMessages(left_join(data.frame(genes_asEntrez=genes, stringsAsFactors = FALSE), df_entrezID %>% mutate(g_id=as.character(g_id)), by=c("genes_asEntrez"="g_id")))
  df_entrezId_joinSymbol = suppressMessages(left_join(data.frame(genes_asSymbol=genes, stringsAsFactors = FALSE), df_entrezID %>% mutate(g_id=as.character(g_id)), by=c("genes_asSymbol"="gene")))

  ret = rep(NA,length(genes))
  ret[which(!is.na(df_entrezId_joinEntrez$gene))] = df_entrezId_joinEntrez$genes_asEntrez[which(!is.na(df_entrezId_joinEntrez$gene))]
  ret[which(!is.na(df_entrezId_joinSymbol$g_id))] = df_entrezId_joinSymbol$g_id[which(!is.na(df_entrezId_joinSymbol$g_id))]

  return(ret)

}

convertEntrezToSymbol = function(genes) {

  df_entrezId_joinEntrez = suppressMessages(left_join(data.frame(genes_asEntrez=genes, stringsAsFactors = FALSE), df_entrezID %>% mutate(g_id=as.character(g_id)), by=c("genes_asEntrez"="g_id")))
  df_entrezId_joinSymbol = suppressMessages(left_join(data.frame(genes_asSymbol=genes, stringsAsFactors = FALSE), df_entrezID %>% mutate(g_id=as.character(g_id)), by=c("genes_asSymbol"="gene")))

  ret = rep(NA,length(genes))
  ret[which(!is.na(df_entrezId_joinEntrez$gene))] = df_entrezId_joinEntrez$gene[which(!is.na(df_entrezId_joinEntrez$gene))]
  ret[which(!is.na(df_entrezId_joinSymbol$g_id))] = df_entrezId_joinSymbol$genes_asSymbol[which(!is.na(df_entrezId_joinSymbol$g_id))]

  return(ret)

}




ranks_Pbates_pipe = function(df) {
  df %>%
    dplyr::select(.data$g_id, .data$zrank_signed_mean, .data$P.bates) %>%
    mutate(rank = ifelse(.data$zrank_signed_mean >= (Ua+Ub)/2, .data$P.bates+1e-50, -.data$P.bates-1e-50)) %>%
    mutate(rank = ifelse(rank<=0, -1-rank, 1-rank)) %>%
    dplyr::select(.data$g_id, rank) %>%
    deframe
}

# Bates CDF
batesCDF = function(x, n, a=Ua, b=Ub) {
  if(n<1 | n>nrow(batesCDFtable))
    return(NA)

  if(x<=a)
    return(0)
  else if(x>=b)
    return(1)

  xs = as.double(colnames(batesCDFtable))
  xi = which.min(abs(xs-x))

  return(as.double(batesCDFtable[n,xi]))
  #pnorm(x, 0.5, 1/sqrt(12*n))
}

bates2sidedPval = function(x,n,a=Ua,b=Ub) {
  if(x<=(a+b)/2)
    return(batesCDF(x, n, a=a, b=b))
  if(x>(a+b)/2)
    return(batesCDF(a+b-x ,n, a=a, b=b))
}

bates2sidedPval_vec = function(xs,ns, a=Ua, b=Ub) {
  return(sapply(1:length(xs), function(i)  bates2sidedPval(xs[i], ns[i], a=a, b=b )))
}


# Test of variance of uniformly distributed [Ua,Ub] random variabels
# Null model: (approx) normal distribution with mean=(b-a)^2/12 and sd = (b-a)/sqrt(12*N).
# (X-0.5)^2: u=1/12, sd^2=(b-a)^4/180 (do the integral)
# One-sided: P<0.05 should be used (don't care about unusually "centered" gene sets)
chi2test_var = function(vals, a=Ua, b=Ub, testM=(b-a)^2/12, testVar=(b-a)^4/(180*length(vals))) {
  m=(a+b)/2
  N = length(vals)
  var_vals = sum((vals-m)^2)/N
  if(var_vals>testM) {
    p = 1-pnorm(var_vals, mean=testM, sd=sqrt(testVar))
  } else {
    p = 1
  }
  return(p)
}


####################################
######## EXTERNAL FUNCTIONS ########
####################################

#' Test pathway gene signature enrichment using SPEED2
#' @importFrom magrittr %>%
#' @import dplyr
#' @import tibble
#'
#' @param genes list of genes in Entrez or gene symbols format to test
#'    for enrichment in SPEED2. Should be at least 10 genes, maximally 500 genes.
#' @param background_genes (optional) list of background genes in Entrez or gene
#'    symbol from which \code{genes} were selected. If not provided, the full
#'    set of SPEED2 genes are used.
#' @param shuffle (optional) shuffle identities of genes in SPEED2, for control
#'     experiments.
#' @param custom_signatures (optional, advanced) user provided custom pathway
#'     gene signatures to use instead of SPEED2, provided as a tibble with the
#'     following columns: p_id (pathway id), g_id (gene id), zrank_signed_mean
#'     (average normalized z score across many experiments, between -1 and +1),
#'     P.bates (p-value associated with zrank_signed_mean). Users are expected
#'     to handle Entrez conversion and background gene list themselves.
#' @return List with four items: \code{df_stattest} a tibble with enrichment
#'     scores, \code{df_rankcoords} coordinates for GSEA plot (see publication),
#'     and two lists with unmatched items in \code{genes} and
#'     \code{background_genes}.
#' @examples
#' \code{ret = speed2run(genes=speed2:::speed2_signatures$g_id[1:50])}
#' @export
speed2run = function(genes, background_genes=c(), shuffle=F, custom_signatures=NULL) {
  genes=unlist(unname(genes))
  unknown_items_backgrund = character(0)
  unknown_items_geneset = character(0)

  if(all(c("p_id", "g_id", "zrank_signed_mean", "P.bates" ) %in% colnames(custom_signatures))) {
    df_signatures = custom_signatures # user provided signatures
    pathways = unique(df_signatures$p_id)
  } else {

    df_signatures = speed2_signatures %>% 
      left_join(df_entrezID,by = "g_id") %>% 
      left_join(speed2_pathways, by = "pathway_number") %>% 
      mutate(g_id=as.character(g_id)) %>% 
      mutate(regulation=if_else(regulation,"UP","DOWN"))
    pathways = unique(df_signatures$p_id)

    are_genes_in_speed = genes %in% as.character(df_entrezID$g_id) |
      genes %in% df_entrezID$gene
    unknown_items_geneset = genes[which(are_genes_in_speed == FALSE)]

    genes = convertSymbolToEntrez(genes)  # (entrez, symbol) -> entrez

    if(length(background_genes)>0) {
      are_genes_in_speed = background_genes %in% as.character(df_entrezID$g_id) |
        background_genes %in% df_entrezID$gene
      unknown_items_backgrund = background_genes[which(are_genes_in_speed == FALSE)]
      background_genes = convertSymbolToEntrez(background_genes)  # (entrez, symbol) -> entrez

      genes = intersect(genes, background_genes)
      df_signatures = df_signatures %>% filter(.data$g_id %in% background_genes)
      if(length(genes)<10 | (nrow(df_signatures)<1000))
        return("number of genes to low (<10) or number of matched background genes to low (<1000) ")
    }

  }

  df_stattest = tibble()
  df_rankcoords = tibble()
  for(pathway_ix in 1:length(pathways)) {
    pathway = pathways[pathway_ix]
    ranks = df_signatures %>% filter(.data$p_id==pathway) %>% ranks_Pbates_pipe
    ranks = ranks[order(ranks)]
    ranks_names = names(ranks)
    ranks = 1:length(ranks)
    names(ranks) = ranks_names
    if(shuffle)
      names(ranks) = sample(ranks_names)
    ranks = scales::rescale(ranks, to=c(Ua,Ub))
    matchRanksGeneset = match(genes, names(ranks))
    geneset_ranks = ranks[matchRanksGeneset[!is.na(matchRanksGeneset)]]
    if(length(geneset_ranks)<1)
      next
    mean_rank = mean(geneset_ranks)
    chi2TestStatistic = sum((geneset_ranks-(Ua+Ub)/2)^2)/length(geneset_ranks)

    #ztest = ztest_mean(geneset_ranks, a=Ua, b=Ub)
    if(length(geneset_ranks)>nrow(batesCDFtable)) {
      print("Warning: too many input genes, truncating.")
      geneset_ranks = geneset_ranks[2:nrow(batesCDFtable)]
      mean_rank = mean(geneset_ranks)
    }
    p_ztest = bates2sidedPval(x=mean(geneset_ranks), n=length(geneset_ranks), a=Ua, b=Ub )
    ztest_reg = ifelse(mean(geneset_ranks)<(Ua+Ub)/2,"DOWN", "UP"  )
    p_chi2test = chi2test_var(geneset_ranks, a=Ua, b=Ub)

    df_stattest = df_stattest %>%
      #rbind(tibble(speed2=pathway, p_ztest=ztest$p, ztest_reg=ztest$reg, p_chi2test=p_chi2test, mean_rank=mean_rank, chi2TestStatistic=chi2TestStatistic ))
      rbind(tibble(speed2=pathway, p_ztest=p_ztest, ztest_reg=ztest_reg, p_chi2test=p_chi2test, mean_rank=mean_rank, chi2TestStatistic=chi2TestStatistic ))

    df_rankcoords = df_rankcoords %>% rbind(tibble(speed2=rep(pathway, length(geneset_ranks)), rank=geneset_ranks  ))

  }
  df_stattest = df_stattest %>%
    arrange(p_ztest, abs(mean_rank-(Ua+Ub)/2))


  return(list(df_stattest=df_stattest, df_rankcoords=df_rankcoords,
              unknown_items_geneset=unknown_items_geneset, unknown_items_backgrund=unknown_items_backgrund))
}

#' SPEED2 enrichment visualization
#' @import ggplot2
#' @import cowplot
#'
#' @param df_stattest output from \code{speed2run}.
#' @param df_rankcoords output from \code{speed2run}.
#' @param test statistical test to use in plot, either 'z' or 'chi2'.
#' @param plot_title (optional) plot title.
#' @param save_plot (optional) boolean.
#' @param plot_filename (optional) file name of plot (if saved).
#' @param version (optional) "fullsize" for big plot, otherwise small plot.
#' @return pathway gene signature enrichment plot similar to publication.
#' @examples
#' \code{ret = speed2run(genes=speed2:::speed2_signatures$g_id[1:50])}
#' \code{speed2plot(ret$df_stattest, ret$df_rankcoords, test="chi2")}
#' @export
speed2plot = function(df_stattest, df_rankcoords, test="z", plot_title="", save_plot=F, plot_filename="plot.png", version="fullsize") {
  df = df_stattest
  df_coords = df_rankcoords
  if(test=="z") {
    df$p = df$p_ztest
  } else if(test=="chi2") {
    df$p = df$p_chi2test
  } else { stop('the parameter test should be either z or chi2'); } 
  minPpos = min(df$p[df$p>0])^2
  df$p[df$p==0] = minPpos
  df$p = p.adjust(df$p, method="BH")

  if(version=="fullsize") {
    fontsize_base = 12
    fontsize_geompoint = 5
    fontsize_geomtext = 4
    barheight = 0.9
  } else {
    fontsize_base = 8
    fontsize_geompoint = 2
    fontsize_geomtext = 2.5
    barheight = 0.7
  }


  df$alpha = ifelse(df$p>0.05, "Opaque", "Solid")
  if(test=="z") {
    df$y = df$mean_rank
    df$speed2 = ordered(df$speed2, levels=df$speed2[order(abs(df$mean_rank-(Ua+Ub)/2) )])
    df_coords$speed2 = ordered(df_coords$speed2,levels=levels(df$speed2))
    #df$y  = 2*(df$y -0.5)
    label_offset = 0.2 # 0.13
  } else if(test=="chi2") {
    df$y = df$chi2TestStatistic
    df$speed2 = ordered(df$speed2, levels=df$speed2[order(abs(df$chi2TestStatistic-(Ub-Ua)^2/12))])
    df_coords$speed2 = ordered(df_coords$speed2,levels=levels(df$speed2))
    df$y = df$y-(Ub-Ua)^2/12 # 1/12 is expected variance for uniform distribution
    label_offset = 0.2#0.1*(max(df$chi2TestStatistic)-min(df$chi2TestStatistic))
  }
  df$label = format(round(df$y, 2), nsmall = 2)

  redgradient_fun = colorRampPalette(c("#FFCCCC", "red"))
  pvalLimit=3
  # ggplot(df, aes(x=speed2, y=-log10(p), fill=reg, alpha=alpha, label=signif(mean_rank, digits=2) )) +
  # label=signif(mean_rank, digits=2) )
  gg = ggplot(df, aes(x=.data$speed2, y=.data$y, fill=-log10(.data$p), alpha=alpha, label=.data$label)) +
    geom_col()+
    geom_text(aes(y = .data$y+sign(.data$y)*label_offset), size=fontsize_geomtext)+ # , nudge_y=label_offset
    coord_flip()+
    theme_minimal()+
    theme(axis.text = element_text(size=fontsize_base),
          axis.title = element_text(size=fontsize_base+2),
          legend.title = element_text(size=fontsize_base),
          legend.text=element_text(size=fontsize_base),
          legend.position="bottom",
          #legend.key.height=unit(3, "cm"),
          axis.ticks.length = unit(.15, "cm"),
          plot.title = element_text(size = fontsize_base+2))+
    #xlab("SPEED 2 ranked lists")+
    scale_alpha_manual(values=c("Opaque" = 0.4, "Solid"=1))+
    xlab("")+
    #ylab(expression(-log[10](Padj)))+
    scale_fill_gradientn(expression(-log[10](P-adj)),
                         limit=c(1, 10), # -log10(min(df$p))
                         #values=rescale(c(0, 1.3-0.01, 1.3, 4 )),
                         breaks=c(1, 1.3, 5, 10),
                         colors=c("gray",redgradient_fun(3)),
                         label=c("", "n.s.", "5", "10" ), # as.character(floor(-log10(min(df$p))))
                         oob = scales::squish)+
    guides(alpha=FALSE, fill=guide_colourbar(barheight=0.7))+
    #scale_fill_manual(drop=FALSE, values=c("UP" = "#F8766D", "DOWN" = "#3498db", "UP+DOWN" = "#2ecc71"))+
    ggtitle(plot_title)

  if(test=="z") {
    gg = gg +
      scale_y_continuous(breaks=c(-1,-0.5,0,0.5,1), limits=c(-1.23,1.23))+
      ylab("Downregulated                   Upregulated\n\nMean rank")
  } else if(test=="chi2") {
    gg = gg +
      scale_y_continuous(breaks=c(-1,-0.5,0,0.5,1), limits=c(-1.23,1.23))+
      ylab("Variance compared to uniform rank distribution")
  }
  gg_coord = ggplot(df_coords, aes(y=.data$speed2, x=rank)) +
    geom_point(alpha=0.15, shape="|", size=fontsize_geompoint)+
    theme_minimal()+
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.title = element_blank(),
          #axis.text.x = element_blank(),
          axis.text = element_text(size=fontsize_base),
          #axis.title = element_text(size=12)
          axis.text.y = element_blank())+
    scale_x_continuous(breaks=c(-1,0,1), limits=c(-1,1))+
    xlab("Rank")

  gg_combined = plot_grid(gg, gg_coord, nrow=1, align="h", axis="tb", rel_widths = c(0.75,0.25))


  if(save_plot)
    save_plot(filename=plot_filename, gg_combined, base_width = 8, base_height = 8)

  return(gg_combined)


  # Base64-encode file
  # plot_file <- tempfile(fileext = ".png")
  # ggsave(filename=plot_file, gg, width = 20, height = 12, units = "cm")
  # txt <- base64Encode(readBin(plot_file, "raw", file.info(plot_file)[1, "size"]), "txt")
  # unlink(plot_file)
  # txt[1]
}
