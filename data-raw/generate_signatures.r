options(stringsAsFactors=FALSE)
library(tidyverse)
library(rjson)
library(speed2)

# INSTRUCTIONS:
# 1) Download the GEO microarray data listed in speed2_experiments.tsv
# 2) Normalize microarray data using the included normalization scripts.
# 2) Run parse_zvalues.R
# 3) Run generate_signatures.R (this file). Q-values need to be computed seperately
#    on a cluster resource.

#############################################
########### HELPER FUNCTIONS ################
#############################################

MINCOUNT = 2
Ua = speed2:::Ua
Ub = speed2:::Ub

LoadToEnvironment <- function(RData, env = new.env()){
  load(RData, env)
  return(env)
}

signZvalue = function(zscore, effect) {
  if(effect %in% c("activating", "activation")) {
    return(zscore)
  } else if(effect %in% c("inhibiting", "inhibition", "inhibting")) {
    return(-zscore)
  }
  return(NA)
}

signZvalue_vec = function(zscore, effect) {
  return(sapply(1:length(zscore), function(ix)  signZvalue(zscore[ix], effect[ix])))
}


regFromZscores = function(zscores) {
  if(length(which(zscores>0)) > length(which(zscores<0)))
    return("UP")
  else
    return("DOWN")
}


bates_pipe = . %>%
  group_by(p_id, gene, g_id) %>%
  summarise(zrank_signed_mean=mean(zrank_signed),
            regulation=regFromZscores(zvalue_signed),
            nexp=n()) %>%
  filter(nexp>=MINCOUNT) %>%
  mutate(P.bates=speed2:::bates2sidedPval_vec(zrank_signed_mean,nexp, a=Ua, b=Ub)) %>%
  arrange(p_id, P.bates)


formatDataFrameForWeb = function(df) {
  df_web = df
  df_web$P.bates = signif(df_web$P.bates, digits=2)
  df_web$qval = signif(df_web$qval, digits=2)
  df_web$zrank_signed_mean = signif(df_web$zrank_signed_mean, digits=2)
  df_web$gene[is.na(df_web$gene)] = '*unnamed*'
  return(df_web)
}

saveGeneSignatureWebTable = function(df, maxGenesPerPathway=300, filename=file.path(SITE_DIR,"speed2_gene_signatures.json")) {
  df_web = formatDataFrameForWeb(df)
  df_web = df_web %>%
    dplyr::select(p_id, gene, regulation, nexp, zrank_signed_mean, P.bates, qval) %>%
    dplyr::rename(Pval=P.bates, Regulation=regulation, Qval=qval, Gene=gene, Nexp=nexp, Pathway=p_id, Zrank=zrank_signed_mean) %>%
    filter(Qval<0.05) %>%
    group_by(Pathway) %>%
    top_n(maxGenesPerPathway, -Pval)
  x <- toJSON(unname(split(df_web, 1:nrow(df_web))))
  fileConn<-file(filename)
  writeLines(x, fileConn)
  close(fileConn)
  #return(cat(x))
}



#############################################
############## PARSE ZVALUES ################
#############################################

# /groups/nils/members/mattias/SPEED2/
df_zvalues = read_tsv(file = "speed2_zvalues.tsv.gz")
df_exp = read_tsv(file = "speed2_experiments.tsv") %>%
  dplyr::select(e_id,effect)

speed2_signatures = left_join(df_zvalues,df_exp) %>%
  #head(n=100000) %>%
  group_by(p_id,e_id) %>%
  filter(expression > median(expression)) %>%
  ungroup %>%
  filter(effect %in% c("activation", "inhibition")) %>% # filter out unknown effect
  group_by(p_id, gene, g_id, e_id, effect) %>%
  summarise(zvalue=zvalue[which.max(expression)]) %>%
  ungroup %>%
  mutate(zvalue_signed = signZvalue_vec(zvalue, effect)) %>%
  group_by(p_id, e_id) %>%
  mutate(zrank_signed=abs(Ub-Ua)*percent_rank(zvalue_signed)-max(Ua,Ub)) %>%
  ungroup

# filter out Wnt activating experiments (unreliable)
speed2_signatures = speed2_signatures %>%
  filter(!(p_id=="Wnt" & effect=="activation"))

df_numExpPerPathway = speed2_signatures %>%
  group_by(p_id, e_id) %>%
  dplyr::count() %>%
  group_by(p_id) %>%
  summarize(nexp_pathway=n())
  #%>% filter(!(p_id %in% BAD_PATHWAYS))

write_tsv(df_numExpPerPathway, path="df_numExpPerPathway.tsv")

z <- gzfile("speed2_signatures_noPval_noQval.tsv.gz")
write_tsv(speed2_signatures, path=z, col_names = T)


#############################################
######## COMPUTE P-VALUES (BATES) ###########
#############################################

speed2_signatures = speed2_signatures %>%
  bates_pipe

write_tsv(speed2_signatures, path="speed2_signatures_noQval.tsv", col_names = T)


#############################################
#### ADD Q-VALUES (SEPERATELY COMPUTED) #####
#############################################

# get p-values for randomly shuffled control (do many times)
#speed2_signatures_randomized = speed2_signatures %>%
#  group_by(p_id, e_id) %>%
#  mutate(zrank_signed = sample(zrank_signed))

#


# Compute q-values function
computeQ = function(p_id, P.bates) {
  p_id = unique(p_id)[1]
  pvals =  P.bates

  # /groups/nils/members/mattias/SPEED2/control/
  env = LoadToEnvironment(paste("control/control_",p_id,".RData", sep=""))
  pvals_control = env$control$P.bates
  pvals_sorted_ix = order(pvals)
  pvals_sorted = pvals[pvals_sorted_ix]
  pvals_control_sorted = pvals_control[order(pvals_control)]

  num_leq_in_control = rep(NA, length(pvals))
  control_ix_prev = 1
  control_ix_curr = 1
  while(pvals_control_sorted[control_ix_curr]<=pvals_sorted[1]) {
    control_ix_curr = control_ix_curr + 1
  }
  num_leq_in_control[1] = control_ix_curr-control_ix_prev

  for(pval_i in 2:length(pvals_sorted)) {
    pval = pvals_sorted[pval_i]
    while((pvals_control_sorted[control_ix_curr]<=pval) & (control_ix_curr<length(pvals_control_sorted))) {
      control_ix_curr = control_ix_curr + 1
    }
    num_leq_in_control[pval_i] = num_leq_in_control[pval_i-1]+(control_ix_curr-control_ix_prev)
    control_ix_prev = control_ix_curr
  }

  len_pvals_control = length(pvals_control_sorted)
  len_pvals = length(pvals_sorted)

  fdrs = sapply(1:length(pvals_sorted), function(i) (max(num_leq_in_control[i], 1)/len_pvals_control)/(i/len_pvals))

  qvals = rep(NA,len_pvals)
  qvals[length(qvals)]=fdrs[length(fdrs)]
  for(i in (length(qvals)-1):1) {
    qvals[i] = min(qvals[i+1], fdrs[i])
  }
  qvals = qvals[match(1:length(qvals),pvals_sorted_ix)]
  return(qvals)
}


speed2_signatures = speed2_signatures %>%
  group_by(p_id) %>%
  mutate(qval=computeQ(p_id, P.bates)) %>%
  ungroup %>%
  mutate(P.bates = abs(P.bates)) %>%
  mutate(P.BH = p.adjust(P.bates, method="BH"))


#save(speed2_signatures_bates_qval, file="speed2_signatures_bates_qval.RData")
speed2_signatures = left_join(speed2_signatures, df_numExpPerPathway)
write_tsv(speed2_signatures, path="speed2_signatures.tsv", col_names = T)

# print pathway statistics pipe
speed2_signatures %>%
  filter(qval<0.05) %>%
  #filter(P.BH<0.05) %>%
  arrange(p_id, P.bates) %>%
  group_by(p_id) %>%
  summarize(n=n()) %>%
  arrange(desc(n))

#save.image("workspace15.RData")
