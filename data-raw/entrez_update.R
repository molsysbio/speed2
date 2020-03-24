options(stringsAsFactors=FALSE)
library(RSQLite)
library(DBI)
library(tidyverse)
library(httr)

# legacy, SQL no longer used
con = dbConnect(RSQLite::SQLite(), dbname=paste(DATABASE_DIR, "SPEED2.db", sep=""))
query <- "SELECT probe_id,g_id FROM Data"
res <- dbGetQuery(con,query)
dbDisconnect(con)

env = LoadToEnvironment("speed2_maxprobe_signed_bates.RData")

gene_ids = unique(speed2_maxprobe_signed_bates$g_id)
BATCH_MAX = 1000
query = 'scopes=entrezgene,retired&q='
updated_ids = list()
for(ix in 1:length(gene_ids)) {
  gene_id = gene_ids[ix]
  if((ix %% BATCH_MAX)==1)
    sep=""
  else
    sep=","

  query = paste(query,as.character(gene_id), sep=sep)

  if((ix %% BATCH_MAX)==0) {
    query = paste(query, "&fields=name,symbol", sep="")
    print(query)

    res <- POST('http://mygene.info/v3/query', body = query,
                encode="json",
                add_headers(`content-type`='application/x-www-form-urlencoded'),
                verbose())
    updated_ids = append(updated_ids, content(res))
    query = 'scopes=entrezgene,retired&q='

  }
}

# POST slack
if((ix %% BATCH_MAX)!=0) {
  query = paste(query, "&fields=name,symbol", sep="")
  print(query)

  res <- POST('http://mygene.info/v3/query', body = query,
              encode="json",
              add_headers(`content-type`='application/x-www-form-urlencoded'),
              verbose())
  updated_ids = append(updated_ids, content(res))

}

save(updated_ids, file="entrez_update.RData")

df_updated_entrezId = do.call(rbind, lapply(updated_ids, function(l) {
  if("notfound" %in% names(l) ) {
    return(data.frame(old_id=l$query, new_id=NA, symbol=NA))
  }
  return(data.frame(old_id=l$query, new_id=l$`_id`, symbol=l$symbol) )
}))

df_updated_entrezId = as_tibble(df_updated_entrezId)

save(df_updated_entrezId, file="df_updated_entrezId.RData")


### GET STATISTICS ###
N_retired = 0
N_updated = 0
N_ok = 0
retired_gene_ids = c()
updated_gene_ids = c()
for(i in 1:length(updated_ids)) {
  gene_obj = updated_ids[[i]]
  if('notfound' %in% names(gene_obj)) {
    print(paste("Retired entrez id: ",gene_obj$query ))
    N_retired = N_retired + 1
    retired_gene_ids = c(retired_gene_ids, gene_obj$query)
  } else if(gene_obj$query != gene_obj[["_id"]]) {
    N_updated = N_updated + 1
    print(paste("Updated entrez id: ",gene_obj$query, " -> ",gene_obj[["_id"]] ,sep=""))
    updated_gene_ids = c(updated_gene_ids, gene_obj$query)
  } else {
    N_ok = N_ok + 1
  }
}

print(paste("N_retired:", N_retired))
print(paste("N_updated:", N_updated))
print(paste("N_ok:", N_ok))
