options(stringsAsFactors=FALSE)
library(tidyverse)
library(httr)

# Microarrays should already be normalized before running this script  
# i.e. each experimnent folder should contain a zvalues.dat file.

DATA_DIR = 'data_combined_23_05_2019'

#############################################
############# HELPER FUNCTIONS ##############
#############################################

clean_line = function(str) {
  ret = gsub("\t+", "\t", str )
  ret = gsub(" +", " ", ret )
  ret = gsub("\n", "", ret )
  return(ret)
}

# Use unified notation for effect + correct misspellings
spell_effect_correctly = function(str) {
  ret = str
  ret = gsub("inhibiting|inhibting", "inhibition", ret)
  ret = gsub("activating", "activation", ret)
  return(ret)
}


pickFirstGeneIds = function(ids) {
  
  pickFirstGeneId = function(str) {
    genes = strsplit( as.character(str), "///")
    if(length(genes) != 1)
      return(NA)
    if(length(genes[[1]])<1)
      return(NA)
    
    return(str_trim(genes[[1]][1]))
  }
  
  return(sapply(1:length(ids), function(i) pickFirstGeneId(ids[i])))
  
}

#############################################
############# PARSE GEO DATA ################
#############################################

pathways = list.dirs(DATA_DIR, recursive = F, full.names=F)
#pathways = "Hippo"

speed2_zvals = tibble()
speed2_experiments = tibble()

for(pathway in pathways) {
  experiments = list.files(file.path(DATA_DIR,pathway), pattern="GSE", recursive = F, full.names=F)
  for(experiment in experiments) {
    
    print(paste("Parsing: ", pathway,"/",experiment,sep=""))
    
    exp_path = file.path(DATA_DIR,pathway, experiment)
    description_file = file.path(DATA_DIR,pathway, experiment, 'Description.txt')        
    gsm_file = file.path(DATA_DIR,pathway, experiment, 'gsm_accession.dat')        
    zvalues_file = file.path(DATA_DIR,pathway, experiment, 'zvalues2.dat')

    if(!file.exists(gsm_file)) {
      print(paste("Missing gsm_accession.dat files for experiment: ", exp_path,sep=""))
      next
    }
    
    # PARSE EXPERIMENT DESCRIPTION
    cell_line = modification = effect = time = comment = ''
    if(file.exists(description_file)) {
      con = file(description_file, "r")
      lines = readLines(con)
      close(con)
      lines = as.vector( sapply(lines, function(str) clean_line(str)) )
      for(line in lines) {
        splitline = strsplit(line, "\t")
        if(length(splitline)!=1)
          next
        if(length(splitline[[1]]) != 2)
          next
        
        descriptor = str_trim(splitline[[1]][1])
        value = str_trim(splitline[[1]][2])
        
        if(descriptor=="cell line")
          cell_line = value
        else if(descriptor=="modification")
          modification = value
        else if(descriptor=="effect")
          effect = spell_effect_correctly(value)
        else if(descriptor=="time")
          time = value
        else if(descriptor=="comment")
          comment = value
        
      }
    }
    
    # PARSE GSM
    stimuli = c()
    controls = c()
    if(file.exists(gsm_file)) {
      con = file(gsm_file, "r")
      lines = readLines(con)
      close(con)
      lines = as.vector( sapply(lines, function(str) clean_line(str)) )
      for(line in lines) {
        splitline = strsplit( str_trim(line), " ")
        if(length(splitline)!=1)
          next
        if(length(splitline[[1]]) != 2)
          next
        
        gsm = str_trim(splitline[[1]][1])
        value = str_trim(splitline[[1]][2])
        
        if(value=="1")
          stimuli = c(stimuli, gsm)
        else if(value=="2")
          controls = c(controls, gsm)
      }
    } else {
      print(paste("Missing gsm_accession.dat files for experiment: ", exp_path,sep=""))
      next
    }
    
    # PARSE ZVALUES
    if(file.exists(zvalues_file)) {
      zvals = read_tsv(zvalues_file, col_names = T)
      zvals = zvals %>% 
        filter(!(is.na(g_id))) %>% 
        mutate(g_id = pickFirstGeneIds(g_id)) %>% 
        mutate(g_id = as.integer(g_id)) %>% 
        filter(!(is.na(g_id))) 
        
      zvals$p_id = pathway
      zvals$e_id = experiment
      
      zvals = zvals %>% 
        dplyr::select(p_id, e_id, g_id, zvalue, expression)
      speed2_zvals = bind_rows(speed2_zvals, zvals)
      
      
    } else {
      print(paste("Missing zvalues2.dat files for experiment: ", exp_path,sep=""))
      next
    }
    
    speed2_experiments = bind_rows(speed2_experiments, tibble(p_id=pathway,
                                                              e_id=experiment,
                                                              cell_line=cell_line,
                                                              modification=modification,
                                                              effect=effect,
                                                              time=time,
                                                              comment=comment,
                                                              stimuli=paste(stimuli, collapse=","),
                                                              controls=paste(controls, collapse=",")))
    
    
    
  }
}

#z <- gzfile("speed2_zvalues_obsoluteEntrez.tsv.gz")
#write_tsv(speed2_zvals, path=z, col_names = T)
write_tsv(speed2_experiments, path="speed2_experiments.tsv", col_names = T)



#############################################
############# UPDATE ENTREZ ID ##############
#############################################

if(file.exists("entrez_updated.tsv")) {
  #env = LoadToEnvironment("df_updated_entrezId.RData")
  #df_updated_entrezId = env$df_updated_entrezId
  df_updated_entrezId = read_tsv("entrez_updated.tsv")
} else {
  
  # fetch updated entrez id from web api
  gene_ids = unique(speed2_zvals$g_id)
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
      updated_ids = append(updated_ids, httr::content(res))
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
    updated_ids = append(updated_ids, httr::content(res))
    
  }
  
  save(updated_ids, file="updated_ids.RData")
  
  df_updated_entrezId = do.call(rbind, lapply(updated_ids, function(l) {
    if("notfound" %in% names(l) ) {
      return(data.frame(oldEntrez=l$query, newEntrez=NA, newSymbol=NA))
    }
    return(data.frame(oldEntrez=l$query, newEntrez=l$`_id`, newSymbol=l$symbol) )
  }))
  
  df_updated_entrezId = as_tibble(df_updated_entrezId)
  df_updated_entrezId$oldEntrez = as.integer(df_updated_entrezId$oldEntrez)
  df_updated_entrezId$newEntrez = as.integer(df_updated_entrezId$newEntrez)
  df_updated_entrezId$newSymbol = as.character(df_updated_entrezId$newSymbol)
  
  write_tsv(df_updated_entrezId, path="entrez_updated.tsv", col_names = T)
  
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
      #print(paste("Retired entrez id: ",gene_obj$query ))
      N_retired = N_retired + 1 
      retired_gene_ids = c(retired_gene_ids, gene_obj$query)
    } else if(gene_obj$query != gene_obj[["_id"]]) {
      N_updated = N_updated + 1
      #print(paste("Updated entrez id: ",gene_obj$query, " -> ",gene_obj[["_id"]] ,sep=""))
      updated_gene_ids = c(updated_gene_ids, gene_obj$query)
    } else {
      N_ok = N_ok + 1
    }
  }
  
  print(paste("N_retired:", N_retired))
  print(paste("N_updated:", N_updated))
  print(paste("N_ok:", N_ok))
  
}

#df_updated_entrezId$oldEntrez = as.integer(df_updated_entrezId$oldEntrez)
#df_updated_entrezId$newEntrez = as.integer(df_updated_entrezId$newEntrez)
#df_updated_entrezId$newSymbol = as.character(df_updated_entrezId$newSymbol)
speed2_zvals_updatedEntrez = left_join(speed2_zvals, df_updated_entrezId, by=c('g_id' = 'oldEntrez'))

speed2_zvals_updatedEntrez = speed2_zvals_updatedEntrez %>% 
  filter(!is.na(newEntrez)) %>% 
  mutate(gene=newSymbol) %>% 
  mutate(g_id = newEntrez) %>% 
  dplyr::select(-newEntrez, -newSymbol) %>% 
  dplyr::select(p_id,e_id,g_id,gene,zvalue,expression)

z <- gzfile("speed2_zvalues.tsv.gz")
write_tsv(speed2_zvals_updatedEntrez, path=z, col_names = T)

