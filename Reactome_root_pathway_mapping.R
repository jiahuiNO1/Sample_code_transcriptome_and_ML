#Identify the root pathway of each Reactome sub-pathway

require(data.table)
library(ggplot2)
library(dplyr)
library(plyr)
library(tidyverse)


Cfolder = "D:/AD_gx/Brain_gx/piano"
if(!dir.exists(Cfolder)){
  dir.create(Cfolder)
}
setwd(Cfolder)

##--> 1. Extract gene sets used in enrichment analysis for each brain region
extract_func = function(x){
test = fread(paste0("20230618_AD_Blood_Brain_piano_reactome_meanTest_summary_",x,".csv"))
summary = data.frame(setnames = test$setnames)

reactome = fread("20220329_AD_brain_piano_reactome_annotation.csv")
colnames(reactome)[3:4] = c("Term","ID")
reactome$Database = "REACTOME"
reactome = unique(reactome[,c("ID","Term","Database","PATHWAY")])

summary$setnames2 = gsub("gsc.","", summary$setnames)

summary_mer = merge(summary, reactome, by.x = 'setnames2', by.y = "PATHWAY" )

return(summary_mer)
}

#
regions = unique(res_all$region)
set_list = lapply(regions, extract_func)
#lapply(set_list, nrow)
names(set_list) = regions

set_df = ldply(set_list)
write.csv(set_df, file = "20240604_AD_piano_Reactome_blood_brain_genesets_in_analysis.csv")


##--> 2. Find root pathways for all the subpathways in reactome annotation file
require(stringr)

reactome = fread("20220329_AD_brain_piano_reactome_annotation.csv")
reactome$PATHNAME = str_trim(reactome$PATHNAME)
reactome = unique(reactome[,c("PATHNAME","REACTOMEID")])
nrow(reactome)
#[1] 2495

reactome_target = reactome

#Reactome pathway hierarchy relationship file
hier = fread("D:/AD_gx/Brain_gx/WGCNA/ReactomePathwaysRelation.txt", header = F)
colnames(hier) = c("parent","child")

#All 29 Reactome root pathways
root_path = c("R-HSA-9612973", "R-HSA-1640170","R-HSA-1500931","R-HSA-8953897","R-HSA-4839726","R-HSA-400253","R-HSA-1266738",
		      "R-HSA-8963743","R-HSA-1643685","R-HSA-73894","R-HSA-69306","R-HSA-9748784","R-HSA-1474244","R-HSA-74160","R-HSA-109582",
		      "R-HSA-168256","R-HSA-1430728","R-HSA-392499","R-HSA-8953854","R-HSA-397014","R-HSA-112316","R-HSA-1852241","R-HSA-5357801",
		      "R-HSA-9609507","R-HSA-1474165","R-HSA-9709957","R-HSA-162582","R-HSA-382551","R-HSA-5653656")

#Map sub-pathway to root pathway
ID_get = reactome_target$REACTOMEID[reactome_target$REACTOMEID %in% c(hier$child)]

root_list = list()
for(i in 1:length(ID_get)){
    
    tar_path = ID_get[i]

	level = hier[hier$child == tar_path,]
	condition1 = sum(level$parent %in% hier$child)
	condition2 = sum(level$parent %in% root_path)
	
	parent_hold = c()

	if(nrow(level) > 1 & condition2 > 0){
    	parent_hold = c(parent_hold, level$parent[level$parent %in% root_path])
    	level = level[!level$parent %in% parent_hold,]
    	condition1 = sum(level$parent %in% hier$child)
	    condition2 = sum(level$parent %in% root_path)
    }

	while(condition1 > 0 & condition2 == 0){

	level = hier[hier$child %in% level$parent,]
	condition1 = sum(level$parent %in% hier$child)
	condition2 = sum(level$parent %in% root_path)
    
    if(nrow(level) > 1 & condition2 > 0){
    	parent_hold = c(parent_hold, level$parent[level$parent %in% root_path])
    	level = level[!level$parent %in% parent_hold,]
    	condition1 = sum(level$parent %in% hier$child)
	    condition2 = sum(level$parent %in% root_path)
    }

	}

	root_list[[i]] = data.frame(tar_path = ID_get[i], root_path = c(parent_hold, unique(level$parent)))	
	
}

root_df = ldply(root_list)

#add 29 rooth pathways
add_df = data.frame(tar_path = c(root_path), root_path = c(root_path))

root_df_add = rbind(root_df, add_df)
root_df_add = root_df_add[!duplicated(root_df_add),]

write.csv(root_df_add, file = "20240604_AD_Reactome_all_subpathway_root.csv")


##--> 3. Map Reactome sub-pathways included in the enrichment analysis to their root pathway 
set_mer = merge(set_df, root_df_add, by.x = "ID", by.y = "tar_path", all.x = T)

set_count = set_mer %>% dplyr::group_by(.id, root_path) %>% dplyr::summarise(count = n()) 
write.csv(set_count, file = "20240605_AD_Reactome_subpathway_count_within_root_across_tissue.csv")

