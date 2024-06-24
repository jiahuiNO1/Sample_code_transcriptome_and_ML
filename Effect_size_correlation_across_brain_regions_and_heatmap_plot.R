#1. Examine the per-gene effect size correlation across brain regions
#2. Heatmap plot of effect size correlation across brain regions

##--> 1. Examine the per-gene effect size correlation across brain regions across the whole transcriptome
#Load imputed brain data
res_df = fread("20240429_BrainGENIE_meta_summary_all_genes_10regions.txt") %>% 
				select(GeneSymbol, Log2FC, region) %>%
                mutate(region = case_when(
                              region == "Cerebellar_Hemisphere" ~ "cerebellum",
                              region == "Frontal_Cortex_BA9" ~ "DLPFC",
                              region == "Hippocampus" ~ "hippocampus",
                              TRUE ~region
                            ))

#Load directly measured brain data
res_brain = fread("D:/AD_gx/Brain_gx/20230617_brain_meta_DEGs_6region_PCgene.txt") %>% select(GeneSymbol, Log2FC, region) %>%
                  filter(region %in% c("cerebellum", "DLPFC", "hippocampus")) %>% mutate(region = paste0("measured_", region))

#Load directly measured blood data
res_blood = fread("D:/AD_gx/Blood_gx/DEG/20220725_blood_meta_all8_geneLOC_pcGene.txt")
res_blood$region = "blood"
res_blood = res_blood %>% select(GeneSymbol, Log2FC, region)

#Combine imputed and measured blood and brain data
sum_get = data.frame(rbind(res_df, res_brain, res_blood))

res_test = sum_get %>% spread(key = region, value = Log2FC)
dim(res_test)
#[1] 21440    15

out <- combn(res_test[,2:ncol(res_test)], 2, FUN = function(x) 
             cor.test(x[[1]], x[[2]], method = "spearman", conf.level = 0.95), simplify = FALSE)

names(out) <- combn(names(res_test[,2:ncol(res_test)]), 2, FUN = paste, collapse='/')

mydata = do.call(rbind, Map(cbind, corgroups = names(out), 
		         unname(lapply(out, function(x)
		         data.frame(cor.value = x$estimate, cor.pvalue = x$p.value)))))

mydata$group1 = ldply(strsplit(mydata$corgroups, "/"))[,1]
mydata$group2 = ldply(strsplit(mydata$corgroups, "/"))[,2]

fwrite(data.table(mydata), file ="20240429_AD_BrainGENIE_DEGs_effect_size_cor_between_10R_blood_brain.csv", col.names = TRUE, row.names = FALSE, sep = ",")


##--> 2. Heatmap plot of effect size correlation across brain regions
mydata_all = fread("20240429_AD_BrainGENIE_DEGs_effect_size_cor_between_10R_blood_brain.csv")

transform_group = function(group_ori) {
  dplyr::case_when(
    group_ori == "Anterior_cingulate_cortex_BA24" ~ "ACC",
    group_ori == "blood" ~ "Blood (measured)",
    group_ori == "Caudate_basal_ganglia" ~ "Caudate",
    group_ori == "cerebellum" ~ "Cerebellum",
    group_ori == "hippocampus" ~ "Hippocampus",
    group_ori == "measured_cerebellum" ~ "Cerebellum (measured)",
    group_ori == "measured_DLPFC" ~ "DLPFC (measured)",
    group_ori == "measured_hippocampus" ~ "Hippocampus (measured)",
    group_ori == "Nucleus_accumbens_basal_ganglia" ~ "Nucleus accumbens",
    group_ori == "Putamen_basal_ganglia" ~ "Putamen",
    group_ori == "Substantia_nigra" ~ "Substantia nigra",
    TRUE ~ group_ori
  )
}

mydata_all =  data.frame(mydata_all) %>%
			  dplyr::rename(group1_ori = group1, group2_ori = group2) %>%
			  dplyr::mutate(group1 = transform_group(group1_ori),
			                group2 = transform_group(group2_ori))
			  
#
mydata = mydata_all %>% select(group1, group2, cor.value)

add_df = mydata[1:13,]
add_df$group1 = add_df$group2
add_df$group2 = "Amygdala"

mydata1 = mydata2 = rbind( mydata, add_df) 

tissues = unique(mydata2$group1)
for(t in 1:length(tissues)){
  add_df = data.frame(group1 = NA, group2 = NA, cor.value = NA)
  add_df$group1 = tissues[t]
  add_df$group2 = tissues[t]
  add_df$cor.value = 1
  mydata2 = rbind(mydata2, add_df)
}

#
require(Matrix)

mydata_get1 = mydata1 %>% pivot_wider(names_from = group2, values_from = cor.value)
mydata_get1 = as.data.frame(mydata_get1)
mydata_get1 = mydata_get1[,c(1, ncol(mydata_get1), 2:(ncol(mydata_get1)-1))]
row.names(mydata_get1) = mydata_get1$group1
mydata_get1 = as.matrix(mydata_get1[,-1])
mydata_mt1 = forceSymmetric(mydata_get1 , uplo = "U")

#
mydata_get2 = mydata2 %>% pivot_wider(names_from = group2, values_from = cor.value)
mydata_get2 = as.data.frame(mydata_get2)
mydata_get2 = mydata_get2[,c(1, ncol(mydata_get2), 2:(ncol(mydata_get2)-1))]
row.names(mydata_get2) = mydata_get2$group1
mydata_get2 = as.matrix(mydata_get2[,-1])
mydata_mt2 = forceSymmetric(mydata_get2 , uplo = "U")

#Plot the heatmap
library("pheatmap")
g1 = pheatmap(
  mydata_mt1, 
  clustering_distance_cols = as.dist(1 - mydata_mt2),
  na_col = "grey90",
  #fontsize = 14,
  clustering_distance_rows = as.dist(1 - mydata_mt2),
  #cutree_rows = 7,
  #cutree_cols = 7
  #border_color = "black"
  #show_colnames = FALSE
  )


#add tissue annotation
cluster_df = data.frame(Cluster = cutree(g1$tree_col, k=2))
cluster_df$study_ID = row.names(cluster_df)
cluster_df = cluster_df[g1$tree_col$order, ]

cluster_df$Data = ifelse(grepl("Blood|measured", cluster_df$study_ID), "Measured","Imputed")
cluster_df = cluster_df[,3, drop = F]


#add p-value annotation in cells
pdata = mydata_all %>% select(group1, group2, cor.pvalue)

add_df = pdata[1:13,]
add_df$group1 = add_df$group2
add_df$group2 = "Amygdala"

pdata1 = data.frame(rbind(pdata, add_df)) 

pdata_get1 = pdata1 %>% pivot_wider(names_from = group2, values_from = cor.pvalue)
pdata_get1 = as.data.frame(pdata_get1)
pdata_get1 = pdata_get1[,c(1, ncol(pdata_get1), 2:(ncol(pdata_get1)-1))]
row.names(pdata_get1) = pdata_get1$group1
pdata_get1 = as.matrix(pdata_get1[,-1])
pdata_mt1 = forceSymmetric(pdata_get1 , uplo = "U")

pdata_df = as.data.frame(as.matrix(pdata_mt1)) %>% mutate_all(~ if_else((. < 0.05 | is.na(.)), "", "NS"))


#
#library(RColorBrewer)
#brewer.pal(7,"Greens")
annoCol<-list(Data_type=c(Measured = "#008000",Imputed = "#DDA0DD"))

require(dichromat)
clrsp <- colorRampPalette(c("skyblue","white", "pink","red","darkred"))   
clrs <- clrsp(182) 

breaks1 <- seq(-0.2, 0.8, length.out = 182)

#Plot the heatmap with annotation
g1 = pheatmap(
  mydata_mt1, 
  clustering_distance_cols = as.dist(1 - mydata_mt2),
  na_col = "grey90",
  clustering_distance_rows = as.dist(1 - mydata_mt2),
  cutree_rows = 2,
  cutree_cols = 2,
  annotation_col = cluster_df,
  annotation_colors = annoCol,
  display_numbers = as.matrix(pdata_df),
  color = clrs,
  breaks = breaks1
  )

png("20240429_AD_BrainGENIE_DEGs_efect_size_cor_between_tissues_heatmap.png",res=300,units="in",height = 9, width = 10)
print(g1)
dev.off()


#what are the optimal number of clusters?
library(factoextra)
library(cluster)

hclusCut <- function(x, k, d.meth = "euclidean", ...)
            list(cluster = cutree(hclust(dist(x, method=d.meth), ...), k=k))

gap_stat <- clusGap(mydata_mt2, FUN = hclusCut, K.max = 7, B = 500)

dat <- data.table(gap_stat$Tab)
#Add a column 'k' with row numbers
dat[, k := .I]

p <- ggplot(dat, aes(k, gap)) + geom_line() + geom_point(size = 3) +
     geom_errorbar(aes(ymax = gap + SE.sim, ymin = gap - SE.sim), width = 0.25) +
     ggtitle("Clustering Results") +
     labs(x = "Number of Clusters", y = "Gap Statistic") +
     theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
           axis.title = element_text(size = 12, face = "bold"))

png("20240429_AD_BrainGENIE_DEGs_efect_size_cor_between_tissues_optimal_cluster_number.png",res=300,units="in",height = 4, width = 5)
print(p)
dev.off()

