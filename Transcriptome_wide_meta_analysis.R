#Transcriptome-wide meta-analysis for hippocampus

require(data.table)
require(plyr)
require(org.Hs.eg.db)
require(AnnotationDbi)
require(EnsDb.Hsapiens.v86)
require(dtangle)
require(dtangle.data)
require(limma)
require(sva)
require(ggplot2)
require(dplyr)
require(corrplot)
require(GEOquery)
require(tidyverse)
require(ggrepel)


##--> 1. Differential gene-expression analysis in each study of hippocampus
affy_study = c("Blalock04","Hokama13","Berchtold14")
affy_dir = list.dirs(path = "D:/AD_gx/Brain_gx/affy" ,full.names = T,recursive = F)

illu_study = c("Miller13")
illu_dir = list.dirs(path = "D:/AD_gx/Brain_gx/illumina" ,full.names = T,recursive = F)


seq_study  = c("ACT","Semick19","Rooji18")
seq_dir = list.dirs(path = "D:/AD_gx/Brain_gx/seq" ,full.names = T,recursive = F) 

study = c(affy_study, illu_study,  seq_study)
all_dir = c(affy_dir,illu_dir, seq_dir)


sva_cor = list()
results = list()
for_meta = list()

for(z in 1:length(study)){

dir = all_dir[grepl(study[z], all_dir)]	

exprs_files = list.files(path = dir, pattern = "_GX_nodup_covarQced_combat_dtangle_noLow_symbolQced_20220712.txt", full.names = TRUE)

exprs_read = fread(exprs_files[grepl("Hippocampus|hippocampus|HIPPO",exprs_files)], header = T)
exprs_read = data.frame(exprs_read)

exprs_df = exprs_read[,!grepl("Sample_ID|FACTOR_",colnames(exprs_read)) ]
row.names(exprs_df) = exprs_read$Sample_ID 

predictors = exprs_read[,grepl("Sample_ID|FACTOR_",colnames(exprs_read)) ]
predictors = data.frame(predictors)
row.names(predictors) = exprs_read$Sample_ID 


##-->
# get_stats = predictors %>% dplyr::group_by(FACTOR_dx) %>% dplyr::mutate(age_mean = mean(FACTOR_age), age_sd = sd(FACTOR_age), sex_portion = sum(FACTOR_sex == "female")/length(FACTOR_sex) )
# get_stats = unique(get_stats[,c("FACTOR_dx","age_mean","age_sd","sex_portion")])
# get_stats$age_stats = paste(format(round(get_stats$age_mean, 2), nsmall = 2), paste0("(", format(round(get_stats$age_sd, 2), nsmall = 2), ")"), sep = " ")
# get_stats$female_stats = format(round(get_stats$sex_portion*100, 1), nsmall = 1)
# get_stats_ad = get_stats[get_stats$FACTOR_dx == "AD", c("age_stats","female_stats")]
# colnames(get_stats_ad) = paste("AD", colnames(get_stats_ad), sep = "_")
# get_stats_ctl = get_stats[get_stats$FACTOR_dx == "CTL", c("age_stats","female_stats")]
# colnames(get_stats_ctl) = paste("CTL", colnames(get_stats_ctl), sep = "_")
# get_stats_df = cbind(get_stats_ad, get_stats_ctl)
# get_stats_df$studyID = study[z]

# results[[z]] = get_stats_df
# }
# results_df = ldply(results)
# fwrite(results_df, file = "D:/AD_gx/Brain_gx/20221103_AD_brain_hippocampus_study_demographic_table.csv")
##<--


for(col in 1:ncol(predictors)){
	if(colnames(predictors[,col, drop = FALSE]) %in% c("FACTOR_dx")){
		predictors[,col] = factor(predictors[,col], levels = c("CTL","AD"))
	} else if(colnames(predictors[,col, drop = FALSE]) %in% c("FACTOR_sex","FACTOR_ethnicity","FACTOR_APOE","FACTOR_batch","FACTOR_Lot","FACTOR_ScanDate","FACTOR_pres","FACTOR_site" )){
		predictors[,col] = factor(predictors[,col])

	} else if(colnames(predictors[,col, drop = FALSE]) %in% c("FACTOR_tissue","Sample_ID","Sample_ID2","Sample_ID3")){
		predictors[,col] = as.character(predictors[,col])

	}else {
		predictors[,col] = as.numeric(predictors[,col])

	}
}

predictors = predictors[,! grepl("batch|APOE|Lot|pres|ScanDate|_site|tissue|Sample|FACTOR_microglia",colnames(predictors))]

predictors_mod = predictors
pheno = predictors_mod[, !colnames(predictors_mod) == "FACTOR_dx"]
mod = model.matrix(~ ., data= predictors_mod) # model with outcome of interest and covariates
mod0 = model.matrix(~.,data= pheno) #  model with covariates for adjustment

##--> Find SV
#column: samples
exprs = t(exprs_df)
n.sv = min(c(num.sv(exprs,mod, method = 'leek'), num.sv(exprs,mod, method = 'be')))

if(n.sv > 0){

cat("\ndetected sv in: ", study[z])
cat("\nnumber of sv: ",  n.sv)

svobj = sva(exprs,mod, mod0, n.sv = n.sv)
svdf = data.frame(svobj$sv)	


if(ncol(svdf) > 0){
colnames(svdf) = paste("FACTOR_SV",1:ncol(svdf), sep = "")

cell = predictors[, colnames(predictors) %in% c("FACTOR_astrocytes","FACTOR_microglia","FACTOR_neurons","FACTOR_oligodendrocytes")]
cors = cor(data.frame(svdf, cell),use='pairwise.complete.obs')

png(paste(dir,"/dtangle_SVDF_correlation_", study[z],"_hippocampus.png", sep =""), res=300,units="in",height=6,width=6)
corrplot::corrplot(cors, tl.col = 'black', number.cex = .5,
                   main = paste(study[z],". dtangle. Tissue: hippocampus", sep = ""),
                   mar=c(0,0,3,0),
                   order='hclust',method='color', addCoef.col = "black",
                   addrect = 3, rect.lwd = 5,
                   rect.col = 'black', outline = T,
                   tl.srt = 45, tl.cex = 0.75)
dev.off()

lm.get = list()
code = ifelse(predictors$FACTOR_dx == "AD", 1, 0)

for(n in 1:ncol(svdf)){
data = 	svdf[,n, drop = FALSE]
lm.cor = glm(code ~., data = data)
lm.sum = data.frame(summary(lm.cor)$coefficients)
lm.get[[n]] = t(lm.sum[!grepl("Intercept",row.names(lm.sum)), colnames(lm.sum) == "Pr...t..", drop = FALSE]) 
}
lm.get_df = data.frame(studyID = study[z], do.call(cbind,lm.get))
sva_cor[[z]] = lm.get_df

predictors = cbind(predictors, svdf)	

}
}

cat("\nall covariates in study ", study[z], " : ", colnames(predictors),"\n")

exprs_com = cbind(exprs_df, predictors)
fwrite(exprs_com, file = paste(dir,"/", study[z],"_for_hippocampus_DEG_addSV_20230617.txt", sep =""),row.names = T, col.names = T, sep= "\t")


##--> Get DEG
#column: samples
exprs = t(exprs_df)

cat("\n", study[z], "- final design covariates: ", colnames(predictors))
mod = model.matrix(~ ., data = predictors)

#rows corresponding to genes and columns to samples
fit = lmFit(exprs,mod)

eb = eBayes(fit)
top = topTable(eb, coef = "FACTOR_dxAD", p.value =1,  number=Inf, sort.by="p",adjust="fdr", confint=TRUE)

results[[z]] = top

mt = match(row.names(fit$coefficients), row.names(top))
top_mt = top[mt,]

coefs_get = fit$coefficients
s2 = sqrt(eb$s2.post)
stdev = fit$stdev.unscaled[,colnames(fit$stdev.unscaled) == "FACTOR_dxAD"]
sdr_get = s2*stdev

for_meta[[z]] = data.frame(coefs_get = coefs_get, sdr_get = sdr_get, n_CN = sum(grepl("CTL",predictors$FACTOR_dx)), n_AD = sum(grepl("AD",predictors$FACTOR_dx)), adjP = top_mt$adj.P.Val )
}


sva_cor_df = do.call(dplyr::bind_rows, sva_cor)
fwrite(sva_cor_df, file = "D:/AD_gx/Brain_gx/20230617_DGEs_hippocampus_SVs_AD_association.txt", row.names = F, col.names = T, sep= "\t")

names(results) = study
results_gene = lapply(results, function(x) cbind(x, geneID = row.names(x)))
results_df = ldply(results_gene)

fwrite(results_df,
           file = "D:/AD_gx/Brain_gx/20230617_DGEs_table_hippocampus_NoAPOEcovar.txt", 
            row.names = F, col.names = T, sep= "\t")


names(for_meta) = study
for_meta_gene = lapply(for_meta, function(x) cbind(x, geneID = row.names(x)))
for_meta_df = ldply(for_meta_gene)

fwrite(for_meta_df,
           file = "D:/AD_gx/Brain_gx/20230617_for_meta_for_hippocampus_NoAPOEcovar.txt", 
            row.names = F, col.names = T, sep= "\t")



##--> 2. Meta-analysis to identify AD-associated genes in each of six brain regions
require(metafor)

#Protein coding genes were retained
pc_gene = fread("D:/AD_gx/Brain_gx/20220727_all_protein_coding_genes_in_studies.txt")

for_meta_files = list.files(path = "D:/AD_gx/Brain_gx", pattern = "20220727_for_meta", full.names = TRUE)
for_meta_files[grepl("hippocampus", for_meta_files)] = "D:/AD_gx/Brain_gx/20230617_for_meta_for_hippocampus_NoAPOEcovar.txt"

res_region = list()
for(r in 1:length(for_meta_files)){

region = strsplit(basename(for_meta_files[[r]]),"_")[[1]][5]

for_meta_df = fread(for_meta_files[r])
for_meta_df = data.frame(for_meta_df)
for_meta_df = for_meta_df[,colnames(for_meta_df) %in% c("coefs_get.FACTOR_dxAD","sdr_get","n_CN","n_AD","geneID",".id","adjP")]
colnames(for_meta_df)[1:3] = c("studyID","estimate","sdr")

print(table(for_meta_df$studyID))

genes = unique(for_meta_df$geneID)
genes_pc = genes[gsub("[.]","-",genes) %in% pc_gene$gene_name] 
#length(genes_pc)
#19815

res_save = list()

for( i in 1:length(genes_pc)){
  
  cat("\rMeta-analysis:",i,"/",length(genes_pc),"\n")
  sub = for_meta_df[for_meta_df$geneID %in% genes_pc[[i]], ]
  sub$N = sub$n_AD + sub$n_CN
  
  sub_combn = data.frame(estimate = c(sub$estimate),
                         SE = c(sub$sdr),
                         N = c(sub$N),
                         N_case = c(sub$n_AD),
                         N_control = c(sub$n_CN),
                         studyID = sub$studyID)
  
  if(nrow(sub) <= 2) next

  direction = sign(sub_combn$estimate)
  direction = ifelse(direction == 1, "+", "-")
  direction = paste(direction, collapse= "", sep = "")
  
  res = NA
  
  res = metafor::rma(yi = sub_combn$estimate,
                                          sei = sub_combn$SE, 
                                          verbose = FALSE,
                                          weighted = TRUE, method = "DL")
  
  res_wei = weights(res)
  
  res = data.frame(GeneSymbol = genes_pc[[i]],
                   Direction = direction, 
                   Nstudy = nrow(sub_combn),
                   Nsample = sum(sub_combn$N),
                   Ncase = sum(sub_combn$N_case),
                   Ncontrol = sum(sub_combn$N_control),
                   Log2FC = res$beta, 
                   SE = res$se, 
                   weights = paste(res_wei, collapse =  ", ", sep = ""),
                   # BetaStd = res.std$beta,
                   # SEstd = res.std$se,
                   P = res$pval,
                   HetQ = res$QE,
                   HetP = res$QEp,
                   study_kept = paste(sub_combn$studyID, collapse =  ", ", sep = ""))
  
  res_save[[i]] = res
}
res_save_df = ldply(res_save)
res_save_df$region = region
res_save_df$FDR = p.adjust(res_save_df$P, "fdr")
res_save_df$BONF = p.adjust(res_save_df$P, "bonferroni")

res_region[[r]] = res_save_df
}


res_df = ldply(res_region)
res_df = res_df[order(res_df$FDR, decreasing = F),]

fwrite(data.table(res_df), 
       sep = "\t",
       file = "D:/AD_gx/Brain_gx/20230617_brain_meta_DEGs_6region_PCgene.txt",
       quote  = F, row.names = F)

table(res_df$region[res_df$FDR<0.05])
# cerebellum       DLPFC  entorhinal     frontal hippocampus    temporal 
#       1521        3901         116          17          30        2967 
