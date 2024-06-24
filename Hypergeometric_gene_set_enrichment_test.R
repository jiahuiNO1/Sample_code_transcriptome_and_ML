#Hypergeometric gene-set enrichment tests for AD-associated WGCNA modules

require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(org.Hs.eg.db)
require(GO.db)
require(reactome.db)
require(data.table)
require(plyr)
require(dplyr)


#load gene names for hg19
genes = genes(TxDb.Hsapiens.UCSC.hg19.knownGene,single.strand.genes.only=FALSE)
genes = as.data.frame(genes)
genename = AnnotationDbi::select(org.Hs.eg.db, keys = as.character(genes$group_name), keytype="ENTREZID", columns=c("GENENAME","SYMBOL", "GO"))

#pull out gene annotations from reactome and gene ontology databases
#for reactome
reactome = AnnotationDbi::select(reactome.db, keys=as.character(unique(genename$ENTREZID)), keytype='ENTREZID', columns=c('PATHNAME', 'REACTOMEID'))
reactome$PATHNAME = gsub("Homo sapiens:", "REACTOME:", reactome$PATHNAME)
conv =  AnnotationDbi::select(org.Hs.eg.db, keys=as.character(unique(genename$ENTREZID)), keytype='ENTREZID', columns='SYMBOL')
reactome = merge(conv, reactome, by='ENTREZID')
colnames(reactome)[3] = 'PATHWAY'

reactome$PATHWAY = gsub("REACTOME: ", "", reactome$PATHWAY)
reactome$PATHWAY = paste(reactome$REACTOMEID, reactome$PATHWAY, sep = ":")

reactome = reactome[,colnames(reactome) %in% c("SYMBOL", "PATHWAY")]
reactome = reactome[!duplicated(reactome), ]

#for gene ontology
go_terms = AnnotationDbi::select(GO.db, keys=as.character(genename$GO), keytype="GOID", columns=c("TERM", "ONTOLOGY"))
go_terms = go_terms[!duplicated(go_terms), ]
table(go_terms$ONTOLOGY)
#BP    CC    MF 
#12396  1751  4419 

colnames(genename)[colnames(genename) == "GO"] = "GOID"
go_terms = merge(go_terms, genename[,colnames(genename) %in% c("GOID", "SYMBOL")], by="GOID")
go_terms = go_terms[!duplicated(go_terms), ]
go_terms$PATHWAY = paste(go_terms$GOID,":", gsub("[,| ]", "_", go_terms$TERM), sep = "")
go_terms = go_terms[!is.na(go_terms$GOID), ]
go_terms = go_terms[,colnames(go_terms) %in% c("SYMBOL", "PATHWAY")]
go_terms = ldply(list(go_terms, reactome))

#keep genes in the WGCNA analysis
datExpr = fread("20220918_AD_Brain_hippo_WGCNA_SamplesGenes_Combat_datExpr.txt")
datExpr = data.frame(datExpr)
row.names(datExpr) = datExpr$V1
datExpr = datExpr[,-1]

go_terms = go_terms[go_terms$SYMBOL %in% colnames(datExpr), ] # only look at genes in the datExpr object

table(grepl("R-HSA", unique(go_terms$PATHWAY)))
#FALSE  TRUE 
#15967  2318  
table(grepl("GO", unique(go_terms$PATHWAY)))
#FALSE  TRUE 
#2319 15966  

#keep gene sets with size between 10 and 100
SetSize = table(go_terms$PATHWAY)
keep_path = names(SetSize[SetSize >= 10 & SetSize <= 100]) 

go_terms = go_terms[go_terms$PATHWAY %in% keep_path,]
View(table(go_terms$PATHWAY))

All_terms = go_terms

##--> Enrichment test using GO gene-sets
go_terms = All_terms[!grepl("R-HSA",All_terms$PATHWAY),]

#for each AD-associated module
module = fread("20220918_AD_Brain_hippo_wgcna_module_membership_ori_pcGene.txt",data.table=F)
sig_mods = gsub("ME","", get_cof$.id[get_cof$bonP < 0.05])

hypStats_save = list()
for( i in 1:length(sig_mods)){
  
  cat("\nRunning hypergeometric test for module:",sig_mods[[i]])
  
  gene_grab = module[module$label %in% sig_mods[[i]], ]
  gene_grab$symbol = as.character(gene_grab$symbol)
  gene_grab = gene_grab[gene_grab$symbol %in% go_terms$SYMBOL, ]
  
  symbolSelect = gene_grab$symbol
  
  msig = go_terms
  msig = msig[!is.na(msig$SYMBOL), ]
  
  msig.split = split(msig, msig$PATHWAY)
  
  #calculate hypergeometric input values
  universe = length(colnames(datExpr)) # total # of genes in transcriptome 
  overlaps = lapply(msig.split, function(x) length(intersect(x$SYMBOL, symbolSelect)))
  geneSetSize = lapply(msig.split, nrow)
  geneSetOut = universe - unlist(geneSetSize)
  
  #hypergeomtric p-values
  hypStats = data.frame(
    Module = sig_mods[[i]],
    PATHWAY = names(msig.split),
    Population = universe, 
    Sample_success = unlist(overlaps),
    Population_success = unlist(geneSetSize),
    Population_fails = unlist(geneSetOut),
    Sample_size = length(unique(symbolSelect))
  )
  
  hypStats = hypStats[hypStats$Sample_success >= 1, ]
  
  #enrichment p-value test
  pvals = phyper(hypStats$Sample_success - 1, 
                 hypStats$Population_success, 
                 hypStats$Population_fails, 
                 hypStats$Sample_size, 
                 lower.tail = FALSE, log.p = FALSE)
  
  hypStats$P = pvals
  
  hypStats = hypStats[order(hypStats$P, decreasing = F), ]
  hypStats$FDR = p.adjust(hypStats$P, "fdr")
  hypStats$BONF = p.adjust(hypStats$P, "bonferroni")
  
  names(hypStats)[names(hypStats) %in% "PATHWAY"] = "SET_NAME"
  hypStats_save[[i]] = hypStats # save sum stats to list obj 
  
}

#combine hypergeometric test results across module eigengenes
hypStats_save_GO = ldply(hypStats_save)

#calculate fold enrichment score
hypStats_save_GO$FES = (hypStats_save_GO$Sample_success/hypStats_save_GO$Sample_size)/(hypStats_save_GO$Population_success/hypStats_save_GO$Population)


table(hypStats_save_GO$Module[hypStats_save_GO$FDR < 0.05])
#    1  17   2   3   4   6 
#  123   1  30   3  69   3 
View(hypStats_save_GO[hypStats_save_GO$FDR < 0.05,])
