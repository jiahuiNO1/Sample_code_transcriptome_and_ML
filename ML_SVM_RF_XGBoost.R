#SVM, RF, XGBoost
#use gene number up tp 200
#gene number the same as risk scoring

require(data.table)
require(dplyr)
require(ggplot2)
require(caret)
require(pROC)
require(randomForest)
require(HGNChelper)
require(genefilter)
require(cytominer)
require(exprso)
require(xgboost)
require(plyr)
require(tidyquant)
require(Rcpp)

#Function to calculate metrics
calc_metrics_roc = function(prob_case = NULL, actual = NULL, plotSkip = F) {
  if(is.null(prob_case) == T) stop("Probability for caseness required!");
  if(is.null(actual) == T) stop("Must supply actual class labels (binary)");
  # Find optimal cutoff based on distance from top-left corner
  p <- ROCR::prediction(prob_case, actual)
  perf <- ROCR::performance(p, measure = "tpr", x.measure = "fpr")
  index <- which.min(sqrt((1 - perf@y.values[[1]])^2 + (0 - perf@x.values[[1]])^2))
  if(!plotSkip) plot(perf, col = rainbow(10))
  if(!plotSkip) graphics::points(perf@x.values[[1]][index], perf@y.values[[1]][index], col = "blue")
  # Get performance for optimal cutoff
  acc <- ROCR::performance(p, "acc")@y.values[[1]][index]
  sens <- ROCR::performance(p, "sens")@y.values[[1]][index]
  spec <- ROCR::performance(p, "spec")@y.values[[1]][index]
  prec <- ROCR::performance(p, "prec")@y.values[[1]][index]
  f1 <- 2 * (prec * sens) / (prec + sens) #Precision-recall F measur
  auc <- ROCR::performance(p, "auc")@y.values[[1]]
  return(data.frame(acc, sens, spec, prec, f1, auc))
}


##--> Load combined and scaled blood gene-expression data
#sample in row
datExpr = fread("D:/AD_gx/Blood_gx/WGCNA/20220926/20220926_AD_Blood_WGCNA_SamplesGenes_Combat_datExpr.txt")
datExpr = data.frame(datExpr)
row.names(datExpr) = datExpr$V1
datExpr = datExpr[,-1]

phenos = fread("D:/AD_gx/Blood_gx/WGCNA/20220926/20220926_AD_Blood_WGCNA_SamplesPhenos_of_datExpr.txt")
phenos = data.frame(phenos)
row.names(phenos) = phenos$V1
phenos = phenos[,-1]

table(phenos$FACTOR_studyID, phenos$FACTOR_dx)            
#               AD CTL
#  adni        114 245
#  Lovestone_1  94 106
#  Lovestone_2 117 130
#  ROSMAP_b1b2  55  87
#  Samsudin16   90  90
table(phenos$FACTOR_studyID)
#        adni Lovestone_1 Lovestone_2 ROSMAP_b1b2  Samsudin16 
#        359         200         247         142         180 


class_pairs = data.frame((unique(phenos$FACTOR_dx)))
top_k_genes = c(2, 5,10, 15, 25, 50, 100, 125, 150, 200)

#use Lovestone_2 as held-out testing data
train_phe = phenos[!grepl("Lovestone_2", phenos$FACTOR_studyID), ]

for(col in 1:ncol(train_phe)){
  if(colnames(train_phe[,col, drop = FALSE]) %in% c("FACTOR_dx")){
    train_phe[,col] = factor(train_phe[,col], levels = c("CTL","AD"))
    
  } else if(colnames(train_phe[,col, drop = FALSE]) %in% c("FACTOR_sex","FACTOR_ethnicity","FACTOR_Plate","FACTOR_batch","FACTOR_site","FACTOR_APOE","FACTOR_studyID")){
    train_phe[,col] = factor(train_phe[,col])
        

  } else {
    train_phe[,col] = as.numeric(train_phe[,col])

  }
}

#
valid_phe = phenos[grepl("Lovestone_2", phenos$FACTOR_studyID), ]

for(col in 1:ncol(valid_phe)){
  if(colnames(valid_phe[,col, drop = FALSE]) %in% c("FACTOR_dx")){
    valid_phe[,col] = factor(valid_phe[,col], levels = c("CTL","AD"))
    
  } else if(colnames(valid_phe[,col, drop = FALSE]) %in% c("FACTOR_sex","FACTOR_ethnicity","FACTOR_Plate","FACTOR_batch","FACTOR_site","FACTOR_APOE","FACTOR_studyID")){
    valid_phe[,col] = factor(valid_phe[,col])
        

  } else {
    valid_phe[,col] = as.numeric(valid_phe[,col])

  }
}

##--> Normalize expression data for known covariates that could influence classification performance
train_edat = datExpr[row.names(datExpr) %in% row.names(train_phe),]

#get gene-expression residuals after adjusting age, sex and estimated cell proportions
mydata = train_phe[, colnames(train_phe) %in% c("FACTOR_age","FACTOR_sex","FACTOR_studyID","FACTOR_B","FACTOR_T","FACTOR_Dendritic","FACTOR_Monocytes", "FACTOR_Mast","FACTOR_NK","FACTOR_Granulocytes") ]
train_edat_adj = resid(lm(as.matrix(train_edat) ~ ., data = mydata))

##--> Adjust testing data based on the scaling of training data
means = colMeans(train_edat_adj)
vars = apply(train_edat_adj, 2, sd)

valid_edat = datExpr[row.names(datExpr) %in% row.names(valid_phe),]
norm = sweep(valid_edat, 2, STATS = means, FUN = `-`)
valid_edat_res = sweep(norm, 2, STATS = vars, FUN=`/`)

train_edat_res = train_edat_adj

##--> Find features that are associated with AD diagnosis in the training data
ttest = genefilter::colttests(x = as.matrix(train_edat_res), fac = train_phe$FACTOR_dx)
ttest = ttest[order(ttest$p.value, decreasing = F), ]
ttest = ttest[ttest$p.value < 0.05, ]

ranked_features = data.frame(ttest, gene = rownames(ttest), classes = paste(class_pairs[,1], sep="", collapse="|"))
ranked_features$rank = 1:nrow(ranked_features)

#
Cfolder = "D:/AD_gx/Blood_gx/MLA/20221030"
if(!dir.exists(Cfolder)){
  dir.create(Cfolder)
}
setwd(Cfolder)


fwrite(ranked_features, file = "20221030_AD_Blood_MLs_ranked_features.txt", row.names = T, col.names = T, sep = "\t")


#use top ranked genes up to 200
ranked_features = ranked_features[1:max(top_k_genes), ]


##--> Model training
k_perfs = list()

for(k in 1:length(top_k_genes)){
  
  # select k features
  features = rownames(ttest)[1:top_k_genes[[k]]]
  
  ##--> Train RF model
  n_trees = c(100,200,500,1000)
  xinput = exprso(train_edat_res[,colnames(train_edat_res) %in% features], y= train_phe$FACTOR_dx)
  tune_rf = plGrid(array.train = xinput, how = 'buildRF', trees=n_trees, aucSkip = F, plCV.acc = 'auc', fold = 5)
  
  # proc 
  model = tune_rf@machs[[which.max(tune_rf@summary$train.plCV)]]@mach
  perf_train = predict(model, type='prob')
  auc_train = pROC::roc(train_phe$FACTOR_dx, perf_train[,2], direction="<")
  ci_train = ci.auc(auc_train)
  
  ##--> Test RF model 
  deploy = predict(model, valid_edat_res[,colnames(valid_edat_res) %in% features], type='prob')
  auc_valid = pROC::roc(valid_phe$FACTOR_dx, deploy[,2], direction="<")
  ci_valid = ci.auc(auc_valid)
  
  train_metrics = calc_metrics_roc(prob_case = perf_train[,2], actual = train_phe$FACTOR_dx)
  valid_metrics = calc_metrics_roc(prob_case = deploy[,2], actual = valid_phe$FACTOR_dx)
  
  ##--> Extract performances metrics
  rf_k_perfs = data.frame(k = top_k_genes[[k]],
                          train.auc = auc_train$auc,
                          train.lower = ci_train[[1]],
                          train.upper = ci_train[[3]],
                          train.sens = train_metrics$sens,
                          train.spec = train_metrics$spec,
                          train.f1 = train_metrics$f1,
                          valid.auc = auc_valid$auc,
                          valid.lower = ci_valid[[1]],
                          valid.upper = ci_valid[[3]],
                          valid.sens = valid_metrics$sens,
                          valid.spec = valid_metrics$spec,
                          valid.f1 = valid_metrics$f1)
  rf_k_perfs = data.frame(build = 'rf', rf_k_perfs)
  
  
  ##--> Train SVM model
  costs = 10^-3:3
  tune_svm = plGrid(array.train = xinput, how = 'buildSVM', kernel='linear', cost = costs, aucSkip = F, plCV.acc = 'auc', fold = 5)
  besttune = which.max(tune_svm@summary$train.plCV)

  # proc
  model = tune_rf@machs[[besttune]]@mach
  perf_train = predict(model, type='prob')
  auc_train = pROC::roc(train_phe$FACTOR_dx, perf_train[,2], direction="<")
  ci_train = ci.auc(auc_train)

  ##--> Test SVM model 
  deploy = predict(model, valid_edat_res[,colnames(valid_edat_res) %in% features], type='prob')
  auc_valid = pROC::roc(valid_phe$FACTOR_dx, deploy[,2], direction="<")
  ci_valid = ci.auc(auc_valid)

  train_metrics = calc_metrics_roc(prob_case = perf_train[,2], actual = train_phe$FACTOR_dx)
  valid_metrics = calc_metrics_roc(prob_case = deploy[,2], actual = valid_phe$FACTOR_dx)

  svm_k_perfs = data.frame(k = top_k_genes[[k]],
                          train.auc = auc_train$auc,
                          train.lower = ci_train[[1]],
                          train.upper = ci_train[[3]],
                          train.sens = train_metrics$sens,
                          train.spec = train_metrics$spec,
                          train.f1 = train_metrics$f1,
                          valid.auc = auc_valid$auc,
                          valid.lower = ci_valid[[1]],
                          valid.upper = ci_valid[[3]],
                          valid.sens = valid_metrics$sens,
                          valid.sepc = valid_metrics$spec,
                          valid.f1 = valid_metrics$f1)
  svm_k_perfs = data.frame(build = 'svmlinear', svm_k_perfs)
  
  
  ##--> Train XGBoost model
  params = expand.grid(eta = seq(0.01, 0.3, length.out = 5), nrounds = c(5, 10, 25))
  input_x = train_edat_res[,colnames(train_edat_res) %in% features]
  inputDat = xgboost::xgb.DMatrix(input_x, info = list(label = as.numeric(train_phe$FACTOR_dx)-1))
  save_res = list()
  for(p in 1:nrow(params)){
    xgb_perf = xgboost::xgb.cv(data = inputDat, nrounds = params$nrounds[p], metrics='logloss',
                               objective='binary:logistic', alpha=0.5,
                               nfold = 5, eta = params$eta[p], subsample = 0.5)
    res = data.frame(xgb_perf$evaluation_log)
    res = res[,colnames(res) %in% c("iter", "train_logloss_mean", "test_logloss_mean")]
    res$iter = as.factor(res$iter)
    res$eta = params$eta[p]
    res$nrounds = params$nrounds[p]
    res = res[nrow(res), ]
    save_res[[p]] = res
  }
  save_res_xgb = ldply(save_res)
  
  bestTune = save_res_xgb[which.min(save_res_xgb$test_logloss_mean),]
  
  model = xgboost::xgboost(data = inputDat, nrounds = bestTune$nrounds[1], 
                           alpha=0.5, metrics='logloss', eta = bestTune$eta[1], 
                           objective='binary:logistic',subsample = 0.5)
  
  
  # proc 
  perf_train = predict(model, input_x, type='prob')
  auc_train = pROC::roc(train_phe$FACTOR_dx, perf_train, direction="<")
  ci_train = ci.auc(auc_train)
  
  ##--> Test XGBoost model
  valid_input = as.matrix(valid_edat_res[,colnames(valid_edat_res) %in% features])
  deploy = predict(model, valid_input, type='prob')
  auc_valid = pROC::roc(as.numeric(valid_phe$FACTOR_dx)-1, deploy, direction = "<")
  ci_valid = ci.auc(auc_valid)
  
  train_metrics = calc_metrics_roc(prob_case = perf_train, actual = train_phe$FACTOR_dx)
  valid_metrics = calc_metrics_roc(prob_case = deploy, actual = valid_phe$FACTOR_dx)
  
  xgb_k_perfs = data.frame(k = top_k_genes[[k]],
                           train.auc = auc_train$auc,
                           train.lower = ci_train[[1]],
                           train.upper = ci_train[[3]],
                           train.sens = train_metrics$sens,
                           train.spec = train_metrics$spec,
                           train.f1 = train_metrics$f1,
                           valid.auc = auc_valid$auc,
                           valid.lower = ci_valid[[1]],
                           valid.upper = ci_valid[[3]],
                           valid.sens = valid_metrics$sens,
                           valid.spec = valid_metrics$spec,
                           valid.f1 = valid_metrics$f1)
  xgb_k_perfs = data.frame(build = 'xgboost', xgb_k_perfs)
  
  # save performances for all models
  k_perfs[[k]] = ldply(list(xgb_k_perfs, rf_k_perfs, svm_k_perfs))
}

# -- concat performance metrics
class_perf = ldply(k_perfs)
class_perf$classes = paste(class_pairs[,1], collapse="|", sep="")
class_perf$n_group_1 = table(train_phe$FACTOR_dx)[[1]]
class_perf$n_group_2 = table(train_phe$FACTOR_dx)[[2]]

# -- concat ranked features for interpretability
all_fs = ranked_features
cast_fs = reshape2::acast(all_fs, gene ~ classes, value.var='rank')


# performance
all_perfs = class_perf
sum_perfs = data.table(all_perfs)


sum_perfs$valid.se = (sum_perfs$valid.upper - sum_perfs$valid.lower)/3.92

##--> One tailed z-test of significance 
sum_perfs$zscore = (sum_perfs$valid.auc - 0.5)/(sum_perfs$valid.se)
sum_perfs$pval = pnorm(sum_perfs$zscore, lower.tail = F)
sum_perfs$build = toupper(sum_perfs$build)
sum_perfs$build = gsub("XGBOOST", "XGBoost", sum_perfs$build)
sum_perfs$build = gsub("SVMLINEAR", "SVM-Linear", sum_perfs$build)

fwrite(sum_perfs, file = "20221030_Summary_Of_MLA_Performance_AD_blood.txt", row.names = T, col.names = T, sep = "\t")
