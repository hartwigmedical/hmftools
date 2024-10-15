options(stringsAsFactors=F)

#========= Paths =========#
base_dir <- list(
   hpc='/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/',
   mnt='/Users/lnguyen/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/',
   local='/Users/lnguyen/Documents/P0013_WGS_patterns_Diagn/'
)

for(i in base_dir){
   if(dir.exists(i)){
      base_dir <- i
      break
   }
}

library(randomForest)
devtools::load_all(paste0(base_dir,'/CHORD/processed/scripts_main/mltoolkit/'))

#========= Load training data =========#
#training_data <- read.delim('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORDv2/training/training_data/mTrain_expandedIndels.txt.gz')

####################################################################################################
# Core functions
####################################################################################################
#' Default randomForest() settings for training CHORD
#'
#' @param df A dataframe containing the features and a response (class label) column
#' @param colname.response Name of the response column
#' @param ... Other arguments that can be passed to randomForest()
#'
#' @return A random forest object
#' @export
#'
randomForestDefault <- function(df, colname.response='response', seed=NULL, ...){
   randomForest(
      x = df[,colnames(df) != colname.response],
      y = as.factor(df[,colname.response]),
      strata = y, importance = T, proximity = F, ...
   )
}

#' Cross-validation wrapper for randomForest.chord()
#'
#' @param df A dataframe containing the features and a response (class label) column
#' @param colname.response Name of the response column
#' @param cv.folds Number of cross-validation folds
#' @param pred.only If TRUE, only returns the predictions from the test sets
#' @param train.function A custom training function that returns a random forest object. Used to
#' nest functions for nested cross-validation
#' @param mc Multicore processing? Will invoke foreach()
#' @param seed An integer specifying the random seed
#' @param verbose Print messages?
#' @param ... Other arguments that can be passed to randomForest.chord()
#'
#' @return A list of random forests and/or the predictions
#' @export
#'
randomForestCv <- function(
   df, colname.response='response',
   cv.folds=10, pred.only=F,
   train.function=NULL, mc=FALSE, seed=1, verbose=T, ...
){

   # df
   # colname.response='response'
   # cv.folds=10
   # pred.only=F
   # train.function=NULL
   # n.cores=1
   # verbose=T

   if(verbose){
      message('\nMaking CV train test sets...')
   }
   
   set.seed(seed)
   folds <- createCvTrainTestSets(df, k=cv.folds, stratify.by.col=colname.response)

   trainTestOneFold <- function(fold, fold.seed, ...){
      set.seed(fold.seed)
      if(is.null(train.function)){
         rf <- randomForestDefault(df=fold$train, ...)
         #rf <- randomForestDefault(df=fold$train)
      } else {
         rf <- train.function(df=fold$train, ...)
      }

      fold$test[is.na(fold$test)] <- 0
      pred <- as.data.frame(predict(
         object = rf,
         newdata = fold$test[,colnames(fold$test) != colname.response],
         type = "prob"
      ))
      pred$response <- fold$test$response

      if(pred.only){
         return(pred)
      } else {
         rf$pred <- pred
         return(rf)
      }
   }

   message('Performing CV...')
   if(!mc){
      lapply(1:length(folds), function(i){
         if(verbose){ message('> Fold: ',i) }
         fold <- folds[[i]]
         fold_seed <- as.integer(paste0(seed,i))
         rf <- trainTestOneFold(fold, fold.seed=fold_seed)
         rf$seed <- fold_seed
         return(rf)
      })
   }

   else {
      library(foreach)
      cl <- parallel::makeForkCluster(detectCores())
      doParallel::registerDoParallel(cl)
      foreach(i=1:length(folds)) %do% {
         if(verbose){ message('> Fold: ',i) }
         fold <- folds[[i]]
         fold_seed <- as.integer(paste0(seed,i))
         rf <- trainTestOneFold(fold, fold.seed=fold_seed)
         rf$seed <- fold_seed
         return(rf)
      }
      parallel::stopCluster(cl)
   }
}

#' CHORD core training procedure
#'
#' @param df A dataframe containing the features and a response (class label) column
#' @param colname.response Name of the response column
#' @param neg.class A character vector specifying the negative class(es)
#' @param pos.class A character vector specifying the positive class(es)
#' @param resample.ratios A list indicating resampling ratios for each class. A grid search is then
#' performed to determine the optimal set of resampling ratios.
#' @param wilcox.max.p For feature selection with wilcox test. Features with p-values higher than
#' this will be removed
#' @param verbose Show messages?
#'
#' @export
#'
randomForestTrain <- function(
   df, colname.response='response',
   neg.class='none', pos.class=c('BRCA2','BRCA1'),
   resample.ratios=list(BRCA2=c(1, 1, 1),BRCA1=c(1, 1.5, 2),none=c(1, 0.5, 0.25)),
   wilcox.max.p=0.01, seed=1, verbose=T
){
   #--------- Init ---------#
   # df=training_data
   # colname.response='response'
   # neg.class='none'
   # wilcox.max.p=1e-10
   # verbose=T
   # pos.class=c('BRCA2','BRCA1')

   x <- as.matrix(df[,colnames(df)!=colname.response])
   y <- df[,colname.response]

   if(!is.factor(y)){ stop('Response column must be a factor') }
   #df <- cbind(as.data.frame(x),response=y)

   if(is.null(pos.class)){
      pos.class <- levels(y)[ !(levels(y) %in% neg.class) ]
   }

   feature_names <- colnames(x)

   #--------- Univariate feature selection ---------#
   if(verbose){ message('Performing feature selection by wilcox test...') }
   ## Remove features that:
   ## - decrease in BRCA1/BRCA2 vs none class
   ## - have zero change
   x_split <- split(as.data.frame(x), y)
   
   wilcox_p <- lapply(pos.class, function(i){
      sapply(feature_names, function(j){
         wilcox.test(
            x_split[[i]][,j], x_split[[neg.class]][,j],
            alternative='greater'
         )$p.value
      })
   })
   names(wilcox_p) <- pos.class
   
   #wilcox_p <- as.data.frame(wilcox_p)
   #wilcox_p[order(wilcox_p$BRCA1,wilcox_p$BRCA1),]
   
   feature_whitelist <- (function(){
      v <- unlist(unname(wilcox_p))
      unique( names(v)[v<wilcox.max.p] )
   })()
   
   ## Select only relevant features
   df <- cbind(
      as.data.frame(x[,feature_whitelist]),
      response=y
   )

   #--------- Class balancing ---------#
   if(verbose){ message('Determining best class resampling ratios by CV...') }
   
   set.seed(seed)
   
   resampleClasses <- function(y, ratios=c(BRCA2=1, BRCA1=2, none=0.5)){
      #ratios <- c(BRCA2=1, BRCA1=2, none=0.5)

      tab <- table(y)
      ratios <- ratios[names(tab)]
      target_n <- round(tab*ratios)

      y_indexed <- structure(1:length(y),names=as.character(y))
      y_indexed_split <- split(y_indexed, y)

      y_resampled <- sapply(names(y_indexed_split), function(i){
         #i='none'
         indexes <- y_indexed_split[[i]]
         ratio <- ratios[[i]]
         n_samples <- target_n[[i]]

         if(ratio==1){ return(indexes) }

         sample(
            indexes, n_samples,
            replace=if(ratio>1){ TRUE } else { FALSE }
         )
      })

      sort(unlist(y_resampled, use.names=F))
   }

   ## Generate resampled datasets and run CV
   search_grid <- unique(expand.grid(resample.ratios))
   
   search_perf <- unlist(lapply(1:nrow(search_grid),function(i){
      #i=1
      resample_ratios <- unlist(search_grid[i,])

      if(verbose){
         message(
            '[',i,'/',nrow(search_grid),']: ',
            paste0(names(resample_ratios),collapse='/'),': ',
            paste0(resample_ratios,collapse='/')
         )
      }

      df_rs <- df[resampleClasses(df$response, resample_ratios),]
      cv_preds <- randomForestCv(df_rs, colname.response, pred.only=T,verbose=F)
      cv_preds_merged <- do.call(rbind, cv_preds)

      ## Calc AUPRC
      confusion <- 
         confusionMatrix(
            cv_preds_merged$BRCA1 + cv_preds_merged$BRCA2,
            ifelse(cv_preds_merged$response!=neg.class,TRUE,FALSE)
         )

      pr <- calcPerfCompound(confusion,'pr', metric.names.as.x.y=T)
      calcAUC(pr[,'x'],pr[,'y'])
   }))

   best_resample_ratios <- unlist(search_grid[which.max(search_perf),])

   #--------- Train final model ---------#
   if(verbose){ message('Training final model...') }
   set.seed(seed)
   df_rs <- df[resampleClasses(df$response, best_resample_ratios),]
   rf <- randomForestDefault(df_rs, colname.response)
   rf$feature_selection <- list(wilcox_pvalues=wilcox_p, whitelist=feature_whitelist, wilcox_max_p=wilcox.max.p)
   rf$resample_ratio <- best_resample_ratios
   
   return(rf)
}

# #========= Debugging =========#
# chordTrain <- function(){
#    pos.class=c('BRCA2','BRCA1')
#    neg.class='none'
#    sample.sel.cv.iters=10
#    verbose=T
#    
#    #--------- Do repeated CV to determine prediction stability ---------#
#    cl <- makeForkCluster(detectCores())
#    doParallel::registerDoParallel(cl)
#    
#    sample_sel_cv <- foreach(i=1:sample.sel.cv.iters) %do% {
#       if(verbose){ message('\n## Repeat: ',i) }
#       inner_cv <- randomForestCv(training_data, train.function=randomForestTrain, seed=i)
#       #saveRDS(inner_cv, paste0(base_dir,'/CHORDv2/training/scripts/train_rf/inner_cv_example.rds'))
#       do.call(rbind,lapply(inner_cv,`[[`,'pred'))
#    }
#    parallel::stopCluster(cl)
#    saveRDS(sample_sel_cv, paste0(base_dir,'/CHORDv2/training/scripts/train_rf/sample_sel_cv_example.rds'))
#    
#    ## Debugging
#    pred <- do.call(rbind,lapply(inner_cv,`[[`,'pred'))
#    sample_sel_cv <- lapply(1:10, function(i){ pred })
#    
#    lapply(sample_sel_cv, function(i){
#       #i=sample_sel_cv[[1]]
#    })
#    
# }
# 
# confusion <- with(cv_preds,{
#    confusionMatrix(
#       BRCA1 + BRCA2,
#       toBinaryResponse(response, pos.class, 1, neg.class, 0)
#    )
# })
# 
# pr <- calcPerfCompound(confusion,'pr', metric.names.as.x.y=T)
# calcAUC(pr[,'x'],pr[,'y'])
# 
# imp <- as.data.frame(importance(inner_cv[[1]]))
# imp[order(imp$MeanDecreaseAccuracy),]

####################################################################################################
# After cross-validation to select samples
####################################################################################################
#' Filter samples based on repeated CV predictions
#'
#' @param preds.dir Path to directory with prediction txt files
#' @param cutoff Classification threshold to consider a sample in the positive class
#' @param pos.class.min.freq Minimum proportion of times a positive class sample must be predicted
#' to be in the positive class to be in the whitelist
#' @param neg.class.max.freq Maximum proportion of times a negative class sample must be predicted
#' to be in the negative class to be in whitelist
#' @param pos.class A character vector specifying the names of the positive class(es)
#' @param neg.class A character vector specifying the names of the negative class(es)
#' @param verbose 
#'
#' @return A dataframe indicated which samples are in the whitelist
#' @export
#'
whitelistSamplesFromCvPreds <- function(
   preds.dir,
   cutoff=0.5, pos.class.min.freq=0.6, neg.class.max.freq=0.4,
   pos.class=c('BRCA1','BRCA2'), neg.class='none',
   verbose=T
){
   # preds.dir=paste0(base_dir,'/CHORDv2/training/models/2.00_expandedIndels/cv_sel_samples/')
   # pos.class=c('BRCA1','BRCA2')
   # neg.class='none'
   # cutoff=0.5
   # pos.class.min.freq=0.6
   # neg.class.max.freq=0.4
   
   ##
   if(verbose){ message('Reading CV preds...') }
   pred_files <- list.files(preds.dir, pattern='pred.*txt.gz',full.names=T)
   preds <- lapply(pred_files, function(i){ 
      df <- read.delim(i,stringsAsFactors=T)
      df <- df[order(rownames(df)),] ## Ensure that samples rows align
      return(df)
   })
   
   ##
   if(verbose){ message('Calculating number of times sample is above cutoff across CV repeats...') }
   m_preds <- do.call(cbind, lapply(preds, function(i){ rowSums(i[pos.class]) }))
   
   cv_summary <- data.frame(
      sample=rownames(m_preds),
      n_times_pos=apply(m_preds,1,function(i){ sum(i>=cutoff) }),
      response=preds[[1]]$response,
      row.names=NULL
   )
   
   if(is.null(pos.class)){
      pos.class <- with(cv_summary,{
         levels(response)[ !(levels(response) %in% neg.class) ]
      })
   }

   ##
   if(verbose){ message('Determining whitelist samples...') }
   ## Calculate absolute number of times for filtering neg/pos class
   n_cv_repeats <- length(preds)
   pos_class_min_times <- pos.class.min.freq * n_cv_repeats
   neg_class_max_times <- neg.class.max.freq * n_cv_repeats
   
   cv_summary_split <- split(cv_summary, cv_summary$response %in% neg.class)
   names(cv_summary_split) <- c('pos_samples','neg_samples')

   cv_summary_split <- within(cv_summary_split,{
      pos_samples$whitelist <- pos_samples$n_times_pos >= pos_class_min_times
      neg_samples$whitelist <- neg_samples$n_times_pos < neg_class_max_times
   })

   out <- do.call(rbind, cv_summary_split); rownames(out) <- NULL
   ##out[out$whitelist==F,]
   return(out)
}


####################################################################################################
# Outer CV
####################################################################################################
spawnOuterCvData <- function(df, out.dir, ...){
   #df=read.delim(paste0(base_dir,'/CHORDv2/training/scripts/train_rf/training_data.txt.gz'))
   #out.dir=paste0(base_dir,'/CHORDv2/training/models/2.00_expandedIndels/outer_cv/')
   
   train_test_sets <- createCvTrainTestSets(df, ...)
   
   sub_dir_names <- formatC(
      1:length(train_test_sets), 
      width=nchar(length(train_test_sets)), 
      format="d", flag="0"
   )
   
   sub_dir_paths <- paste0(out.dir,'/',sub_dir_names)
   
   for(i in 1:length(train_test_sets)){
      dir.create(sub_dir_paths[[i]], showWarnings=F, recursive=T)
      
      write.table(
         train_test_sets[[i]]$train,
         gzfile(paste0(sub_dir_paths[[i]],'/training_data.txt.gz')),
         sep='\t',quote=F
      )
      
      write.table(
         train_test_sets[[i]]$test,
         gzfile(paste0(sub_dir_paths[[i]],'/test_data.txt.gz')),
         sep='\t',quote=F
      )
   }
}















