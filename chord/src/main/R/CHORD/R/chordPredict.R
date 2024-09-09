#' Predict the probability of homogolous recombination deficiency using mutational signatures
#'
#' @description A wrapper for predict.randomForest() from the randomForest package
#'
#' @param features The output of extractSigsChord(), which is a dataframe containing the SNV, indel
#' and SV context counts.
#' @param rf.model The random forest model. Defaults to CHORD.
#' @param hrd.cutoff Default=0.5. Samples greater or equal to this cutoff will be marked as HRD 
#' (is_hrd==TRUE).
#' @param trans.func Function used to transform raw features. Raw features should be in the format:
#' list(snv=matrix(), indel=matrix(), sv=matrix())
#' @param min.indel.load Default=50. The minimum number of indels required to make an accurate HRD
#' prediction. Samples with fewer indels than this value will be marked as is_hrd==NA (HR status 
#' could not be confidently determined).
#' @param min.sv.load Default=30. The minimum number of SVs required to make an accurate prediction
#' of BRCA1-type vs. BRCA2-type HRD. Samples with fewer SVs than this value will be marked as 
#' hrd_type==NA (HRD type could not be confidently determined).
#' @param min.msi.indel.rep Default=14000 (changing this value is not advised). Samples with more 
#' indels within repeats than this threshold will be considered to have microsatellite instability.
#' @param do.bootstrap Test the stability of prediction probabilities? NOTE: this is computationally
#' expensive. Resamples the feature vector for each sample (number of times provided by
#' bootstrap.iters) and calculates HRD probabilities for each iteration. Returns the probabilities
#' at the quantiles specifying in bootstrap.quantiles
#' @param bootstrap.iters Number of resampling iterations for determining the confidence intervals
#' @param bootstrap.quantiles A numeric vector of length 2 specifying the quantiles used to calculate the 
#' confidence intervals
#' @param detailed.remarks If TRUE, shows min.indel.load and min.sv.load numbers in the remarks columns
#' @param show.features If TRUE, appends features to output
#' @param verbose Show messages?
#'
#' @return A dataframe containing the HRD probabilities, bootstrap probabilities, and input features
#' @export
#' 
#' @examples
#' ## Extract mutation contexts
#' vcf_dir <- '/path_to_vcfs/'
#' vcf_snv <- paste0(vcf_dir,'SampleX_post_processed_v2.2.vcf.gz')
#' vcf_indel <- paste0(vcf_dir,'SampleX_post_processed_v2.2.vcf.gz')
#' vcf_sv <- paste0(vcf_dir,'SampleX_somaticSV_bpi.vcf.gz')
#' contexts <- extractSigsChord(vcf_snv, vcf_indel, vcf_sv, sample.name='SampleX')
#' 
#' ## Predict HRD probability with CHORD
#' chordPredict(contexts)

chordPredict <- function(
  features, rf.model=CHORD, hrd.cutoff=0.5, 
  trans.func=NULL,
  
  ## QC thresholds
  min.indel.load=100, min.sv.load=30, min.msi.indel.rep=14000,
  
  ## Confidence interval estimation
  do.bootstrap=F, bootstrap.iters=20, bootstrap.quantiles=c(0.05, 0.5, 0.95),
  
  ## Other
  detailed.remarks=T, show.features=F, verbose=T
){
  
  # features=read.delim('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/scripts_main/CHORD/example/output/merged_contexts.txt', check.names=F)
  # features=readRDS('/Users/lnguyen/hpc/cog_bioinf/cuppen/project_data/Luan_projects/CHORD/HMF_DR010_DR047/matrices/merged_contexts.rds')
  # features=do.call(cbind, unname(features))
  # features=features[1:20,]
  
  if(verbose){ message('Preprocessing features...') }
  ## Converts the raw signature counts from extractSigsChord() to features used by CHORD
  if(is.list(features) & !is.data.frame(features)){
    features_split <- features
  } else {
    features_split <- mutSigExtractor::splitDfRegex(features, c(snv='>',indel='[a-z]{3}[.]',sv='[A-Z]{3}'))
  }
  #features_split <- features_split[c('snv','indel','sv')] ## Remove other mut types if they exist
  
  ## Default feature transformation function
  if(is.null(trans.func)){
    
    trans.func <- function(features_split){
      ## Merge del mh bimh 2-5
      indels <- features_split$indel
      indel_types <- c('del.rep','ins.rep','del.mh','ins.mh','del.none','ins.none')
      indels_split <- mutSigExtractor::splitDfRegex(indels, indel_types)
      names(indels_split) <- indel_types
      
      counter <- 0
      indels_custom <- lapply(indels_split, function(i){
        counter <<- counter + 1
        df <- as.data.frame(rowSums(i))
        colnames(df) <- indel_types[counter]
        return(df)
      })
      
      indels_custom$del.mh <- with(indels_split,{
        cbind(
          del.mh['del.mh.bimh.1'],
          'del.mh.bimh.2.5'=rowSums(del.mh[!(colnames(del.mh) %in% 'del.mh.bimh.1')])
        )
      })
      
      ## Add new indels back to list
      l <- features_split[c('snv','indel','sv')]
      l$indel <- do.call(cbind, unname(indels_custom))
      
      ## Simplify SNVs to the 6 types of base substitions (C>A, C>G, C>T, T>A, T>C, T>G) by ignoring
      ## the immediate flanking nucleotides
      ## Convert absolute counts to relative counts. Calculated per variant type.
      m <- mutSigExtractor::transformContexts(l, simplify.types='snv', rel.types='all')
      return(m)
    }
  }
  
  features_processed <- trans.func(features_split)
  
  #--------- Prediction ---------#
  if(verbose){ message('Calculating HRD probabilities...') }
  doPredict <- function(m){
    df <- as.data.frame(
      randomForest:::predict.randomForest(rf.model, m, type='prob'),
      stringsAsFactors=F
    )
    df <- df[,c('none','BRCA1','BRCA2')]
    colnames(df) <- paste0('p_',colnames(df))
    df$p_hrd <- df$p_BRCA1 + df$p_BRCA2
    #df <- cbind(sample=rownames(df), df); rownames(df) <- NULL
    
    return(df)
  }
  
  df <- doPredict(features_processed)
  
  #--------- Bootstrap predictions ---------#
  if(do.bootstrap){
    if(verbose){ message('Calculating bootstrap confidence intervals...') }
    
    resampleFeatureMatrix <- function(m, repeats=bootstrap.iters){
      #m=features_split$sv
      
      feature_ids <- as.factor(1:ncol(m)) ## Replace feature names by integers to save memory
      
      resampleFeatureVector <- function(v){
        #v=unlist(m[229,])
        total_counts <- sum(v)
        if(total_counts==0){ ## Prevent divide by 0 errors
          return(
            matrix(0, nrow=repeats, ncol=ncol(m))
          )
        }
        
        prob <- v/total_counts
        t(replicate(
          repeats,
          table( sample(feature_ids, total_counts, replace=T, prob=prob) )
        ))
      }
      
      if(verbose){ pb <- txtProgressBar(max=nrow(m), style=3) }
      out <- lapply(1:nrow(m), function(i){
        #i=1
        #print(i)
        if(verbose){ setTxtProgressBar(pb, i) }
        m_resampled <- resampleFeatureVector(m[i,])
        colnames(m_resampled) <- colnames(m)
        
        return(m_resampled)
      })
      message('\n')
      names(out) <- rownames(m)
      
      return(out)
    }
    
    ## list structure: mut type -> sample name
    features_split_resampled <- lapply(names(features_split), function(i){
      if(verbose){ message('> Resampling ',i,' matrix...') }
      resampleFeatureMatrix(features_split[[i]])
    })
    names(features_split_resampled) <- names(features_split)
    
    if(verbose){ message('> Performing bootstrap predictions...') }
    bootstrap_pred <- lapply(1:nrow(features), function(i){
      #i=229
      #print(i)
      ## convert to: sample name --> mut type
      l <- list(
        snv=features_split_resampled$snv[[i]],
        indel=features_split_resampled$indel[[i]],
        sv=features_split_resampled$sv[[i]]
      )
      
      ## Transform features
      m <- trans.func(l)
      
      ## Predict
      pred <- doPredict(m)
      
      ## Calculate quantiles
      unlist(lapply(pred,function(pred.class){ ## Unlist automatically prepends pred class names
        quantile(pred.class, bootstrap.quantiles)
      }))
    })
    bootstrap_pred <- do.call(rbind, bootstrap_pred)
    rownames(bootstrap_pred) <- rownames(features)
  }
  
  #--------- QC ---------#
  if(verbose){ message('Performing QC checks...') }
  qc <- list()
  qc$has_msi <- with(features_split,{
    rowSums(indel[,grep('rep',colnames(indel)),drop=F]) > min.msi.indel.rep
  })
  
  qc$low_indel_load <- rowSums(features_split$indel) < min.indel.load
  qc$low_sv_load <- rowSums(features_split$sv) < min.sv.load
  qc <- as.data.frame(qc)
  
  failed_qc <- with(qc,{ has_msi | low_indel_load | low_sv_load })
  
  ## Informative QC tags
  if(!detailed.remarks){
    failed_qc_strings <- colnames(qc)
    names(failed_qc_strings) <- colnames(qc)
  } else {
    ## Note to self: make sure order of qc names is correct upon editing
    failed_qc_strings <- c(
      has_msi=paste0('Has MSI (>',min.msi.indel.rep,' indel.rep)'),
      low_indel_load=paste0('<',min.indel.load,' indels'),
      low_sv_load=paste0('<',min.sv.load,' SVs')
    )
  }
  
  ## is_hrd qc tags
  df_qc_is_hrd <- data.frame(
    has_msi=ifelse(qc$has_msi,failed_qc_strings['has_msi'],''),
    low_indel_load=ifelse(qc$low_indel_load,failed_qc_strings['low_indel_load'],'')
  )
  
  qc_is_hrd <- with(df_qc_is_hrd,{
    paste(has_msi, low_indel_load,sep=';')
  })
  qc_is_hrd <- gsub('^;|;$','',qc_is_hrd)
  qc_is_hrd[nchar(qc_is_hrd)==0] <- ''
  
  ## hrd_type qc tags
  qc_hrd_type <- ifelse(
    qc$low_sv_load,
    failed_qc_strings['low_sv_load'],
    ''
  )
  
  qc_out <- data.frame(
    remarks_hr_status=qc_is_hrd, 
    remarks_hrd_type=qc_hrd_type
  )
  
  #--------- Determine HR status ---------#
  if(verbose){ message('Finalizing output...') }
  
  ## Is HRD?
  df$hr_status <- ifelse(
    df$p_hrd >= hrd.cutoff,
    'HR_deficient','HR_proficient'
  )
  
  ## Which HRD subtype?
  df$hrd_type <- unlist(
    Map(function(p_BRCA1, p_BRCA2, p_hrd){
      if(p_hrd>=hrd.cutoff){
        c('BRCA1_type','BRCA2_type')[ which.max(c(p_BRCA1,p_BRCA2)) ]
      } else {
        'none'
      }
    }, df$p_BRCA1, df$p_BRCA2, df$p_hrd)
  )
  
  ## Add qc
  df <- cbind(df, qc_out)
  
  ## Clean up tags
  df$hr_status[ qc$low_indel_load | qc$has_msi ] <- 'cannot_be_determined'
  
  df$hrd_type[ qc$low_sv_load ] <- 'cannot_be_determined'
  df$hrd_type[ df$hr_status %in% c('HR_proficient','cannot_be_determined') ] <- 'none'
  df$remarks_hrd_type[ df$hr_status=='HR_proficient' ] <- ''
  
  #df <- df[,!grepl('^p_none',colnames(df))] ## Remove redunant 'none' class
  
  #--------- Gather outputs ---------#
  out <- cbind(sample=rownames(df), df)
  if(do.bootstrap){ out <- cbind(out, bootstrap_pred) }
  if(show.features){ out <- cbind(out, features_processed) }
  out <- out[,!grepl('^p_none',colnames(out))] ## Remove redunant 'none' class
  rownames(out) <- NULL

  return(out)
}




