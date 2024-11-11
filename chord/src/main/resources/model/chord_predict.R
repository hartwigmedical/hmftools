suppressPackageStartupMessages(library('randomForest'))

options(stringsAsFactors=F) # to avoid invalid factor level warning

## Command line args ================================
args <- commandArgs(TRUE)

CHORD_MODEL_PATH <- args[1]
MUTATION_CONTEXTS_PATH <- args[2]
OUTPUT_FILE_PATH <- args[3]
GLOBAL_LOG_LEVEL <- args[4]

## Logging ================================
GLOBAL_LOG_LEVEL <- toupper(GLOBAL_LOG_LEVEL)

LOG_LEVEL <- list(
   TRACE = list(ordinal = 1, name = "TRACE", display_name = "TRACE"),
   DEBUG = list(ordinal = 2, name = "DEBUG", display_name = "DEBUG"),
   INFO  = list(ordinal = 3, name = "INFO",  display_name = "INFO "),
   WARN  = list(ordinal = 4, name = "WARN",  display_name = "WARN "),
   ERROR = list(ordinal = 5, name = "ERROR", display_name = "ERROR"),
   FATAL = list(ordinal = 6, name = "FATAL", display_name = "FATAL")
)

logMessage <- function(log_level, string){
   
   current_time <- format(Sys.time(), "%H:%H:%OS3")
   
   log_message <- sprintf("%s [R] [%s] %s", current_time, log_level$display_name, string)
   
   if(log_level$ordinal >= LOG_LEVEL[[LOG_LEVEL$ERROR$name]]$ordinal)
      stop(log_message)
   
   if(log_level$ordinal >= LOG_LEVEL[[GLOBAL_LOG_LEVEL]]$ordinal)
      message(log_message)
}

## Parse args ================================
logMessage(LOG_LEVEL$INFO, sprintf("Starting CHORD predict"))

logMessage(LOG_LEVEL$INFO, sprintf("Loading CHORD model: %s", CHORD_MODEL_PATH))
CHORD_MODEL <- readRDS(CHORD_MODEL_PATH)

logMessage(LOG_LEVEL$INFO, sprintf("Loading mutation contexts: %s", MUTATION_CONTEXTS_PATH))
MUTATION_CONTEXTS <- read.delim(MUTATION_CONTEXTS_PATH, check.names = FALSE, row.names = 1)

## Transform features ================================
SNV_SUBSTITUTIONS <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")

INDEL_NEW_CONTEXTS <- c(
   del.rep         = "del.rep", 
   ins.rep         = "ins.rep", 
   del.mh.bimh.1   = "del.mh.bimh.1", 
   del.mh.bimh.2.5 = "del.mh.bimh.[2,3,4,5]", 
   ins.mh          = "ins.mh", 
   del.none        = "del.none", 
   ins.none        = "ins.none"
)

MUTATION_TYPE_REGEX <- c(
   snv =   "^\\w\\[\\w>\\w\\]\\w$",
   indel = "^[a-z]{3}[.]",
   sv =    '^[A-Z]{3}'
)

getColumns <- function(df, regex, fixed=FALSE){
   column_indexes <- grep(regex, colnames(df), fixed=fixed)
   df[,column_indexes, drop=FALSE]
}

groupColumns <- function(df, regex_list, fixed=FALSE){

   df_list <- lapply(regex_list, function(regex){
      getColumns(df, regex=regex, fixed=fixed)
   })
   
   names(df_list) <- names(regex_list)
   
   old_ncol <- ncol(df)
   new_ncol <- sum(sapply(df_list, ncol))
   if(old_ncol != new_ncol){
      logMessage(
         LOG_LEVEL$ERROR,
         sprintf("Output no. of columns (%s) != original no. of columns (%s)", new_ncol, old_ncol),
      )
   }
   
   return(df_list)
}

mergeSnvContexts <- function(snv_contexts){
   
   contexts_grouped <- groupColumns(df = snv_contexts, regex_list = SNV_SUBSTITUTIONS, fixed = TRUE)
   contexts_merged <- do.call(cbind, lapply(contexts_grouped, rowSums))
   colnames(contexts_merged) <- sub(">", ".", SNV_SUBSTITUTIONS)
   contexts_merged <- as.data.frame(contexts_merged)
   
   return(contexts_merged)
}

mergeIndelContexts <- function(indel_contexts){
   
   contexts_grouped <- groupColumns(df = indel_contexts, regex_list = INDEL_NEW_CONTEXTS)
   contexts_merged <- do.call(cbind, lapply(contexts_grouped, rowSums))
   colnames(contexts_merged) <- names(INDEL_NEW_CONTEXTS)
   contexts_merged <- as.data.frame(contexts_merged)
   
   return(contexts_merged)
}

mergeMutationContextsByType <- function(mutation_contexts){
   
   logMessage(LOG_LEVEL$INFO, "Transforming CHORD input features")
   
   contexts_by_type <- groupColumns(df = mutation_contexts, regex_list = MUTATION_TYPE_REGEX)
   
   logMessage(LOG_LEVEL$DEBUG, "Merging SNV contexts")
   contexts_by_type$snv <- mergeSnvContexts(contexts_by_type$snv)
   
   logMessage(LOG_LEVEL$DEBUG, "Merging INDEL contexts")
   contexts_by_type$indel <- mergeIndelContexts(contexts_by_type$indel)
   
   contexts_by_type$sv <- contexts_by_type$sv
   
   return(contexts_by_type)
}

MUTATION_CONTEXTS_BY_TYPE <- mergeMutationContextsByType(MUTATION_CONTEXTS)

## Predict ================================
HRD_PROBABILITY_CUTOFF <- 0.5

HR_STATUS <- list(
   DEFICIENT  = 'HR_deficient',
   PROFICIENT = 'HR_proficient',
   UNCLEAR =    'cannot_be_determined'
)

HRD_TYPE <- list(
   BRCA1 = "BRCA1_type",
   BRCA2 = "BRCA2_type",
   none  = "none" 
)

MSI_MIN_INDEL_REP <- 14000
LOW_INDEL_LOAD_THRES <- 100
LOW_SV_LOAD_THRES <- 30

chordPredict <- function(mutation_contexts_by_type){
   
   logMessage(LOG_LEVEL$INFO, "Calculating mutation context relative counts")
   context_rel_counts <- lapply(mutation_contexts_by_type, function(contexts){
      contexts <- contexts/rowSums(contexts)
      contexts[is.na(contexts)] <- 0
      return(contexts)
   })
   
   context_rel_counts <- do.call(cbind, unname(context_rel_counts))
   
   
   logMessage(LOG_LEVEL$INFO, "Predicting HRD probabilities from random forest")
   probs <- randomForest:::predict.randomForest(object=CHORD_MODEL, newdata=context_rel_counts, type='prob')
   probs <- data.frame(probs)
   
   
   logMessage(LOG_LEVEL$DEBUG, "Determining HR status")
   probs_brca_type <- probs[c("BRCA1","BRCA2")]
   prob_hrd <- rowSums(probs_brca_type)
   
   hr_status <- ifelse(prob_hrd >= HRD_PROBABILITY_CUTOFF, HR_STATUS$DEFICIENT, HR_STATUS$PROFICIENT)
   
   hrd_type <- colnames(probs_brca_type)[ max.col(probs_brca_type) ]
   hrd_type <- unname(unlist(HRD_TYPE[hrd_type]))
   
   predictions <- data.frame(
      sample = rownames(probs),
      p_BRCA1 = probs$BRCA1,
      p_BRCA2 = probs$BRCA2,
      p_hrd = prob_hrd,
      hr_status = hr_status,
      hrd_type = hrd_type
   )
   
   
   logMessage(LOG_LEVEL$INFO, "Performing QC checks")
   indel_rep_counts <- rowSums(cbind(
      mutation_contexts_by_type$indel$del.rep,
      mutation_contexts_by_type$indel$ins.rep
   ))
   
   indel_load <- rowSums(mutation_contexts_by_type$indel)
   sv_load <- rowSums(mutation_contexts_by_type$sv)
   
   has_msi <- indel_rep_counts > MSI_MIN_INDEL_REP
   has_low_indel_load <- indel_load < LOW_INDEL_LOAD_THRES
   has_low_sv_load <- sv_load < LOW_SV_LOAD_THRES
   
   logMessage(LOG_LEVEL$DEBUG, "Modifying predictions based on QC status")
   predictions$hr_status[ has_low_indel_load | has_msi ] <- HR_STATUS$UNCLEAR
   predictions$hrd_type[ predictions$hr_status %in% c(HR_STATUS$PROFICIENT, HR_STATUS$UNCLEAR) ] <- HRD_TYPE$none
   
   
   logMessage(LOG_LEVEL$DEBUG, "Adding QC remarks")
   remarks_hr_status <- paste(
      ifelse(has_msi, sprintf("Has MSI (>%s indel.rep)", MSI_MIN_INDEL_REP), ""),
      ifelse(has_low_indel_load, sprintf("<%s indels", LOW_INDEL_LOAD_THRES), ""),
      sep=";"
   )
   remarks_hr_status <- gsub("^;|;$","",remarks_hr_status)
   
   remarks_hrd_type <- ifelse(has_low_sv_load, sprintf("<%s SVs", LOW_SV_LOAD_THRES), "")
   remarks_hrd_type[ predictions$hr_status == HR_STATUS$PROFICIENT ] <- ""
   
   predictions$remarks_hr_status <- remarks_hr_status
   predictions$remarks_hrd_type <- remarks_hrd_type
   
   return(predictions)
}

PREDICTIONS <- chordPredict(MUTATION_CONTEXTS_BY_TYPE)

## Output ================================
logMessage(LOG_LEVEL$INFO, sprintf("Writing CHORD predictions to: %s", OUTPUT_FILE_PATH))
write.table(PREDICTIONS, OUTPUT_FILE_PATH, sep="\t", quote=FALSE, row.names=FALSE)

logMessage(LOG_LEVEL$INFO, sprintf("Completed CHORD predict"))