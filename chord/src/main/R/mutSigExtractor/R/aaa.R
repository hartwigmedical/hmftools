## Global variables
.onLoad <- function(libname, pkgname){

   ## hg19 ref genome --------------------------------
   DEFAULT_GENOME <- tryCatch({
      BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
   }, error=function(e){ return() })

   if(is.null(DEFAULT_GENOME)){
      warning("
   No reference genome loaded. Please install and load a BSgenome.
   For example:
      install.packages('BiocManager')
      BiocManager::install('BSgenome.Hsapiens.UCSC.hg19')
      library('BSgenome.Hsapiens.UCSC.hg19')

   Then specify the BSgenome to the ref.genome arguemnts to the relevant functions.
   For example:
      extractSigsSnv(..., ref.genome=BSgenome.Hsapiens.UCSC.hg19)
")
   }

   assign('DEFAULT_GENOME', DEFAULT_GENOME, envir=parent.env(environment()))

   assign(
      'SIG_METADATA_PATH',
      system.file('/sigs_v3.2_metadata.txt',package='mutSigExtractor'),
      envir=parent.env(environment())
   )
}


## SNV context types --------------------------------
SUBSTITUTIONS <- c('C>A','C>G','C>T','T>A','T>C','T>G')

C_TRIPLETS <- c(
   "ACA", "ACC", "ACG", "ACT",
   "CCA", "CCC", "CCG", "CCT",
   "GCA", "GCC", "GCG", "GCT",
   "TCA", "TCC", "TCG", "TCT")

T_TRIPLETS <- c(
   "ATA", "ATC", "ATG", "ATT",
   "CTA", "CTC", "CTG", "CTT",
   "GTA", "GTC", "GTG", "GTT",
   "TTA", "TTC", "TTG", "TTT")

CONTEXTS_96 <- c(rep(C_TRIPLETS, 3), rep(T_TRIPLETS, 3))
SUBSTITUTIONS_96 <- rep(SUBSTITUTIONS, each=16)
SUBS_CONTEXTS_96 <- paste0(substr(CONTEXTS_96,1,1), "[", SUBSTITUTIONS_96, "]", substr(CONTEXTS_96,3,3))


## Indel context types --------------------------------
INDEL_CONTEXTS <- (function(){
   rep_len_del <- c(1:5,'6+')
   rep_len_ins <- c(0:4,'5+')
   big_indel_len <- c(2:4,'5+')

   homopolymer_contexts <- c(
      paste0('del.1.C.',rep_len_del),
      paste0('del.1.T.',rep_len_del),
      paste0('ins.1.C.',rep_len_ins),
      paste0('ins.1.T.',rep_len_ins)
   )

   del_rep_contexts <- paste0(
      'del.',rep(big_indel_len,each=length(rep_len_del)),'.rep.',rep(rep_len_del,length(big_indel_len))
   )
   ins_rep_contexts <- paste0(
      'ins.',rep(big_indel_len,each=length(rep_len_del)),'.rep.',rep(rep_len_ins,length(big_indel_len))
   )

   mh_len <- c(1,1:2,1:3, 1:4,'5+')
   mh_del_len <- c(2,3,3,4,4,4,'5+','5+','5+','5+','5+')
   del_mh_contexts <- paste0('del.',mh_del_len,'.mh.',mh_len)

   c(homopolymer_contexts,del_rep_contexts,ins_rep_contexts,del_mh_contexts)
})()

## DBS context types --------------------------------
DBS_TYPES <- as.data.frame(
   matrix(
      c(
         'AC','AC>CA','GT>TG',
         'AC','AC>CG','GT>CG',
         'AC','AC>CT','GT>AG',
         'AC','AC>GA','GT>TC',
         'AC','AC>GG','GT>CC',
         'AC','AC>GT','GT>AC',
         'AC','AC>TA','GT>TA',
         'AC','AC>TG','GT>CA',
         'AC','AC>TT','GT>AA',
         'AT','AT>CA','AT>TG',
         'AT','AT>CC','AT>GG',
         'AT','AT>CG','AT>CG',
         'AT','AT>GA','AT>TC',
         'AT','AT>GC','AT>GC',
         'AT','AT>TA','AT>TA',
         'CC','CC>AA','GG>TT',
         'CC','CC>AG','GG>CT',
         'CC','CC>AT','GG>AT',
         'CC','CC>GA','GG>TC',
         'CC','CC>GG','GG>CC',
         'CC','CC>GT','GG>AC',
         'CC','CC>TA','GG>TA',
         'CC','CC>TG','GG>CA',
         'CC','CC>TT','GG>AA',
         'CG','CG>AT','CG>AT',
         'CG','CG>GC','CG>GC',
         'CG','CG>GT','CG>AC',
         'CG','CG>TA','CG>TA',
         'CG','CG>TC','CG>GA',
         'CG','CG>TT','CG>AA',
         'CT','CT>AA','AG>TT',
         'CT','CT>AC','AG>GT',
         'CT','CT>AG','AG>CT',
         'CT','CT>GA','AG>TC',
         'CT','CT>GC','AG>GC',
         'CT','CT>GG','AG>CC',
         'CT','CT>TA','AG>TA',
         'CT','CT>TC','AG>GA',
         'CT','CT>TG','AG>CA',
         'GC','GC>AA','GC>TT',
         'GC','GC>AG','GC>CT',
         'GC','GC>AT','GC>AT',
         'GC','GC>CA','GC>TG',
         'GC','GC>CG','GC>CG',
         'GC','GC>TA','GC>TA',
         'TA','TA>AT','TA>AT',
         'TA','TA>CG','TA>CG',
         'TA','TA>CT','TA>AG',
         'TA','TA>GC','TA>GC',
         'TA','TA>GG','TA>CC',
         'TA','TA>GT','TA>AC',
         'TC','TC>AA','GA>TT',
         'TC','TC>AG','GA>CT',
         'TC','TC>AT','GA>AT',
         'TC','TC>CA','GA>TG',
         'TC','TC>CG','GA>CG',
         'TC','TC>CT','GA>AG',
         'TC','TC>GA','GA>TC',
         'TC','TC>GG','GA>CC',
         'TC','TC>GT','GA>AC',
         'TG','TG>AA','CA>TT',
         'TG','TG>AC','CA>GT',
         'TG','TG>AT','CA>AT',
         'TG','TG>CA','CA>TG',
         'TG','TG>CC','CA>GG',
         'TG','TG>CT','CA>AG',
         'TG','TG>GA','CA>TC',
         'TG','TG>GC','CA>GC',
         'TG','TG>GT','CA>AC',
         'TT','TT>AA','AA>TT',
         'TT','TT>AC','AA>GT',
         'TT','TT>AG','AA>CT',
         'TT','TT>CA','AA>TG',
         'TT','TT>CC','AA>GG',
         'TT','TT>CG','AA>CG',
         'TT','TT>GA','AA>TC',
         'TT','TT>GC','AA>GC',
         'TT','TT>GG','AA>CC'
      ),
      ncol=3, byrow=T, dimnames=list(NULL, c('ref','context','context_rev_comp'))
   )
)
