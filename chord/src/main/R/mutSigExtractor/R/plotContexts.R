MUT_SUBTYPE_COLORS <- list(
   snv=c(
      "C>A"="#06BAEB",
      "C>G"="#737373",
      "C>T"="#E12825",
      "T>A"="grey",
      "T>C"="#A1CD63",
      "T>G"="#EDC4C5"
   ),

   dbs=c(
      "AC>NN"="#03BCEC",
      "AT>NN"="#0165CB",
      "CC>NN"="#9FCC62",
      "CG>NN"="#006401",
      "CT>NN"="#FE9796",
      "GC>NN"="#E22823",
      "TA>NN"="#FCB066",
      "TC>NN"="#FC8100",
      "TG>NN"="#CC98FD",
      "TT>NN"="#4C0199"
   ),

   indel=c(
      "del.1.C.#"="#F3BF7B",
      "del.1.T.#"="#EE8633",
      "ins.1.C.#"="#B9DA94",
      "ins.1.T.#"="#569D40",
      "del.2.rep.#"="#F5CBB7",
      "del.3.rep.#"="#EE8E72",
      "del.4.rep.#"="#DE523E",
      "del.5+.rep.#"="#AD2D25",
      "ins.2.rep.#"="#D2E0EF",
      "ins.3.rep.#"="#9DC3DC",
      "ins.4.rep.#"="#5D97C4",
      "ins.5+.rep.#"="#2E65A5",
      "del.2.mh.#"="#E1E1ED",
      "del.3.mh.#"="#B5B6D5",
      "del.4.mh.#"="#8585B8",
      "del.5+.mh.#"="#5D4494",
      "del.#.mh.#"="#5D4494"
   ),

   sv=c(
      "DEL"="#9DD1C7",
      "DUP"="#FFFCBB",
      "INV"="#BDBAD7",
      "TRA"="#EB8677"
   )
)

####################################################################################################
#' Plot contexts
#'
#' @param x A matrix or vector of context counts
#' @param mut.type A string specifying the mutation type Can be 'snv', 'dbs', 'indel'
#' @param group A vector indicating the which group each sample belongs to. If provided, a plot will
#' be produced for each group. If unspecified, the function will assume that all rows of `x` belong
#' to the same group. If NULL, all rows will be considered as separate samples
#' @param y.axis.var.scale If TRUE, variable y-axis scale for each sample group will be used
#' @param horizontal.group.labels If TRUE, horizontal text for the sample group labels will be used
#' @param force.group.labels If TRUE, sample group label will be shown even if there is only one
#' group
#'
#' @return A ggplot2 object
#' @export
#'
plotContexts <- function(
   x, mut.type='auto', group='group1',
   y.axis.var.scale=T, horizontal.group.labels=F, force.group.labels=F
){

   ## Checks --------------------------------
   require(ggplot2)

   if(length(mut.type)!=1){ stop("`mut.type` must be a single string") }
   if(mut.type=='sbs'){ mut.type <- 'snv' }
   if(!(mut.type %in% c('snv','indel','dbs','auto','all'))){
      stop("`mut.type` must be one of the following: 'snv','indel', 'dbs', 'auto','all")
   }

   if(is.vector(x)){ x <- t(x) }
   if(length(colnames(x))==0){ stop("`x` must have colnames if a matrix, or names if a vector") }

   MUT_CONTEXTS_VALID <- list(
      snv=SUBS_CONTEXTS_96,
      dbs=DBS_TYPES$context,
      indel=INDEL_CONTEXTS
   )

   if(mut.type=='auto'){
      ## Infer mut type based on column names of matrix
      mut_type_tmp <- NULL
      for(i in names(MUT_CONTEXTS_VALID)){
         #i='snv'
         #x=cbind(x, test=NA)
         if(all(MUT_CONTEXTS_VALID[[i]] %in% colnames(x))){
            mut_type_tmp <- i
            break
         }
      }
      if(is.null(mut_type_tmp)){
         stop("`mut.type` is set to 'auto' but could not be inferred from colnames of `x`")
      }
      mut.type <- mut_type_tmp
      rm(mut_type_tmp)

   } else {
      if(mut.type=='all'){
         ## Join all mut types into one vector
         mut_contexts <- unlist(MUT_CONTEXTS_VALID, use.names=F)
      } else {
         mut_contexts <- MUT_CONTEXTS_VALID[[mut.type]]
      }

      invalid_mut_contexts <- colnames(x)[!all(mut_contexts %in% colnames(x))]
      if(length(invalid_mut_contexts)>0){
         stop('names of `x` contain invalid or missing mutation types:\n', paste(invalid_mut_contexts, collapse=', '))
      }
      rm(mut_contexts, invalid_mut_contexts)
   }

   ## Prep data --------------------------------
   ## Subset for relevant contexts
   df <- as.data.frame(x)
   if(mut.type=='all'){
      df <- df[,unlist(MUT_CONTEXTS_VALID, use.names=F)]
   } else {
      df <- df[, MUT_CONTEXTS_VALID[[mut.type]] ]
   }

   ## Add sample name column
   df$sample <- if(length(rownames(x))==0){
      paste0('sample',1:nrow(x))
   } else {
      rownames(x)
   }
   rownames(df) <- NULL

   ## Add group name column
   df$group <- if(is.null(group)){
      ## Plot all rows separately
      factor(1:nrow(x), 1:nrow(x))
   } else {
      factor(group, unique(group))
   }

   ## Convert to long form dataframe
   df <- reshape2::melt(
      df,
      id.vars=colnames(df)[!(colnames(df) %in% colnames(x))] ## Get metadata columns by removing context names
   )

   ## Summary stats
   if(nrow(x)!=1){
      agg <- aggregate(
         df$value,
         list(group=df$group, context=df$variable),
         function(x){ c( median(x), quantile(x, c(0.25,0.75)) ) }
      )
      agg <- cbind(
         subset(agg, select=-x),
         structure(as.data.frame(agg$x), names=c('median','q1','q3'))
      )
   } else {
      agg <- structure(df, names=c('sample','group','context','median'))
   }
   rm(df)

   ## Get mut subtype --------------------------------
   ## Make a vector of mut types, with names being the context
   mut_type_lookup <- structure(
      rep(names(MUT_CONTEXTS_VALID), sapply(MUT_CONTEXTS_VALID, length)),
      names=unlist(MUT_CONTEXTS_VALID, use.names=F)
   )
   mut_type_lookup <- mut_type_lookup[names(mut_type_lookup) %in% colnames(x)]

   ## Use regex to get mut subtype
   MUT_SUBTYPE_PARSERS <- list(
      snv=function(chr){ gsub('(^\\w\\[)|(\\]\\w$)','',chr) },
      indel=function(chr){ sub('\\d+\\+*$','#',chr) },
      dbs=function(chr){ sub('\\w{2}$','NN',chr) },
      sv=function(chr){ grep('^[[:upper:]]{3}',chr,value=T) }
   )

   mut_subtype_lookup <- sapply(names(mut_type_lookup), function(i){
      #i=names(mut_type_lookup)[[1]]
      MUT_SUBTYPE_PARSERS[[ mut_type_lookup[[i]] ]] (i)
   })

   ## Assign mut types and subtypes
   agg$context <- as.character(agg$context)
   agg$mut_type <- unname(mut_type_lookup[agg$context])
   #agg$mut_subtype_pre <- unname(mut_subtype_lookup[agg$context])
   agg$context <- factor(agg$context, unique(agg$context))

   ## Merge delmh contexts into one group
   agg$mut_subtype <- unname(mut_subtype_lookup[agg$context])
   agg$mut_subtype <- factor(agg$mut_subtype, unique(agg$mut_subtype))
   agg$mut_subtype_color <- agg$mut_subtype ## Keep original del mh context for bar coloring
   levels(agg$mut_subtype)[grep('del.*mh',levels(agg$mut_subtype))] <- 'del.#.mh.#'

   ## Get mut subtype colors
   colors <- unlist(unname(MUT_SUBTYPE_COLORS))
   colors <- colors[names(colors) %in% levels(agg$mut_subtype_color)]

   ## Plot --------------------------------
   ## Simplify context names for x axis
   MUT_SUBTYPE_PARSERS_2 <- list(
      snv=function(chr){ sub('\\[.+\\]','.',chr) },
      indel=function(chr){ paste( strsplit(chr,'[.]')[[1]][c(2,4)], collapse='_' ) },
      dbs=function(chr){ sub('^\\w+>','',chr) },
      sv=function(chr){ sub('^[[:upper:]]{3}_',chr) }
   )

   context_simple_names <- sapply(names(mut_type_lookup), function(i){
      #i=names(mut_type_lookup)[[1]]
      MUT_SUBTYPE_PARSERS_2[[ mut_type_lookup[[i]] ]] (i)
   })

   ## Init
   p <- ggplot(agg, aes(x=context, y=median))

   ## Facetting
   scales <- if(y.axis.var.scale){ 'free' } else { 'free_x' }
   if(force.group.labels || is.null(group) || length(unique(group))>1){
      p <- p + facet_grid(group~mut_subtype, scales=scales, space='free_x')
   } else {
      p <- p + facet_grid(~mut_subtype, scales=scales, space='free_x')
   }

   ## Main
   p <- p +
      geom_bar(aes(fill=mut_subtype_color),stat='identity') +
      scale_fill_manual(values=colors, name='Mutation type', guide='none') +
      scale_x_discrete(labels=context_simple_names)

   ## Error bars
   if(!is.null(group)){
      if(any(table(group)>1)){
         p <- p + geom_linerange(aes(ymin=q1, ymax=q3))
      }
   }

   ## Formatting
   p <- p +
      ylab('Contribution') +
      xlab('Mutation context') +
      theme_bw() +
      theme(
         panel.grid.minor.y=element_blank(),
         #panel.grid.minor.x=element_blank(),
         #panel.grid.major.x=element_blank(),
         axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
         panel.spacing.x=unit(0, "lines"),
         legend.position='left'
      )

   if(horizontal.group.labels){
      p <- p + theme(
         strip.text.y=element_text(angle=0, hjust=0)
      )
   }

   return(p)
}
