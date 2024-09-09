#' Cosine similarity
#'
#' @description Calculate cosine similarity between two vectors of the same length, or cosine
#' similarity between two matrices
#'
#' @param x A numeric vector or matrix
#' @param y A numeric vector or matrix
#' @param all.vs.all (When x and y are matrices) If TRUE, all rows of x are compared with all rows
#' of y. If FALSE, rows of x are compared with rows of y 'side-by-side'
#' @param by.row (When x and y are matrices) If TRUE, compare rows. If FALSE, compare columns
#'
#' @return A numeric vector or matrix
#' @rdname cosSim
#' @export
#'
#' @examples
#' set.seed(1)
#' x <- matrix(runif(200), nrow=20)
#' y <- matrix(runif(200), nrow=20)
#'
#' ## Works with matrices
#' cosSim(x, y)
#'
#' ## ...or vectors
#' cosSim(x[1,], y[1,])
#'
cosSim <- function(x, ...) {
   UseMethod("cosSim", x)
}

#' @rdname cosSim
#' @method cosSim default
#' @export
cosSim.default <- function(x, y){
   if(length(x)!=length(y)){ stop('`x` and `y` must be the same length') }

   ## Convert to numeric to prevent integer overflow
   x <- as.numeric(x)
   y <- as.numeric(y)

   sum(x*y)/sqrt(sum(x^2)*sum(y^2))
}

#' @rdname cosSim
#' @method cosSim matrix
#' @export
cosSim.matrix <- function(x, y, all.vs.all=FALSE, by.row=TRUE) {

   ## Convert dataframe to matrices
   ## Convert to numeric to prevent integer overflow
   x <- matrix(as.numeric(as.matrix(x)), nrow=nrow(x), dimnames=dimnames(x))
   y <- matrix(as.numeric(as.matrix(y)), nrow=nrow(y), dimnames=dimnames(y))

   if(!by.row){
      x <- t(x)
      y <- t(y)
   }

   ## Compare rows of x with y 'side-by-side' --------------------------------
   if(!all.vs.all){
      if(!identical(dim(x), dim(y))){ stop('`x` and `y` must have the same dimensions when `all.vs.all`=FALSE') }
      return(
         rowSums(x*y) / sqrt( rowSums(x^2) * rowSums(y^2) )
      )
   }

   ## Compare every row of x with every row of y --------------------------------
   if(ncol(x)!=ncol(y)){
      stop(sprintf(
         '`x` and `y` must have the same number of %s when `all.vs.all`=TRUE',
         if(by.row){ 'columns' } else { 'rows' }
      ))
   }

   ## Slow version
   # apply(x, 1, function(i){
   #    apply(y, 1, function(j){
   #       sum(i*j)/sqrt(sum(i^2)*sum(j^2))
   #    })
   # })

   ## Fast version
   apply(x,1,function(i){
      #i <- x[1,]
      numer <- rowSums(sweep(y,2,i,'*'))
      denom <- sqrt( sum(i^2) * rowSums(y^2) )
      numer/denom
   })
}

#' @rdname cosSim
#' @method cosSim data.frame
#' @export
cosSim.data.frame <- function(...){
   cosSim.matrix(...)
}

