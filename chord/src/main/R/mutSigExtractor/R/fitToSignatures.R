## Slow R implementation ##########################################################################
## Machine epsilon (floating-point relative accuracy). Adapted from eps() function from pracma R package
EPS <- (function(x=1.0) {
   x <- max(abs(x))

   if (x < .Machine$double.xmin) {
      .Machine$double.xmin
   } else {
      2^floor(log2(x)) * .Machine$double.eps
   }
})()

#' Non-negative least squares fitting
#'
#' @rdname lsqnonneg
#'
#' @description Ripped least squares fitting function from pracma
#'
#' @param mut.context.counts A numeric vector of mutation context counts
#' @param signature.profiles A matrix containing the mutational signature profiles, where rows are
#' the mutation contexts and the columns are the mutational signatures.
#'
#' @return A vector returned containing the the absolute contribution of each signature
#'
lsqnonneg <- function (x, ...) {
   UseMethod("lsqnonneg", x)
}

#' @rdname lsqnonneg
#' @method lsqnonneg default
#' @export
lsqnonneg.default <- function(mut.context.counts, signature.profiles){
   m <- nrow(signature.profiles)
   n <- ncol(signature.profiles)

   if (m != length(mut.context.counts)){
      stop(
         "No. of contexts in mut.context.counts (",length(mut.context.counts),")",
         "does not match no. of contexts in signature.profiles (",m,")"
      )
   }

   tol <- 10 * EPS * norm(signature.profiles, type = "2") * (max(n, m) + 1)

   x <- rep(0, n)             # initial point
   P <- logical(n); Z <- !P   # non-active / active columns

   resid <- mut.context.counts - signature.profiles %*% x
   w <- t(signature.profiles) %*% resid
   wz <- numeric(n)

   # iteration parameters
   outeriter <- 0; it <- 0
   itmax <- 3 * n; exitflag <- 1

   while (any(Z) && any(w[Z] > tol)) {
      outeriter <- outeriter + 1
      z <- numeric(n)
      wz <- rep(-Inf, n)
      wz[Z] <- w[Z]
      im <- which.max(wz)
      P[im] <- TRUE; Z[im] <- FALSE
      z[P] <- qr.solve(signature.profiles[, P], mut.context.counts)

      while (any(z[P] <= 0)) {
         it <- it + 1
         if (it > itmax) stop("Iteration count exceeded.")

         Q <- (z <= 0) & P
         alpha <- min(x[Q] / (x[Q] - z[Q]))
         x <- x + alpha*(z - x)
         Z <- ((abs(x) < tol) & P) | Z
         P <- !Z
         z <- numeric(n)
         z[P] <- qr.solve(signature.profiles[, P], mut.context.counts)
      }
      x <- z
      resid <- mut.context.counts - signature.profiles %*% x
      w <- t(signature.profiles) %*% resid
   }

   names(x) <- colnames(signature.profiles)

   return(x)
}

#' @rdname lsqnonneg
#' @method lsqnonneg matrix
#' @export
lsqnonneg.matrix <- function(mut.context.counts, signature.profiles, verbose=F){

   mut.context.counts <- as.matrix(mut.context.counts)

   if(verbose){
      message('Fitting input matrix to signature profiles...')
      counter <- 0
      pb <- txtProgressBar(max=nrow(mut.context.counts), style=2)
   }
   out <- apply(mut.context.counts, 1, function(i){
      if(verbose){ counter <<- counter + 1; setTxtProgressBar(pb, counter) }
      #if(verbose){ counter <<- counter + 1; print(rownames(mut.context.counts)[counter]) }
      lsqnonneg.default(i, signature.profiles)
   })
   if(verbose){ message('\n') }
   out <- t(out)
   colnames(out) <- colnames(signature.profiles)
   rownames(out) <- rownames(mut.context.counts)

   return(out)
}

#' @rdname lsqnonneg
#' @method lsqnonneg data.frame
#' @export
lsqnonneg.data.frame <- lsqnonneg.matrix

####################################################################################################
if( !('NNLM' %in% installed.packages()) ){
   USE_LSQ_R <- TRUE
   message(
      'NNLM package not installed. A (slow) R implementation will be used for fitting mutation contexts to signatures. ',
      '\nThe fast C++ implementation from NNLM is recommended when fitting context matrices with many samples.'
   )

} else {
   USE_LSQ_R <- FALSE
}

#' Linear least-squares (non-negative) fitting
#'
#' @rdname fitToSignatures
#'
#' @description Calculate the contribution of each mutational signature in a sample, given a vector
#' of mutation contexts. This is done by least-squares fitting, i.e. find a linear non-negative
#' combination of mutation signatures that optimally reconstructs the mutation matrix. This is
#' performed by `fitToSignatures()`
#'
#' Only performing least-squares fitting is however potentially prone to overfitting.
#' `fitToSignaturesStrict()` tries to solve this problem by performing an initial fit, then
#' iteratively removing the signature with the lowest contribution and with fitting then repeated.
#' Each time the cosine distance between the original and reconstructed mutation context profile is
#' calculated. Iterations are stopped when cosine distance > max.delta. The second-last set of
#' signatures is then returned.
#'
#' `fitToSignaturesFast()` and `fitToSignaturesFastStrict()` are wrappers for
#' `fitToSignatures(..., use.lsq.r=F)` and `fitToSignaturesStrict(..., use.lsq.r=F)`
#' for backwards compatibility
#'
#' @param signature.profiles A matrix containing the mutational signature profiles, where rows are
#' the mutation contexts and the columns are the mutational signatures.
#' @param mut.contexts A vector of mutation contexts, or a matrix where rows are samples and columns
#' are mutation contexts
#' @param max.delta See description.
#' @param detailed.output Only for `fitToSignaturesStrict()`. Also return results from the
#' iterative fitting? Includes: cosine similarity between the original and reconstructed mutation
#' context; signatures removed at each iteration. Useful for plotting fitting performance.
#' @param use.lsq.r If TRUE, least squares fitting will be performed using the slow R
#' function `lsqnonneg()`. If FALSE, the much faster NNLM::nnlm() will be used (written in C++)
#' @param scale.contrib If TRUE, the signature contributions will be scaled so that the total
#' signature contribution is equal to the total number of mutations. The reason for this is that
#' when performing least squares fitting, there are unaccounted for mutations (aka residual; not
#' necessarily an integer amount of mutations).
#' @param verbose Show progress messages?
#'
#' @return If vector is provided to mut.contexts, a vector returned containing the the absolute
#' contribution of each signature (i.e. the number of mutations contributing to each mutational
#' signature). If a matrix is provided, a matrix of absolute contributions is returned
#' @export
#'
fitToSignatures <- function(
   mut.context.counts, signature.profiles,
   use.lsq.r=USE_LSQ_R, scale.contrib=T, verbose=F
){

   ## Checks --------------------------------
   context_names <- if(is.vector(mut.context.counts)){
      names(mut.context.counts)
   } else {
      mut.context.counts <- as.matrix(mut.context.counts)
      colnames(mut.context.counts)
   }
   signature.profiles <- as.matrix(signature.profiles)

   if( length(context_names)!=nrow(signature.profiles) ){
      stop(
         "No. of contexts in mut.context.counts (",length(context_names),") ",
         "does not match no. of contexts in signature.profiles (",nrow(signature.profiles),")"
      )
   }

   if( !(identical(context_names, rownames(signature.profiles))) ){
      warning("Context names of mut.context.counts and signature.profiles do not match. Fitting may not be correct.")
   }

   ## Least squares fitting --------------------------------
   if(use.lsq.r){
      sig_contrib <- lsqnonneg(mut.context.counts, signature.profiles)

   } else {
      if(is.vector(mut.context.counts)){
         sig_contrib <- NNLM::nnlm(
            signature.profiles,
            matrix(mut.context.counts, ncol=1)
         )$coefficients[,1]
      } else {
         sig_contrib <- t(
            NNLM::nnlm(
               signature.profiles,
               t(mut.context.counts)
            )$coefficients
         )
      }
   }

   if(!scale.contrib){ return(sig_contrib) }

   ## Scale to total mutational load --------------------------------
   if(is.vector(mut.context.counts)){
      sig_contrib <- sig_contrib * (sum(mut.context.counts) / sum(sig_contrib))
   } else {
      sig_contrib <- sig_contrib * rowSums(mut.context.counts) / rowSums(sig_contrib)
   }

   return(sig_contrib)
}

#' @rdname fitToSignatures
fitToSignaturesFast <- function(...){
   fitToSignatures(..., use.lsq.r=F)
}

#' @rdname fitToSignatures
fitToSignaturesStrict <- function(
   mut.context.counts, signature.profiles, max.delta=0.004,
   detailed.output=F, use.lsq.r=USE_LSQ_R, scale.contrib=T, verbose=F
){

   ## Checks --------------------------------
   context_names <- if(is.vector(mut.context.counts)){
      names(mut.context.counts)
   } else {
      mut.context.counts <- as.matrix(mut.context.counts)
      colnames(mut.context.counts)
   }
   signature.profiles <- as.matrix(signature.profiles)

   if( length(context_names)!=nrow(signature.profiles) ){
      stop(
         "No. of contexts in mut.context.counts (",length(context_names),") ",
         "does not match no. of contexts in signature.profiles (",nrow(signature.profiles),")"
      )
   }

   if( !(identical(context_names, rownames(signature.profiles))) ){
      warning("Context names of mut.context.counts and signature.profiles do not match. Fitting may not be correct.")
   }

   ## Choose R or C++ implementation of lsq fitting --------------------------------
   if(is.vector(mut.context.counts)){
      mut.context.counts <- matrix(mut.context.counts, ncol=1, dimnames=list(names(mut.context.counts), NULL))
   } else {
      mut.context.counts <- t(mut.context.counts)
   }

   if(use.lsq.r){
      if(verbose){ message('Using R implementation for least squares fitting...') }
      if(ncol(mut.context.counts)==1){
         f_fit <- function(mut.context.counts, signature.profiles){
            matrix(
               lsqnonneg.default(mut.context.counts[,1], signature.profiles),
               ncol=1, dimnames=list(colnames(signature.profiles), NULL)
            )
         }
      } else {
         f_fit <- function(mut.context.counts, signature.profiles){
            t(lsqnonneg.matrix(t(mut.context.counts), signature.profiles))
         }
      }
   } else {
      if(verbose){ message('Using C++ implementation for least squares fitting...') }
      require(NNLM)
      f_fit <- function(mut.context.counts, signature.profiles){
         nnlm(signature.profiles, mut.context.counts)$coefficients
      }
   }

   ## --------------------------------
   if(verbose){ message('Performing first fit...') }

   #mut.context.counts <- t(mut.context.counts)
   fit_init <- f_fit(mut.context.counts, signature.profiles)
   #print(fit_init)

   sig_pres <- rowSums(fit_init) != 0
   my_signatures_total <- signature.profiles[, sig_pres, drop=FALSE]

   ## --------------------------------
   if(verbose){ message('Performing signature selection per sample...') }

   n_sigs <- ncol(my_signatures_total)
   n_samples <- ncol(mut.context.counts)

   l_contribs <- list()
   l_sims <- list()
   l_removed_sigs <- list()

   if(sum(sig_pres)==1){
      l_contribs[[1]] <- fit_init[,1]
   } else {
      if(verbose==2){ pb <- txtProgressBar(max=n_samples, style=3) }

      for (i in 1:n_samples) {
         #i=1

         if(verbose==1){ message(' [',i,'/',n_samples,'] ', colnames(mut.context.counts)[i]) }
         if(verbose==2){ setTxtProgressBar(pb, i) }

         my_signatures <- my_signatures_total
         mut_mat_sample <- mut.context.counts[, i, drop=FALSE]

         ## Initialize output vectors
         ## Keep track of the cosine similarity and which signatures are removed.
         sims <- rep(NA, n_sigs)
         removed_sigs <- rep(NA, n_sigs)

         contrib <- structure(
            rep(0, ncol(signature.profiles)),
            names=colnames(signature.profiles)
         )
         if(sum(mut_mat_sample[,1])==0){
            l_contribs[[i]] <- contrib
            l_sims[[i]] <- sims
            l_removed_sigs[[i]] <- removed_sigs
            next
         }

         ## Fit again
         fit_res <- list()
         #fit_res$contribution <- NNLM::nnlm(my_signatures, mut_mat_sample)$coefficients
         fit_res$contribution <- f_fit(mut_mat_sample, my_signatures)
         fit_res$reconstructed <- my_signatures %*% fit_res$contribution

         sim <- cosSim.default(fit_res$reconstructed[,1], mut_mat_sample[,1])
         sims[[1]] <- sim

         ## Sequentially remove the signature with the lowest contribution
         for (j in 2:n_sigs) {
            #j=2

            # Remove signature with the weakest relative contribution
            contri_order <- order(fit_res$contribution[,1] / sum(fit_res$contribution[,1]))
            weakest_sig_index <- contri_order[1]
            weakest_sig <- colnames(my_signatures)[weakest_sig_index]
            removed_sigs[[j]] <- weakest_sig
            signatures_sel <- my_signatures[, -weakest_sig_index, drop=FALSE]

            # Fit with new signature selection
            fit_res <- list()
            #fit_res$contribution <- NNLM::nnlm(signatures_sel, mut_mat_sample)$coefficients
            fit_res$contribution <- f_fit(mut_mat_sample, signatures_sel)
            fit_res$reconstructed <- signatures_sel %*% fit_res$contribution

            sim_new <- cosSim.default(fit_res$reconstructed, mut_mat_sample)

            if(is.na(sim_new)){
               sim_new <- 0
               # if(verbose){
               #    warning("New similarity between the original and the reconstructed
               #             spectra after the removal of a signature was NaN.
               #             It has been converted into a 0.
               #             This happened with the following fit_res:")
               #    print(fit_res)
               # }
            }
            sims[[j]] <- sim_new

            # Check if the loss in cosine similarity between the original vs reconstructed after removing the signature is below the cutoff.
            if(sim-sim_new <= max.delta){
               my_signatures <- signatures_sel
               sim <- sim_new
            } else {
               break
            }
         }

         ## Fit with the final set of signatures
         ## Fill in 0 for absent signatures
         contrib_pre <- f_fit(mut_mat_sample, my_signatures)[,1]
         contrib[names(contrib_pre)] <- contrib_pre

         l_contribs[[i]] <- contrib
         l_sims[[i]] <- sims
         l_removed_sigs[[i]] <- removed_sigs
      }
   }

   ## --------------------------------
   if(verbose==2){ message('\n') }
   if(verbose){ message('Returning output...') }

   ## Signature contributions
   m_contribs <- do.call(rbind, l_contribs)
   rownames(m_contribs) <- colnames(mut.context.counts)

   if(scale.contrib){
      m_contribs <- m_contribs * colSums(mut.context.counts) / rowSums(m_contribs)
   }

   if(!detailed.output){ return(m_contribs) }

   ## Stats from removing ref signatures
   m_sim <- do.call(rbind, l_sims)
   rownames(m_sim) <- colnames(mut.context.counts)

   m_removed_sigs <- do.call(rbind, l_removed_sigs)
   rownames(m_removed_sigs) <- colnames(mut.context.counts)

   list(
      sig_contrib=m_contribs,
      removed_sigs_common=names(sig_pres)[ !sig_pres ],
      removed_sigs_iter=m_removed_sigs,
      cos_sims=m_sim
   )
}

#' @rdname fitToSignatures
fitToSignaturesFastStrict <- function(...){
   fitToSignaturesStrict(..., use.lsq.r=F)
}
