#' @title type
#' @description For a data.frame of mutations,provide mutation type.
#'  meaning:substitution, insertion or deletion
#' @param df a data.frame of mutations
#'
#' @return a function to be applied on such data.frame
#' @export
#' @examples
#' vcffile <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
#' vcf <- VariantAnnotation::readVcf(vcffile, "hg38")
#' refBase <- VariantAnnotation::ref(vcf)
#' refBase <- unlist(lapply(refBase, as.character))
#' altBase <- VariantAnnotation::alt(vcf)
#' altBase <- unlist(lapply(altBase, as.character))
#' startRef <- stats::start(vcf)
#' endRef <- stats::end(vcf)
#' chr <- GenomeInfoDb::seqnames(vcf)
#' Basedf <- data.frame(ref = refBase, alt=altBase, startRef= startRef, endRef= endRef, chr = chr)
#' Basedf$type <- apply(Basedf,1,type)
#' head(Basedf)
#'
  type <- function(df){
    ref <- df[1]
    alt <- df[2]
    if (nchar(ref)> nchar(alt)){
      df[6] <-'del'
    }
    else if (nchar(ref)< nchar(alt)){
        df[6] <-'ins'
    }
    else{
        df[6] <- 'sub'}
        }
