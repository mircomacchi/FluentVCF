#' mut_summary
#'
#' @param df a data.frame of mutations with mutation types
#'
#' @return a summary of counts per mutation
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
#' mut_summary <- mut_summary(Basedf)
#' mut_summary

  mut_summary <- function(df){
    df <- df %>% dplyr::count(type)}
