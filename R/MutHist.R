#' @title MutHist
#' @description plot an histogram of mutation type counts
#' @param counts a data.frame of mutations counts per type
#'
#' @return a graphical output: histogram of counts per mutation type
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
#' MutHist(mut_summary)
#'

MutHist <- function(counts){
  n<- NULL
  mycolors <- brewer.pal(3, "Set2")
  mycolors <- colorRampPalette(mycolors)(3)
  ggplot(counts, aes(x = type, y = n, fill=n)) +
    geom_bar(stat = "identity", fill = mycolors) +
    coord_flip() +
    geom_text(aes(label = n), hjust =-0.25) +
    ylim(0, 11550) + xlab(NULL) + ylab("count")
}
