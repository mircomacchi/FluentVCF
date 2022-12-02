#' @title mut_type
#' @description determines for each mutation
#'     the corresponding mutation type in VCF file
#' @param mutation_list an input VCF file
#' @param ref.genome reference genome
#' @param context_length distance from mutation in ref.genome
#'
#' @return a vector of mutation signatures
#' @export
#' @examples
#' vcffile <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
#' vcf <- VariantAnnotation::readVcf(vcffile, "hg38")
#' Hsapiens <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
#' vcfDF <- as(MatrixGenerics::rowRanges(vcf[1:100]), "DataFrame")
#' muts <- mut_type(vcfDF, ref.genome = Hsapiens)
#' head(muts)
#'
  mut_type <- function(mutation_list,ref.genome,context_length=3){
    if (context_length < 3){
      stop("context_length must be >= 3!")
    }
    mutations <- NULL
    chr <- mutation_list[[1]]@seqnames@values
    for (i in seq_len(mutation_list@nrows)){
      start_coord <- mutation_list[[1]]@ranges@start[i]
      gr <- GRanges(paste0("chr",as.vector(chr),":", as.character(start_coord),":*"), seqinfo=seqinfo(ref.genome))
      segment <- getSeq(ref.genome, gr+context_length)

      signature_end <- 2
      mutation <- (paste0(substr(as.character(segment),context_length,context_length),
                          "[",as.character(mutation_list$REF[i]),">",as.character(mutation_list$ALT[[i]]),"]",
                          substr(as.character(segment),start = context_length+2,context_length+signature_end)))

      while (nchar(mutation)-4 < context_length){
        mutation <- (paste0(substr(as.character(segment),context_length,context_length),
                            "[",as.character(mutation_list$REF[i]),">",as.character(mutation_list$ALT[[i]]),"]",
                            substr(as.character(segment),context_length+2,context_length+signature_end)))
        signature_end <- signature_end +1
      }

      rev_compl_segment <- reverseComplement(segment)
      rev_compl <- paste0(substr(rev_compl_segment,context_length,context_length),
                          "[",reverseComplement(mutation_list$REF[i]),">",reverseComplement(mutation_list$ALT[[i]]),"]",
                          substr(rev_compl_segment,context_length+2,context_length+signature_end-1))

      if (mutation %in% mutations | rev_compl %in% mutations){
        next
      }
      else{
        mutations <- c(mutations,mutation)
      }
    }
    return(mutations)
  }

