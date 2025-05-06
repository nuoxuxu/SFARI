#' Read a VCF into a tidy tibble
#'
#' @param file  Path to a VCF file.
#' @param info  Character vector of INFO‑field keys to pull out into columns.
#'
#' @return A tibble with at least
#'         chrom, start, end, ref, alt and one column per element of `info`.
#' @examples
#' df <- read_vcf("example.vcf", info = c("CLNVC", "GENEINFO"))
read_vcf <- function(file, info = c("CLNVC", "GENEINFO")) {
  library(readr)     # read_tsv()
  library(dplyr)     # mutate(), select()
  library(stringr)   # str_match()

  ## 1.  read the VCF body ----------------------------------------------------
  vcf <- read_tsv(
    file,
    comment   = "#",                       # skip metadata / header lines
    col_names = c("chrom", "start", "id",
                  "ref", "alt", "qual",
                  "filter", "info"),
    col_types = cols(chrom = col_character())
  )

  ## 2.  coordinates ----------------------------------------------------------
  vcf <- vcf %>%
    mutate(end = start + nchar(ref) - 1)   # 1‑based, inclusive

  ## 3.  pull requested INFO attributes --------------------------------------
  for (attribute in info) {
    pattern <- paste0(attribute, "=([^;]+)")
    vcf <- vcf %>%
      mutate(
        !!attribute := str_match(.data$info, pattern)[, 2]
      )
  }

  ## 4.  tidy output ----------------------------------------------------------
  vcf %>%
    select(chrom, start, end, ref, alt, everything(), all_of(info))
}
