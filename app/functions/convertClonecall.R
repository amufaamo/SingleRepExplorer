.convertClonecall <- function(x) {
  
  clonecall_dictionary <- hash::hash(
    "gene" = "CTgene",
    "genes" = "CTgene",
    "ctgene" = "CTgene",
    "ctstrict" = "CTstrict",
    "nt" = "CTnt",
    "nucleotide" = "CTnt",
    "nucleotides" = "CTnt",
    "ctnt" = "CTnt",
    "aa" = "CTaa",
    "amino" = "CTaa",
    "ctaa" = "CTaa",
    "gene+nt" = "CTstrict",
    "strict" = "CTstrict",
    "ctstrict" = "CTstrict",
    "clonotypeid" = "raw_clonotype_id",
    "subclonotypeid" = "exact_subclonotype_id"
  )
  
  if (!is.null(clonecall_dictionary[[tolower(x)]])) {
    return(clonecall_dictionary[[tolower(x)]])
  } else {
    warning("A custom variable ", x, " will be used to call clones")
    return(x)
  }
}
