#' Title
#'
#' Description
#'
#' Details
#'
#' @param snp_data To be completed
#' @param kept_markers To be completed
#' @param center_pos To be completed
#'
#' @return To be completed
#'
#' @examples
#' NULL
#'
haplotypes <- function(snp_data, kept_markers, center_pos = NULL) {

  # Haplotype assignment will be done from markers that are in perfect LD
  #  with a least one of the final set of markers ("Selected_clusters")

  # Assigning linkage disequilibrium data to a local object and setting NA to 0
  LD_data <- snp_data$LD
  LD_data[is.na(LD_data)] <- 0
  # We select columns corresponding
  LD_data <- LD_data[ , as.character(kept_markers$rs)]
  # Identifying markers that are in perfect (1) LD with one of the selected clusters
  which_markers <- rownames(LD_data)[apply(LD_data, 1, function(x) any(x == 1))]
  # The set of markers for haplotype assignment comprises selected markers and their
  #  perfect linkage disequilibrium mates
  which_markers <- unique(c(which_markers, kept_markers$rs))

  # Extracting genotype data from the object
  genotypes <- snp_data[["Genotypes"]]@.Data
  # We keep the markers of interest based on column names
  genotypes <- genotypes[ , as.character(which_markers)]
  # Keeping track of individual names throughout the function
  ind_names <- rownames(genotypes)

  # Coercing the genotype data to a character data.frame and setting "00" to NA
  genotypes <- as.data.frame(apply(genotypes, 2, as.character), stringsAsFactors = FALSE)
  genotypes <- as.data.frame(lapply(genotypes, function(x) ifelse(x == "00", NA, x)), stringsAsFactors = FALSE)

  # Initializing the haplotypes data.frame
  haplotypes <- genotypes[0,]
  # Reordering the individual by increasing number of missing data
  #  This allows individuals without missing data to be processed first
  n_missing <- apply(genotypes, 1, function(x) sum(is.na(x)))
  genotypes <- genotypes[order(n_missing), ]
  ind_names <- ind_names[order(n_missing)] # Keeping track of IDs

  # Initializing a named vector of the haplotypes assigned
  hap_assignment <- numeric(nrow(genotypes))
  names(hap_assignment) <- ind_names

  # Iterating over all the individuals
  for(i in 1:nrow(genotypes)) {

    #gen_string is the "genotype string" of individual i
    gen_string <- as.character(genotypes[i, ])
    #a vector indicating which haplotypes individual is compatible with
    compatible <- logical(nrow(haplotypes))

      if(nrow(haplotypes)) { # The first genotype is automatically a haplotype
        for(j in 1:nrow(haplotypes)) {
          # hap_string is the "haplotype string" of haplotype j
          hap_string <- as.character(haplotypes[j, ])
          if(sum(gen_string == hap_string, na.rm = TRUE) ==
             sum(!is.na(gen_string) & !is.na(hap_string))) {
            # Checking if genotype i is compatible with haploype j
            #  i.e. whether there are no conflicts for any marker
            compatible[j] <- TRUE
          }
        }
      }

      if(!any(compatible)) {
        # New haplotype if it is not compatible with any
        haplotypes <- rbind(haplotypes, genotypes[i, , drop = FALSE])
        rownames(haplotypes)[nrow(haplotypes)] <- as.character(nrow(haplotypes))
        # Haplotype number of individual i is the last haplotype added
        hap_assignment[i] <- nrow(haplotypes)
      } else if (sum(compatible) == 1) {
        # Assigning a pre-existing haplotype number
        hap_assignment[i] <- which(compatible)
      } else if (sum(compatible) > 1) {
        # Haplotype is NA if it cannot be unambiguously assigned
        hap_assignment[i] <- NA
      }
  }

  # Haplotypes get reordered as a function of the genotypes of markers
  #  close to gene center (if center_pos is provided)
  if(!is.null(center_pos)) {
    # Ordering marker names from closest to furthest from gene center
    marker_order <- snp_data$Markers$rs[order(abs(snp_data$Markers$pos - center_pos))]
    # Keeping only those markers that are of interest in this function
    marker_order <- marker_order[marker_order %in% which_markers]
    # Generating a list (actually a data.frame) used to reorder haplotypes
    marker_order <- haplotypes[ , marker_order]
    # Effectively reordering the haplotypes
    haplotypes <- haplotypes[do.call("order", marker_order), ]
  }

  # Here we are going to generate a vector of haplotype names
  #  The first 26 haplotypes get single letters, then two letters, and so on

  # Determining the maximum number of letters a haplotype name can have
  n_letters <- ceiling(nrow(haplotypes) / 26)
  # Initializing a vector of haplotype names
  hapnames <- character()

  # Filling the vector with single letters, then doubles, and so on
  for(i in 1:n_letters) {
    hapnames <- c(hapnames, do.call("paste0", replicate(i, LETTERS, simplify = FALSE)))
  }

  # Cutting the vector at the desired length
  hapnames <- hapnames[1:nrow(haplotypes)]

  # This step renames the haplotypes in hap_assignment
  names(hapnames) <- as.character(rownames(haplotypes))
  hap_assignment <- hapnames[as.character(hap_assignment)]

  # Preparing the assignment data.frame for output
  ind_h <- data.frame(ind = ind_names,
                      haplotype = hap_assignment,
                      stringsAsFactors = FALSE)

  # Preparing the haplotype data.frame for output
  haplotypes <- suppressWarnings(apply(haplotypes, 2, as.raw))
  rownames(haplotypes) <- hapnames
  haplotypes <- haplotypes[ , as.character(kept_markers$rs)]

  # Preparing the marker data.frame for output
  output_markers <- snp_data$Markers
  marker_pos <- match(kept_markers$rs, output_markers$rs)
  output_markers <- output_markers[marker_pos, ]

  # Building and returning the output
  return(list(Genotypes  = new("SnpMatrix", .Data = haplotypes),
              Markers    = output_markers,
              Assignment = ind_h))
}
