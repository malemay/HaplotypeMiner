#' Title
#'
#' Description
#'
#' Details
#'
#' @param snp_data To be completed
#' @param gene_pos To be completed
#' @param kept_markers To be completed
#' @param assignment To be completed
#' @param name_order To be completed
#'
#' @return To be completed
#' @export
#'
#' @examples
#' NULL
#'
genotype_plot <- function(snp_data, gene_pos = NULL,
                          kept_markers = NULL,
                          assignment = NULL,
                          name_order = FALSE) {

  # Extracting data from snp_data for easier manipulation
  marker_data <- snp_data$Markers
  genotype_data   <- snp_data$Genotypes@.Data

  # Generates a data.frame with one row per marker and one column per individual
  genotype_data <- as.data.frame(apply(genotype_data, 1, as.character),
                                 stringsAsFactors = FALSE)

  # Creating a named (ind. names) character vector of haplotypes
  if(!is.null(assignment)) {
    # If haplotypes have been assigned, individuals are given their haplotype name
    haplotypes <- assignment$haplotype
    # The "name" of a non-unambiguously assigned haplotype is "-"
    haplotypes[is.na(haplotypes)] <- "-"
    names(haplotypes) <- assignment$ind
  } else {
    # Else individuals are assigned a "haplotype string" concatenated from their genotypes
    haplotypes <- apply(genotype_data, 2, function(x) paste0(x, collapse = ""))
    names(haplotypes) <- names(genotype_data)
  }

  # Adding some metadata columns to genotype_data
  genotype_data$alleles <- marker_data$alleles
  genotype_data$x_pos <- order(marker_data$pos)
  genotype_data$rs <- marker_data$rs

  # Generating a long-form data.frame, better suited for plotting with ggplot2
  #  This data.frame will contain one row for each tile to be plotted.
  #  (This corresponds to the number of individuals * the number of markers)
  genotype_data <- reshape2::melt(genotype_data,
                                  id.vars = c("alleles", "rs", "x_pos"),
                                  variable.name = "ind",
                                  value.name = "genotype")

  # Haplotypes are retrieved from object "haplotypes" via individual names
  genotype_data[["haplotypes"]] <- haplotypes[as.character(genotype_data[["ind"]])]

  # Individuals are optionally ordered according to their names or their haplotype strings
  #  Ordering of the data.frame must occur prior to assigning the y-positions
  if(name_order) {
    genotype_data <- genotype_data[order(genotype_data$ind, decreasing = TRUE), ]
  } else {
    genotype_data <- genotype_data[order(genotype_data$haplotypes, decreasing = TRUE), ]
  }

  # Adding some metadata columns for plotting purposes
  genotype_data[["y_pos"]]  <- as.integer(factor(genotype_data[["ind"]], levels = unique(genotype_data[["ind"]])))
  genotype_data[["ref"]]    <- substring(genotype_data[["alleles"]], 1, 1)
  genotype_data[["alt"]]    <- substring(genotype_data[["alleles"]], 3, 3)
  # There should be an option in case the marker is heterozygous
  genotype_data[["allele"]] <- ifelse(genotype_data$genotype == "01", genotype_data[["ref"]],
                                      ifelse(genotype_data$genotype == "03", genotype_data[["alt"]],
                                             ifelse(genotype_data$genotype == "02", "HET", "N")))


  # Adding asterisks to selected markers if these are provided as input
  if(!is.null(kept_markers)) {
    genotype_data[["rs"]] <- ifelse(genotype_data[["rs"]] %in% kept_markers[["rs"]],
                                paste0(genotype_data[["rs"]], "**"),
                                genotype_data[["rs"]])
  }

  # Setting character size depending on the number of individuals or markers
  #  There might be a better default, but this temporarily does the job
  text_size <- ifelse(length(unique(genotype_data[["ind"]])) >= nrow(marker_data),
                      90 / length(unique(genotype_data[["ind"]])),
                      90 / nrow(marker_data))

  text_size <- min(text_size, 6)

  lab_size <- ifelse(length(unique(genotype_data[["ind"]])) >= nrow(marker_data),
                     450 / length(unique(genotype_data[["ind"]])),
                     450 / nrow(marker_data))

  lab_size <- min(lab_size, 12)

  # Creating the graphical object
  hplot <-
    # Initializing the plot with genotype_data
    ggplot2::ggplot(data = genotype_data,
                    mapping = ggplot2::aes_string(x = "x_pos", y = "y_pos")) +
    # Every genotype is plotted as a tile according to the marker (x) and individual (y)
    ggplot2::geom_tile(mapping = ggplot2::aes_string(fill = "genotype"),
                       height = 0.9, col = "black") +
    # The nucleotide at every position is added at the center of every tile
    ggplot2::geom_text(mapping = ggplot2::aes_string(label = "allele"),
                       size = text_size, vjust = "center") +
    # Controlling the display of marker names
    ggplot2::scale_x_continuous(name = "Marker",
                                breaks = unique(genotype_data$x_pos),
                                labels = unique(genotype_data$rs),
                                position = "top") +
    # Controlling the display of individual names
    ggplot2::scale_y_continuous(name = "Individual",
                                breaks = unique(genotype_data$y_pos),
                                labels = unique(genotype_data$ind)) +
    # Controlling the color of tiles
    ggplot2::scale_fill_manual(name = "Genotype",
                               values = c("00" = "red", "01" = "cornflowerblue",
                                          "02" = "pink", "03" = "darkorange"),
                               labels = c("00" = "Missing", "01" = "Major allele",
                                          "02" = "Heterozygote", "03" = "Minor allele")) +
    # Setting a baseline default theme
    ggplot2::theme_bw() +
    # Making some more adjustements to the theme
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, size = lab_size),
                   axis.text.y = ggplot2::element_text(size = lab_size), #temporarily removed
                   panel.grid = ggplot2::element_blank())

  # LAST ADDITIONS TO THE GRAPHICAL OBJECT

  # A red line indicating the central position of the gene is optionally added
  if(!is.null(gene_pos)) {
    # Rescaling the gene position for plotting the line at the right place
    rescaled_pos <- sum(gene_pos > marker_data[["pos"]]) + 0.5
    hplot <- hplot + ggplot2::geom_vline(xintercept = rescaled_pos, col = "red")
  }

  # The name of the haplotype to which the individual is assigned is optionally added
  # TO CORRECT : add every label only once per row
  if(!is.null(assignment)) {
    hplot <-
      hplot + ggplot2::geom_text(data = genotype_data[!duplicated(genotype_data$ind), ],
                                 x = max(genotype_data$x_pos) + 0.75,
                                 mapping = ggplot2::aes_string(label = "haplotypes"),
                                 size = text_size, vjust = "center")
  }

  return(hplot)
}
