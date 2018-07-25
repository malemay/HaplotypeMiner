#' Generates a tile plot of the LD of the given snp_data.
#'
#' Description
#'
#' Details
#'
#' @param snp_data The snp_data to be plotted.
#' @param center_pos The central position on either side of which markers must be in LD.
#' @param kept_markers To be completed.
#'
#' @return To be completed
#'
#' @export
#' @examples
#' NULL
#'
ld_plot <- function(snp_data, center_pos, kept_markers = NULL) {

  melt_df <- reshape2::melt(as.matrix(Matrix::forceSymmetric(snp_data$LD)),
                            varnames = c("SNP_1", "SNP_2"), value.name = "Rsquared")

  melt_df$SNP_1 <- factor(melt_df$SNP_1, levels = snp_data$Markers$rs)
  melt_df$SNP_2 <- factor(melt_df$SNP_2, levels = rev(snp_data$Markers$rs))

  if(!is.null(kept_markers)) {
    SNP_1_levels <- ifelse(levels(melt_df$SNP_1) %in% kept_markers$rs,
                           paste0("**", levels(melt_df$SNP_1)),
                           levels(melt_df$SNP_1))

    SNP_2_levels <- ifelse(levels(melt_df$SNP_2) %in% kept_markers$rs,
                           paste0("**", levels(melt_df$SNP_2)),
                           levels(melt_df$SNP_2))

    melt_df$SNP_1 <- factor(ifelse(melt_df$SNP_1 %in% kept_markers$rs,
                                   paste0("**", melt_df$SNP_1),
                                   as.character(melt_df$SNP_1)),
                            levels = SNP_1_levels)

    melt_df$SNP_2 <- factor(ifelse(melt_df$SNP_2 %in% kept_markers$rs,
                                   paste0("**", melt_df$SNP_2),
                                   as.character(melt_df$SNP_2)),
                            levels = SNP_2_levels)
  }

  #Get middle position
  middle_x <- sum(snp_data$Markers$pos <= center_pos) + 0.5
  middle_y <- nrow(snp_data$Markers) - middle_x + 1
  melt_df$Rsquared[melt_df$SNP_1 == melt_df$SNP_2] <- NA

  ld_plot <-
    ggplot2::ggplot(melt_df) +
    ggplot2::geom_tile(mapping = ggplot2::aes_string(x = "SNP_1", y = "SNP_2", fill = "Rsquared"),
                       color = "grey") +
    ggplot2::geom_hline(yintercept = middle_y, size = 1) +
    ggplot2::geom_vline(xintercept = middle_x, size = 1) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ggplot2::scale_fill_distiller(name = r2expr(snp_data[["LD_measure"]]),
                                  type = "seq", palette = "YlOrRd",
                                  direction = 1, limits = c(0,1.001)) +
    ggplot2::xlab("Marker 1") + ggplot2::ylab("Marker 2")

  return(ld_plot)
}
