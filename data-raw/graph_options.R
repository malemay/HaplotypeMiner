# Generating objects used internally by package HaplotypeMiner
.null_graphs <- list("All_markers" = character(0),
                      "Filtered_markers" = character(0),
                      "Clustered_markers" = character(0),
                      "Selected_clusters" = character(0),
                      "Selected_markers" = character(0),
                      "Haplotypes" = character(0))

.default_graphs <- list("All_markers" = character(0),
                        "Filtered_markers" = character(0),
                        "Clustered_markers" = character(0),
                        "Selected_clusters" = character(0),
                        "Selected_markers" = character(0),
                        "Haplotypes" = character(0))

.input_options <- c("All_markers", "Filtered_markers", "Clustered_markers",
                    "Selected_clusters", "Selected_markers",
                    "Haplotypes")

.output_options <- c("density", "matrix", "distance", "genotypes")

devtools::use_data(.null_graphs, .default_graphs,
                   .input_options, .output_options,
                   internal = TRUE, overwrite = TRUE)
