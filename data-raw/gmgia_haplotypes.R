# Required packages
library(snpStats)
library(HaplotypeMiner)

# Storing the parameters of the analysis in an object
params <- haplo_params(
  input_file = system.file("extdata", "SNP_data.hmp.txt",
                           package = "HaplotypeMiner"),
  structure_file = system.file("extdata", "structure.txt",
                               package = "HaplotypeMiner"),
  kinship_file = system.file("extdata", "kinship.txt",
                             package = "HaplotypeMiner"),
  gene_db_file = system.file("extdata", "gene_db.txt",
                             package = "HaplotypeMiner"),
  chr_db_file = system.file("extdata", "gmax_chr_sizes.txt",
                            package = "HaplotypeMiner"),
  gene_name = "GmGia",
  R2_measure = "r2vs",
  cluster_R2 = "r2vs",
  max_missing_threshold = 0.6,
  max_het_threshold = 0.05,
  min_alt_threshold = 0.05,
  min_allele_count = 4,
  cluster_threshold = 0.8,
  max_marker_to_gene_distance = 250000,
  max_flanking_pair_distance = 250000,
  marker_independence_threshold = 0.5)

# Performing the computation
gmgia_haplotypes <- haplo_selection(params, verbose = TRUE)

# Saving as a dataset in package HaplotypeMiner
devtools::use_data(gmgia_haplotypes)
