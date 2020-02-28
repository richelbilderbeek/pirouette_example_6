#
# Standard
#

library(pirouette)
suppressMessages(library(ggplot2))

################################################################################
# Constants
################################################################################
is_testing <- is_on_travis()
example_no <- 6
rng_seed <- 314
folder_name <- paste0("example_", example_no, "_", rng_seed)

################################################################################
# Create phylogeny
################################################################################
set.seed(rng_seed)
phylogeny <- create_yule_tree(n_taxa = 6, crown_age = 10)

################################################################################
# Setup pirouette
################################################################################
pir_params <- create_std_pir_params(folder_name = folder_name)
if (is_testing) {
  pir_params <- shorten_pir_params(pir_params)
}

pir_out <- pir_run(
  phylogeny,
  pir_params = pir_params
)

pir_save(
  phylogeny = phylogeny,
  pir_out = pir_out,
  folder_name = folder_name
)
