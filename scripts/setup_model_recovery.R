## -----------------------------------------------------------------------------
## Author: Francisco Garre-Frutos
## -----------------------------------------------------------------------------

# Packages and CmdStan needed to run scripts/model_recovery.R on a new machine

options(timeout = 1800) # CmdStan download (~45 MB)

if (identical(unname(getOption("repos")["CRAN"]), "@CRAN@")) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

packages <- c("pacman", "posterior", "loo", "dplyr", "tidyr", "purrr", "tibble",
              "ggplot2", "ggdist", "ggthemes", "future", "furrr", "Rmisc")
missing <- setdiff(packages, rownames(installed.packages()))
if (length(missing) > 0) {
  install.packages(missing)
}

# cmdstanr is not on CRAN
if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
}

# Same CmdStan version as the main fits
cmdstanr::check_cmdstan_toolchain(fix = TRUE)
cmdstanr::install_cmdstan(version = "2.40.0", cores = parallel::detectCores())

cat("cmdstanr", as.character(packageVersion("cmdstanr")), "/ CmdStan",
    cmdstanr::cmdstan_version(), "\n")
