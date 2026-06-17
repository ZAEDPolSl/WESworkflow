conda_lib <- file.path(Sys.getenv("CONDA_PREFIX"), "lib/R/library")
if (!dir.exists(conda_lib)) stop("Conda R library not found: ", conda_lib)

.libPaths(conda_lib)
options(repos = c(CRAN = "https://cloud.r-project.org"))

cat("Using R library paths:\n")
print(.libPaths())

cran_packages <- c("dpGMM", "rnndescent")

for (pkg in cran_packages) {
	if (!requireNamespace(pkg, quietly = TRUE)) {
		install.packages(pkg, lib = conda_lib)
	}
}

if (!requireNamespace("remotes", quietly = TRUE)) {
	stop("Package 'remotes' is missing. It should be installed from envs/wes-env.yml.")
}

if (!requireNamespace("kBET", quietly = TRUE)) {
	remotes::install_github("theislab/kBET", lib = conda_lib, upgrade = "never", dependencies = FALSE)
}

cat("Additional R package installation finished.\n")
