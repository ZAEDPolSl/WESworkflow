#!/usr/bin/env bash
set -euo pipefail

CONFIG="${1:-config/example_config.yaml}"

if [[ ! -f "$CONFIG" ]]; then
	echo "ERROR: Config file not found: $CONFIG"
	exit 1
fi

echo "Checking command-line tools..."

check_cmd() {
	local cmd="$1"
	if command -v "$cmd" >/dev/null 2>&1; then
		echo "OK: $cmd"
	else
		echo "ERROR: Missing command: $cmd"
		exit 1
	fi
}

for cmd in fastqc multiqc trimmomatic bwa samtools bcftools bedtools picard CrossMap perl python Rscript; do
	check_cmd "$cmd"
done

echo
echo "Checking Python packages..."

python - <<'PY'
import importlib.util

packages = ["argparse", "re", "numpy", "pandas", "sklearn", "cyvcf2", "parc", "yaml"]
missing = [pkg for pkg in packages if importlib.util.find_spec(pkg) is None]

if missing:
	raise SystemExit("ERROR: Missing Python packages: " + ", ".join(missing))

print("OK: Python packages")
PY

echo
echo "Checking R packages..."

R_LIBS= R_LIBS_USER= R_LIBS_SITE= R_PROFILE_USER=/dev/null Rscript --vanilla - <<'RS'
.libPaths(file.path(Sys.getenv("CONDA_PREFIX"), "lib/R/library"))

packages <- c(
	"circlize", "colorspace", "ComplexHeatmap", "data.table", "dplyr",
	"future", "future.apply", "ggplot2", "ggrepel", "kBET", "lisi",
	"patchwork", "RANN", "RColorBrewer", "reshape2", "reticulate",
	"rnndescent", "stringr", "tibble", "tidyr", "uwot", "dpGMM", "remotes"
)

installed <- sapply(packages, requireNamespace, quietly = TRUE)

if (!all(installed)) {
	stop("Missing R packages: ", paste(packages[!installed], collapse = ", "))
}

cat("OK: R packages\n")
RS

echo
echo "Checking YAML config syntax..."

python - "$CONFIG" <<'PY'
import sys
import yaml

with open(sys.argv[1]) as f:
	config = yaml.safe_load(f)

required_sections = ["project", "environment", "tools", "resources", "parameters", "directories"]
missing = [section for section in required_sections if section not in config]

if missing:
	raise SystemExit("ERROR: Missing config sections: " + ", ".join(missing))

print("OK: YAML config syntax")
PY

echo
echo "Checking external paths from config..."

python - "$CONFIG" <<'PY'
import os
import sys
import yaml

config_path = sys.argv[1]

with open(config_path) as f:
	config = yaml.safe_load(f)

def flatten(d, prefix=""):
	out = {}
	for k, v in d.items():
		key = f"{prefix}.{k}" if prefix else k
		if isinstance(v, dict):
			out.update(flatten(v, key))
		else:
			out[key] = v
	return out

paths = flatten({k: config.get(k, {}) for k in ["tools", "resources", "directories"]})
errors = []
warnings = []

dir_keys = (
	"annovar_dir", "annovar_humandb", "reference_panel_dir",
	"raw_fastq", "input_bam", "output_dir", "temporary_dir"
)

for key, path in paths.items():
	if not isinstance(path, str):
		continue

	if path.startswith("/path/to"):
		warnings.append((key, path))
		continue

	if key.endswith(dir_keys):
		if not os.path.isdir(path):
			errors.append((key, path))
	else:
		if not os.path.exists(path):
			errors.append((key, path))

for key, path in warnings:
	print(f"WARNING: Placeholder path not configured: {key} = {path}")

if errors:
	for key, path in errors:
		print(f"ERROR: Missing path: {key} = {path}")
	raise SystemExit(1)

print("OK: External paths checked")
PY

echo
echo "Installation check finished."
