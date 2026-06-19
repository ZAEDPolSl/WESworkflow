#!/usr/bin/env bash
set -euo pipefail

CONFIG="${1:-config/local_config.yaml}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ -f "$SCRIPT_DIR/scripts/read_config.py" ]]; then
	REPO_DIR="$SCRIPT_DIR"
elif [[ -f "$SCRIPT_DIR/../scripts/read_config.py" ]]; then
	REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
else
	REPO_DIR="$(pwd)"
fi

READ_CONFIG="$REPO_DIR/scripts/read_config.py"

if [[ "$CONFIG" != /* ]]; then
	CONFIG="$REPO_DIR/$CONFIG"
fi

if [[ ! -f "$CONFIG" ]]; then
	echo "ERROR: User config file not found: $CONFIG"
	echo "Create it with: cp config/example_config.yaml config/local_config.yaml"
	exit 1
fi

if [[ ! -f "$READ_CONFIG" ]]; then
	echo "ERROR: read_config.py not found: $READ_CONFIG"
	exit 1
fi

resolve_path() {
	local path="$1"

	if [[ "$path" == /* ]]; then
		echo "$path"
	else
		echo "$REPO_DIR/$path"
	fi
}

read_cfg() {
	python "$READ_CONFIG" "$CONFIG" "$1"
}

check_cmd() {
	local cmd="$1"

	if command -v "$cmd" >/dev/null 2>&1; then
		echo "OK: $cmd"
	else
		echo "ERROR: Missing command: $cmd"
		exit 1
	fi
}

check_file_from_config() {
	local key="$1"
	local path

	path="$(read_cfg "$key")"

	if [[ "$path" == /path/to* || -z "$path" || "$path" == "None" ]]; then
		echo "WARNING: Placeholder or empty config value: $key = $path"
		return 0
	fi

	path="$(resolve_path "$path")"

	if [[ -f "$path" ]]; then
		echo "OK: $key = $path"
	else
		echo "ERROR: Missing file configured in $key: $path"
		exit 1
	fi
}

check_dir_from_config() {
	local key="$1"
	local path

	path="$(read_cfg "$key")"

	if [[ "$path" == /path/to* || -z "$path" || "$path" == "None" ]]; then
		echo "WARNING: Placeholder or empty config value: $key = $path"
		return 0
	fi

	path="$(resolve_path "$path")"

	if [[ -d "$path" ]]; then
		echo "OK: $key = $path"
	else
		echo "ERROR: Missing directory configured in $key: $path"
		exit 1
	fi
}

echo "Checking command-line tools..."

for cmd in fastqc multiqc trimmomatic bwa samtools bcftools bedtools bgzip tabix picard CrossMap perl python Rscript java docker parallel glnexus_cli; do
	check_cmd "$cmd"
done

echo
echo "Checking Python packages..."

python - <<'PY'
import importlib.util

packages = ["numpy", "pandas", "sklearn", "cyvcf2", "parc", "yaml"]
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
echo "Installation check finished."