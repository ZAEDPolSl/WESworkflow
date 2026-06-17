#!/usr/bin/env bash
set -euo pipefail

CONFIG="${1:-config/example_downstream_config.yaml}"
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

cd "$REPO_DIR"

RESULTS_DIR="$(python scripts/read_config.py "$CONFIG" directories.results_dir)"

if [[ "$RESULTS_DIR" != /* ]]; then
	RESULTS_DIR="$REPO_DIR/$RESULTS_DIR"
fi

FEATURE_SOURCE="$REPO_DIR/Data/example/features"
FEATURE_TARGET="$RESULTS_DIR/Features"
RESULTS_SUBDIR="$RESULTS_DIR/Results"
METADATA="$REPO_DIR/Data/example/sample_path_map_example.tsv"

if [[ ! -d "$FEATURE_SOURCE" ]]; then
	echo "ERROR: Example feature directory not found: $FEATURE_SOURCE"
	exit 1
fi

if [[ ! -f "$METADATA" ]]; then
	echo "ERROR: Example metadata file not found: $METADATA"
	exit 1
fi

rm -rf "$RESULTS_DIR"
mkdir -p "$FEATURE_TARGET" "$RESULTS_SUBDIR"

cp -r "$FEATURE_SOURCE"/* "$FEATURE_TARGET"/

echo "Checking example feature files..."

RESULTS_DIR="$RESULTS_DIR" python - <<'PY'
from pathlib import Path
import csv
import os

results_dir = Path(os.environ["RESULTS_DIR"])
feature_dir = results_dir / "Features"
results_subdir = results_dir / "Results"
metadata_file = Path("Data/example/sample_path_map_example.tsv")

features = sorted(feature_dir.glob("*/*.feature.txt"))
if not features:
	raise SystemExit("ERROR: No .feature.txt files found.")

with metadata_file.open() as f:
	reader = csv.DictReader(f, delimiter="\t")
	if not {"Sample", "Dataset", "Kit"}.issubset(reader.fieldnames or []):
		raise SystemExit("ERROR: Metadata must contain Sample, Dataset and Kit columns.")
	metadata_samples = {row["Sample"] for row in reader}

feature_samples = set()
feature_rows = []

for path in features:
	sample = path.name.replace(".feature.txt", "")
	chromosome = path.parent.name
	feature_samples.add(sample)

	with path.open() as f:
		reader = csv.DictReader(f, delimiter="\t")
		required = {"Gene", "CADD_weighted_avg_AF"}
		if not required.issubset(reader.fieldnames or []):
			raise SystemExit(f"ERROR: Missing required columns in {path}")
		n_rows = sum(1 for _ in reader)

	feature_rows.append((chromosome, sample, str(path), n_rows))

missing_in_metadata = feature_samples - metadata_samples
if missing_in_metadata:
	raise SystemExit(
		"ERROR: Samples missing in metadata: " + ", ".join(sorted(missing_in_metadata))
	)

report_file = results_subdir / "example_downstream_check.txt"
table_file = results_subdir / "example_feature_files.tsv"

with report_file.open("w") as f:
	f.write("Example downstream dry run report\n")
	f.write("=================================\n\n")
	f.write(f"Feature files checked: {len(features)}\n")
	f.write(f"Samples found in feature files: {len(feature_samples)}\n")
	f.write(f"Samples found in metadata: {len(metadata_samples)}\n")
	f.write("Required feature columns: Gene, CADD_weighted_avg_AF\n")
	f.write("Required metadata columns: Sample, Dataset, Kit\n")
	f.write("\nStatus: PASSED\n")

with table_file.open("w") as f:
	f.write("Chromosome\tSample\tFile\tRows\n")
	for chromosome, sample, path, n_rows in feature_rows:
		f.write(f"{chromosome}\t{sample}\t{path}\t{n_rows}\n")

print(f"OK: {len(features)} feature files checked.")
print(f"OK: {len(feature_samples)} samples found.")
print(f"Report written to: {report_file}")
print(f"Feature summary written to: {table_file}")
PY

echo "Example downstream dry run completed successfully."
echo "Example results directory: $RESULTS_DIR"
