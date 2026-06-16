#!/usr/bin/env python3
import sys
import yaml

if len(sys.argv) != 3:
	raise SystemExit("Usage: read_config.py <config.yaml> <dot.separated.key>")

with open(sys.argv[1]) as f:
	config = yaml.safe_load(f)

value = config
for key in sys.argv[2].split("."):
	value = value[key]

print(value)
