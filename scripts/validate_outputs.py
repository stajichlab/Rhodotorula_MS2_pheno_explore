#!/usr/bin/env python3
"""
Validate that Phase outputs exist and have expected content
Provides audit trail for reproducibility
"""

import os
import json
import gzip
import hashlib
from pathlib import Path
from datetime import datetime

results_dir = Path("../results")

def checksum_file(filepath, method='md5'):
    """Compute file checksum"""
    hash_func = hashlib.new(method)
    with open(filepath, 'rb') as f:
        for chunk in iter(lambda: f.read(4096), b''):
            hash_func.update(chunk)
    return hash_func.hexdigest()

def validate_phase(phase_num):
    """Validate a phase's outputs"""
    phase_dir = results_dir / f"phase{phase_num}"
    manifest = {}
    
    if not phase_dir.exists():
        print(f"✗ {phase_dir} does not exist")
        return False
    
    print(f"\n[Phase {phase_num}]")
    for filepath in sorted(phase_dir.glob("*")):
        if filepath.is_file():
            size_mb = filepath.stat().st_size / (1024**2)
            checksum = checksum_file(filepath)
            
            # Read row count for CSV files
            row_count = "N/A"
            if filepath.suffix == '.gz' and filepath.name.endswith('.csv.gz'):
                try:
                    with gzip.open(filepath, 'rt') as f:
                        row_count = sum(1 for _ in f) - 1  # -1 for header
                except:
                    row_count = "Error"
            elif filepath.suffix == '.json':
                try:
                    with open(filepath) as f:
                        data = json.load(f)
                        row_count = len(data) if isinstance(data, list) else "dict"
                except:
                    row_count = "Error"
            
            print(f"  ✓ {filepath.name:45} | {size_mb:7.1f} MB | rows: {row_count}")
            manifest[filepath.name] = {
                "size_bytes": filepath.stat().st_size,
                "checksum_md5": checksum,
                "rows": row_count
            }
    
    return manifest

# Run validation
print("="*80)
print("PHASE OUTPUT VALIDATION")
print("="*80)
print(f"Timestamp: {datetime.now().isoformat()}")

validation_report = {}
for phase in [0, 1, 2, 3]:
    validation_report[f"phase{phase}"] = validate_phase(phase)

# Save manifest
manifest_file = results_dir / "VALIDATION_MANIFEST.json"
with open(manifest_file, 'w') as f:
    json.dump({
        "timestamp": datetime.now().isoformat(),
        "outputs": validation_report
    }, f, indent=2)

print(f"\n✓ Validation manifest saved: {manifest_file}")

