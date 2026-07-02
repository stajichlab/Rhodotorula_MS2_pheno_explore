# Pipeline Logs

This directory contains execution logs from pipeline runs.

## File Naming Convention

- `pipeline_YYYYMMDD_HHMMSS.log` — Main pipeline execution log
- `phase0_YYYYMMDD_HHMMSS.log` — Phase 0 detailed log
- `phase1_YYYYMMDD_HHMMSS.log` — Phase 1 detailed log
- `phase2_YYYYMMDD_HHMMSS.log` — Phase 2 detailed log
- `phase3_YYYYMMDD_HHMMSS.log` — Phase 3 detailed log
- `pipeline_manifest_YYYYMMDD_HHMMSS.json` — Execution manifest (git commit, timestamps)

## Running the Pipeline with Logging

```bash
bash scripts/run_pipeline.sh
```

This creates timestamped logs for each run and generates a validation manifest.

## Log Inspection

View the main pipeline log:
```bash
tail -f logs/pipeline_*.log
```

Check a specific phase:
```bash
cat logs/phase2_*.log | tail -50
```

View validation manifest:
```bash
cat logs/pipeline_manifest_*.json
```

## Cleaning Up Old Logs

Logs are retained for reproducibility but can be archived:
```bash
# Compress logs older than 30 days
find logs -name "*.log" -mtime +30 -exec gzip {} \;

# Remove logs older than 1 year
find logs -name "*.log.gz" -mtime +365 -delete
```

