#!/bin/bash
#
# Master pipeline orchestration script
# Runs all phases in sequence with logging and validation
#

set -e  # Exit on error

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SCRIPTS_DIR="$PROJECT_DIR/scripts"
LOGS_DIR="$PROJECT_DIR/logs"
RESULTS_DIR="$PROJECT_DIR/results"

# Create logs directory if not exists
mkdir -p "$LOGS_DIR"

# Timestamp for this run
RUN_ID=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="$LOGS_DIR/pipeline_${RUN_ID}.log"
MANIFEST_FILE="$LOGS_DIR/pipeline_manifest_${RUN_ID}.json"

echo "============================================" | tee "$LOG_FILE"
echo "Rhodotorula MS2 Analysis Pipeline" | tee -a "$LOG_FILE"
echo "Start time: $(date)" | tee -a "$LOG_FILE"
echo "Run ID: $RUN_ID" | tee -a "$LOG_FILE"
echo "============================================" | tee -a "$LOG_FILE"

# Save pipeline config
cat > "$MANIFEST_FILE" << JSON
{
  "run_id": "$RUN_ID",
  "start_time": "$(date -Is)",
  "project_dir": "$PROJECT_DIR",
  "git_commit": "$(cd $PROJECT_DIR && git rev-parse HEAD 2>/dev/null || echo 'unknown')",
  "phases": {}
}
JSON

# Function to run a phase
run_phase() {
  local phase_num=$1
  local script=$2
  
  echo "" | tee -a "$LOG_FILE"
  echo "======== Phase $phase_num: $(date) ========" | tee -a "$LOG_FILE"
  
  local phase_log="$LOGS_DIR/phase${phase_num}_${RUN_ID}.log"
  
  if [ -f "$SCRIPTS_DIR/$script" ]; then
    echo "Running: python3 $script" | tee -a "$LOG_FILE"
    
    if python3 "$SCRIPTS_DIR/$script" > "$phase_log" 2>&1; then
      echo "✓ Phase $phase_num completed" | tee -a "$LOG_FILE"
      tail -5 "$phase_log" | tee -a "$LOG_FILE"
    else
      echo "✗ Phase $phase_num FAILED" | tee -a "$LOG_FILE"
      tail -20 "$phase_log" | tee -a "$LOG_FILE"
      exit 1
    fi
  else
    echo "✗ Script not found: $script" | tee -a "$LOG_FILE"
    exit 1
  fi
}

# Run all phases
run_phase 0 "phase0_batch_assessment.py"
run_phase 1 "phase1_feature_filtering.py"
run_phase 2 "phase2_correlation_analysis.py"
run_phase 3 "phase3_stratified_species_analysis.py"

echo "" | tee -a "$LOG_FILE"
echo "======== Validation: $(date) ========" | tee -a "$LOG_FILE"
python3 "$SCRIPTS_DIR/validate_outputs.py" | tee -a "$LOG_FILE"

echo "" | tee -a "$LOG_FILE"
echo "======== Generating visualizations ========" | tee -a "$LOG_FILE"
if python3 "$SCRIPTS_DIR/phase3_visualizations.py" >> "$LOG_FILE" 2>&1; then
  echo "✓ Visualizations generated" | tee -a "$LOG_FILE"
else
  echo "✗ Visualization generation failed (non-blocking)" | tee -a "$LOG_FILE"
fi

echo "" | tee -a "$LOG_FILE"
echo "============================================" | tee -a "$LOG_FILE"
echo "✓ Pipeline complete: $(date)" | tee -a "$LOG_FILE"
echo "============================================" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"
echo "Logs saved to: $LOGS_DIR/" | tee -a "$LOG_FILE"
echo "Results in: $RESULTS_DIR/" | tee -a "$LOG_FILE"

