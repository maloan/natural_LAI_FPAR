#!/usr/bin/env bash
# Build the water/ice-only LAI and FPAR trend products at 0.25 degrees.

set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo_dir=$(dirname "$script_dir")
cd "$repo_dir"

for command in gdal_translate cdo Rscript; do
  command -v "$command" >/dev/null || {
    echo "Missing dependency: $command" >&2
    exit 1
  }
done

mkdir -p "$repo_dir/analysis/unmasked/0p25"

for var in LAI FPAR; do
  input_dir="$repo_dir/output/nonvegetated_only_0p25/$var"
  output_dir="$repo_dir/analysis/unmasked/0p25"
  work_dir="$repo_dir/analysis/unmasked/tmp_$var"
  log_file="$repo_dir/analysis/unmasked/build_${var}.log"

  rm -rf "$work_dir"
  mkdir -p "$work_dir"

  mapfile -t input_files < <(
    find "$input_dir" -maxdepth 1 -type f \
      -name "${var}_*_0p25_masked_nonvegetated.tif" -print | sort
  )

  if [[ ${#input_files[@]} -ne 516 ]]; then
    echo "Expected 516 $var inputs, found ${#input_files[@]}" >&2
    exit 1
  fi

  : > "$log_file"
  echo "[$(date '+%F %T')] Converting $var" | tee -a "$log_file"

  monthly_files=()
  for input_file in "${input_files[@]}"; do
    yyyymm=$(basename "$input_file" | grep -oE '[0-9]{6}' | head -1)
    raw_nc="$work_dir/${var}_${yyyymm}_raw.nc"
    dated_nc="$work_dir/${var}_${yyyymm}_dated.nc"

    gdal_translate -of NetCDF "$input_file" "$raw_nc" >>"$log_file" 2>&1
    cdo -O settaxis,"${yyyymm:0:4}-${yyyymm:4:2}-01",00:00:00,1month \
      "$raw_nc" "$dated_nc" >>"$log_file" 2>&1
    rm -f "$raw_nc"
    monthly_files+=("$dated_nc")
  done

  monthly="$output_dir/${var}_georef_monthly_0p25.nc"
  cdo -O mergetime "${monthly_files[@]}" "$monthly" >>"$log_file" 2>&1

  n_time=$(cdo -s ntime "$monthly")
  if [[ "$n_time" -ne 516 ]]; then
    echo "Expected 516 monthly timesteps for $var, found $n_time" >&2
    exit 1
  fi

  for metric in yearmean yearmax yearmin; do
    cdo -O "$metric" "$monthly" \
      "$output_dir/${var}_georef_${metric}_0p25.nc" >>"$log_file" 2>&1
  done

  cdo -O sub \
    "$output_dir/${var}_georef_yearmax_0p25.nc" \
    "$output_dir/${var}_georef_yearmin_0p25.nc" \
    "$output_dir/${var}_georef_yearamp_0p25.nc" >>"$log_file" 2>&1

  for metric in yearmean yearmax yearmin yearamp; do
    annual="$output_dir/${var}_georef_${metric}_0p25.nc"
    intercept="$output_dir/${var}_georef_${metric}_trend_intercept_0p25.nc"
    slope="$output_dir/${var}_georef_${metric}_trend_slope_peryear_0p25.nc"
    cdo -O trend "$annual" "$intercept" "$slope" >>"$log_file" 2>&1
  done

  mk_pids=()
  for metric in yearmean yearmax yearmin yearamp; do
    Rscript "$repo_dir/trends/compute_mk_pval.R" unmasked "$var" "$metric" \
      >>"$log_file" 2>&1 &
    mk_pids+=("$!")
  done

  mk_failed=0
  for pid in "${mk_pids[@]}"; do
    wait "$pid" || mk_failed=$((mk_failed + 1))
  done
  if [[ "$mk_failed" -gt 0 ]]; then
    echo "Mann-Kendall failed for $mk_failed $var metric(s)" >&2
    exit 1
  fi

  eps=0.02
  [[ "$var" == LAI ]] && eps=0.05

  for metric in yearmean yearmax yearmin yearamp; do
    annual="$output_dir/${var}_georef_${metric}_0p25.nc"
    slope="$output_dir/${var}_georef_${metric}_trend_slope_peryear_0p25.nc"
    mean_file="$work_dir/${var}_${metric}_mean.nc"
    valid_file="$work_dir/${var}_${metric}_valid.nc"
    relative="$output_dir/${var}_georef_${metric}_trend_relative_percent_peryear_0p25.nc"
    metric_eps="$eps"
    [[ "$metric" == yearamp ]] && metric_eps=0.01

    cdo -O timmean "$annual" "$mean_file" >>"$log_file" 2>&1
    cdo -O gec,"$metric_eps" "$mean_file" "$valid_file" >>"$log_file" 2>&1
    cdo -O setname,trend_relative_percent \
      -ifthen "$valid_file" -mulc,100 -div "$slope" "$mean_file" \
      "$relative" >>"$log_file" 2>&1
  done

  rm -rf "$work_dir"
  echo "[$(date '+%F %T')] Completed $var" | tee -a "$log_file"
done

echo "Water/ice-only trend products completed."
