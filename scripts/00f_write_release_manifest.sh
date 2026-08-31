#!/usr/bin/env bash
set -euo pipefail
REL="/mnt/longGWAS_disk_100GB/long_gwas/11.fg3_Proteomics/08.genewiz_batch3/01.QCed_Batch03_Release_Aug2026"
cd "$REL"
V="v.04"
OUT="MANIFEST_${V}.txt"
{
  echo "FG3 Batch 03 Olink proteomics — release package manifest"
  echo "Version:   $V"
  echo "Generated: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "Path:      $REL"
  echo
  echo "Every file that constitutes the delivery. Anything present in this directory and absent"
  echo "from this list is not part of the release."
  echo
  printf "%-72s %10s %34s\n" "FILE" "BYTES" "MD5"
  printf "%-72s %10s %34s\n" "$(printf '%.0s-' {1..72})" "----------" "$(printf '%.0s-' {1..34})"
  # find rather than ls, so that the release-note figures in figures/ are covered too
  for f in $(find . -type f -not -name 'MANIFEST_*' -printf '%P\n' | sort); do
    printf "%-72s %10d %34s\n" "$f" "$(stat -c%s "$f")" "$(md5sum "$f" | cut -d' ' -f1)"
  done
  echo
  stale=$(ls -1 MANIFEST_* 2>/dev/null | grep -v "^${OUT}$" || true)
  if [ -n "$stale" ]; then
    echo "WARNING - superseded manifest(s) present. Only ${OUT} describes this package;"
    echo "the following are from earlier versions and must be removed before handover:"
    echo "$stale" | sed 's/^/  /'
    echo
  fi
  hid=$(find . -name '.*' -not -name '.' -printf '%P\n' || true)
  if [ -n "$hid" ]; then
    echo "WARNING - hidden entries present in the release directory. These are NOT part of the"
    echo "delivery and must be removed before handover:"
    echo "$hid" | sed 's/^/  /'
    echo
  fi
  echo "Row and column counts for the tabular artefacts:"
  for f in FG3_batch03_delivery_metadata_fg3_batch_03.tsv \
           qc_annotated_metadata_all_6191_samples_fg3_batch_03.tsv \
           comprehensive_outliers_list_fg3_batch_03.tsv \
           FG3_batch03_sample_metadata_schema_fg3_batch_03.tsv; do
    [ -f "$f" ] || continue
    r=$(( $(wc -l < "$f") - 1 ))
    c=$(head -1 "$f" | tr '\t' '\n' | wc -l)
    printf "  %-64s %6d rows x %4d cols\n" "$f" "$r" "$c"
  done
  echo
  echo "NPX matrices (rows = samples, columns = assays):"
  for f in npx_matrix_*.tsv; do
    [ -f "$f" ] || continue
    c=$(head -1 "$f" | tr '\t' '\n' | wc -l)
    r=$(( $(wc -l < "$f") - 1 ))
    printf "  %-64s %6d rows x %5d cols\n" "$f" "$r" "$c"
  done
} > "$OUT"
echo "wrote $OUT"
grep -c . "$OUT"
