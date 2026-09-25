#!/usr/bin/env bash
set -euo pipefail
REL="/mnt/longGWAS_disk_100GB/long_gwas/11.fg3_Proteomics/08.genewiz_batch3/01.QCed_Batch03_Release_Aug2026"
cd "$REL"
V="v.07"
OUT="MANIFEST_${V}.txt"
PKG="FG3_Batch_3_QCed_release/${V}"
# Files present in the release directory that are BUILD INPUTS, not delivered artefacts, and are
# therefore excluded from the distribution to the release bucket. The manifest must describe the
# distributed package only: listing a build input makes every checksum verification against the
# bucket report a miss for a file that was never meant to be there.
NOT_DISTRIBUTED='^figures/|\.tex$'
{
  echo "FG3 Batch 03 Olink proteomics — release package manifest"
  echo "Version:   $V"
  echo "Generated: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "Package:   $PKG"
  echo
  echo "Every file distributed as part of this release, with its size, its md5 checksum and, for"
  echo "the tabular artefacts, its row and column counts. Anything absent from this list is not"
  echo "part of the release. The list is complete against the release bucket: the manifest itself"
  echo "is the only object in the package that it does not describe."
  echo
  printf "%-72s %10s %34s\n" "FILE" "BYTES" "MD5"
  printf "%-72s %10s %34s\n" "$(printf '%.0s-' {1..72})" "----------" "$(printf '%.0s-' {1..34})"
  n_files=0
  while IFS= read -r f; do
    printf "%-72s %10d %34s\n" "$f" "$(stat -c%s "$f")" "$(md5sum "$f" | cut -d' ' -f1)"
    n_files=$((n_files + 1))
  done < <(find . -type f -not -name 'MANIFEST_*' -printf '%P\n' | grep -Ev "$NOT_DISTRIBUTED" | sort)
  echo
  echo "Distributed files: ${n_files}, plus this manifest."
  echo
  # Record the exclusions explicitly, so that a recipient who sees them referenced elsewhere knows
  # they were withheld by design rather than lost in transit.
  excl=$(find . -type f -not -name 'MANIFEST_*' -printf '%P\n' | grep -E "$NOT_DISTRIBUTED" | sort || true)
  if [ -n "$excl" ]; then
    echo "Present in the build directory and deliberately NOT distributed (source and build inputs"
    echo "for the release note; the note itself ships as .md and .pdf):"
    echo "$excl" | sed 's/^/  /'
    echo
  fi
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
