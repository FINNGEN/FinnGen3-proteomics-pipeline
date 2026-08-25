#!/usr/bin/env bash
# Build a FG3 release note: markdown -> LaTeX -> PDF, using the house preamble in
# docs/fg3_release_note.latex (matching the Batch 02 December-2025 release note conventions).
#
# Usage: build_release_note.sh <input.md> <title> <footer-left>
set -euo pipefail
MD="$1"; TITLE="$2"; FOOT="$3"
DIR="$(cd "$(dirname "$MD")" && pwd)"; STEM="$(basename "$MD" .md)"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TPL="$REPO/docs/fg3_release_note.latex"

# The markdown carries its own "# Title" and cover block; strip them so the LaTeX title page owns them.
awk 'NR==1 && /^# /{next} {print}' "$MD" > "$DIR/.$STEM.body.md"

pandoc "$DIR/.$STEM.body.md" \
  --from=markdown+pipe_tables+raw_tex \
  --to=latex \
  --template="$TPL" \
  --toc-depth=3 \
  --metadata title="$TITLE" \
  --metadata footer-left="$FOOT" \
  --output "$DIR/$STEM.tex"
echo "wrote $DIR/$STEM.tex"

( cd "$DIR" && tectonic -X compile "$STEM.tex" --outdir . >/dev/null 2>&1 ) \
  && echo "wrote $DIR/$STEM.pdf" \
  || { echo "tectonic failed; retrying with diagnostics"; ( cd "$DIR" && tectonic -X compile "$STEM.tex" --outdir . 2>&1 | tail -25 ); }
