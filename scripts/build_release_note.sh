#!/usr/bin/env bash
# Build a FG3 release note: markdown -> LaTeX -> PDF, using the house preamble and cover page in
# docs/fg3_release_note.latex. The output matches the conventions of the FG3 cross-batch harmonised
# release note, the CKD progression report and the ATC residual-variance report: a centred cover page
# with a right-aligned metadata tabular and a grey scope box, an author/Broad footer on every page,
# and teal URLs.
#
# The markdown is the single source and must open like this:
#
#     # <document title>                     <- \LARGE\bfseries, " | " splits lines
#     <blank>
#     *<subtitle>*                           <- optional, \Large, " | " splits lines
#     <blank>
#     _<version line>_                       <- optional, \normalsize
#     <blank>
#     **Release Date**: August 2026          <- cover metadata rows, one bold key per line
#     **Platform**: ...
#     <blank>
#     > Scope. ...                           <- optional blockquote, becomes the grey tcolorbox
#     <blank>
#     ---                                    <- first horizontal rule closes the cover block
#
# Everything after that rule is the body. Body headings are shifted up one level so the markdown's
# "##" sections become LaTeX \section, matching the reference documents.
#
# Usage: build_release_note.sh <input.md> <title-fallback> <unused, kept for call compatibility>
set -euo pipefail
MD="$1"; TITLE_FALLBACK="${2:-}"
DIR="$(cd "$(dirname "$MD")" && pwd)"; STEM="$(basename "$MD" .md)"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TPL="$REPO/docs/fg3_release_note.latex"

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT           # a real temp dir, so nothing is left beside the release

# md2tex <file>  : markdown to LaTeX, flattened to one line (single-line cover fields)
# mdblock <file>  : markdown to LaTeX, paragraph breaks preserved (the scope box)
md2tex() { pandoc "$1" --from=markdown --to=latex | sed -e 's/[[:space:]]*$//' | tr '\n' ' ' | sed -e 's/[[:space:]]*$//'; }
mdblock() { pandoc "$1" --from=markdown --to=latex; }
str2tex() { printf '%s\n' "$1" > "$WORK/frag.md"; md2tex "$WORK/frag.md"; }

# Guard: a non-ASCII glyph with no \newunicodechar declaration is silently mangled or dropped by
# the T1 text fonts (± came out as "s-acute", rho vanished). Fail loudly instead.
UNDECLARED="$(python3 - "$MD" "$TPL" <<'PY'
import io, re, sys
md  = io.open(sys.argv[1], encoding="utf-8").read()
tpl = io.open(sys.argv[2], encoding="utf-8").read()
declared = set(re.findall(r'\\newunicodechar\{(.)\}', tpl))
bad = sorted({c for c in md if ord(c) > 127} - declared - set("\u2018\u2019\u201c\u201d\u2013\u2014\u2026\u00a0"))
print(" ".join(f"{c}(U+{ord(c):04X})" for c in bad))
PY
)"
if [ -n "$UNDECLARED" ]; then
  echo "ERROR: glyphs used in $MD with no \\newunicodechar in $TPL:" >&2
  echo "       $UNDECLARED" >&2
  echo "       Declare them in the template, or they will render wrongly." >&2
  exit 1
fi

# Guard: a release note is read outside the FinnGen Sandbox, so no participant or tube identifier
# belongs in it. Deny-by-default; a token that must stay goes in docs/identifier_allowlist.txt with
# a recorded reason. This is a hard gate -- it runs before any LaTeX is produced.
python3 "$(dirname "${BASH_SOURCE[0]}")/check_no_identifiers.py" "$MD" || {
  echo "ERROR: identifier guard failed for $MD -- refusing to build." >&2
  exit 1
}

COVER_MD="$WORK/cover.md"; BODY_MD="$WORK/body.md"; SCOPE_MD="$WORK/scope.md"
: > "$COVER_MD"; : > "$BODY_MD"; : > "$SCOPE_MD"

TITLE=""; SUBTITLE=""; VERSION=""
ROWS="$WORK/rows.tex"; : > "$ROWS"

seen_rule=0
while IFS= read -r line || [ -n "$line" ]; do
  if [ "$seen_rule" -eq 1 ]; then printf '%s\n' "$line" >> "$BODY_MD"; continue; fi
  case "$line" in
    '# '*)                       TITLE="${line#\# }" ;;
    '---'|'---'[[:space:]]*)     seen_rule=1 ;;
    '> '*)                       printf '%s\n' "${line#> }" >> "$SCOPE_MD" ;;
    '>')                         printf '\n' >> "$SCOPE_MD" ;;
    '**'*'**:'*)                 printf '%s\n' "$line" >> "$COVER_MD" ;;
    '*'*'*')                     if [ -z "$SUBTITLE" ]; then SUBTITLE="${line#\*}"; SUBTITLE="${SUBTITLE%\*}"; fi ;;
    '_'*'_')                     if [ -z "$VERSION" ]; then VERSION="${line#_}"; VERSION="${VERSION%_}"; fi ;;
    *)                           : ;;   # blank lines and stray prose in the cover block are dropped
  esac
done < "$MD"

[ -n "$TITLE" ] || TITLE="$TITLE_FALLBACK"

# The cover metadata rows: "**Key**: value" becomes "\textbf{Key} & value \\[4pt]".
while IFS= read -r line; do
  key="${line#\*\*}"; key="${key%%\*\*:*}"
  val="${line#*\*\*: }"
  printf '\\textbf{%s} & %s \\\\[4pt]\n' "$(str2tex "$key")" "$(str2tex "$val")" >> "$ROWS"
done < "$COVER_MD"
# Drop the trailing row separator so the tabular does not end with an empty ruled row
sed -i '$ s/ \\\\\[4pt\]$//' "$ROWS"

# A literal " | " in the title or subtitle is a line break on the cover page. It is written that way
# in the markdown rather than as raw LaTeX, so that pandoc does not escape the backslashes.
COVER_TITLE="$(str2tex "$TITLE" | sed 's/ \\textbar{} /\\\\[0.3em] /g')"
COVER_SUBTITLE=""; [ -n "$SUBTITLE" ] && COVER_SUBTITLE="$(str2tex "$SUBTITLE" | sed 's/ \\textbar{} /\\\\[0.2em] /g')"
COVER_VERSION="";  [ -n "$VERSION" ]  && COVER_VERSION="$(str2tex "$VERSION")"
COVER_SCOPE="";    [ -s "$SCOPE_MD" ] && COVER_SCOPE="$(mdblock "$SCOPE_MD")"
COVER_ROWS="$(cat "$ROWS")"

pandoc "$BODY_MD" \
  --from=markdown+pipe_tables+raw_tex \
  --to=latex \
  --template="$TPL" \
  --shift-heading-level-by=-1 \
  --toc-depth=2 \
  --metadata title="$TITLE" \
  --variable covertitle="$COVER_TITLE" \
  --variable coversubtitle="$COVER_SUBTITLE" \
  --variable coverversion="$COVER_VERSION" \
  --variable coverrows="$COVER_ROWS" \
  --variable coverscope="$COVER_SCOPE" \
  --output "$DIR/$STEM.tex"
echo "wrote $DIR/$STEM.tex"

( cd "$DIR" && tectonic -X compile "$STEM.tex" --outdir . >/dev/null 2>&1 ) \
  && echo "wrote $DIR/$STEM.pdf" \
  || { echo "tectonic failed; retrying with diagnostics"; ( cd "$DIR" && tectonic -X compile "$STEM.tex" --outdir . 2>&1 | tail -30 ); }
