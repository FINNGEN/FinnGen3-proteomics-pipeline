#!/usr/bin/env python3
# ==============================================================================================
# Pre-release identifier guard.
#
# Fails the build if a participant- or sample-level identifier reaches a document that leaves the
# FinnGen Sandbox. A release note and its README are read by people who are not inside the Sandbox
# perimeter, so no identifier token belongs in either, however pseudonymous it is.
#
# The guard is deliberately a DENY-BY-DEFAULT list with an explicit allowlist. A token that must
# stay has to be written into docs/identifier_allowlist.txt together with the reason it is there.
# That converts "nobody noticed" into "somebody signed for it", which is the failure this guard
# exists to prevent.
#
# Usage: check_no_identifiers.py <file> [<file> ...]
#        check_no_identifiers.py --allowlist <path> <file> ...
# ==============================================================================================

import io
import os
import re
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_ALLOWLIST = os.path.join(REPO, "docs", "identifier_allowlist.txt")

# Token families. Each is (name, compiled regex, what it identifies).
FAMILIES = [
    ("FINNGENID",
     re.compile(r"\bFG[A-Z0-9]{8}\b"),
     "FinnGen participant pseudonym"),
    ("SAMPLE_ID",
     re.compile(r"\bP[0-9]{20}(?:-rep)?\b"),
     "tube-level sample pseudonym"),
    ("BROAD_COLLECTION_ID",
     re.compile(r"\bP1516_[0-9]+\b"),
     "vendor sub-project collection identifier"),
    ("PLATE_ID",
     re.compile(r"\b90-[0-9]{10}[_-][AB][0-9]+[_-][A-Z]+\b"),
     "laboratory plate identifier"),
    ("CONTAINER_BARCODE",
     re.compile(r"\b1000113[0-9]{5}\b"),
     "plasma container barcode"),
]

# Language that describes exactly one participant. A regex cannot catch a disclosure made by
# arithmetic -- v.06 never stated one participant's sex, it printed two classifier probabilities
# next to the scale that decodes them -- so this is a WARNING that asks for a human read, not a
# hard failure.
SINGLING = re.compile(
    r"(?:\bthe (?:one|single|sole) (?:participant|individual|sample|tube|aliquot)\b"
    r"|\bthis participant\b"
    r"|\bone of the two aliquots\b"
    r"|\bpredicts (?:male|female)\b"
    r"|\bgenetic sex\b)",
    re.I,
)


def load_allowlist(path):
    """token -> reason. Blank lines and # comments ignored; format is  TOKEN<whitespace>reason."""
    allowed = {}
    if not os.path.exists(path):
        return allowed
    for raw in io.open(path, encoding="utf-8"):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(None, 1)
        allowed[parts[0]] = parts[1].strip() if len(parts) > 1 else "(no reason recorded)"
    return allowed


def scan(path, allowed):
    text = io.open(path, encoding="utf-8", errors="replace").read()
    lines = text.splitlines()
    violations, allowed_hits, warnings = [], [], []

    for name, rx, desc in FAMILIES:
        for n, line in enumerate(lines, 1):
            for tok in rx.findall(line):
                if tok in allowed:
                    allowed_hits.append((name, tok, n, allowed[tok]))
                else:
                    violations.append((name, tok, n, desc))

    for n, line in enumerate(lines, 1):
        if SINGLING.search(line):
            warnings.append((n, line.strip()[:110]))

    return violations, allowed_hits, warnings


def main(argv):
    allowlist_path = DEFAULT_ALLOWLIST
    files = []
    i = 0
    while i < len(argv):
        if argv[i] == "--allowlist":
            allowlist_path = argv[i + 1]
            i += 2
        else:
            files.append(argv[i])
            i += 1

    if not files:
        print("usage: check_no_identifiers.py [--allowlist PATH] <file> [<file> ...]", file=sys.stderr)
        return 2

    allowed = load_allowlist(allowlist_path)
    total_violations = 0

    for f in files:
        if not os.path.exists(f):
            print(f"ERROR: {f} does not exist", file=sys.stderr)
            return 2
        violations, allowed_hits, warnings = scan(f, allowed)
        base = os.path.basename(f)

        if violations:
            total_violations += len(violations)
            print(f"ERROR: participant or sample identifiers found in {base}:", file=sys.stderr)
            seen = set()
            for name, tok, n, desc in violations:
                if (name, tok) in seen:
                    continue
                seen.add((name, tok))
                print(f"       line {n:>5}  [{name}] {tok}   <- {desc}", file=sys.stderr)

        for name, tok, n, reason in {(a[0], a[1], a[2], a[3]) for a in allowed_hits}:
            print(f"       allowlisted: [{name}] {tok} (line {n}) -- {reason}", file=sys.stderr)

        if warnings:
            print(f"       NOTE: {len(warnings)} passage(s) in {base} describe a single participant.", file=sys.stderr)
            print("             A regex cannot detect an attribute disclosed by arithmetic, so these", file=sys.stderr)
            print("             need a human disclosure read before the document is circulated:", file=sys.stderr)
            for n, snippet in warnings[:8]:
                print(f"               line {n:>5}  {snippet}", file=sys.stderr)
            if len(warnings) > 8:
                print(f"               ... and {len(warnings) - 8} more", file=sys.stderr)

    if total_violations:
        print("", file=sys.stderr)
        print("       A document that leaves the Sandbox must not name a participant or a tube.", file=sys.stderr)
        print("       Either remove the token, or add it to", file=sys.stderr)
        print(f"         {allowlist_path}", file=sys.stderr)
        print("       with the reason it has to stay. Do not silence this guard any other way.", file=sys.stderr)
        return 1

    print(f"identifier guard: clean ({len(files)} file(s) scanned)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
