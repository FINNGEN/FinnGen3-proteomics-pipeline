#!/bin/bash
# sync_version.sh - Synchronize version across project files
#
# This script updates the version number in:
# - VERSION file
# - package.json (version field)
# - Dockerfile (LABEL version and comment)
# - README.md (all version references)
# - Documentation files (docs/*.md and docs/*.tex) - pipeline version references only
#
# Usage: ./scripts/sync_version.sh <version>
# Example: ./scripts/sync_version.sh 1.2.2

set -euo pipefail

VERSION="${1:-}"
if [ -z "$VERSION" ]; then
    echo "Error: Version argument required" >&2
    echo "Usage: $0 <version>" >&2
    exit 1
fi

# Remove 'v' prefix if present (semantic-release provides version without 'v')
VERSION_CLEAN="${VERSION#v}"

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

cd "$PROJECT_ROOT"

echo "Syncing version to ${VERSION_CLEAN}..."

# 1. Update VERSION file
echo "${VERSION_CLEAN}" > VERSION
echo "✓ Updated VERSION file"

# 2. Update package.json
if [ -f "package.json" ]; then
    # Check if jq is available for safer JSON manipulation
    if command -v jq &> /dev/null; then
        # Use jq for robust JSON manipulation
        jq --arg ver "${VERSION_CLEAN}" '.version = $ver' package.json > package.json.tmp && mv package.json.tmp package.json
    else
        # Fallback to sed if jq is not available
        sed -i "s/\"version\": \"[^\"]*\"/\"version\": \"${VERSION_CLEAN}\"/" package.json
    fi
    
    echo "✓ Updated package.json"
fi

# 3. Update Dockerfile
if [ -f "Dockerfile" ]; then
    # Update version comment
    sed -i "s/^# Version: .*/# Version: ${VERSION_CLEAN}/" Dockerfile
    
    # Update LABEL version
    sed -i "s/LABEL version=\".*\"/LABEL version=\"${VERSION_CLEAN}\"/" Dockerfile
    
    echo "✓ Updated Dockerfile"
fi

# 4. Update README.md
if [ -f "README.md" ]; then
    # Update version in **Version**: line (multiple patterns)
    sed -i "s/\*\*Version\*\*: [0-9.]*/**Version**: ${VERSION_CLEAN}/g" README.md
    # Update version in - **Version**: line
    sed -i "s/- \*\*Version\*\*: [0-9.]*/- **Version**: ${VERSION_CLEAN}/g" README.md
    # Update version in example comments (e.g., v1.2.1)
    sed -i "s/v[0-9]\+\.[0-9]\+\.[0-9]\+/v${VERSION_CLEAN}/g" README.md
    
    echo "✓ Updated README.md"
fi

# 5. Update documentation files (docs/) - Pipeline version only
if [ -d "docs" ]; then
    # Update pipeline version in markdown files
    find docs -name "*.md" -type f | while read -r file; do
        # Update "Pipeline version: vX.X.X" patterns
        sed -i "s/Pipeline version: v[0-9.]*/Pipeline version: v${VERSION_CLEAN}/g" "$file"
        # Update "**Pipeline Version**: vX.X.X" patterns
        sed -i "s/\*\*Pipeline Version\*\*: v[0-9.]*/\*\*Pipeline Version\*\*: v${VERSION_CLEAN}/g" "$file"
    done
    
    # Update pipeline version in LaTeX files
    find docs -name "*.tex" -type f | while read -r file; do
        # Update "Pipeline version: vX.X.X" patterns
        sed -i "s/Pipeline version: v[0-9.]*/Pipeline version: v${VERSION_CLEAN}/g" "$file"
        # Update "\textbf{Pipeline Version}: vX.X.X" patterns
        sed -i "s/\\\\textbf{Pipeline Version}: v[0-9.]*/\\\\textbf{Pipeline Version}: v${VERSION_CLEAN}/g" "$file"
    done
    
    echo "✓ Updated documentation files"
fi

echo "Version synchronization complete: ${VERSION_CLEAN}"
