### Issue Description

The GitHub Actions job `61216901782` is failing due to an invalid Docker image tag produced in the workflow file `.github/workflows/docker-build.yml`. The image tag sometimes ends up empty for PR events, which results in an 'invalid reference format' error.

### Steps to Reproduce:
1. Trigger the workflow on a pull request.

### Expected Behavior:
- The image tag should have a valid format.

### Actual Behavior:
- The workflow fails with an invalid Docker reference format error.

### Suggested Solution:
- Update the `docker/metadata-action` tags in the workflow to ensure that tags use PR number or non-empty branch/SHA values for PRs. For example:
  - Use the pattern `pr-${{ github.event.pull_request.number }}` for PR tags.
  - Prefix `sha-` for SHA tags.

### Current Date:
- 2026-01-23 00:18:08

### User:
- RezaJF