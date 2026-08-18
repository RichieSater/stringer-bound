#!/usr/bin/env bash
set -euo pipefail

# GitHub's Ubuntu package mirrors occasionally leave apt waiting indefinitely.
# Bound every network/package attempt and retry on a fresh apt invocation so a
# transient mirror failure cannot consume the entire algebra-job time budget.
for attempt in 1 2 3; do
  echo "Singular installation attempt ${attempt}/3"
  if sudo timeout --signal=TERM --kill-after=30s 600s \
      env DEBIAN_FRONTEND=noninteractive \
      apt-get \
        -o Acquire::Retries=3 \
        -o Acquire::http::Timeout=30 \
        -o Acquire::https::Timeout=30 \
        -o DPkg::Lock::Timeout=120 \
        update \
    && sudo timeout --signal=TERM --kill-after=30s 600s \
      env DEBIAN_FRONTEND=noninteractive \
      apt-get \
        -o Acquire::Retries=3 \
        -o Acquire::http::Timeout=30 \
        -o Acquire::https::Timeout=30 \
        -o DPkg::Lock::Timeout=120 \
        install --no-install-recommends -y singular; then
    Singular --version
    exit 0
  fi
  if [[ ${attempt} -lt 3 ]]; then
    sleep $((15 * attempt))
  fi
done

echo "Singular installation failed after three bounded attempts" >&2
exit 1
