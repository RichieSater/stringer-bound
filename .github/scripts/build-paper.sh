#!/usr/bin/env bash
set -euo pipefail

# Tectonic obtains its TeX bundle lazily.  Preserve its cache across attempts
# and retry transient relay throttling without weakening the warning check in
# the Makefile target.
for attempt in 1 2 3; do
  echo "Manuscript build attempt ${attempt}/3"
  if make paper; then
    exit 0
  fi
  if [[ ${attempt} -lt 3 ]]; then
    sleep $((30 * attempt))
  fi
done

echo "Manuscript build failed after three attempts" >&2
exit 1
