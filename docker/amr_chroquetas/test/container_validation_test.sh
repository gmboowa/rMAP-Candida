#!/usr/bin/env bash
set -euo pipefail

IMAGE="${1:-gmboowa/rmap-myc-candida-amr:2026.07-chroquetas-v9}"

echo "Validating AMR container: ${IMAGE}"

docker run --rm "${IMAGE}" \
  bash -lc 'set -euo pipefail; echo "Checking required executables..."; which run_fungamr_scan; which ChroQueTas.sh; which miniprot; echo "Checking miniprot version..."; miniprot --version; echo "AMR container validation completed successfully."'
