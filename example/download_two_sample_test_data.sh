#!/usr/bin/env bash
set -euo pipefail

OUTDIR="${1:-example/fastq}"
mkdir -p "$OUTDIR"
ACCESSIONS=(ERR263534 ERR331060)
API="https://www.ebi.ac.uk/ena/portal/api/filereport"

command -v curl >/dev/null 2>&1 || { echo "ERROR: curl is required" >&2; exit 1; }

for ACC in "${ACCESSIONS[@]}"; do
  echo "Querying ENA for ${ACC}..."
  REPORT=$(curl -fsSL "${API}?accession=${ACC}&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes&format=tsv")
  URLS=$(printf '%s\n' "$REPORT" | awk -F '\t' 'NR==2 {print $2}')
  if [[ -z "$URLS" || "$URLS" == "" ]]; then
    echo "ERROR: no FASTQ URLs found for ${ACC}" >&2
    exit 1
  fi
  IFS=';' read -r -a URL_ARRAY <<< "$URLS"
  for URL in "${URL_ARRAY[@]}"; do
    [[ -z "$URL" ]] && continue
    FILE=$(basename "$URL")
    DEST="${OUTDIR}/${FILE}"
    if [[ -s "$DEST" ]]; then
      echo "Already present: $DEST"
    else
      echo "Downloading $FILE"
      curl -L --retry 5 --retry-delay 10 -o "$DEST" "https://${URL}"
    fi
  done
done

echo "FASTQ files available in ${OUTDIR}:"
ls -lh "$OUTDIR"/*.fastq.gz
