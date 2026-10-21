#!/usr/bin/env bash
#SBATCH --account=maize_cpu
#SBATCH --partition=xfer
#SBATCH --qos=xfer
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=2-00:00:00
#SBATCH --output=/share/maize/rrellan/MaizeExpressC/logs/sra_%x_%j.out
# Download one study's .sra files from NCBI's AWS open-data bucket (xfer node).
#   sbatch --job-name=sra_<study> download_sra.sh <study> [parallel, default 4]
# Each file is resumed until its size equals the bucket's Content-Length; integrity is
# then checked by vdb-validate in convert_sra.sbatch (SRA files carry their own checksums).
# Re-submitting is safe: complete files are skipped. Writes <sra_dir>/download.log.
set -uo pipefail
STUDY="${1:?study_id}"
JOBS="${2:-4}"
PIPE=/rsstu/users/r/rrellan/sara/MaizeExpressC/pipeline
SRA_DIR=/share/maize/rrellan/MaizeExpressC/sra/$STUDY
mkdir -p "$SRA_DIR"

fetch() {
  local run="$1" url="$2" dest="$SRA_DIR/$1.sra" want got try
  want=$(curl -sI --max-time 60 "$url" | tr -d '\r' | awk 'tolower($1) == "content-length:" {print $2}')
  [[ -z "$want" ]] && { echo -e "$(date +%F_%T)\t$run\tNO_SIZE" >> "$SRA_DIR/download.log"; return 1; }
  # a finished file whose FASTQ already exists, or a complete .sra, is skipped
  if [[ -s "$dest" && "$(stat -c %s "$dest")" == "$want" ]]; then
    echo -e "$(date +%F_%T)\t$run\tSKIP" >> "$SRA_DIR/download.log"; return 0
  fi
  for try in $(seq 1 30); do
    curl -fsSL --connect-timeout 60 --speed-limit 100000 --speed-time 120 -C - -o "$dest.part" "$url"
    got=$(stat -c %s "$dest.part" 2>/dev/null || echo 0)
    [[ "$got" == "$want" ]] && { mv "$dest.part" "$dest"; echo -e "$(date +%F_%T)\t$run\tOK\t$want" >> "$SRA_DIR/download.log"; return 0; }
    (( got > want )) && rm -f "$dest.part"    # overshoot: start again
    sleep $(( try < 6 ? try * 10 : 60 ))
  done
  echo -e "$(date +%F_%T)\t$run\tINCOMPLETE\t$got/$want" >> "$SRA_DIR/download.log"; return 1
}
export -f fetch; export SRA_DIR

status=0
tail -n +2 "$PIPE/manifests/$STUDY.sra.tsv" | while IFS=$'\t' read -r run url out1 out2 n; do
  # runs already converted (both FASTQ present) need no .sra
  [[ -s "$out1" && ( -z "$out2" || -s "$out2" ) ]] && continue
  printf '%s\t%s\n' "$run" "$url"
done | xargs -P "$JOBS" -L 1 bash -c 'fetch "$0" "$1"' || status=1
grep -c -P "\tOK\t|\tSKIP" "$SRA_DIR/download.log" 2>/dev/null
exit $status
