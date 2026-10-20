#!/usr/bin/env bash
# Split single-file paired runs into R1 / R2 (run by download_study.sbatch after the download).
#   bash split_halves.sh manifests/<study>.split.tsv
# Columns: run, src, n_records, out1, out2.
# Some paired runs are served by ENA as one FASTQ that holds all R1 reads followed by all
# R2 reads in the same order (Zhou2024, Friero2023: SRA loaded the mates as separate spots;
# checked 20 Oct 2026: spot i and spot i + n/2 overlap as mates, shifted pairs never do).
# The first n/2 records go to out1, the rest to out2, in one pass; the record count must
# equal n_records (ENA read_count) or nothing is kept. The source file is removed after a
# successful split. Re-running skips runs whose outputs exist.
set -uo pipefail
LIST="${1:?manifest .split.tsv}"
status=0
while IFS=$'\t' read -r run src n out1 out2; do
  [[ "$run" == "run" ]] && continue
  if [[ -s "$out1" && -s "$out2" && ! -e "$src" ]]; then
    echo "$(date +%F_%T) $run SKIP (already split)"; continue
  fi
  if [[ ! -s "$src" ]]; then
    echo "$(date +%F_%T) $run MISSING source $src"; status=1; continue
  fi
  half=$(( n / 2 ))
  if (( n % 2 )); then echo "$(date +%F_%T) $run ODD record count $n"; status=1; continue; fi
  got=$(zcat "$src" | awk -v h=$(( half * 4 )) -v o1="gzip -1 > '$out1.part'" -v o2="gzip -1 > '$out2.part'" \
        '{ if (NR <= h) print | o1; else print | o2 } END { close(o1); close(o2); print NR }')
  if [[ "$got" == "$(( n * 4 ))" ]]; then
    mv "$out1.part" "$out1" && mv "$out2.part" "$out2" && rm -f "$src"
    echo "$(date +%F_%T) $run SPLIT $half + $half records"
  else
    rm -f "$out1.part" "$out2.part"
    echo "$(date +%F_%T) $run COUNT_MISMATCH lines=$got expected=$(( n * 4 ))"; status=1
  fi
done < "$LIST"
exit $status
