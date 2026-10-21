#!/usr/bin/env bash
#SBATCH --account=maize_cpu
#SBATCH --partition=compute
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --output=/share/maize/rrellan/MaizeExpressC/logs/convert_%x_%j.out
# Convert one study's .sra files to gzipped FASTQ on a compute node (no internet needed).
#   sbatch --job-name=cv_<study> --dependency=afterok:<sra job> convert_sra.sbatch <study>
# Per run: vdb-validate (SRA checksums), fasterq-dump --split-files, record counts must equal
# reads_per_file in manifests/<study>.sra.tsv in each output file, pigz, then the .sra is
# removed. 4 runs at a time, 4 threads each. Re-running skips runs whose FASTQ exist.
set -uo pipefail
STUDY="${1:?study_id}"
PIPE=/rsstu/users/r/rrellan/sara/MaizeExpressC/pipeline
SCR=/share/maize/rrellan/MaizeExpressC
SIF=$SCR/singularity/sra-tools_3.2.1.sif
export TMPDIR=$SCR/tmp/convert_$STUDY HOME=$SCR/tmp
mkdir -p "$TMPDIR"
source /etc/profile.d/modules.sh
module load apptainer >/dev/null 2>&1

one() {
  local run="$1" out1="$2" out2="$3" n="$4"
  local sra="$SCR/sra/$STUDY/$run.sra" dir; dir=$(dirname "$out1"); mkdir -p "$dir"
  local log="$dir/convert.log" t="$TMPDIR/$run"
  if [[ -s "$out1" && ( -z "$out2" || -s "$out2" ) ]]; then echo -e "$(date +%F_%T)\t$run\tSKIP" >> "$log"; return 0; fi
  [[ -s "$sra" ]] || { echo -e "$(date +%F_%T)\t$run\tMISSING_SRA" >> "$log"; return 1; }
  singularity exec -B /share -B /gpfs_common "$SIF" vdb-validate "$sra" > "$t.validate" 2>&1 \
    || { echo -e "$(date +%F_%T)\t$run\tVALIDATE_FAILED" >> "$log"; return 1; }
  mkdir -p "$t"
  singularity exec -B /share -B /gpfs_common "$SIF" fasterq-dump --split-files -e 4 -t "$t" -O "$t" "$sra" > "$t.dump" 2>&1 \
    || { echo -e "$(date +%F_%T)\t$run\tDUMP_FAILED" >> "$log"; rm -rf "$t"; return 1; }
  local f c ok=1
  for f in "$t/${run}_1.fastq" ${out2:+"$t/${run}_2.fastq"}; do
    c=$(( $(wc -l < "$f") / 4 ))
    [[ "$c" == "$n" ]] || { ok=0; echo -e "$(date +%F_%T)\t$run\tCOUNT_MISMATCH\t$(basename "$f") $c != $n" >> "$log"; }
  done
  if [[ $ok -eq 1 ]]; then
    pigz -p 4 -c "$t/${run}_1.fastq" > "$out1.part" && mv "$out1.part" "$out1"
    [[ -n "$out2" ]] && pigz -p 4 -c "$t/${run}_2.fastq" > "$out2.part" && mv "$out2.part" "$out2"
    rm -rf "$t" "$t.validate" "$t.dump" "$sra"
    echo -e "$(date +%F_%T)\t$run\tOK\t$n reads per file" >> "$log"; return 0
  fi
  rm -rf "$t"; return 1
}
export -f one; export SCR STUDY SIF TMPDIR

tail -n +2 "$PIPE/manifests/$STUDY.sra.tsv" | cut -f1,3,4,5 | xargs -P 4 -L 1 bash -c 'one "$0" "$1" "$2" "$3"'
st=$?
echo "converted: $(grep -c -P "\tOK\t|\tSKIP" "$(dirname "$(awk -F'\t' 'NR==2{print $3}' "$PIPE/manifests/$STUDY.sra.tsv")")/convert.log") of $(( $(wc -l < "$PIPE/manifests/$STUDY.sra.tsv") - 1 ))"
exit $st
