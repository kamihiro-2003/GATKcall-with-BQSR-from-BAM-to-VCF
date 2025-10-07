#!/usr/bin/env bash
#SBATCH -J map_preproc
#SBATCH -p epyc
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH -t 2-00:00:00
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err
set -euo pipefail

############################################################
# 要編集パラメータ
############################################################
genome="${genome:-201123_Psic.polished.scaffold.masked.10k}"
REF_DIR="${REF_DIR:-$HOME/working/Pipeline_SNP_call/reference}"
ref_fasta="${REF_DIR}/${genome}.fa"
MAP_DIR="${MAP_DIR:-$HOME/working/map}"
filename=("f5A.10k")      # ← 複数サンプルならここに列挙
chr_list="${chr_list:-$HOME/working/Pipeline_SNP_call/reference/chromosomes.list}"

THREADS="16"
JAVA_MEM_LOCAL="32g"       # SLURMメモリ(40G)の80%程度

GATK_SIF="${GATK_SIF:-/usr/local/biotools/g/gatk4:4.6.1.0--py310hdfd78af_0}"
SINGULARITY_BIND="${SINGULARITY_BIND:--B /lustre10/home,/home}"
GATK_CMD="${GATK_CMD:-singularity exec ${SINGULARITY_BIND} ${GATK_SIF} gatk}"
ANALYZER_CMD="$GATK_CMD"
BCFTOOLS="${BCFTOOLS:-bcftools}"

TMPDIR_ROOT="${TMPDIR_ROOT:-$HOME/working/tmp_gatk}"
SCRIPTS_DIR="${SCRIPTS_DIR:-$(pwd)/scripts}"
mkdir -p "${TMPDIR_ROOT}" "${SCRIPTS_DIR}" \
  logs 09final_gvcf 08qual_check 07bqsr 04rmdup 03fixmate 02addRG 01map_sort 05Haplotypecaller_temp 06GenomicDB

############################################################
# 前提チェック
############################################################
if [ ! -f "${ref_fasta}" ]; then echo "[ERR] ref_fasta not found: ${ref_fasta}" >&2; exit 1; fi
if [ ! -f "${chr_list}" ]; then echo "[ERR] chr_list not found: ${chr_list}" >&2; exit 1; fi
for K in "${filename[@]}"; do
  if [ ! -f "${MAP_DIR}/${K}.bam" ]; then echo "[ERR] input BAM not found: ${MAP_DIR}/${K}.bam" >&2; exit 1; fi
done

if [ ! -f "${ref_fasta}.fai" ]; then samtools faidx "${ref_fasta}"; fi
DICT_PATH="${REF_DIR}/${genome}.dict"
if [ ! -f "${DICT_PATH}" ]; then ${GATK_CMD} CreateSequenceDictionary -R "${ref_fasta}" -O "${DICT_PATH}"; fi

############################################################
# 実行モード（配列ジョブ: 1〜8 / 共同遺伝型: 9〜11）
############################################################
RUN_MODE="${RUN_MODE:-per_sample}"  # per_sample | joint

# 配列ジョブのインデックスからサンプルを決定（未指定時は0）
IDX="${SLURM_ARRAY_TASK_ID:-0}"
if (( IDX < 0 || IDX >= ${#filename[@]} )); then IDX=0; fi
K="${filename[$IDX]}"

if [[ "${RUN_MODE}" == "per_sample" ]]; then
  ##########################################################
  # 1) BAM sort / index
  ##########################################################
  samtools sort -@ "${THREADS}" -o "./01map_sort/${K}.sort.bam" "${MAP_DIR}/${K}.bam"
  samtools index -@ "${THREADS}" "./01map_sort/${K}.sort.bam"

  ##########################################################
  # 2) AddOrReplaceReadGroups
  ##########################################################
  ${GATK_CMD} --java-options "-Xmx${JAVA_MEM_LOCAL} -Djava.io.tmpdir=${TMPDIR_ROOT}" AddOrReplaceReadGroups \
    --INPUT "./01map_sort/${K}.sort.bam" \
    --OUTPUT "./02addRG/${K}.sort.addRG.bam" \
    --RGID "FLOWCELLID" \
    --RGLB "N${K}_library1" \
    --RGPU "H0164ALXX140820.2" \
    --RGPL "illumina" \
    --RGSM "${K}" \
    --VALIDATION_STRINGENCY LENIENT

  ##########################################################
  # 3) FixMateInformation
  ##########################################################
  ${GATK_CMD} --java-options "-Xmx${JAVA_MEM_LOCAL} -Djava.io.tmpdir=${TMPDIR_ROOT}" FixMateInformation \
    --I "./02addRG/${K}.sort.addRG.bam" \
    --O "./03fixmate/${K}.sort.addRG.fixmate.bam" \
    --ADD_MATE_CIGAR true \
    --VALIDATION_STRINGENCY LENIENT \
    --SO coordinate \
    --ASSUME_SORTED true \
    --CREATE_INDEX true

  ##########################################################
  # 4) MarkDuplicatesSpark
  ##########################################################
  ${GATK_CMD} --java-options "-Xmx${JAVA_MEM_LOCAL} -Djava.io.tmpdir=${TMPDIR_ROOT}" MarkDuplicatesSpark \
    -R "${ref_fasta}" \
    -I "./03fixmate/${K}.sort.addRG.fixmate.bam" \
    -O "./04rmdup/${K}.sort.addRG.fixmate.markdups.bam" \
    -M "./04rmdup/${K}.marked_dup_metrics.txt"

  ##########################################################
  # 5) 暫定SNP生成（BQSR用）
  ##########################################################
  ${GATK_CMD} --java-options "-Xmx${JAVA_MEM_LOCAL} -Djava.io.tmpdir=${TMPDIR_ROOT}" HaplotypeCallerSpark \
    -R "${ref_fasta}" \
    --emit-ref-confidence GVCF \
    -I "./04rmdup/${K}.sort.addRG.fixmate.markdups.bam" \
    --heterozygosity 0.001 \
    --indel-heterozygosity 0.001 \
    -O "./05Haplotypecaller_temp/${K}.temp.g.vcf"

  ${GATK_CMD} IndexFeatureFile -I "./05Haplotypecaller_temp/${K}.temp.g.vcf"
  ${GATK_CMD} VcfToIntervalList -I "./05Haplotypecaller_temp/${K}.temp.g.vcf" -O "./05Haplotypecaller_temp/${K}.interval_list"

  ${GATK_CMD} GenomicsDBImport \
    -R "${ref_fasta}" \
    -V "./05Haplotypecaller_temp/${K}.temp.g.vcf" \
    -L "./05Haplotypecaller_temp/${K}.interval_list" \
    --genomicsdb-workspace-path "${K}_db" \
    --tmp-dir "${TMPDIR_ROOT}" \
    --batch-size 20 --reader-threads 2

  ${GATK_CMD} GenotypeGVCFs -R "${ref_fasta}" -V "gendb://${K}_db" -O "./06GenomicDB/${K}.temp.vcf"
  ${GATK_CMD} SelectVariants -R "${ref_fasta}" -V "./06GenomicDB/${K}.temp.vcf" --select-type-to-include SNP -O "./06GenomicDB/${K}.temp.snp.vcf"
  ${GATK_CMD} VariantFiltration -R "${ref_fasta}" -V "./06GenomicDB/${K}.temp.snp.vcf" \
    -O "./06GenomicDB/${K}.temp.snp.filtered.vcf" \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 4.0" --filter-name "SOR4" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"

  ${GATK_CMD} BaseRecalibrator -R "${ref_fasta}" -I "./04rmdup/${K}.sort.addRG.fixmate.markdups.bam" \
    --known-sites "./06GenomicDB/${K}.temp.snp.filtered.vcf" -O "./06GenomicDB/${K}_recal_data.table"

  ##########################################################
  # 6) ApplyBQSR + index
  ##########################################################
  ${GATK_CMD} ApplyBQSR -R "${ref_fasta}" \
    -I "./04rmdup/${K}.sort.addRG.fixmate.markdups.bqsr.bam" || true # in case of missing, regen below
  ${GATK_CMD} ApplyBQSR -R "${ref_fasta}" \
    -I "./04rmdup/${K}.sort.addRG.fixmate.markdups.bam" \
    -bqsr "./06GenomicDB/${K}_recal_data.table" \
    -O "./07bqsr/${K}.sort.addRG.fixmate.markdups.bqsr.bam"
  samtools index "./07bqsr/${K}.sort.addRG.fixmate.markdups.bqsr.bam"

  ##########################################################
  # 7) AnalyzeCovariates CSV（任意）
  ##########################################################
  ${ANALYZER_CMD} BaseRecalibrator -R "${ref_fasta}" \
    -I "./07bqsr/${K}.sort.addRG.fixmate.markdups.bqsr.bam" \
    --known-sites "./06GenomicDB/${K}.temp.snp.filtered.vcf" \
    -O "./07bqsr/${K}_recal_data.table.2" || true

  ${ANALYZER_CMD} AnalyzeCovariates \
    -before "./06GenomicDB/${K}_recal_data.table" \
    -after  "./07bqsr/${K}_recal_data.table.2" \
    -csv    "./08qual_check/${K}_recalibration.csv" || true

  ##########################################################
  # 8) Final HC（BQSR後BAM → per-sample GVCF 全ゲノム）
  ##########################################################
  FINAL_HC_MEM="${FINAL_HC_MEM:-32g}"
  FINAL_HC_THREADS="${FINAL_HC_THREADS:-8}"

  echo "[INFO] Final HC start: ${K}"
  ${GATK_CMD} --java-options "-Xmx${FINAL_HC_MEM} -Djava.io.tmpdir=${TMPDIR_ROOT}" \
    HaplotypeCaller \
    -R "${ref_fasta}" \
    --emit-ref-confidence GVCF \
    -I "./07bqsr/${K}.sort.addRG.fixmate.markdups.bqsr.bam" \
    --native-pair-hmm-threads "${FINAL_HC_THREADS}" \
    -O "./09final_gvcf/${K}.bqsr.g.vcf.gz"

  ${GATK_CMD} IndexFeatureFile -I "./09final_gvcf/${K}.bqsr.g.vcf.gz"
  echo "[INFO] Final HC done: ${K}"
  exit 0
fi

############################################################
# 9) 共同遺伝型（全ゲノム一括）
############################################################
DIR="$(pwd)/geno_${genome}_$(date +%Y%m%d)_1"
while [ -e "${DIR}" ]; do DIR="${DIR%_*}_$(( ${DIR##*_} + 1 ))"; done
RUN_ROOT="${DIR}"
PERCHROM_DIR="${DIR}/per_chrom"
mkdir -p "${RUN_ROOT}" "${PERCHROM_DIR}" "${DIR}/scripts"

# 染色体ごとのディレクトリを作成
while read -r chrom; do
  CHRDIR="${PERCHROM_DIR}/${chrom}_dir"
  mkdir -p "${CHRDIR}/gatkHC/scripts" "${CHRDIR}/scripts"
done < "${chr_list}"

# samplemap（per-sample gVCF の場所は 09final_gvcf に統一）
samplemap="${DIR}/samplemap.txt"; : > "${samplemap}"
for K in "${filename[@]}"; do
  echo -e "${K}\t$(pwd)/09final_gvcf/${K}.bqsr.g.vcf.gz" >> "${samplemap}"
done

export GATK_CMD TMPDIR_ROOT ref_fasta

# 染色体ごとに sbatch → 依存IDを収集
COMBINE_JIDS_ALL=""
while IFS= read -r chrom; do
  [[ -z "${chrom}" ]] && continue
  CHRDIR="${PERCHROM_DIR}/${chrom}_dir"
  mkdir -p "${CHRDIR}/scripts"

  GENDB="${CHRDIR}/gendb"
  script2="${CHRDIR}/scripts/${chrom}_combine.sh"

  cat > "${script2}" <<EOF
#!/usr/bin/env bash
#SBATCH -t 0-24:00:00
#SBATCH -p epyc
#SBATCH --mem=16G
#SBATCH -J ${chrom}_combine
#SBATCH -o ${script2}.log
#SBATCH -e ${script2}.log
set -euo pipefail

rm -rf "${GENDB}"

${GATK_CMD} GenomicsDBImport \
  --genomicsdb-workspace-path "${GENDB}" \
  --batch-size 20 \
  -L ${chrom} \
  --sample-name-map "${samplemap}" \
  --overwrite-existing-genomicsdb-workspace true \
  --reader-threads 2 \
  --tmp-dir "${TMPDIR_ROOT}"

${GATK_CMD} GenotypeGVCFs \
  -R "${ref_fasta}" \
  -V "gendb://${GENDB}" \
  -L ${chrom} \
  --include-non-variant-sites \
  -O "${CHRDIR}/${chrom}.raw.vcf.gz"

tabix -p vcf "${CHRDIR}/${chrom}.raw.vcf.gz" || true
echo "[INFO] Genotype(${chrom}) → ${CHRDIR}/${chrom}.raw.vcf.gz 完了"
EOF
  chmod +x "${script2}"

  jid=$(sbatch --parsable "${script2}")
  COMBINE_JIDS_ALL="${COMBINE_JIDS_ALL:+${COMBINE_JIDS_ALL}:}${jid}"
done < "${chr_list}"

echo "[INFO] submitted per-chrom GVCF jobs: ${COMBINE_JIDS_ALL}"

############################################################
# 10) concat VCFs
############################################################
final_script="${SCRIPTS_DIR}/concat_final.sh"
cat > "${final_script}" <<'EOS'
#!/usr/bin/env bash
set -euo pipefail

CHR_LIST_FILE="${CHR_LIST_FILE}"
PERCHROM_DIR="${PERCHROM_DIR}"
OUT_DIR="${OUT_DIR}"
BCFTOOLS="${BCFTOOLS}"
REPAIR_MISSING="${REPAIR_MISSING}"

LIST="${OUT_DIR}/.vcf_filelist.txt"; : > "${LIST}"
missing=0

while IFS= read -r chrom; do
  chrom="${chrom%$'\r'}"
  [[ -z "${chrom}" ]] && continue
  f="${PERCHROM_DIR}/${chrom}_dir/${chrom}.raw.vcf.gz"
  if [ -s "${f}" ]; then
    echo "${f}" >> "${LIST}"
  else
    echo "[ERR] missing: ${f}" >&2
    missing=1
  fi
done < "${CHR_LIST_FILE}"

find "${PERCHROM_DIR}" \( -name '*.raw.vcf.gz.tbi' -o -name '*.raw.vcf.gz.csi' \) -print -delete || true
xargs -r -a "${LIST}" -n1 -P8 "${BCFTOOLS}" index -f -c

if [ "${REPAIR_MISSING}" = "true" ] && [ $missing -ne 0 ]; then
  echo "[INFO] Repairing missing contigs by creating empty VCFs..."
  hdr_src="$(head -n1 "${LIST}" || true)"
  if [[ -z "${hdr_src}" || ! -s "${hdr_src}" ]]; then
    hdr_src="$(find "${PERCHROM_DIR}" -type f -name '*.raw.vcf.gz' | head -n1 || true)"
  fi
  if [[ -z "${hdr_src}" || ! -s "${hdr_src}" ]]; then
    echo "[ERR] cannot find any VCF to borrow header from. Abort."; exit 4
  fi

  : > "${OUT_DIR}/.missing.txt"
  while IFS= read -r chrom; do
    chrom="${chrom%$'\r'}"
    [[ -z "${chrom}" ]] && continue
    f="${PERCHROM_DIR}/${chrom}_dir/${chrom}.raw.vcf.gz"
    [[ -s "${f}" ]] || echo "${chrom}" >> "${OUT_DIR}/.missing.txt"
  done < "${CHR_LIST_FILE}"

  while IFS= read -r chrom; do
    [[ -z "${chrom}" ]] && continue
    outf="${PERCHROM_DIR}/${chrom}_dir/${chrom}.raw.vcf.gz"
    mkdir -p "$(dirname "${outf}")"
    (
      "${BCFTOOLS}" view -h "${hdr_src}"
      echo "##contig=<ID=${chrom}>"
      echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    ) | bgzip > "${outf}"
    "${BCFTOOLS}" index -f -c "${outf}"
  done < "${OUT_DIR}/.missing.txt"

  : > "${LIST}"
  while IFS= read -r chrom; do
    chrom="${chrom%$'\r'}"
    [[ -z "${chrom}" ]] && continue
    echo "${PERCHROM_DIR}/${chrom}_dir/${chrom}.raw.vcf.gz" >> "${LIST}"
  done < "${CHR_LIST_FILE}"

elif [ $missing -ne 0 ]; then
  echo "[ERR] Missing VCF files. Set REPAIR_MISSING=true or regenerate missing per-chrom VCFs."; exit 3
fi

mkdir -p "${OUT_DIR}/tmp_sort"
"${BCFTOOLS}" concat -f "${LIST}" -O u \
| "${BCFTOOLS}" sort -O z -T "${OUT_DIR}/tmp_sort" -o "${OUT_DIR}/raw.vcf.gz"

"${BCFTOOLS}" index -f -c "${OUT_DIR}/raw.vcf.gz"
echo "[INFO] Merge → ${OUT_DIR}/raw.vcf.gz 完了"
EOS
chmod +x "${final_script}"

REPAIR_MISSING="${REPAIR_MISSING:-false}"
jid_concat=$(sbatch --parsable \
  -J vcf_concat \
  --dependency=afterok:${COMBINE_JIDS_ALL} \
  -p epyc -t 0-24:00:00 --mem=8G \
  -o "${RUN_ROOT}/concat.log" -e "${RUN_ROOT}/concat.log" \
  --export=ALL,CHR_LIST_FILE="${chr_list}",PERCHROM_DIR="${PERCHROM_DIR}",OUT_DIR="${RUN_ROOT}",BCFTOOLS="${BCFTOOLS}",REPAIR_MISSING="${REPAIR_MISSING}" \
  "${final_script}")

echo "[INFO] submit: vcf_concat ${jid_concat}"

############################################################
# 11) フィルタ（GATK）
############################################################
filter_script="${SCRIPTS_DIR}/vcf_filter.sh"
cat > "${filter_script}" <<'EOS'
#!/usr/bin/env bash
set -euo pipefail

REF="${REF}"
BCFTOOLS="${BCFTOOLS}"
IN_VCF="${IN_VCF}"
OUT_DIR="${OUT_DIR}"
GATK_CMD="${GATK_CMD}"

GATK_MARKED_VCF="${OUT_DIR}/gatk_marked.vcf.gz"
GATK_FILTERED_VCF="${OUT_DIR}/gatk_filtered.vcf.gz"
OUT_VCF="${OUT_DIR}/mq40_q30_qd2_sor4_fs60_miss0.2.vcf.gz"

# 既存 index を一掃してから “tabix(.tbi)” を作る
rm -f "${IN_VCF}.tbi" "${IN_VCF}.csi" 2>/dev/null || true
"${BCFTOOLS}" index -f -t "${IN_VCF}"

# GATK: 低品質サイトにフィルタフラグを付与
${GATK_CMD} VariantFiltration \
  -R "${REF}" -V "${IN_VCF}" -O "${GATK_MARKED_VCF}" \
  -filter "QD < 2.0"              --filter-name "QD2" \
  -filter "QUAL < 30.0"           --filter-name "QUAL30" \
  -filter "SOR > 4.0"             --filter-name "SOR4" \
  -filter "FS > 60.0"             --filter-name "FS60" \
  -filter "MQ < 40.0"             --filter-name "MQ40" \
  -filter "MQRankSum < -12.5"     --filter-name "MQRankSum-12.5" \
  -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"

# 出力も “tabix(.tbi)” で索引
"${BCFTOOLS}" index -f -t "${GATK_MARKED_VCF}"

# GATK: フィルタ済みサイトの除外＋NO_VARIATION を残す
${GATK_CMD} SelectVariants \
  -R "${REF}" -V "${GATK_MARKED_VCF}" -O "${GATK_FILTERED_VCF}" \
  --exclude-filtered --set-filtered-gt-to-nocall \
  --select-type-to-include SNP --select-type-to-include NO_VARIATION

# 出力も “tabix(.tbi)” で索引
"${BCFTOOLS}" index -f -t "${GATK_FILTERED_VCF}"

# 欠損率 20% 超のサイトを除去
if ${BCFTOOLS} --plugins 2>/dev/null | grep -q fill-tags; then
  ${BCFTOOLS} +fill-tags "${GATK_FILTERED_VCF}" -O u -- -t F_MISSING \
  | ${BCFTOOLS} view -e 'INFO/F_MISSING>0.2' -O z -o "${OUT_VCF}"
else
  ( ${BCFTOOLS} view -h "${GATK_FILTERED_VCF}";
    ${BCFTOOLS} view -H "${GATK_FILTERED_VCF}" \
    | awk -F'\t' '{
        n=0; miss=0;
        for(i=10;i<=NF;i++){ n++; if($i ~ /^\.\/\./) miss++ }
        if(n==0 || miss/n <= 0.2) print $0;
      }'
  ) | bgzip > "${OUT_VCF}"
fi

# 最終出力も “tabix(.tbi)” で索引
"${BCFTOOLS}" index -f -t "${OUT_VCF}"

echo "[INFO] Filtered VCF → ${OUT_VCF} 完了"
EOS
chmod +x "${filter_script}"

jid_filter=$(sbatch --parsable \
  -J vcf_filter \
  --dependency=afterok:${jid_concat} \
  -p epyc -t 2-00:00:00 --mem=12G \
  -o "${RUN_ROOT}/vcf_filter_%j.log" -e "${RUN_ROOT}/vcf_filter_%j.log" \
  --export=ALL,REF="${ref_fasta}",BCFTOOLS="${BCFTOOLS}",IN_VCF="${RUN_ROOT}/raw.vcf.gz",OUT_DIR="${RUN_ROOT}",GATK_CMD="${GATK_CMD}" \
  "${filter_script}")

echo "[INFO] submit: vcf_filter ${jid_filter}"
echo "[INFO] 進捗確認: squeue -u $USER | egrep '(_combine|vcf_concat|vcf_filter)' || true"
