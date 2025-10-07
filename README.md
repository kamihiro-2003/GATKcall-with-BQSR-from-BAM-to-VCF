# GATKcall-with-BQSR-from-BAM-to-VCF
遺伝研スパコン上でBQSR（ベースクオリティを向上させることでsnpcallの精度を向上させる方法）を採用したGATKcallのパイプラインを掲載する
BQSRについて：https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR

## 1) ジョブヘッダと安全設定
```
#!/usr/bin/env bash
#SBATCH -J map_preproc
#SBATCH -p epyc
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH -t 2-00:00:00
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err
```


この部分では、スパコンにジョブを送る際の実行条件を指定している。

-J はジョブ名を設定し、-p は使用するパーティション（計算ノードのグループ）を選ぶ。

--cpus-per-task は1ジョブあたりに割り当てるCPUコア数、--mem は使用メモリ量、-t は最大実行時間を表す。

最後の -o と -e はそれぞれ標準出力と標準エラーのログ保存先を指定しており、%x はジョブ名、%A はジョブID、%a はジョブ配列のインデックスに自動的に置き換えられる。

```
set -euo pipefail
```
この1行は、スクリプトの安全性と再現性を高めるための設定である。

-e は途中のコマンドが1つでも失敗した場合にスクリプト全体を即終了させる。

-u は未定義の変数を使用した際にエラーを出す。

-o pipefail はパイプでつないだコマンドのどれかが失敗した場合でも、正しくエラーを検知できるようにする。


## 2) 要編集パラメータ（入力・リソース・実行系）
### 2-1) 入力
```
genome="${genome:-201123_Psic.polished.scaffold.masked.10k}"
REF_DIR="${REF_DIR:-$HOME/working/Pipeline_SNP_call/reference}"
ref_fasta="${REF_DIR}/${genome}.fa"
MAP_DIR="${MAP_DIR:-$HOME/working/map}"
filename=("f5A.10k")
chr_list="${chr_list:-$HOME/working/Pipeline_SNP_call/reference/chromosomes.list"
```

このセクションでは、解析で使用する入力データや参照ファイルのパスを定義する。

スクリプトを別環境で使う場合や異なるデータセットを扱う場合は、ここを変更するだけで再利用できる。

### 参照系設定

genome は参照ゲノムのベース名を指定する。

REF_DIR は参照ファイル（FASTAやインデックス）を置くディレクトリ。

ref_fasta は実際に使用されるFASTAファイルのフルパスを自動的に組み立てる。

### 入力ファイル設定

MAP_DIR はBAMファイル（マッピング済みリード）の保存先ディレクトリ。

この中のファイルを処理対象として読み込む。

### サンプル指定

filename 配列に処理対象サンプル名を列挙する。


### 染色体リスト

chr_list は染色体（またはスキャフォルド）のリストファイルへのパス。

joint呼び出し（全サンプル統合解析）時に、染色体単位で分割処理を行う際に使用される。

参考：https://github.com/Zoshoku-GH/NIG-SuperComputer/tree/main/Variant_Calling/GATK_Calling
### 2-2) リソース
```
THREADS="16"
JAVA_MEM_LOCAL="32g"
JAVA_MEM_SLURM="24g"
SLURM_HC_PARTITION="short"
SLURM_HC_TIME="0-1:00:00"
SLURM_DB_PARTITION="epyc"
SLURM_DB_TIME="0-24:00:00"
```

このセクションでは、スレッド数やJavaヒープメモリ、そして各ジョブのSLURM実行設定をまとめている。

主に GATK や Spark ベースのツールを安定して動かすためのリソース指定に関わる。

### スレッド数とメモリ設定
THREADS は並列実行に使用するCPUコア数。

JAVA_MEM_LOCAL はローカル処理（単一ノード内でのJava実行）に割り当てるヒープメモリ量を指定し、

JAVA_MEM_SLURM はSLURMジョブ内でのGATK実行時に使うメモリ上限を示す。


### SLURMジョブ設定
SLURM_HC_PARTITION および SLURM_HC_TIME は、HaplotypeCallerなど短時間で終わる per-chromosome ジョブの実行パーティションと最大実行時間を指定する。

一方で、SLURM_DB_PARTITION と SLURM_DB_TIME は、GenomicsDBImportやVCF統合など、より長時間かかる結合処理（concat/filterジョブ）用の設定。

### 2-3) ツール
```
GATK_SIF="${GATK_SIF:-/usr/local/biotools/g/gatk4:4.6.1.0--py310hdfd78af_0}"
SINGULARITY_BIND="${SINGULARITY_BIND:--B /lustre10/home,/home}"
GATK_CMD="${GATK_CMD:-singularity exec ${SINGULARITY_BIND} ${GATK_SIF} gatk}"
ANALYZER_CMD="$GATK_CMD"
BCFTOOLS="${BCFTOOLS:-bcftools}"
```
### Singularityコンテナ設定
GATK_SIF は使用するGATKコンテナのイメージパスを指定する。

SINGULARITY_BIND はコンテナ内からアクセス可能にするディレクトリを -B オプションでバインド指定している（ここでは /lustre10/home および /home）。

GATK_CMD は実際にGATKを起動するためのコマンドテンプレートであり、以降のスクリプト内で共通的に利用する。

### 解析コマンド
ANALYZER_CMD は GATK の AnalyzeCovariates など、GATKベースの補助ツール実行時に用いるコマンドで、GATK_CMD と同一。

### bcftools設定
BCFTOOLS は bcftools 実行コマンドのパスを指定する。

環境にbcftoolsがPATH上で通っている場合は明示的に指定する必要はない。
### 2-4) ディレクトリの作成
```
TMPDIR_ROOT="${TMPDIR_ROOT:-$HOME/working/tmp_gatk}"
SCRIPTS_DIR="${SCRIPTS_DIR:-$(pwd)/scripts}"
mkdir -p "${TMPDIR_ROOT}" "${SCRIPTS_DIR}" \
  logs 09final_gvcf 08qual_check 07bqsr 04rmdup 03fixmate 02addRG 01map_sort 05Haplotypecaller_temp 06GenomicDB
```

一時ディレクトリと出力系の全ディレクトリを作成。途中で落ちても再実行しやすい構成。

## 3) 前提チェック & 参照インデックス自動生成
```
if [ ! -f "${ref_fasta}" ]; then ...; fi
if [ ! -f "${chr_list}" ]; then ...; fi
for K in "${filename[@]}"; do
  if [ ! -f "${MAP_DIR}/${K}.bam" ]; then ...; fi
done
```

入力の存在チェック。欠けていれば即時終了→早期発見。

```
if [ ! -f "${ref_fasta}.fai" ]; then samtools faidx "${ref_fasta}"; fi
DICT_PATH="${REF_DIR}/${genome}.dict"
if [ ! -f "${DICT_PATH}" ]; then ${GATK_CMD} CreateSequenceDictionary -R "${ref_fasta}" -O "${DICT_PATH}"; fi
```

*.fai / .dict 自動生成。GATK実行に必須。

## 4) 実行モード切替と配列ジョブの対象決定
```
RUN_MODE="${RUN_MODE:-per_sample}"  # per_sample | joint
IDX="${SLURM_ARRAY_TASK_ID:-0}"
if (( IDX < 0 || IDX >= ${#filename[@]} )); then IDX=0; fi
K="${filename[$IDX]}"
```

既定は per-sample。sbatch --array=0-... でサンプルごとに実行(最終GVCFまでの作成はper-sampleで実行)。
### ジョブ配列とサンプル選択

IDX は SLURM の配列ジョブ番号（SLURM_ARRAY_TASK_ID）を受け取る変数で、filename 配列のインデックスとして使用される。

例えば sbatch --array=0-4 として実行すれば、0〜4番目のサンプルを並列処理できる。

また、IDX が範囲外だった場合（例：配列外の番号が指定されたとき）でも、安全に先頭サンプル（filename[0]）にフォールバックするようにしており、想定外のジョブエラーを防いでいる。

最終的に K に現在処理対象のサンプル名が代入され、以降の処理で参照される。

## 5) per-sample モード（1〜8）
```
if [[ "${RUN_MODE}" == "per_sample" ]]; then ... exit 0; fi
```
→ per-sample の完了後に早期 exit 。jointブロックには入らない。

ここまででBQSR処理済みのgVCFの作成が完了する

## 5.1) ソート＆インデックス（01）
```
samtools sort -@ "${THREADS}" -o "./01map_sort/${K}.sort.bam" "${MAP_DIR}/${K}.bam"
samtools index -@ "${THREADS}" "./01map_sort/${K}.sort.bam"
```

座標ソートBAMを作成。以降のGATK処理の前提。

## 5.2) Read Group 付与（02）
```
${GATK_CMD} ... AddOrReplaceReadGroups \
  --RGID "FLOWCELLID" --RGLB "N${K}_library1" --RGPU "H0164ALXX140820.2" \
  --RGPL "illumina" --RGSM "${K}" --VALIDATION_STRINGENCY LENIENT
```

BQSRやjointで必須のRGメタ情報を付与。LENIENTで多少のフォーマット逸脱を許容。

## 5.3) FixMateInformation（03）
```
${GATK_CMD} ... FixMateInformation \
  --ADD_MATE_CIGAR true --SO coordinate --ASSUME_SORTED true --CREATE_INDEX true
```

ペアリード整合性修復。CIGAR追加や座標ソート前提で index も生成。

## 5.4) MarkDuplicatesSpark（04）
```
${GATK_CMD} ... MarkDuplicatesSpark \
  -R "${ref_fasta}" -I "./03fixmate/...bam" -O "./04rmdup/...markdups.bam" -M "./04rmdup/...metrics.txt"
```

PCR重複をマーク。Spark版で高速処理。メトリクスも出力。

## 5.5) 暫定 gVCF → GenomicsDB → Genotype → SNP抽出 → しきい値フィルタ → BQSRテーブル作成（05–06）
```
# 暫定 gVCF
${GATK_CMD} HaplotypeCallerSpark --emit-ref-confidence GVCF ... -O "./05Haplotypecaller_temp/${K}.temp.g.vcf"
${GATK_CMD} IndexFeatureFile -I "./05.../${K}.temp.g.vcf"
${GATK_CMD} VcfToIntervalList -I "./05.../${K}.temp.g.vcf" -O "./05.../${K}.interval_list"
# GenomicsDBImport（単サンプルだがBQSRモデル作りのため）
${GATK_CMD} GenomicsDBImport -V "./05.../${K}.temp.g.vcf" -L "./05.../${K}.interval_list" --genomicsdb-workspace-path "${K}_db" ...
# Genotype → SNP選択 → VariantFiltration
${GATK_CMD} GenotypeGVCFs -V "gendb://${K}_db" -O "./06GenomicDB/${K}.temp.vcf"
${GATK_CMD} SelectVariants --select-type-to-include SNP -V "./06.../${K}.temp.vcf" -O "./06.../${K}.temp.snp.vcf"
${GATK_CMD} VariantFiltration -V "./06.../${K}.temp.snp.vcf" -O "./06.../${K}.temp.snp.filtered.vcf" \
  -filter "QD < 2.0" --filter-name "QD2" \
  -filter "QUAL < 30.0" --filter-name "QUAL30" \
  -filter "SOR > 4.0" --filter-name "SOR4" \
  -filter "FS > 60.0" --filter-name "FS60" \
  -filter "MQ < 40.0" --filter-name "MQ40" \
  -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
  -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"
# BQSRテーブル
${GATK_CMD} BaseRecalibrator -I "./04rmdup/...markdups.bam" --known-sites "./06.../${K}.temp.snp.filtered.vcf" \
  -O "./06GenomicDB/${K}_recal_data.table"
```

目的: 参照変異セットが無い前提で、自前で既知サイト相当を作りBQSRへ。

フィルタ式は GATK Best Practices に準拠した一般的なしきい値。

## 5.6) ApplyBQSR ＋ index（07）
```
${GATK_CMD} ApplyBQSR -I "./04rmdup/...markdups.bam" -bqsr "./06.../${K}_recal_data.table" \
  -O "./07bqsr/${K}.sort.addRG.fixmate.markdups.bqsr.bam"
samtools index "./07bqsr/${K}....bqsr.bam"
```

BQSR適用BAMを作成。以降の最終HCの入力へ。

## 5.7) AnalyzeCovariate(08)
```
${ANALYZER_CMD} BaseRecalibrator ... -O "./07bqsr/${K}_recal_data.table.2" || true
${ANALYZER_CMD} AnalyzeCovariates -before "./06.../${K}_recal_data.table" \
  -after "./07bqsr/${K}_recal_data.table.2" -csv "./08qual_check/${K}_recalibration.csv" || true
```
BQSR前と後でのBAMファイルのベースクオリティをCSVで作成し比較できるようにしている。

## 5.8) 最終 HaplotypeCaller（BQSR後BAM → per-sample gVCF）（09）
```
FINAL_HC_MEM="32g"; FINAL_HC_THREADS="8"
${GATK_CMD} --java-options "-Xmx${FINAL_HC_MEM} ..." HaplotypeCaller \
  --emit-ref-confidence GVCF -I "./07bqsr/...bqsr.bam" --native-pair-hmm-threads "${FINAL_HC_THREADS}" \
  -O "./09final_gvcf/${K}.bqsr.g.vcf.gz"
${GATK_CMD} IndexFeatureFile -I "./09final_gvcf/${K}.bqsr.g.vcf.gz"
```

joint入力となる高品質 gVCF を出力。

※ per-sample モードはここで exit 0。jointは下へ。

## 6) joint モード：全サンプル統合(最終GVCF→VCF)
## 6.1) 出力ルートと per-chrom ディレクトリ
```
DIR="$(pwd)/geno_${genome}_$(date +%Y%m%d)_1"
while [ -e "${DIR}" ]; do DIR="${DIR%_*}_$(( ${DIR##*_} + 1 ))"; done
RUN_ROOT="${DIR}"
PERCHROM_DIR="${DIR}/per_chrom"
mkdir -p "${RUN_ROOT}" "${PERCHROM_DIR}" "${DIR}/scripts"
```

同日複数回実行しても衝突しない連番DIRを自動採番。
```
while read -r chrom; do
  CHRDIR="${PERCHROM_DIR}/${chrom}_dir"
  mkdir -p "${CHRDIR}/gatkHC/scripts" "${CHRDIR}/scripts"
done < "${chr_list}"
```

染色体ごとの作業DIRを事前作成。

## 6.2) samplemap.txt の生成
```
samplemap="${DIR}/samplemap.txt"; : > "${samplemap}"
for K in "${filename[@]}"; do
  echo -e "${K}\t$(pwd)/09final_gvcf/${K}.bqsr.g.vcf.gz" >> "${samplemap}"
done
```

gVCFの場所を GATK 用の sample-name-map 形式で列挙。
フルパス（$(pwd)）で解決するため依存ジョブのCWDに影響されない。

## 6.3) 依存ジョブ：染色体ごとに GenomicsDBImport → GenotypeGVCFs
```
export GATK_CMD TMPDIR_ROOT ref_fasta
COMBINE_JIDS_ALL=""
while IFS= read -r chrom; do
  CHRDIR="${PERCHROM_DIR}/${chrom}_dir"
  GENDB="${CHRDIR}/gendb"
  script2="${CHRDIR}/scripts/${chrom}_combine.sh"

  cat > "${script2}" <<EOF
#!/usr/bin/env bash
#SBATCH -t 0-24:00:00
#SBATCH -p ${SLURM_DB_PARTITION}
#SBATCH --mem=16G
#SBATCH -J ${chrom}_combine
...
${GATK_CMD} GenomicsDBImport --genomicsdb-workspace-path "${GENDB}" -L ${chrom} --sample-name-map "${samplemap}" ...
${GATK_CMD} GenotypeGVCFs -R "${ref_fasta}" -V "gendb://${GENDB}" -L ${chrom} --include-non-variant-sites -O "${CHRDIR}/${chrom}.raw.vcf.gz"
tabix -p vcf "${CHRDIR}/${chrom}.raw.vcf.gz" || true
EOF
chmod +x "${script2}"
jid=$(sbatch --parsable "${script2}")
COMBINE_JIDS_ALL="${COMBINE_JIDS_ALL:+${COMBINE_JIDS_ALL}:}${jid}"
done < "${chr_list}"
```

染色体単位で独立ジョブを大量投入し、jid をコロン連結で集約。

--include-non-variant-sites により NO_VARIATION サイトも保持（後段の一貫性に寄与）。

## 6.4) 依存：結合 concat_final.sh の自動生成と提出（afterok）
```
final_script="${SCRIPTS_DIR}/concat_final.sh"
cat > "${final_script}" <<'EOS'
#!/usr/bin/env bash
...
CHR_LIST_FILE="${CHR_LIST_FILE}"
PERCHROM_DIR="${PERCHROM_DIR}"
OUT_DIR="${OUT_DIR}"
BCFTOOLS="${BCFTOOLS}"
REPAIR_MISSING="${REPAIR_MISSING}"
...
# 1) 存在ファイルを LIST に収集。欠損検出。
# 2) 既存 index を掃除 → すべて再 index（tbi/csi 混在対策）
# 3) 欠損があり、REPAIR_MISSING=true の場合：ヘッダを流用して空VCFを自動生成し補完
# 4) 最終 concat → sort → raw.vcf.gz 作成
EOS
chmod +x "${final_script}"
REPAIR_MISSING="${REPAIR_MISSING:-false}"
jid_concat=$(sbatch --parsable \
  -J vcf_concat \
  --dependency=afterok:${COMBINE_JIDS_ALL} \
  -p "${SLURM_DB_PARTITION}" -t "${SLURM_DB_TIME}" --mem=8G \
  -o "${RUN_ROOT}/concat.log" -e "${RUN_ROOT}/concat.log" \
  --export=ALL,CHR_LIST_FILE="${chr_list}",PERCHROM_DIR="${PERCHROM_DIR}",OUT_DIR="${RUN_ROOT}",BCFTOOLS="${BCFTOOLS}",REPAIR_MISSING="${REPAIR_MISSING}" \
  "${final_script}")
```

REPAIR_MISSING: 欠損染色体があっても、ヘッダを借りて空VCFを自動合成→結合継続。

--export=ALL,... で環境変数を子ジョブに明示伝搬。

## 6.5) 依存：最終フィルタ vcf_filter.sh 自動生成と提出（afterok）
```
filter_script="${SCRIPTS_DIR}/vcf_filter.sh"
cat > "${filter_script}" <<'EOS'
#!/usr/bin/env bash
set -euo pipefail
REF="${REF}"
BCFTOOLS="${BCFTOOLS}"
IN_VCF="${IN_VCF}"
OUT_DIR="${OUT_DIR}"
GATK_CMD="${GATK_CMD}"

# 既存 index の掃除 → tabix索引再作成（csi/tbi混在対策）
rm -f "${IN_VCF}.tbi" "${IN_VCF}.csi" 2>/dev/null || true
"${BCFTOOLS}" index -f -t "${IN_VCF}"

# GATK: フィルタマーク付け（SNP/NO_VARIATION混在でも可）
${GATK_CMD} VariantFiltration ... → "${OUT_DIR}/gatk_marked.vcf.gz"
"${BCFTOOLS}" index -f -t "${GATK_MARKED_VCF}"

# GATK: 除外＆NO_VARIATIONは残す（下流の整合性維持）
${GATK_CMD} SelectVariants --exclude-filtered --set-filtered-gt-to-nocall \
  --select-type-to-include SNP --select-type-to-include NO_VARIATION \
  -V "${GATK_MARKED_VCF}" -O "${GATK_FILTERED_VCF}"
"${BCFTOOLS}" index -f -t "${GATK_FILTERED_VCF}"

# 欠損率 >0.2 を除去（bcftools plugin が無い場合は awk フォールバック）
if ${BCFTOOLS} --plugins 2>/dev/null | grep -q fill-tags; then
  ${BCFTOOLS} +fill-tags ... -t F_MISSING | ${BCFTOOLS} view -e 'INFO/F_MISSING>0.2' -O z -o "${OUT_VCF}"
else
  ( ${BCFTOOLS} view -h ...; ${BCFTOOLS} view -H ... \
    | awk -F'\t' '{ n=0; miss=0; for(i=10;i<=NF;i++){ n++; if($i ~ /^\.\/\./) miss++ } if(n==0 || miss/n <= 0.2) print $0; }'
  ) | bgzip > "${OUT_VCF}"
fi
"${BCFTOOLS}" index -f -t "${OUT_VCF}"
EOS
chmod +x "${filter_script}"

jid_filter=$(sbatch --parsable \
  -J vcf_filter \
  --dependency=afterok:${jid_concat} \
  -p "${SLURM_DB_PARTITION}" -t 2-:00:00 --mem=12G \
  -o "${RUN_ROOT}/vcf_filter_%j.log" -e "${RUN_ROOT}/vcf_filter_%j.log" \
  --export=ALL,REF="${ref_fasta}",BCFTOOLS="${BCFTOOLS}",IN_VCF="${RUN_ROOT}/raw.vcf.gz",OUT_DIR="${RUN_ROOT}",GATK_CMD="${GATK_CMD}" \
  "${filter_script}")
```

フィルタ段階は2層：

VariantFiltration で品質フラグ付与 → SelectVariants で除外（SNP/NO_VARIATION維持）

bcftools +fill-tags があれば F_MISSING > 0.2 を除去（無ければ awk 代替）

すべてtabix再索引で最終整合。

## 7) 進捗確認
```
echo "[INFO] submit: vcf_concat ${jid_concat}"
echo "[INFO] submit: vcf_filter ${jid_filter}"
echo "[INFO] 進捗確認: squeue -u $USER | egrep '(_combine|vcf_concat|vcf_filter)' || true"
```

一連のジョブIDを出力し、簡易ウォッチコマンドも表示。

## 8) 出力物まとめ
```
per-sample
09final_gvcf/<SAMPLE>.bqsr.g.vcf.gz（index付）
07bqsr/...bqsr.bam（index付）
08qual_check/<SAMPLE>_recalibration.csv（任意）
joint（geno_${genome}_YYYYMMDD_N/）
per_chrom/<chrom>_dir/<chrom>.raw.vcf.gz（染色体ごと）
raw.vcf.gz（concat+sort後の全染色体結合）
gatk_marked.vcf.gz / gatk_filtered.vcf.gz
mq40_q30_qd2_sor4_fs60_miss0.2.vcf.gz（最終VCF）
```

## 10) 推奨実行例
### per-sample（1サンプル）
```
# filenameを特定の個体名に変更したうえで
sbatch --array=0-0 map_preproc.sh
```
### per-sample（複数サンプル）
```
# スクリプト先頭近くで:
# filename=("S1" "S2" "S3" ...)
sbatch --array=0-2 map_preproc.sh
```

### joint（結合・フィルタのみ）
```
# 事前に 09final_gvcf/*.g.vcf.gz が揃っていること
sbatch --export=RUN_MODE=joint map_preproc.sh
# 欠損染色体があるかも…という時は
REPAIR_MISSING=true sbatch --export=RUN_MODE=joint,REPAIR_MISSING=true map_preproc.sh
```
