fastqc *.fastq

#!/bin/bash
# Run Trimmomatic on all single-end FASTQ files in the current directory
# Handles files like SRR29454795.fastq, SRR29454795.fastq.gz, *.fq, *.fq.gz

set -euo pipefail

# === Config (edit if needed) ===
THREADS=4
TRIMMOMATIC="trimmomatic-0.39.jar"   # path to Trimmomatic jar
ADAPTERS="TruSeq3-SE.fa"             # adapter file (single-end version)
OUTDIR="01_trim"
MINLEN=50
# ===============================

mkdir -p "$OUTDIR"

# Allow globs that don't match to expand to nothing
shopt -s nullglob

# Patterns to consider (order matters for extension stripping)
patterns=( *.fq.gz *.fastq )

# If no files found, exit gracefully
files_found=0

for f in "${patterns[@]}"; do
  for fq in $f; do
    # skip files already produced by this script (to avoid re-processing)
    if [[ "$fq" == *"_trimmed.fastq.gz" ]] || [[ "$fq" == *"_trimmed.fastq" ]] || [[ "$fq" == *"_trimmed.fq.gz" ]] || [[ "$fq" == *"_trimmed.fq" ]]; then
      echo "Skipping already-trimmed file: $fq"
      continue
    fi

    files_found=1

    # Build sample name by removing known extensions
    SAMPLE=$(basename "$fq")
    SAMPLE=${SAMPLE%.fastq.gz}
    SAMPLE=${SAMPLE%.fq.gz}
    SAMPLE=${SAMPLE%.fastq}
    SAMPLE=${SAMPLE%.fq}

    echo "Processing sample: $SAMPLE (input: $fq)"

    OUT="$OUTDIR/${SAMPLE}_trimmed.fastq.gz"

    # Check trimmomatic jar exists
    if [[ ! -f "$TRIMMOMATIC" ]]; then
      echo "ERROR: Trimmomatic jar not found at '$TRIMMOMATIC'. Update TRIMMOMATIC variable." >&2
      exit 1
    fi

    # Warn if adapter file missing
    if [[ ! -f "$ADAPTERS" ]]; then
      echo "Warning: Adapter file '$ADAPTERS' not found in current directory. Make sure path is correct or Trimmomatic will fail." >&2
    fi

    # Run Trimmomatic (SE mode, explicitly using -phred33)
    java -jar "$TRIMMOMATIC" SE -threads "$THREADS" -phred33 \
      "$fq" "$OUT" \
      ILLUMINACLIP:"$ADAPTERS":2:30:10 \
      LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:"$MINLEN"

    if [[ $? -eq 0 ]]; then
      echo "Done: $SAMPLE -> $OUT"
    else
      echo "Trimmomatic failed for $SAMPLE" >&2
    fi

    echo "----------------------------------------"
  done
done

if [[ $files_found -eq 0 ]]; then
  echo "No FASTQ files found with patterns: ${patterns[*]}. Exiting."
  exit 0
fi

echo "All samples processed. Trimmed reads are in $OUTDIR/"

bwa index GCF_000005845.2_ASM584v2_genomic.fna

#!/bin/bash
set -euo pipefail
shopt -s nullglob

# Number of threads
THREADS=4

# Reference genome
REF="GCF_000005845.2_ASM584v2_genomic.fna"

# Output directory for BAM files
OUTDIR="../02_bam"
mkdir -p "$OUTDIR"

# Sanity check: make sure BWA index files exist
for ext in amb ann bwt pac sa; do
  if [[ ! -f "${REF}.${ext}" ]]; then
    echo "ERROR: Missing index file ${REF}.${ext} in $(pwd). Run: bwa index ${REF}"
    exit 1
  fi
done

# Grab all trimmed FASTQ files
files=(*_trimmed.fastq *.fastq.gz *.fq.gz *.fq)
if ((${#files[@]}==0)); then
  echo "ERROR: No FASTQ files found (expected *_trimmed.fastq or compressed fastq)."
  exit 1
fi

for fq in "${files[@]}"; do
  # Strip extension to get sample name
  SAMPLE=$(basename "$fq")
  SAMPLE=${SAMPLE%.fastq.gz}
  SAMPLE=${SAMPLE%.fq.gz}
  SAMPLE=${SAMPLE%.fastq}
  SAMPLE=${SAMPLE%.fq}

  echo "Aligning sample: $SAMPLE (input: $fq)"

  # Optional read group
  RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA"

  # Alignment + sorting
  bwa mem -t "$THREADS" -R "$RG" "$REF" "$fq" \
    | samtools sort -@ "$THREADS" -o "${OUTDIR}/${SAMPLE}.sorted.bam"

  # Index BAM
  samtools index "${OUTDIR}/${SAMPLE}.sorted.bam"

  # Alignment stats
  samtools flagstat "${OUTDIR}/${SAMPLE}.sorted.bam" > "${OUTDIR}/${SAMPLE}.flagstat.txt"

  echo "Done: ${OUTDIR}/${SAMPLE}.sorted.bam"
  echo "----------------------------------------"
done

echo "All alignments complete. Results are in $OUTDIR/"

bcftools mpileup -Ou -f GCF_000005845.2_ASM584v2_genomic.fna *.sorted.bam \
  | bcftools call -mv -Oz -o variants.raw.vcf.gz
bcftools view -v snps variants.raw.vcf.gz -Oz -o snps.vcf.gz
bcftools view -v indels variants.raw.vcf.gz -Oz -o indel.vcf.gz
# ==================================================================================

#!/usr/bin/env bash
# WSL-safe script to extract variant burden per KEGG category
# Requires: bcftools, bedtools, awk

set -euo pipefail

# ---------------- CONFIG ----------------
VCF="variants.raw.vcf"   # your actual VCF file
GFF="GCF_000005845.2_ASM584v2_genomic.gff"   # matching E. coli GFF
KEGG_CAT="function_vol2_category.csv"        # your KEGG category mapping
OUTDIR="category_burden"
mkdir -p "$OUTDIR"
# ----------------------------------------

# Check required tools
for tool in bcftools bedtools awk; do
    command -v $tool >/dev/null || { echo "ERROR: $tool not found in PATH"; exit 1; }
done

# 1) Get sample names
bcftools query -l "$VCF" > "$OUTDIR/samples.txt"
echo "Samples detected:"
cat "$OUTDIR/samples.txt"

# 2) Extract genotypes in matrix format (GT only)
bcftools query -f '%CHROM:%POS:%REF:%ALT[\t%GT]\n' "$VCF" > "$OUTDIR/genotypes.raw"

# 3) Convert GT (0/0,0/1,1/1,./.) → allele counts (0,1,2,NA)
awk -v OFS="\t" '
{
    id=$1; printf "%s", id;
    for(i=2; i<=NF; i++){
        gt=$i;
        if(gt=="." || gt=="./." || gt=="./") {printf OFS "NA"; continue;}
        gsub("|","/",gt); split(gt,a,"/");
        if(a[1]=="." || a[2]=="."){printf OFS "NA"; continue;}
        cnt=a[1]+a[2];
        printf OFS cnt;
    }
    print "";
}' "$OUTDIR/genotypes.raw" > "$OUTDIR/genotype_counts_wide.tmp"

# 4) Add header
{
    echo -n "variant"
    while read -r s; do echo -n -e "\t$s"; done < "$OUTDIR/samples.txt"
    echo
    cat "$OUTDIR/genotype_counts_wide.tmp"
} > "$OUTDIR/genotype_counts_wide.tsv"
rm "$OUTDIR/genotype_counts_wide.tmp" "$OUTDIR/genotypes.raw"

# 5) Create BED for variants (variant id format CHR:POS:REF:ALT)
cut -f1 "$OUTDIR/genotype_counts_wide.tsv" | tail -n +2 | \
awk -F: -v OFS="\t" '{chr=$1;pos=$2;ref=$3;alt=$4;start=pos-1;end=pos-1+length(ref);print chr,start,end,$0}' \
> "$OUTDIR/variants.bed"

# 6) Extract genes from GFF
awk '$3=="gene" || $3=="CDS" {
  attr=$9; gid="."
  if (match(attr,/locus_tag=([^;]+)/,m)) gid=m[1]
  else if (match(attr,/gene=([^;]+)/,m)) gid=m[1]
  else if (match(attr,/Name=([^;]+)/,m)) gid=m[1]
  else if (match(attr,/ID=([^;]+)/,m)) gid=m[1]
  print $1"\t"($4-1)"\t"$5"\t"gid
}' "$GFF" > "$OUTDIR/genes.bed"

# 7) Intersect variants with genes
bedtools intersect -a "$OUTDIR/variants.bed" -b "$OUTDIR/genes.bed" -wa -wb > "$OUTDIR/variant_gene_hits.tsv"

# 8) Build variant->gene map
awk -F"\t" '{print $4"\t"$8}' "$OUTDIR/variant_gene_hits.tsv" | sort -u > "$OUTDIR/variant_gene_map.csv"

echo "✅ Done. Outputs in $OUTDIR:"
echo " - genotype_counts_wide.tsv : matrix of variants × samples (0/1/2/NA)"
echo " - variants.bed              : BED file of variants"
echo " - genes.bed                 : BED file of genes"
echo " - variant_gene_map.tsv      : mapping of variant → gene"

#!/bin/bash
# Usage: bash variant_region_classification.sh category_burden/variants.bed GCF_000005845.2_ASM584v2_genomic.gff

# ======================
# Input check
# ======================
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <variants.bed> <annotation.gff3>"
    exit 1
fi

VARIANTS=$1
GFF=$2

# ======================
# Step 1. Extract feature regions
# ======================
echo "[1/5] Extracting feature regions from GFF3..."

awk '$3=="exon" {print $1"\t"$4-1"\t"$5"\t"$9}' $GFF > exons.bed
awk '$3=="five_prime_UTR" {print $1"\t"$4-1"\t"$5"\t"$9}' $GFF > utr5.bed
awk '$3=="three_prime_UTR" {print $1"\t"$4-1"\t"$5"\t"$9}' $GFF > utr3.bed
awk '$3=="gene" {print $1"\t"$4-1"\t"$5"\t"$9"\t.\t"$7}' $GFF > genes.bed

# ======================
# Step 2. Generate introns and flanking regions
# ======================
echo "[2/5] Deriving introns, upstream and downstream regions..."

bedtools subtract -a genes.bed -b exons.bed > introns.bed

awk '$3=="gene" {
    if ($7=="+") print $1"\t"($4-1001)"\t"($4-1)"\t"$9"\t.\t"$7;
    else print $1"\t"($5+1)"\t"($5+1001)"\t"$9"\t.\t"$7
}' $GFF > upstream1kb.bed

awk '$3=="gene" {
    if ($7=="+") print $1"\t"($5+1)"\t"($5+1001)"\t"$9"\t.\t"$7;
    else print $1"\t"($4-1001)"\t"($4-1)"\t"$9"\t.\t"$7
}' $GFF > downstream1kb.bed

# ======================
# Step 3. Classify variants
# ======================
echo "[3/5] Classifying variants by region..."

declare -A regions=( \
["Exonic"]="exons.bed" \
["Intronic"]="introns.bed" \
["UTR5"]="utr5.bed" \
["UTR3"]="utr3.bed" \
["Upstream_1kb"]="upstream1kb.bed" \
["Downstream_1kb"]="downstream1kb.bed" )

echo -e "Category\tVariant_Count" > region_counts.tsv
total=$(wc -l < $VARIANTS)

for cat in "${!regions[@]}"; do
    count=$(bedtools intersect -a $VARIANTS -b ${regions[$cat]} -wa -u | wc -l)
    echo -e "$cat\t$count" >> region_counts.tsv
done

# Intergenic = total - all annotated
annotated=$(awk 'NR>1{sum+=$2}END{print sum}' region_counts.tsv)
intergenic=$((total - annotated))
echo -e "Intergenic\t$intergenic" >> region_counts.tsv

# ======================
# Step 4. Output summary
# ======================
echo "[4/5] Writing summary table..."
cat region_counts.tsv | column -t

# ======================
# Step 5. Done!
# ======================
echo "[5/5] Completed successfully!"
echo "Results saved in region_counts.tsv"

./variant_region_classification.sh category_burden/variants.bed GCF_000005845.2_ASM584v2_genomic.gff
# ==========================================================================================

# Calculate transition/transversion ratio
bcftools stats variants.raw.vcf > control_stats.txt
grep "^TSTV" control_stats.txt

# Extract SNP types for spectrum plot
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' variants.raw.vcf | \
awk 'length($3)==1 && length($4)==1 {print $3"->"$4}' | \
sort | uniq -c > control_snp_spectrum.txt

# View results
cat control_snp_spectrum.txt

# Using the variant_gene_map.csv from vcf_to_category_burden.sh
# Count variants per gene
cut -f2 category_burden/variant_gene_map.csv | sort | uniq -c | \
sort -rn > control_genes_variant_count.txt

head control_genes_variant_count.txt

# Create a summary file
echo "=== Healthy Control Samples Summary ===" > control_summary.txt
echo "Date: $(date)" >> control_summary.txt
echo "" >> control_summary.txt

# Total variants
echo "Total variants: $(bcftools view -H variants.raw.vcf | wc -l)" >> control_summary.txt

# SNPs vs Indels
echo "SNPs: $(bcftools view -v snps variants.raw.vcf | grep -v "^#" | wc -l)" >> control_summary.txt
echo "Indels: $(bcftools view -v indels variants.raw.vcf | grep -v "^#" | wc -l)" >> control_summary.txt

# Region distribution
cat control_region_counts.tsv >> control_summary.txt

# Genes affected
echo "Genes affected: $(cut -f2 category_burden/variant_gene_map.csv | sort -u | wc -l)" >> control_summary.txt

# View summary
cat control_summary.txt

# Create SNP types table
bcftools query -f '%REF\t%ALT\n' variants.raw.vcf | \
awk 'length($1)==1 && length($2)==1 {print $1"->"$2}' | \
sort | uniq -c | awk '{print $2"\t"$1}' | \
sort -t'>' -k1,1 -k2,2 > control_snp_types.txt

# View to verify format (should match your CRC format)
cat control_snp_types.txt
Expected output format:

text
A->C    3800
A->G    18000
A->T    3700
C->A    4000
C->G    2800
C->T    19500
G->A    19800
G->C    2700
G->T    3800
T->A    3700
T->C    17800
T->G    3800
2. Generate Indel Length Distribution (like indel_length_distribution.png)
bash
# Create indel length table
bcftools view -v indels variants.raw.vcf | \
awk '{
    ref=length($4); 
    alt=length($5); 
    diff=alt-ref;
    if(diff != 0) print diff
}' | sort -n | uniq -c | awk '{print $2"\t"$1}' > control_indel_lengths.txt

# View to verify format
cat control_indel_lengths.txt
Expected output format:

text
-8      2
-7      5
-6      6
-5      7
-4      8
-3      9
-2      10
-1      11
1       12
2       13
3       14
4       15
5       16
6       17
7       18
3. (Optional) Create the PNG files
If you want to generate the actual PNG files like before, you can use gnuplot:

bash
# For SNP plot
gnuplot -persist <<-EOF
    set terminal png size 800,600
    set output 'control_snp_types.png'
    set title "Control Cohort: SNP Types"
    set style data histogram
    set style fill solid
    set xtics rotate by -45
    plot 'control_snp_types.txt' using 2:xtic(1) notitle
EOF

# For Indel plot
gnuplot -persist <<-EOF
    set terminal png size 800,600
    set output 'control_indel_lengths.png'
    set title "Control Cohort: Indel Length Distribution"
    set style data histogram
    set style fill solid
    set xlabel "Indel Length"
    set ylabel "Count"
    plot 'control_indel_lengths.txt' using 2:xtic(1) notitle
EOF
# ==============================================================================

# ---- Input files ----
vcf_file  <- "E:/variants.raw.vcf"                  # your VCF file
gff_file  <- "E:/GCF_000005845.2_ASM584v2_genomic.gff"   # genome annotation (E. coli)
kegg_file <- "E:/enrich.csv"             # KEGG category/enrichment results

# ---- Step 1. Read VCF and GFF ----
cat("Reading VCF...\n")
vcf <- readVcf(vcf_file)

cat("Reading GFF...\n")
gff <- import.gff(gff_file)

# keep coding sequences
cds <- gff[gff$type == "CDS"]

# variants as genomic ranges
gr <- rowRanges(vcf)

# ---- Step 2. Map variants to genes ----
cat("Mapping variants to genes...\n")
hits <- findOverlaps(gr, cds)
head(hits)
variant_gene_map <- data.frame(
  variant   = names(gr)[queryHits(hits)],
  geneID    = cds$locus_tag[subjectHits(hits)],
  gene_name = cds$Name[subjectHits(hits)]
)

# ---- Step 3. Load KEGG categories ----
cat("Reading KEGG categories...\n")
kegg <- read.csv(kegg_file, sep=",", header=TRUE)

# Ensure geneID is text
kegg$geneID <- as.character(kegg$geneID)

# Expand KEGG gene lists (split "b0001/b0002/..." into multiple rows)
kegg_long <- kegg %>%
  separate_rows(geneID, sep="/")

# ---- Step 4. Join variants with KEGG ----
cat("Joining variant-gene map with KEGG...\n")
final_map <- variant_gene_map %>%
  inner_join(kegg_long, by="geneID")

# ---- Step 5. Save results ----
out_file <- "E:/variant_kegg_map.tsv"
write.table(final_map, out_file, sep="\t", quote=FALSE, row.names=FALSE)

cat("✅ Done. Results saved to", out_file, "\n")

# ---- Step 6. Optional summary ----
summary_file <- "variant_kegg_summary.tsv"
summary_table <- final_map %>%
  group_by(category, subcategory, Description) %>%
  summarise(variant_count = n_distinct(variant),
            gene_count    = n_distinct(geneID),
            .groups="drop")

write.table(summary_table, summary_file, sep="\t", quote=FALSE, row.names=FALSE)
cat("📊 Summary saved to", summary_file, "\n")

map <- read_tsv("E:/variant_gene_hits.tsv") 
colnames(map) <- c("variants", "geneID")
library(clusterProfiler)
library(org.EcK12.eg.db)
# Perform KEGG enrichment analysis
kegg_enrich_vol2 <- enrichKEGG(gene = map$geneID, 
                          organism = "eco")
head(kegg_enrich_vol2)

write.csv(kegg_category, "E:/mutation_related.csv")
dotplot(kegg_enrich_vol2, showCategory = 20)
dotplot(kegg_category, showCategory = 20)
emapplot(kegg_enrich_vol2)
cnetplot(kegg_enrich_vol2, showCategory = 10, circular = TRUE, colorEdge = TRUE)

library(ggplot2)
library(dplyr)

summary_table <- map %>%
  group_by(Description) %>%
  summarise(variant_count = n_distinct(variant),
            gene_count = n_distinct(geneID),
            .groups = "drop") %>%
  arrange(desc(variant_count))

# Plot top 15 affected pathways
ggplot(head(summary_table, 15),
       aes(x = reorder(Description, variant_count), 
           y = variant_count)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Top KEGG Pathways Affected by Mutations",
       x = "KEGG Pathway",
       y = "Number of Variants")

########################################################################

cat_summary <- final_map %>%
  group_by(category) %>%
  summarise(variant_count = n_distinct(variant),
            gene_count = n_distinct(geneID),
            .groups = "drop")

ggplot(cat_summary, aes(x = category, y = variant_count, fill = category)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 15, hjust = 1)) +
  labs(title = "Variants per KEGG Category", x = "KEGG Category", y = "Variant Count")

###############################################################################


library(reshape2)

heat_data <- final_map %>%
  group_by(Description, variant) %>%
  summarise(n = n(), .groups="drop") %>%
  acast(Description ~ variant, value.var="n", fun.aggregate=length)

heatmap(heat_data, scale="row", Colv=NA)

####################################################################


pathway_matrix <- final_map %>%
  group_by(Description, variant) %>%
  summarise(n=1, .groups="drop") %>%
  tidyr::pivot_wider(names_from=variant, values_from=n, values_fill=0)

rownames(pathway_matrix) <- pathway_matrix$Description
pathway_matrix <- pathway_matrix[,-1]

pca_res <- prcomp(pathway_matrix, scale.=TRUE)

plot(pca_res$x[,1:2], col="blue", pch=19,
     xlab="PC1", ylab="PC2",
     main="PCA of KEGG Pathways Based on Mutations")


##########################################################################

# Load required libraries
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tibble)
library(tidyr)

# Step 1: Load your final map file
# (replace with your actual file name if different)
final_map <- read.delim("E:/variant_kegg_map.tsv", header = TRUE, sep = "\t")

# Step 2: Build a matrix of pathways × variants (or pathways × genes)
mat <- table(final_map$Description, final_map$variant)  # pathways × variants

# Step 3: Run PCA
pca_res <- prcomp(mat, scale. = TRUE)

# Step 4: Prepare dataframe with pathway labels
pca_df <- as.data.frame(pca_res$x) %>%
  tibble::rownames_to_column("Pathway")

# Step 5: Plot PCA with pathway names
ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(color = "steelblue", size = 3) +
  geom_text_repel(aes(label = Pathway), size = 3, max.overlaps = 20) +
  theme_minimal() +
  labs(title = "PCA of KEGG Pathways Based on Mutations",
       x = "Principal Component 1",
       y = "Principal Component 2")
