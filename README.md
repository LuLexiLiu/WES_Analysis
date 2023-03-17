# WES_Analysis
Current practice for somatic mutation calling of WES data using Mutect2

Overview of the general steps:
1.	Preprocessing: This step involves quality control and filtering of raw sequencing data to remove low-quality reads, adapter sequences, and other artifacts. The tool to perform reads preprocessing is Cutadapt (https://github.com/marcelm/cutadapt). 
2.	Alignment: The preprocessed reads are then aligned to the reference genome using the mapping tool bwa-mem2 (https://github.com/bwa-mem2/bwa-mem2) and the mapped reads are sorted using Samtools (https://github.com/samtools/samtools).
3.	Post-alignment processing: This step involves marking duplicates and base recalibration to improve the accuracy of the alignment. The tools include GATK toolkit (https://github.com/broadinstitute/gatk) and Picard tools (https://broadinstitute.github.io/picard/).
4.	Variant calling: Somatic variants are identified by comparing the aligned reads from tumor and normal samples. We use popular variant calling tool GATK4 Mutect2 (https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2) to perform somatic mutations calling.
5.	Filtering: The somatic variants are filtered based on various criteria such as quality scores, allele frequency, sequence context artifacts, sample contamination and functional impact to reduce false positives.
6.	Annotation: The somatic variants are then annotated with functional information such as gene name, variant effect, and frequency in population databases using Funcotator in GATK4 toolkit (https://gatk.broadinstitute.org/hc/en-us/articles/360037224432-Funcotator).

Example script:

#!/bin/bash

# Set the path to the references and input/output files
export TumorIN=/path/to/tumor/fastq/input
export NormalIN=/path/to/normal/fastq/input
export REF=/path/to/reference_genome.fa
export OUTPUT_DIR=/path/to/output

# Trim adaptors of sequencing reads
cutadapt \
  --minimum-length 20 \
  --overlap 1 \
  -q 20 \
  -a CTGTCTCTTATACACATCT \
  -A CTGTCTCTTATACACATCT \
  -o $TumorIN/Tumor_sample.trimmed.R1.fastq.gz \
  -p $TumorIN/Tumor_sample.trimmed.R2.fastq.gz \
  $TumorIN/Tumor_sample.untrimmed.R1.fastq.gz \
  $TumorIN/Tumor_sample.untrimmed.R2.fastq.gz
cutadapt \
  --minimum-length 20 \
  --overlap 1 \
  -q 20 \
  -a CTGTCTCTTATACACATCT \
  -A CTGTCTCTTATACACATCT \
  -o $NormalIN/Normal_sample.trimmed.R1.fastq.gz \
  -p $NormalIN/Normal_sample.trimmed.R2.fastq.gz \
  $NormalIN/Normal_sample.untrimmed.R1.fastq.gz \
  $NormalIN/Normal_sample.untrimmed.R2.fastq.gz

# Perform alignment and preprocessing
bwa-mem2 mem \
  -t 36 \
  $REF \
  $TumorIN/Tumor_sample.trimmed.R1.fastq.gz \
  $TumorIN/Tumor_sample.trimmed.R2.fastq.gz \
  |samtools sort \
  -o $OUTPUT_DIR/Tumor_sample_bwa_sorted.bam
bwa-mem2 mem \
  -t 36 \
  $REF \
  $NormalIN/Normal_sample.trimmed.R1.fastq.gz \
  $NormalIN/Normal_sample.trimmed.R2.fastq.gz \
  |samtools sort \
  -o $OUTPUT_DIR/Normal_sample_bwa_sorted.bam

picard AddOrReplaceReadGroups \
  -I $OUTPUT_DIR/Tumor_sample_bwa_sorted.bam \
  -O $OUTPUT_DIR/Tumor_sample_bwa_RG_sorted.bam \
  -ID name \
  -LB lib \
  -PL ILLUMINA \
  -SM sample_name \
  -PU barcodes
picard AddOrReplaceReadGroups \
  -I $OUTPUT_DIR/Normal_sample_bwa_sorted.bam \
  -O $OUTPUT_DIR/Normal_sample_bwa_RG_sorted.bam \
  -ID name \
  -LB lib \
  -PL ILLUMINA \
  -SM sample_name \
  -PU barcodes

samtools view \
  -@10 -b -q 15 \
  $OUTPUT_DIR/Tumor_sample_bwa_RG_sorted.bam \
  |samtools sort \
  -@10 -n - \
  -o $OUTPUT_DIR/Tumor_sample_bwa_sorted_Name.bam
samtools view \
  -@10 -b -q 15 \
  $OUTPUT_DIR/Normal_sample_bwa_RG_sorted.bam \
  |samtools sort \
  -@10 -n - \
  -o $OUTPUT_DIR/Normal_sample_bwa_sorted_Name.bam

samtools fixmate \
  -m \
  $OUTPUT_DIR/Tumor_sample_bwa_sorted_Name.bam -\
  |samtools sort \
  -@10 - \
  -o $OUTPUT_DIR/Tumor_sample_bwa_sorted_Name_cor.bam
samtools fixmate \
  -m \
  $OUTPUT_DIR/Normal_sample_bwa_sorted_Name.bam -\
  |samtools sort \
  -@10 - \
  -o $OUTPUT_DIR/Normal_sample_bwa_sorted_Name_cor.bam

# Mark duplicates
samtools markdup \
  -r -s \
  -f $OUTPUT_DIR/Tumor_dup.txt \
  $OUTPUT_DIR/Tumor_bwa_sorted_Name_cor.bam \
  $OUTPUT_DIR/Tumor_bwa_sorted_q15_rmdup.bam
samtools markdup \
  -r -s \
  -f $OUTPUT_DIR/Normal_dup.txt \
  $OUTPUT_DIR/Normal_bwa_sorted_Name_cor.bam \
  $OUTPUT_DIR/Normal_bwa_sorted_q15_rmdup.bam

# Base recalibration
gatk BaseRecalibrator \
  -I $OUTPUT_DIR/Tumor_bwa_sorted_q15_rmdup.bam \
  -R $REF \
  --known-sites /path/to/GATK/dbsnp_146.hg38.vcf.gz \
  -O $OUTPUT_DIR/Tumor_recal_data.table
gatk BaseRecalibrator \
  -I $OUTPUT_DIR/Normal_bwa_sorted_q15_rmdup.bam \
  -R $REF \
  --known-sites /path/to/GATK/dbsnp_146.hg38.vcf.gz \
  -O $OUTPUT_DIR/Normal_recal_data.table

gatk ApplyBQSR \
  -R $REF \
  -I $OUTPUT_DIR/Tumor_bwa_sorted_RG.bam \
  --bqsr-recal-file $OUTPUT_DIR/Tumor_recal_data.table \
  -O $OUTPUT_DIR/Tumor_final.bam
gatk ApplyBQSR \
  -R $REF \
  -I $OUTPUT_DIR/Normal_bwa_sorted_RG.bam \
  --bqsr-recal-file $OUTPUT_DIR/Normal_recal_data.table \
  -O $OUTPUT_DIR/Normal_final.bam

# Remove temporary bam files
rm $OUTPUT_DIR/Tumor_sample_bwa_RG_sorted.bam \
  $OUTPUT_DIR/Normal_sample_bwa_RG_sorted.bam \
  $OUTPUT_DIR/Tumor_sample_bwa_sorted_Name.bam \
  $OUTPUT_DIR/Normal_sample_bwa_sorted_Name.bam \
  $OUTPUT_DIR/Tumor_sample_bwa_sorted_Name_cor.bam \
  $OUTPUT_DIR/Normal_sample_bwa_sorted_Name_cor.bam \
  $OUTPUT_DIR/Tumor_bwa_sorted_q15_rmdup.bam \
  $OUTPUT_DIR/Normal_bwa_sorted_q15_rmdup.bam 

# Call somatic variants using Mutect2
gatk Mutect2 \
  -R $REF \
  -I $OUTPUT_DIR/Tumor_final.bam \
  -I $OUTPUT_DIR/Normal_final.bam \
  -normal Normal \
  -O $OUTPUT_DIR/somatic_variants.vcf

# Filter somatic variants
gatk GetPileupSummaries \
  -I $OUTPUT_DIR/Tumor_final.bam \
  -V /path/to/gnomAD/small_exac_common_3.hg38.vcf.gz \
  -L /path/to/dbsnp/common_all_20180418.vcf.gz \
  -O $OUTPUT_DIR/Tumor_pileups.table
gatk GetPileupSummaries \
  -I $OUTPUT_DIR/Normal_final.bam \
  -V /path/to/gnomAD/small_exac_common_3.hg38.vcf.gz \
  -L /path/to/dbsnp/common_all_20180418.vcf.gz \
  -O $OUTPUT_DIR/Normal_pileups.table

gatk CalculateContamination \
  -I $OUTPUT_DIR/Tumor_pileups.table \
  -matched $OUTPUT_DIR/Normal_pileups.table \
  -O $OUTPUT_DIR/Tumor_contamination.table

gatk CollectF1R2Counts \
  -I $OUTPUT_DIR/Tumor_final.bam \
  -R $REF \
  -O $OUTPUT_DIR/Tumor_f1r2.tar.gz

gatk LearnReadOrientationModel \
  -I $OUTPUT_DIR/Tumor_f1r2.tar.gz \
  -O $OUTPUT_DIR/Tumor_artifact-prior.tar.gz

gatk FilterMutectCalls \
  -R $REF \
  -V $OUTPUT_DIR/Tumor_somatic.vcf.gz \
  --contamination-table \
  $OUTPUT_DIR/Tumor_contamination.table \
  -ob-priors $OUTPUT_DIR/Tumor_artifact-prior.tar.gz \
  -O $OUTPUT_DIR/Tumor_filtered_interval.vcf.gz \

# Annotate the passing-filter somatic variants
zcat $OUTPUT_DIR/Tumor_iltered_interval.vcf.gz\
  |grep '^#' > $OUTPUT_DIR/Tumor.vcf.head

zcat $OUTPUT_DIR/Tumor_filtered_interval.vcf.gz|\
  awk '$7=="PASS"' >\
  $OUTPUT_DIR/Tumor_PASS_interval.vcf.clean

cat $OUTPUT_DIR/Tumor.vcf.head \
  $OUTPUT_DIR/Tumor_PASS_interval.vcf.clean \
  > $OUTPUT_DIR/Tumor_PASS_interval.vcf

rm $OUTPUT_DIR/Tumor.vcf.head \
  $OUTPUT_DIR/Tumor_PASS_interval.vcf.clean

gatk Funcotator \
  -R $REF \
  -V $OUTPUT_DIR/Tumor_PASS_interval.vcf \
  -O $OUTPUT_DIR/Tumor_PASS_interval.table \
  --output-file-format MAF \
  --data-sources-path /path/to/FuncotatorDataSource \
  --ref-version hg38

Downstream analysis after SNP calling:
1.	Further filtering: After annotation, SNPs could be further filtered based on various criteria, such as frequency, quality, and functional impact. For example, SNPs with a minor allele frequency (MAF) higher than a certain threshold in population databases, or those with low read depth or quality scores, may be filtered out.
2.	Comparison with public databases: Comparison with public databases such as COSMIC (https://cancer.sanger.ac.uk/cosmic) and dbSNP (https://www.ncbi.nlm.nih.gov/snp) can help determine the frequency of the somatic SNPs in other cancer types and populations.
3.	Association with clinical outcomes: Somatic SNPs can be associated with clinical outcomes such as survival or response to treatment to determine their prognostic or predictive value (https://translational-medicine.biomedcentral.com/articles/10.1186/s12967-020-02273-4, https://www.nature.com/articles/s41467-022-31780-9).
4.	Functional prediction: Once the somatic SNPs have been validated, their functional impact can be analyzed using tools like PolyPhen (http://genetics.bwh.harvard.edu/pph2/), SIFT (https://sift.bii.a-star.edu.sg/), or PROVEAN (https://www.jcvi.org/research/provean#overview). These tools predict the effect of a variant on the protein structure and function.
5.	Pathway analysis: Somatic SNPs can be further analyzed by looking at the pathways they affect. This can be done using pathway analysis tools like Gene Set Enrichment Analysis (GSEA, https://www.gsea-msigdb.org/gsea/index.jsp), Kyoto Encyclopedia of Genes Genomes (KEGG, https://www.genome.jp/kegg/) or Ingenuity Pathway Analysis (IPA, https://www.qiagenbioinformatics.com/products/ingenuity-pathway-analysis/).
6.	Clustering analysis: If multiple samples have been analyzed, clustering analysis can be performed to identify subgroups of samples with similar somatic SNP profiles. This can be done using unsupervised clustering algorithms like hierarchical clustering (https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html) and principal component analysis (PCA, https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html) or dimensionality reduction techniques like t-SNE (https://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html) and UMAP (https://umap-learn.readthedocs.io/en/latest/).
7.	Clonality analysis: Somatic SNPs can be used to infer tumor clonality and evolution. This can help identify driver mutations and potential targets for therapy. Popular tools to perform clonality analysis include ABSOLUTE (https://software.broadinstitute.org/cancer/cga/absolute/), PyClone (https://github.com/Roth-Lab/pyclone), SciClone (https://github.com/genome/sciclone), ClonEvol (https://github.com/hdng/clonevol) and CITUP (https://github.com/amcpherson/citup).
8.	Signature analysis: Cosmic signature analysis of somatic SNPs involves identifying the specific mutational signatures that are present in a tumor sample. Mutational signatures are patterns of nucleotide changes that occur due to various mutational processes such as DNA damage and repair mechanisms (https://cancer.sanger.ac.uk/signatures/). These signatures can provide insights into the underlying biological mechanisms of cancer development and progression. Several software tools are available for performing cosmic signature analysis including SigProfilerExtractor (https://github.com/AlexandrovLab/SigProfilerExtractor), SigProfilerExtractorR (https://github.com/AlexandrovLab/SigProfilerExtractorR), deconstructSigs (https://github.com/raerose01/deconstructSigs) and SomaticSignatures (https://bioconductor.org/packages/release/bioc/html/SomaticSignatures.html).
9.	Driver mutation identification: Driver mutations are genetic alterations that confer a selective growth advantage to tumor cells, allowing them to proliferate and evade normal cellular processes. Identifying driver mutations can help to understand the underlying biological mechanisms driving cancer development and progression. These mutations are usually found in genes that are involved in key cellular processes such as cell cycle regulation, DNA damage repair, and signal transduction pathways. To identify driver mutation, tools include dNdScv (https://github.com/im3sanger/dndscv), MutSigCV (https://software.broadinstitute.org/cancer/cga/mutsig), oncodriveFML (http://bbglab.irbbarcelona.org/oncodrivefml/home)and OncodriveCLUST (http://bbglab.irbbarcelona.org/oncodriveclustl/home).

