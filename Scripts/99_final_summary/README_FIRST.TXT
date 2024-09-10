Summary of Key Output Files:

1. 02_Read_Strandness/Job_Summary/Read_Strandness_Stats.txt which
contains information how reads were stranded for strand-specific
RNA-seq data.

2. 03_FASTQC/Job_Summary contains information about quality control
for each separate sample

4. 09d_DE_1_RefSeqLncRNA76k contains information about Up/Down
regulated genes -
Output_DiffExp_1h_featureCounts_RefSeqLncRNA76k_FullGeneBody -
Output_DiffExp_1i_featureCounts_RefSeqLncRNA76k_ExonCollapsed

3. 14_final_summary/output contains summary for all previous steps:

- multiqc_report - html which aggregates all information about
  mapping, quality control, counting

- Segex_09d - information about differential expression for 2
        different counting methods 1. FullGeneBody - counts reads for
        introns and exons for each gene 2. ExonCollapsed - counts
        reads, which were mapped only to exons

- 13d_RefSeqLncRNA76k_ExonCollapsed - PCA, Pearson and Spearman
  correlation plots for groups of samples. They represent information
  on how groups correlate with each other.
