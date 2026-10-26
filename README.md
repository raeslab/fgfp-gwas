# rep-cookbook
Replication Cookbook for microbiome GWAS

### Microbiome processing: Sample inference from amplicon data
Perform sample inference from high-throughput amplicon data using [dada2](https://github.com/benjjneb/dada2).

The dada2 pipeline performs the following steps:
1. Quality control: Perform filtering and trimming of reads.
2. Perform dereplication: Group amplicon reads with same sequence into unique sequences.
3. Sample inference: Using the core sequence-variant inference algorithm to the dereplicated data.
4. Remove chimeras: Identify sequences that are exact bimeras.
5. Taxonomical classification.

Perform these analysis in R using the following script:
```
Seq_dada2.R
```
Transform the dada2 output into count tables and perform rarefaction using the following script:
```
dada2_to_taxtables.R
```

### Genotype QC and imputation

1. FGFP samples were processed on two different arrays, Human Core Exome v1.0 (576 samples) and Human Core Exome v1.1 (2112 samples), resulting in 545, 535 markers in merged raw data set. 
2. Allele calling was performed by GenomeStudio.
3. SNP QC: cross check with sex information, remove Unmapped (chr0) Variants, update to RS ids and remove duplicates; remove variants with >5% missing, remove monomorphic alleles, remove HWE p<1e-5 and remove AT/CG ambiguous sites. 
4. Sample QC: remove sample with >5% SNPs missing, remove samples with heterozygosity <3SD of mean value, remove samples with relatedness >0.025. 
5. Merge with 1KG samples and remove those with outstanding population structure. 
6. Imputation with HRCv1.1 server and EUR population as reference. 
7. Post imputation QC: remove SNPs that are monomorphic, INFO < 0.5 & MAF < 1% .

The scratch code in the following file will help to prepare the genotype data and prepare the files for the next step:
```
set_and_run_GWAS.txt
```

### Association tests
#### 1 - Select overlapping DADA2 between FGFP and Kiel cohorts.
There are a total of 92 taxa and 3 alpha diversity measures. DADA2 list can be found here:
```
analyzable_taxa.txt
```

#### 2 - Transformations and association analysis.
We performed:
   - Linear regression in SNPTEST.
     - Log2 transformed data
     - Rank normal transformed data
    - Generalized linear regression on binary traits in SNPTEST.
      - This is performed on all potential zero-inflated taxa in a hurdle step analysis
        a. Hurdle taxa = those taxa with 5% of more individuals with zero counts
        b. Zero samples set to “0”, non-zero samples set to “1” to for a test of presence vs absence
    - Lude Franke meta analysis pipeline.

The following Rscript is provided to assist in transformation of the taxa abundance and alpha-diversity data.
```
DADA2_GWAS_Replication_Prep.Rmd 
```
It runs the following steps:
   - Identify the taxa with >=5% of samples with zero-values. 
     - While the above step would be the typical procedure we would like for you define the same taxa that we did as the hurdle taxa. These taxa are listed in the file “TaxaNamesHurdleOnly.txt”.
     - Define these as hurdle taxa.
     - create a binary variable for these taxa.
       - 0 = 0
       - !0 =1
     - create zero-truncated abundance data for these taxa.
       - all zero values are turned to NA and all non-zero values are treated as normal.
   - All remaining abundance data (both non-hurdle and hurdle-truncated abundance data) are rank normal transformed and log2 transformed
     - Rank normal transformation: completed with the rntransfrom() function from the GenABEL package
     - Log2 transformation: all 0 values are turned into NA, effectively performing a truncation of zero values just as done in the hurdle step.
   - Residualize RNT and LOG2 data
     - Fit a linear model in R, or your favorite programing language, setting the abundance values as the response and the following variables as explanatory variables. 
       - Top 10 PCs
       - Sex, age
       - Study specific batch variables (exttype, extractyear, aliquotyear, aliquotby, plate)
     - Extract the residuals of the fitted model
     - Take care to return NAs for those samples that were already NAs. This is done to maintain order and structure in the data file for the association analysis.

#### 3 - Perform the association in SNPTEST.
  - Construct the .sample file for SNPTEST that includes all of the 
    - Sample ids
    - Covariates
    - LOG2 residual data
    - RNT residual data
    - Hurdle Binary data
  - Extract Variants to test
    - The ~25K variant positions to test are in the file “Sites_2_Test.txt”
    - This file contains HRCv1.1 SNPs and 1000G indels, as such you may not identify all of the variants in your data. This is OK! Proceed with all of those that you can identify.
      - We find 23,098 significant variants in the HRCv1.1 data.
    - Extract and prepare genotype data with qctools. Below is some code that could be used to do that:

```
for i in ~/YOURDIRECTORYPATH/*.bgen; do echo ${i}; qctool -g ${i} -incl-rsids SNPids2test.txt -ofiletype gen -og sigsites_chr${i# YOURDIRECTORYPATH /data_chr}.gen; done 

wc -l sigsites_chr*.gen  # 23,098 sites identified

## Concatenate all the sigsites together
for i in sigsites_chr*.gen; do echo ${i}; cat ${i} >> all_sigsites2text.gen; done
wc -l all_sigsites2text.gen # 23,098 sites 

## remove chromosome data
rm sigsites_chr*.gen
```
  - Perform association analysis for each taxa transformation

```
## RUN SNPTEST for Hurdle Binary Analysis
VAR1=C_Alphaproteobacteria_HB 
snptest_v2.5 -data all_sigsites2text.gen TransResCovar_MYSTUDY.sample -o $VAR1.out -log $VAR1.log -frequentist 1 -method score -cov_names PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 batchvar1 batchvar2 sex age -pheno $VAR1 -use_raw_phenotypes

##RUN SNPTEST for linear RNT and LOG analysis
VAR1=G_Blautia_RNT_Residuals
snptest_v2.5 -data all_sigsites2text.gen TransResCovar_MYSTUDY.sample -o $VAR1.out -log $VAR1.log -frequentist 1 -method score -pheno $VAR1 -use_raw_phenotypes
```

### Manhattan Plots for the association tests.

The following scripts will help to create the Manhattan plots displaying the significance of each SNP associated with the microbial taxa tested.

Pull out the required information from tables with all the results. For creating the Manhattan plots, only the SNP names ("SNP"), chromosomes ("CHR"), base pair coordinates ("BP") and p-values ("P") are needed (one p-value for each statistical model used to perform the MGWAS analysis).

```
# In a directory containing all the results for each taxa:
cat data_chr*out | grep -v "#" | grep -v "rs" | awk '{if ($19 >= 0.01 && $9 >= 0.3 && $22 >= 0.3) print
$2,$21}' | sed 's/:/\t/g' | sed 's/_.* /\t/g' > rnt.pvalues
```
Then use the following scripts to help you obtaining the manhattan plots:
```
mhplots.R
```

### Phylogenetic trees 
Produce a phylogenetic tree of the taxa under study with overlaying piecharts containing information on loci counts (using different methods) and branches coloured according to heritability traits.

#### Construct a phylogenetic tree using all the taxa in the study
Get all the taxa to include (genera, families, classes…) and load them in phyloT (http://phylot.biobyte.de/). This program generates phylogenetic trees based on the NCBI taxonomy. Choose options “expanded” and “polytomy yes”, otherwise random bifurcated structure is generated. 

**note**: Some names have to be corrected to NCBI taxonomy.

For example:
```
## From this list of taxa:
Bacteroidetes, Barnesiella, Bifidobacterium, Clostridia, Dorea, Prevotella

## phyloT generates this tree:
(((Prevotella,Barnesiella)Bacteroidales,(Dorea,Bifidobacterium)Terrabacteria_group)Bacteria);


                        /-Prevotella
           /Bacteroidales
          |             \-Barnesiella
-- /Bacteria
          |                   /-Dorea
           \Terrabacteria_group
                              \-Bifidobacterium


```

Some clades may present many polytomies. These can be resolved constructing trees based on the 16S rRNA gene sequences from RDP database.

#### Retrieve information on the loci counts and heritability

These can be obtained from the SigLocusCount table.
```
for tax in `cat all_taxa.txt`
do
   if grep -q $tax"_HurdleTruncRNT_Residuals" SNPbasedHeritability_SigLocusCount_v1_BOLTLMM.txt
      then grep $tax"_HurdleTruncRNT_Residuals" SNPbasedHeritability_SigLocusCount_v1_BOLTLMM.txt | cut -f 1,4
   else 
      echo -e $tax" \t0"
   fi
done > loci_counts.txt
```

#### Upload data into iTOL and export tree
Load your .tree and .txt files with loci counts and heritability to https://itol.embl.de/ and format your tree [as desired](https://itol.embl.de/help.cgi)

## Companion website and supplementary data

Companion website and supplementary data for:

> Hughes DA, Bacigalupe R, Wang J, Rühlemann MC, Tito RY, Falony G, Joossens M, Vieira-Silva S, Henckaerts L, Rymenans L, Verspecht C, Ring S, Franke A, Wade KH, Timpson NJ, Raes J. **Genome-wide associations of human gut microbiome variation and implications for causal inference analyses**. *Nature Microbiology* 5(9):1079–1087 (2020). [doi:10.1038/s41564-020-0743-8](https://doi.org/10.1038/s41564-020-0743-8)

It contains the main figures (PNG previews and PDF), the supplementary tables (xlsx), the supplementary information (PDF), and the FGFP data access procedure that is referred to in the article.

**Website:** https://raeslab.github.io/fgfp-gwas/ (GitHub Pages, served from `docs/` on the `master` branch)

This material was previously hosted on the Raes Lab website at `raeslab.org/companion/fgfp-gwas/` and has been migrated to this repository.

### Contents

| Path | Size | Description |
|------|------|-------------|
| `supplementary/FGFP_SI.pdf` | 3.5 MB | Supplementary Information (methods, supplementary figures and tables 6, 15 and 16). |
| `supplementary/figures/Figure_1.pdf` | 81 KB | Figure 1 – MT-associated loci and heritability (PDF) |
| `supplementary/figures/Figure_2.pdf` | 3.1 MB | Figure 2 (i) – Genomic variants associated with microbial traits (PDF) |
| `supplementary/figures/Figure_2.pdf` | 3.1 MB | Figure 2 (ii) – Genomic variants associated with microbial traits (PDF) |
| `supplementary/tables/Supplementary_Table1.xlsx` | 10 KB | Summary of analytical methods as published in previous human, 16S rRNA microbiome GWAS studies. The table includes publication details of the studies (the title, first author, website url), their sample size, 16S rRNA gene variable region amplified and sequenced, software or pipeline used to cluster 16S reads into taxonomic and/or operational taxonomic units (OUT). In addition, the number of microbial traits (MT) analyzed in the mGWAS, transformation applied to the MT abundance data, software or test applied to perform association between genotype and MT, and methodologies used to account for population structure and for multiple testing burden. |
| `supplementary/tables/Supplementary_Table11.xlsx` | 162 KB | Results derived from PhenoScanner V2 analysis. Including all previously reported GWAS associations for the top 14 variants at a p-value threshold of 0.01 (PhenoScanner_GWAS); all previously reported expression quantitative loci associations, (PhenoScanner_eQTL) and all previously reported methylation quantitative loci associations (PhenoScanner_methQTL). |
| `supplementary/tables/Supplementary_Table12.xlsx` | 1012 KB | Meta-analysis GENE2FUNC enrichment. All enrichment results derived from GENE2FUNC (http://fuma.ctglab.nl), containing all results “GENE2FUNC_All_Enrich_Result” and only those from the GWAS Catalog and Canonical Pathways “GENE2FUNC_GWAS_&_CanonicalPath” that we used for the study. In addition, the gene expression (GTEx v7) analysis for general tissues (“GTEx_v7_GeneralTissue_DEG”) and specific tissues (“GTEx_v7_SpecificTissue_DEG”) differential expression groups are provided. Table includes the category identifier, the gene set or specific pathways, the number of genes in the gene set, the number of genes associated to meta-supported variants present in the gene set, the p-value, the Benjamin-Hochberg correct p-value, and the HUGO or ENSEMBL ids for all of the meta-supported genes in said gene set. |
| `supplementary/tables/Supplementary_Table13.xlsx` | 87 KB | Bi-directional Mendelian randomization results. Mendelian randomization (MR) results for MT of disease phenotypes and disease of MTs are provided. For MT of disease only a single SNP was used as the instrument, so Wald Ratio results are presented for those MT in Table 1 and for Bifidobacterium AB separately. Tables include the exposure (MT) outcome (disease) the number of SNPs used as instruments, the model, and parameter estimates beta (presented as log odd ratios for binary outcomes and standard deviation units of change for continuous), standard error and p-value. Disease-on-MT MR analysis include results for inverse variance weighted, weighted median, weighted mode, and MR Egger methods, as well as additional results tables for Heterogeneity and Pleiotropy sensitivity analysis. |
| `supplementary/tables/Supplementary_Table14.xlsx` | 33 KB | SNP rids for all genetic variants used to proxy disease phenotypes in MR analyses. |
| `supplementary/tables/Supplementary_Table17.xlsx` | 24 KB | FGFP cohort summary statistics for analyzed taxa. Table includes taxon name (taxa), mean abundance (mean), standard deviation, Shapiro’s W-statistic for the raw, log2, and rank normal data transformed distributions, and accompanying p-values for the Shapiro normality test. In addition, the proportion of individuals with zero counts (equivalent to 1-prevelance) and a flag (0 = no, 1 = yes) declaring if a taxon should be processed through the 2-step hurdle binary (HB) model, including a P/A and zero-truncated abundance trait. Finally, for HB taxon there is also an estimate of the zero-truncated mean abundance and its standard deviation. |
| `supplementary/tables/Supplementary_Table18.xlsx` | 244 KB | Cohort microbiome summary statistics. The taxa in FGFP, FOCUS and PopGen are included, with information on their mean abundances, standard deviation, Shapiro’s W-statistic for the raw, log2, and rank normal data transformed distributions, and accompanying p-values for the Shapiro normality test. In addition, the proportion of individuals with zero counts (equivalent to 1-prevelance) and a flag (0 = no, 1 = yes) declaring if a taxon should be processed through the 2-step hurdle binary (HB) model, including a P/A and zero-truncated abundance trait. Finally, for HB taxon there is also an estimate of the zero-truncated mean abundance and its standard deviation. |
| `supplementary/tables/Supplementary_Table19.xlsx` | 44 KB | DEPICT gene enrichment results for individuals MTs. This table includes GO, KEGG, Reactome, and Ensembl gene subnetworks. The gene set ID and functional category source along with description of the functional category gene set are included. |
| `supplementary/tables/Supplementary_Table2.xlsx` | 7.8 MB | Microbial trait data used in this study. The table includes individual level DADA2 integer count data for all phylogenetic levels (highlighted in green); the (zero-truncated)-abundance data rank normal transformed and regressed (RNTRes) against technical batch covariates, PC1-10, sex, and age (residuals extracted) (in blue); the derived presence (1) or absence(0) (P/A) hurdle binary (HB) data for zero-inflated counts (yellow); and the diversity metrics, untransformed or regressed and enterotype calls. |
| `supplementary/tables/Supplementary_Table20.xlsx` | 74 KB | Tissue enrichment results for individuals MTs. Tissue specific DEPICT enrichment, MeSH refers to Medical Subject Headings. |
| `supplementary/tables/Supplementary_Table21.xlsx` | 11 KB | DEPICT gene priority enrichment results for individuals MTs. All variants with a p-value less than 1x10-5 in the initial FGFP GWAS screen are included. |
| `supplementary/tables/Supplementary_Table22.xlsx` | 322 KB | Garfield enrichment results for MTs. Tissue specific epigenetic regulatory function enrichment analysis derived from GARFIELD using FGFP variants at an association smaller than p-value of 1x10-5. The table includes MT identifier, odds ratio, p-value, Benjamini-Hochberg correct p-value, along with the specific regulatory annotation, cell type, tissue and epigenetic category of the association. |
| `supplementary/tables/Supplementary_Table3.xlsx` | 49 KB | GWAS taxa (92) result summaries. The table contains the 92 taxon used in the mGWAS with their corresponding phylogenetic description (phylum, class, order, family, genus), prevalence and abundance in FGFP, the number of sites and loci associated with each taxon derived from the meta-analysis, and from FGFP alone. In addition, the heritability point estimates, standard errors, and p-values for presence/absence (P/A), rank normal transformed abundances (AB), as well as log and box-cox transformations of the integer count data for cross study comparison. |
| `supplementary/tables/Supplementary_Table5.xlsx` | 1.9 MB | Results and annotation for all 3,321 meta-supported associations. Meta-supported associations refer to those with smaller p-value in the meta-analysis compared to FGFP alone. The table includes the microbial traits, rsids, snpids, a flag (TagLocus) indicating if a variant is a locus tag or top result in a linkage disequilibrium block, parameter estimates for the meta-analysis, FGFP, FoCUS and PopGen including, beta, se, association p-value, sample size, variant info score, minor allele frequency, effect allele frequency, and annotation of variants with DEPICT, and ENSEMBL v75 via biomaRt. |
| `supplementary/tables/Supplementary_Table7.xlsx` | 153 KB | Catalog of MT-SNP matched hits previously reported in mGWAS associations. The table provides the estimates for the like or a higher order taxonomic level microbial trait, when available, from the FGFP (-method score) alongside the publicly available data from previous studies. |
| `supplementary/tables/Supplementary_Table8.xlsx` | 133 KB | Catalog of top hits previously reported in mGWAS associations. Estimates for the STONGEST associated (smallest p-value, REGARDLESS OF THE TRAIT) microbial trait, when available, from the FGFP (-method score) is provided alongside the publicly available data from previous studies. |
| `supplementary/tables/Supplementary_Table9.xlsx` | 23 KB | The effect estimates for all FGFP microbial traits at the MCM6 rs4988235 variant. The table includes MT as defined by FGFP, the rsid in FGFP, the SNPID in FGFP, defined as a compound identifier of chromosome:basepair_allele1_allele2, chromosome as reported by FGFP, base position as reported by FGFP (hg19), reference allele in FGFP, effect allele in FGFP, minor allele frequency in FGFP, imputation information score in FGFP, beta or effect estimate for the MT, standard error in FGFP, p-value in FGFP. |

Total supplementary data: 21.9 MB.

#### Data access

Summary statistics (≈110 GB) are available at <https://doi.org/10.5523/bris.22bqn399f9i432q56gt3wfhzlc>.

The process of enquiry for data access is as follows: upon data request by email to jeroen.raes@kuleuven.be, the
FGFP data access committee will evaluate access permission, which will be granted upon signature of a data use
agreement between the governing legal entities.

#### Notes

* Supplementary Tables 6, 15 and 16 are part of the Supplementary Information PDF.
* Supplementary Tables 4 and 10 were listed on the original website but the files were not present in the archived site.

### Repository layout

```
docs/                 Static, self-contained website (GitHub Pages source; no external CSS/JS/fonts)
  index.html
  assets/             Stylesheet, script, logo, figure previews
supplementary/        Downloadable data files
```

#### Data files

The files in `supplementary/` are stored as normal git files (≈ 22 MB, largest 8 MB). Git LFS is not used because
GitHub does not accept LFS uploads to public forks, and this repository is a fork of `kul-fgfpgwas/rep-cookbook`.
The website's download buttons point to `https://github.com/raeslab/fgfp-gwas/raw/master/supplementary/...`;
when `docs/index.html` is opened from a local checkout the links switch to the relative `../supplementary/` paths.

#### Publishing

In the repository settings, enable GitHub Pages with source *Deploy from a branch* → `master` / `/docs`.

### How to cite

If you use these data, please cite:

> Hughes DA, Bacigalupe R, Wang J, Rühlemann MC, Tito RY, Falony G, Joossens M, Vieira-Silva S, Henckaerts L, Rymenans L, Verspecht C, Ring S, Franke A, Wade KH, Timpson NJ, Raes J. **Genome-wide associations of human gut microbiome variation and implications for causal inference analyses**. *Nature Microbiology* 5(9):1079–1087 (2020). [doi:10.1038/s41564-020-0743-8](https://doi.org/10.1038/s41564-020-0743-8)

```bibtex
@article{hughes2020genome,
  title   = {{Genome-wide associations of human gut microbiome variation and implications for causal inference analyses}},
  author  = {Hughes, D. A. and Bacigalupe, R. and Wang, J. and Rühlemann, M. C. and Tito, R. Y. and Falony, G. and Joossens, M. and Vieira-Silva, S. and Henckaerts, L. and Rymenans, L. and Verspecht, C. and Ring, S. and Franke, A. and Wade, K. H. and Timpson, N. J. and Raes, J.},
  journal = {Nature Microbiology},
  year    = {2020},
  volume  = {5},
  number  = {9},
  pages   = {1079--1087},
  doi     = {10.1038/s41564-020-0743-8}
}
```

### Contact

Raes Lab, VIB–KU Leuven Center for Microbiology – Prof. Jeroen Raes (jeroen.raes@kuleuven.be).
