# Evaluation of copy number variation in *Lottia gigantea* for the investigation of their northward range expansion
Project Managed by: Julianna Porter
PI: Rachael Bay

## Overview of Methods
I called copy number variants using two programs: Delly and Manta. Afterwards, I combined the results using the program SURVIVOR and only retained variants called by both callers with the same genotype. Additionally, variants were filtered for a maximum no call fraction of 20% and a minor allele frequency of 0.05. This was done so CNV results were comparable to SNP results. CNV results were then used for PCA, Structure, and Genotype x Environment Associations. 
**Resources:**
Delly [Documentation and Gtihub]
 (https://github.com/dellytools/delly)
Manta [Github]
 (https://github.com/Illumina/manta)
Manta [Documentation]
 (https://github.com/Illumina/manta/blob/master/docs/userGuide/README.md)
SURVIVOR [Github]
 (https://github.com/fritzsedlazeck/SURVIVOR)
SURVIVOR [Documentation]
 (https://github.com/fritzsedlazeck/SURVIVOR/wiki)
## Pipeline for Calling Copy Number Variants
### Step 0: Test Manta on all samples
**See scripts:**
    manta.prb.chd.sh
**Description:**
When running manta, a select few BAM files gave an error that caused Manta to exit (11/559 samples). The issue was with a single read in the BAM file; I removed the problematic read and Manta ran just fine. I ran Manta on individual BAM files, and for all the files that gave an exit code of 1, I extracted the problematic read from the error log and removed it from the BAM file using picard FilterSamReads. 
**Excerpt from error log:**
    Failed on command: ‘makeLocusGraph’
    FATAL_ERROR: 2023-Mar-28 18:27:50 /ocean/projects/deb200006p/jporter/install/manta-1.6.0.release_src/src/c++/lib/manta/SVLocusScanner.cpp(1204): Throw in function void getSVLociImpl
    Dynamic exception type: boost::exception_detail::clone_impl<illumina::common::GeneralException>
    std::exception::what: Unexpected breakend pattern proposed from bam record.
    local_breakend: Breakend: GenomeInterval: 114:[0,-90) COMPLEX
    SVBreakendLowResEvidence: pair: 0 local_pair: 0 cigar: 0 softclip: 0 semialign: 1 shadow: 0 split_align: 0 unknown: 0
    remote_breakend: Breakend: GenomeInterval: 114:[0,-90) UNKNOWN
    SVBreakendLowResEvidence: pair: 0 local_pair: 0 cigar: 0 softclip: 0 semialign: 0 shadow: 0 split_align: 0 unknown: 0
    bam_record: E100048427L1C017R04003413207/2 tid:pos:strand 114:20:+ cigar: 20I130S templSize: 144 mate_tid:pos:strand 114:20:-
Key Info: The bam_record name. Copy that name (without the /2) into a file that will be used to exclude that read
### Step 1: Split genome into segments
**See scripts:**
    0.split.genome.regions.sh
    1.fasta2bed.sh
**Description:**
I split the genome into multiple regions to reduce computational load when calling copy number variants on all samples. Delly wants a file with regions to exclude on each line (created using script 0). I extracted scaffold names from the reference index (file ending in fa.fai), split them into ten separate genome regions, and then for each region I used the comm command to create a file that contained all scaffolds EXCEPT the ones for a given region. Manta wants a sorted BED file of regions to include that is bg-zipped and tabix indexed (created using script 1). I converted the reference fasta to a BED file using faidx and split that BED file into ten regions. 
### Step 2: Run CNV callers on all samples by region
**See scripts:**
    2.delly.ALL.sh
    3.manta.ALL.sh
**Description:**
I called copy number variants using two different programs: Delly and Manta. Both took BAM files and the reference genome as input. For Delly only INS, DEL and DUP were analyzed. Both output VCF files. 
### Step 3: Process CNV calls for merging regions
**See script:**
    4.preprocess.sh
**Description:**
For Delly calls:
I used bcftools norm to ensure that the reference allele in the VCF file matches the fasta file. The VCF is then bg-zipped and tabix indexed.
For Manta calls:
Files are already bg-zipped so they are only tabix-indexed. I also renamed the files to include region number and moved them to a common directory. Reason: all files generated by Manta have the same output names located in the working directory created during Manta’s configuration.
### Step 4: Merge all regions of the respective callers
**See script:**
    5.merge.regions.sh
**Description:**
I used bcftools concat to combine all genome regions for delly and manta, respectively. 
### Note: At this step I renamed the samples in the VCF using bcftools reheader. See script: 6.rename.samples.sh
This was to make following steps that extract sample names easier.
### Step 5: Split each vcf by sample for merging with SURIVOR
Reason for this approach: when merging Delly and Manta with multiple samples, the resulting VCF did not retain all samples. Instead, it would have two sample columns for the first individual (one from each caller). So the solution was to split each VCF by individual, and then merge individual Delly and Manta results. 
**See script:**
    7.split.vcf.sh
**Description:**
Use bcftools +split plugin for splitting VCF by sample. The output option creates a new directory where the files are named as {sample}.vcf 
### Step 6: Merge Delly and Manta results by sample using SURVIVOR
Note: the following steps were run using Snakemake. See the snakefile: rules/cnv.processing.smk for the commands. 
**See rules:**
    gunzip
    surv_callset
    survivor
**Description**
VCFs following the per sample split were bg-zipped and tabix indexed. SURVIVOR requires uncompressed VCF files, so all VCF files were unzipped using bgzip -d. SURVIVOR takes as input a file that contains the paths to the input files to be merged on new lines, which was generated in the rule surv_callset. For SURVIVOR, I used options '1000 2 1 0 0 50'. This meant SURVIVOR allowed a maximum distance of 1kb between variants merged, required variants retained to have support from 2 callers and agree on type (1 in third position, otherwise 0), did not need variants to agree on strand (0 in fourth position), and had a minimum variant length of 50 base pairs. Each sample then had a vcf with two sample columns: one with Delly results, and another with Manta results. 
### Step 7: Filter for matching genotype between Manta and Delly calls
**See rules:**
    match_gt
    filt_prep
    filt_gt
**Description:**
Using awk scripting, I extracted the scaffold, position, and the genotype from either the Delly or Manta column from each VCF separately. Using grep I got the names of loci that had the same genotype for Manta and Delly for a given sample. I then formatted these loci into a tab-delimited file that could be used to filter the per sample VCFs that contained both Delly and Manta results. With bcftools view, I filtered the per sample VCFs to only keep loci with the same genotype for Manta and Delly, and removed the Manta sample column, leaving only one sample column per VCF (the one from Delly).
### Step 8: Merge all per-sample VCFs into a single VCF
**See rules:**
    merge_prep
    check_ref
    merge_geno
**Description:**
In order to merge all VCFs using bcftools I had to first bgzip and tabix index all VCFs. I then confirmed that the reference allele in the VCFs matched that in the reference genome using bcftools norm. I then merged all the per-sample VCFs into a single BCF using bcftools merge.
### Step 9: Filter VCF of all CNVs
**See rules:**
    max_no_call_loose
    max_no_call_strict
    maf
**Description:**
Using GATK SelectVariants I filtered loci based on a maximum no call fraction. I tried a max no call fraction of 0.20 and 0.50. I also filtered out monomorphic sites by requiring a minor allele frequency of 0.05. I filtered based on minor allele frequency using vcftools. 
## Analyses carried out with Copy Number Variants
**Note** For the following analyses the bash scripts were run on Bridges2 however the R scripts were run locally on my computer. For that reason the file paths will differ between the two.
### PCA
**See scripts:**
    ld.prune.sh
    pca.sh
**Description:**

### Structure
**See scripts**
    install.structure.sh
    structure.sh
    hierarch.struc.sh
**Description**

### Genotype by Environment Analysis
**See scripts**
    extract.gt.sh
**Description:**

### Pairwise Fst for all populations
**See scripts**

**Description:**
