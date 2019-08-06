# DKFZ/ODCF Reference File Setup Scripts

Code to set up reference genomes and associated reference data for the DKFZ/ODCF workflows. This includes assemblies used in the OTP platform.

> Former eilslabs reference data, including much of what is used in the OTP platform, is notoriously undocumented. In many cases, no scripts or any other information are available to set up these assemblies.

## General Remarks on Human Genome Assemblies

### Sequence Content 
  
In general, there are a number of special FASTA entries in downloaded reference genome files. The actual chromosomes (with centromers and everything) are just 1 to 22, X, and Y. Beyond these you can find:

  * **"random" contigs** can be located to a specific chromosome but not fitted in at a specific place
  * **"unplaced" contigs** can not even be assigned to a chromosome
  * **"alt" contigs** are sequences from alternative haplotypes show some degree of variation among humans
  * **human leukocyte antigen (HLA)** sequences are highly variable regions
  * **phiX** is a bacteriophage genome that is frequently used for color calibration in the Illumina sequencers ([1](https://support.illumina.com/bulletins/2017/02/what-is-the-phix-control-v3-library-and-what-is-its-function-in-.html), [2](http://www.support.illumina.com/content/dam/illumina-marketing/documents/products/technotes/hiseq-phix-control-v3-technical-note.pdf), [3](https://support.illumina.com/sequencing/sequencing_software/igenome.html)). We used the RTA build.
  * **lambda** bacteriophage genome ...
  * **Human herpesvirus 4 (EBV)** us usually a decoy sequence
  * **other human viruses** also constitute decoy sequences (e.g. in the GDC reference genome)
  * other **decoy sequences** include for instance "hs37d5" (for hg19) or "hs38d1" (for hg38)

### Assembly Identifiers

In general there is the following relation between human assemblies: 
  * `hg19 | decoys = hs37d5 = 1KGRef`. We will refer to `hg19` as the human base assembly `hg19` without any decoys. If decoys are added,  we use either `hs37d5` or `1KGRef` but will prefer `hs37d5`.
  * `hg19 = GRCh37`
  * `hg38 = GRCh38`
 
## Legacy Reference Data

The identifiers of assemblies used in the logging output and emails of OTP contain suffixes and prefixes

  * bwa06,  bwa07,  etc.: BWA version used to build the index. Note that many patch versions (last number) are compatible. Changes happened at arbitrary version numbers. In particular, a new index format was introduced at 0.6.9.  
  * methylCtools: Bisulfite-converted assembly used for WGBS data analysis 
  * suffixes for specific sequences included in the FASTA, e.g. phiX or lambda
  * hD3A: ???

Unfortunately, no systematic naming format was used in the past. Therefore, the following table with general information about the assemblies can be helpful:

| ngs_share identifier | identifier | chromosomes  | description  |
|--------|---------------|--------------|--------------|
| TBD            | bwa06_1KGRef | 1-22, X, Y, MT, GL*, NC_007605, hs37d5 | TBD |
| TBD            | bwa06_1KGRef_PhiX | TBD | TBD |
| TBD            | bwa06_hg19_chr | TBD | TBD |
| TBD            | bwa06_methylCtools_hs37d5_PhiX_Lambda | TBD | TBD |
| | | | |
| TBD            | bwa06_hg38 | 1-22, X, Y, MT, KI*, GL* | TBD |
| legacy.hg_GRCh38      | bwa06_hg38_CGA_000001405.15-no_alt_analysis_set | **with** chr-prefix: 1-22, X, Y, MT, KI*, GL*, EBV | Originally created for the DEEP project by Barbara Hutter. |
| TBD            | bwa07_hg38_PhiX | 1-22, X, Y, MT, KI*, GL*, phiX174 | TBD |
| TBD            | legacy.hg_GRCh38 | 1-22, X, Y, MT, KI*, GL* | TBD |
| TBD            | bwa06_methylCtools_hg38_PhiX_Lambda | TBD | TBD |
| TBD            | bwa06_methylCtools_hg38_PhiX_Lambda_Benchmark | TBD | TBD |
| | | | |
| TBD            | bwa06_hs37d5_Bovine_Phix | TBD | Why is _Bos taurus_ in here? |
| | | | |
| TBD            | bwa06_hs37d5_GRCm38mm | TBD | For xenograft mapping. Only mouse chromosomes with chr-prefix. |
| TBD            | bwa06_hs37d5_GRCm38mm_PhiX | TBD | For xenograft mapping. Only mouse chromosomes with chr-prefix. |
| TBD            | bwa06_hg19_GRCh37|mm10 | TBD | For xenograft mapping. Only mouse chromosomes with chr-prefix. |
| TBD            | bwa06_methylCtools_hs37d5_GRCm38mm10_PhiX_Lambda | TBD | For xenograft mapping. Only mouse chromosomes with chr-prefix. |
| | | | |
| TBD            | bwa06_GRCm38mm10 | TBD | TBD |
| TBD            | bwa06_GRCm38mm10_PhiX | TBD | TBD |
| TBD            | bwa06_GRCm38mm10_PhiX_hD3A | TBD | TBD |
| TBD            | bwa06_methylCtools_GRCm38mm10_PhiX_Lambda | TBD | TBD |
| TBD            | bwa06_methylCtools_GRCm38mm10_PhiX_Lambda_hD3A | TBD | TBD |
| TBD            | bwa06_methylCtools_mm10_UCSC_PhiX_Lambda | TBD | TBD |
| TBD            | bwa06_methylCtools_mm9_PhiX_Lambda | TBD | TBD |
| TBD            | bwa06_mm9_PhiX_Lambda_old | TBD | TBD |
| TBD            | hg_GRCm38 | TBD | TBD |

Currently, OTP assemblies are installed in `/icgc/dkfzlsdf/project/ODCF/reference_genomes/` and "ngs_share" assemblies in `/icgc/ngs_share/assemblies/` or its subdirectories.

## Newer Reference Data

For newer setups of reference data the directory structure was change. For instance, index-directories now contain full version numbers, e.g. `bwa-0.7.15`. Furthermore, assemblies are grouped by their [primary assembly](https://www.ncbi.nlm.nih.gov/grc/help/definitions/) -- the actual chromosomes of the organism plus the unplaced and unlocalized sequences representing a non-redundant haploid genome.

| ngs_share identifier | identifier | chromosomes  | description  |
|--------|---------------|--------------|--------------|
| - | GRCh38/GRCh38_decoy_ebv_phiX  | 1-22, X, Y, M, random, unplaced, phix | "chr" prefixes were dropped. |
| - | GRCh38/GRCh38_decoy_ebv_phiX_alt_hla | 1-22, X, Y, M, random, unplaced, alt, hla, phix | "chr" prefixes were dropped. Without phix this is the same assembly as the one used ICGC-ARGO (checked via Picard NormalizeFasta and md5sum) and by the 1000 Genomes Project, BROAD and KidsFirst project at CHOPS ([according to personal communication by Junjun Zhang](https://object.cancercollaboratory.org:9080/swift/v1/genomics-public-data/reference-genome/GRCh38_hla_decoy_ebv/GRCh38_hla_decoy_ebv.fa). |
| - | GRCh38/GDC_GRCh38_d1_vd1_phiX | 1-22, X, Y, M, random, unplaced, viruses, phix | see [here](https://gdc.cancer.gov/about-data/data-harmonization-and-generation/gdc-reference-files), phiX was added. "chr" prefixes were dropped. | 

## What does it all mean?

Please refer to one of the following sides for more information:

  * https://software.broadinstitute.org/gatk/blog?id=8180
  * https://software.broadinstitute.org/gatk/documentation/article?id=8017
  * https://software.broadinstitute.org/gatk/documentation/article?id=8017
  * http://sourceforge.net/p/bio-bwa/mailman/message/32600693/
  
# Installation and Usage

First, install the Conda environment:

```bash
conda env create -n setup-reference-data -f "$repoRoot/conda.yml"
```

You can then install a specific assembly and associated reference data by calling the `prepare.sh` script in the corresponding directory in the repository. The currently maintained assemblies are in the directory named like the primary assembly (e.g. `GRCh38`) and a variant subdirectory (e.g. `GRCh38_decoy_ebv_phiX`).

Note that the scripts try to reduce downloads by saving caching the downloaded files in a `cache/` directory, in particular if you want to build reference files for multiple related assemblies. Obviously, this directory is not automatically removed after running the scripts.
 
The scripts for the legacy assemblies should be in the `src/legacy` directory as soon as they are written. Until then only general information or protocol information are put into these directories (if possible).