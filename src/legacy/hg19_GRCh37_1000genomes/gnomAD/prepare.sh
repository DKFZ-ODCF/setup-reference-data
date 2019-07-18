#!/usr/bin/env bash
# Original by Barbara Hutter, 2013

# To download
wget -c https://data.broadinstitute.org/gnomAD/release-170228/gnomAD.release.170228.tar

# Untar
tar -xcf gnomAD.release.170228.tar

# For concatenated WGS file

(seq 1 22 ; echo X ) | while read a; do ls /icgc/dkfzlsdf/analysis/UniExome/GermlineAnnotation/annotationInfo/gnomAD/web/macarthurlab-distribution/gnomAD/release-170228/genomes/vcf/gnomad.genomes.r2.0.1.sites.$a.vcf.gz ;done | tr '\n' ' '
vcf-concat-0.1.10 /icgc/dkfzlsdf/analysis/UniExome/GermlineAnnotation/annotationInfo/gnomAD/web/macarthurlab-distribution/gnomAD/release-170228/genomes/vcf/gnomad.genomes.r2.0.1.sites.1.vcf.gz /icgc/dkfzlsdf/analysis/UniExome/GermlineAnnotation/annotationInfo/gnomAD/web/macarthurlab-distribution/gnomAD/release-170228/genomes/vcf/gnomad.genomes.r2.0.1.sites.2.vcf.gz /icgc/dkfzlsdf/analysis/UniExome/GermlineAnnotation/annotationInfo/gnomAD/web/macarthurlab-distribution/gnomAD/release-170228/genomes/vcf/gnomad.genomes.r2.0.1.sites.3.vcf.gz /icgc/dkfzlsdf/analysis/UniExome/GermlineAnnotation/annotationInfo/gnomAD/web/macarthurlab-distribution/gnomAD/release-170228/genomes/vcf/gnomad.genomes.r2.0.1.sites.4.vcf.gz /icgc/dkfzlsdf/analysis/UniExome/GermlineAnnotation/annotationInfo/gnomAD/web/macarthurlab-distribution/gnomAD/release-170228/genomes/vcf/gnomad.genomes.r2.0.1.sites.5.vcf.gz /icgc/dkfzlsdf/analysis/UniExome/GermlineAnnotation/annotationInfo/gnomAD/web/macarthurlab-distribution/gnomAD/release-170228/genomes/vcf/gnomad.genomes.r2.0.1.sites.6.vcf.gz /icgc/dkfzlsdf/analysis/UniExome/GermlineAnnotation/annotationInfo/gnomAD/web/macarthurlab-distribution/gnomAD/release-170228/genomes/vcf/gnomad.genomes.r2.0.1.sites.7.vcf.gz /icgc/dkfzlsdf/analysis/UniExome/GermlineAnnotation/annotationInfo/gnomAD/web/macarthurlab-distribution/gnomAD/release-170228/genomes/vcf/gnomad.genomes.r2.0.1.sites.8.vcf.gz /icgc/dkfzlsdf/analysis/UniExome/GermlineAnnotation/annotationInfo/gnomAD/web/macarthurlab-distribution/gnomAD/release-170228/genomes/vcf/gnomad.genomes.r2.0.1.sites.9.vcf.gz /icgc/dkfzlsdf/analysis/UniExome/GermlineAnnotation/annotationInfo/gnomAD/web/macarthurlab-distribution/gnomAD/release-170228/genomes/vcf/gnomad.genomes.r2.0.1.sites.10.vcf.gz /icgc/dkfzlsdf/analysis/UniExome/GermlineAnnotation/annotationInfo/gnomAD/web/macarthurlab-distribution/gnomAD/release-170228/genomes/vcf/gnomad.genomes.r2.0.1.sites.11.vcf.gz /icgc/dkfzlsdf/analysis/UniExome/GermlineAnnotation/annotationInfo/gnomAD/web/macarthurlab-distribution/gnomAD/release-170228/genomes/vcf/gnomad.genomes.r2.0.1.sites.12.vcf.gz /icgc/dkfzlsdf/analysis/UniExome/GermlineAnnotation/annotationInfo/gnomAD/web/macarthurlab-distribution/gnomAD/release-170228/genomes/vcf/gnomad.genomes.r2.0.1.sites.13.vcf.gz /icgc/dkfzlsdf/analysis/UniExome/GermlineAnnotation/annotationInfo/gnomAD/web/macarthurlab-distribution/gnomAD/release-170228/genomes/vcf/gnomad.genomes.r2.0.1.sites.14.vcf.gz /icgc/dkfzlsdf/analysis/UniExome/GermlineAnnotation/annotationInfo/gnomAD/web/macarthurlab-distribution/gnomAD/release-170228/genomes/vcf/gnomad.genomes.r2.0.1.sites.15.vcf.gz /icgc/dkfzlsdf/analysis/UniExome/GermlineAnnotation/annotationInfo/gnomAD/web/macarthurlab-distribution/gnomAD/release-170228/genomes/vcf/gnomad.genomes.r2.0.1.sites.16.vcf.gz /icgc/dkfzlsdf/analysis/UniExome/GermlineAnnotation/annotationInfo/gnomAD/web/macarthurlab-distribution/gnomAD/release-170228/genomes/vcf/gnomad.genomes.r2.0.1.sites.17.vcf.gz /icgc/dkfzlsdf/analysis/UniExome/GermlineAnnotation/annotationInfo/gnomAD/web/macarthurlab-distribution/gnomAD/release-170228/genomes/vcf/gnomad.genomes.r2.0.1.sites.18.vcf.gz /icgc/dkfzlsdf/analysis/UniExome/GermlineAnnotation/annotationInfo/gnomAD/web/macarthurlab-distribution/gnomAD/release-170228/genomes/vcf/gnomad.genomes.r2.0.1.sites.19.vcf.gz /icgc/dkfzlsdf/analysis/UniExome/GermlineAnnotation/annotationInfo/gnomAD/web/macarthurlab-distribution/gnomAD/release-170228/genomes/vcf/gnomad.genomes.r2.0.1.sites.20.vcf.gz /icgc/dkfzlsdf/analysis/UniExome/GermlineAnnotation/annotationInfo/gnomAD/web/macarthurlab-distribution/gnomAD/release-170228/genomes/vcf/gnomad.genomes.r2.0.1.sites.21.vcf.gz /icgc/dkfzlsdf/analysis/UniExome/GermlineAnnotation/annotationInfo/gnomAD/web/macarthurlab-distribution/gnomAD/release-170228/genomes/vcf/gnomad.genomes.r2.0.1.sites.22.vcf.gz /icgc/dkfzlsdf/analysis/UniExome/GermlineAnnotation/annotationInfo/gnomAD/web/macarthurlab-distribution/gnomAD/release-170228/genomes/vcf/gnomad.genomes.r2.0.1.sites.X.vcf.gz > /icgc/dkfzlsdf/analysis/UniExome/GermlineAnnotation/annotationInfo/gnomAD/web/macarthurlab-distribution/gnomAD/release-170228/genomes/vcf/gnomad.genomes.r2.0.1.sites.1-22-X.vcf

bgzip /icgc/dkfzlsdf/analysis/UniExome/GermlineAnnotation/annotationInfo/gnomAD/web/macarthurlab-distribution/gnomAD/release-170228/genomes/vcf/gnomad.genomes.r2.0.1.sites.1-22-X.vcf
tabix -p vcf /icgc/dkfzlsdf/analysis/UniExome/GermlineAnnotation/annotationInfo/gnomAD/web/macarthurlab-distribution/gnomAD/release-170228/genomes/vcf/gnomad.genomes.r2.0.1.sites.1-22-X.vcf.gz


# FileProcessing
python /home/paramasi/projects/repository/ngs2/trunk/pipelines/TRIO_analysis_pipeline/dataFileProcessing/gnomAD_processing.toPipeFormat.py

## File details
## WES files processed for AF
#gnomad.exomes.r2.0.1.sites.INDELs.vcf.gz
#gnomad.exomes.r2.0.1.sites.SNVs.vcf.gz
#gnomad.exomes.r2.0.1.sites.vcf.gz # Origin file
#
## WGS files processed for AF and also file only with common alleles (MAF > 0.001)
#gnomad.genomes.r2.0.1.sites.1-22-X.INDELs.Common.vcf.gz
#gnomad.genomes.r2.0.1.sites.1-22-X.INDELs.vcf.gz
#gnomad.genomes.r2.0.1.sites.1-22-X.SNVs.Common.vcf.gz
#gnomad.genomes.r2.0.1.sites.1-22-X.SNVs.vcf.gz
#gnomad.genomes.r2.0.1.sites.1-22-X.vcf.gz # Origin file


