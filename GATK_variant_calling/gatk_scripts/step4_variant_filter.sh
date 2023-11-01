#---------------------------------------------------
# Exclude variants marked by flags other than "PASS"
#---------------------------------------------------
../gatk-4.1.9.0/gatk SelectVariants \
    -V ../wgs/cohort.unfiltered.gt.vcf.gz \
    --exclude-filtered true\
    -O ../wgs/cohort.filtered.raw.gt.vcf.gz
num_variants_step1=$(bcftools view -H ../wgs/cohort.filtered.raw.gt.vcf.gz | wc -l)

#------------------------------------------
# select variants: SNPs, INDELs, Mixed-type
#------------------------------------------
../gatk-4.1.9.0/gatk SelectVariants \
    -V ../wgs/cohort.filtered.raw.gt.vcf.gz \
    -select-type SNP \
    --exclude-filtered true\
    -O ../wgs/cohort.filtered.raw.gt.snp.vcf.gz
num_snps_step1=$(bcftools view -H ../wgs/cohort.filtered.raw.gt.snp.vcf.gz | wc -l)

../gatk-4.1.9.0/gatk SelectVariants \
    -V ../wgs/cohort.filtered.raw.gt.vcf.gz \
    -select-type INDEL \
    --exclude-filtered true\
    -O ../wgs/cohort.filtered.raw.gt.indel.vcf.gz
num_indels_step1=$(bcftools view -H ../wgs/cohort.filtered.raw.gt.indel.vcf.gz | wc -l)

../gatk-4.1.9.0/gatk SelectVariants \
    -V ../wgs/cohort.filtered.raw.gt.vcf.gz \
    -select-type MIXED \
    --exclude-filtered true\
    -O ../wgs/cohort.filtered.raw.gt.mixed.vcf.gz
num_mixed_step1=$(bcftools view -H ../wgs/cohort.filtered.raw.gt.mixed.vcf.gz | wc -l)

#------------------------------------------
# Hard-filtering using GATK recommendations
# For SNP, label SNP clusters and filter
#------------------------------------------
../gatk-4.1.9.0/gatk VariantFiltration \
  -V ../wgs/cohort.filtered.raw.gt.snp.vcf.gz \
  -filter "QD < 2.00" --filter-name "QD2" \
  -filter "QUAL < 30.0" --filter-name "QUAL30" \
  -filter "SOR > 3.00" --filter-name "SOR3" \
  -filter "FS > 60.00" --filter-name "FS60" \
  -filter "MQ < 40.0" --filter-name "MQ40" \
  -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
  -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
  -O ../wgs/cohort.filtered.hf.gt.snp.vcf.gz

../gatk-4.1.9.0/gatk VariantFiltration \
  -V ../wgs/cohort.filtered.raw.gt.indel.vcf.gz \
  -filter "QD < 2.00" --filter-name "QD2" \
  -filter "QUAL < 30.0" --filter-name "QUAL30" \
  -filter "FS > 200.0" --filter-name "FS200" \
  -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
  -O ../wgs/cohort.filtered.hf.gt.indel.vcf.gz

# mixed variants are evaluated with the indel model
../gatk-4.1.9.0/gatk VariantFiltration \
  -V ../wgs/cohort.filtered.raw.gt.mixed.vcf.gz \
  -filter "QD < 2.00" --filter-name "QD2" \
  -filter "QUAL < 30.0" --filter-name "QUAL30" \
  -filter "FS > 200.0" --filter-name "FS200" \
  -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
  -O ../wgs/cohort.filtered.hf.gt.mixed.vcf.gz

#-------------------------
# Merge hard filtered VCFs
#-------------------------
java -jar ../picard.jar MergeVcfs \
  I=../wgs/cohort.filtered.hf.gt.snp.vcf.gz \
  I=../wgs/cohort.filtered.hf.gt.indel.vcf.gz \
  I=../wgs/cohort.filtered.hf.gt.mixed.vcf.gz \
  O=../wgs/cohort.filtered.hf.tmp.gt.vcf.gz
../gatk-4.1.9.0/gatk SelectVariants \
     -V ../wgs/cohort.filtered.hf.tmp.gt.vcf.gz \
     --exclude-filtered true\
     -O ../wgs/cohort.filtered.hf.gt.vcf.gz
num_variants_step2=$(bcftools view -H ../wgs/cohort.filtered.hf.gt.vcf.gz | wc -l)

#-------------------------------------------------------------------------------------------------------------------------------------------
# Further filtering 1: minimum read depth (DP) and genotype quality (GQ), and ref-to-alt AD ratio
# Note that the -S option changes genotype of samples that do not fulfill the requirement to ./. but does not remove the entire variant site
# & and | apply multiple filters to the same sample simultaneously, while && and || apply to different samples independently
#-------------------------------------------------------------------------------------------------------------------------------------------
bcftools filter -S . -e 'FMT/DP<10 | FMT/GQ<20' -O z -o ../wgs/cohort.filtered.hf.DP10.GQ20.gt.vcf.gz ../wgs/cohort.filtered.hf.gt.vcf.gz
num_variants_step3=$(bcftools view -H ../wgs/cohort.filtered.hf.DP10.GQ20.gt.vcf.gz | wc -l)

#---------------------------------------------------------------------------------------------------------------------------------------------
# Further filtering 2: remove monomorphic SNPs/INDELs, multiallelic SNPs and indels, SNPs in the close proximity of INDELS, and INDEL clusters
# & and | apply multiple filters to the same sample simultaneously, while && and || apply to different samples independently
#---------------------------------------------------------------------------------------------------------------------------------------------
# AC==0: no variants (all the same as the reference)
# AC==AN: only alternative alleles are called. Do not include AC==AN if a site where only alternative alleles are called is considered as a variant.
# IndelGap: filter clusters of indels separated by INT or fewer base pairs allowing only one to pass
# SnpGap: filter SNPs within INT base pairs of an indel or other variant type
# bcftools view: -m2 means at least 2 alleles, -M2 means at most 2 alleles, -O z means output compressed VCF
bcftools filter -e 'AC==0 || AC==AN' --IndelGap 5 --SnpGap 10 ../wgs/cohort.filtered.hf.DP10.GQ20.gt.vcf.gz | bcftools view -m2 -M2 -O z -o ../wgs/cohort.filtered.hf.DP10.GQ20.allele.gt.vcf.gz
#bcftools filter -e 'AC==0' --IndelGap 5 --SnpGap 10 ../wgs/cohort.filtered.hf.DP10.GQ20.gt.vcf.gz | bcftools view -m2 -M2 -O z -o ../wgs/cohort.filtered.hf.DP10.GQ20.allele.gt.vcf.gz
num_variants_step4=$(bcftools view -H ../wgs/cohort.filtered.hf.DP10.GQ20.allele.gt.vcf.gz | wc -l)

#--------------------------------------------------------------
# Further filtering 3: remove variants in repetitive regions
# Note that tabix would fail if gzip, instead of bgzip, is used
#--------------------------------------------------------------
bedtools subtract -header -a ../wgs/cohort.filtered.hf.DP10.GQ20.allele.gt.vcf.gz \
    -b ../RepeatMasker/output/CD630.fasta.bed \
    > ../wgs/cohort.filtered.hf.DP10.GQ20.allele.repmask.gt.vcf
bedtools subtract -header -a ../wgs/cohort.filtered.hf.DP10.GQ20.allele.repmask.gt.vcf \
    -b ../TRF/CD630.fasta.2.5.7.80.10.50.2000.bed \
    > ../wgs/cohort.filtered.hf.DP10.GQ20.allele.repmask.trf.gt.vcf
bgzip ../wgs/cohort.filtered.hf.DP10.GQ20.allele.repmask.trf.gt.vcf
tabix -fp vcf ../wgs/cohort.filtered.hf.DP10.GQ20.allele.repmask.trf.gt.vcf.gz
num_variants_step5=$(bcftools view -H ../wgs/cohort.filtered.hf.DP10.GQ20.allele.repmask.trf.gt.vcf.gz | wc -l)

#----------------------------------------
# Further filtering 4: remove SNPclusters
#----------------------------------------
# A SNP cluster is defined as 3 SNPs within a window of 10 bases
../gatk-4.1.9.0/gatk VariantFiltration \
    -V ../wgs/cohort.filtered.hf.DP10.GQ20.allele.repmask.trf.gt.vcf.gz\
    -cluster 3 -window 10\
    -O ../wgs/cohort.filtered.hf.DP10.GQ20.allele.repmask.trf.sc.tmp.gt.vcf.gz
../gatk-4.1.9.0/gatk SelectVariants \
    -V ../wgs/cohort.filtered.hf.DP10.GQ20.allele.repmask.trf.sc.tmp.gt.vcf.gz \
    --exclude-filtered true\
    -select "FILTER == SnpCluster" --invertSelect \
    -O ../wgs/cohort.filtered.hf.DP10.GQ20.allele.repmask.trf.sc.gt.vcf.gz
num_variants_step6=$(bcftools view -H ../wgs/cohort.filtered.hf.DP10.GQ20.allele.repmask.trf.sc.gt.vcf.gz | wc -l)

#-----------------------------------------------------------
# Create vcf that includes only variants with filter PASS
#-----------------------------------------------------------
bcftools view -f PASS ../wgs/cohort.filtered.hf.DP10.GQ20.allele.repmask.trf.sc.gt.vcf.gz -O z -o ../wgs/cohort.filtered.gt.vcf.gz
tabix -fp vcf ../wgs/cohort.filtered.gt.vcf.gz
num_variants_step7=$(bcftools view -H ../wgs/cohort.filtered.gt.vcf.gz | wc -l)
num_snps_step7=$(bcftools view -H --types snps ../wgs/cohort.filtered.gt.vcf.gz | wc -l)
num_indels_step7=$(bcftools view -H --types indels ../wgs/cohort.filtered.gt.vcf.gz | wc -l)
num_mixed_step7=$(bcftools view -H --types other ../wgs/cohort.filtered.gt.vcf.gz | wc -l)

# print number of variants at each filtering step
echo "Number of variants at each filtering step:" > ../wgs/filtering_stats.txt
echo "Number of unfiltered variants: #Total=$num_variants_step1(#SNP=$num_snps_step1,#INDEL=$num_indels_step1,#MIXED=$num_mixed_step1)" >> ../wgs/filtering_stats.txt
echo "Number of variants after hard filtering: #Total=$num_variants_step2" >> ../wgs/filtering_stats.txt
echo "Number of variants after filtering by DP and GQ: #Total=$num_variants_step3" >> ../wgs/filtering_stats.txt
echo "Number of variants after filtering by number of alleles: #Total=$num_variants_step4" >> ../wgs/filtering_stats.txt
echo "Number of variants after filtering by repetitive regions: #Total=$num_variants_step5" >> ../wgs/filtering_stats.txt
echo "Number of variants after removing SNP clusters: #Total=$num_variants_step6" >> ../wgs/filtering_stats.txt
echo "Number of filtered variants: #Total=$num_variants_step7(#SNP=$num_snps_step7,#INDEL=$num_indels_step7,#MIXED=$num_mixed_step7)" >> ../wgs/filtering_stats.txt

# convert to table
../gatk-4.1.9.0/gatk VariantsToTable \
     -V ../wgs/cohort.filtered.gt.vcf.gz \
     -F CHROM -F POS -F ID -F TYPE -F REF -F ALT -F HET -F HOM-REF -F HOM-VAR -F NO-CALL -F VAR -F NSAMPLES -F NCALLED -F MULTI-ALLELIC -GF AD -GF GT \
     --error-if-missing-data \
     -O ../wgs/cohort.filtered.gt.txt

# remove temporary files
rm ../wgs/cohort.filtered.raw.gt.vcf.gz
rm ../wgs/cohort.filtered.raw.gt.snp.vcf.gz
rm ../wgs/cohort.filtered.raw.gt.indel.vcf.gz
rm ../wgs/cohort.filtered.raw.gt.mixed.vcf.gz
rm ../wgs/cohort.filtered.hf.gt.snp.vcf.gz
rm ../wgs/cohort.filtered.hf.gt.indel.vcf.gz
rm ../wgs/cohort.filtered.hf.gt.mixed.vcf.gz
rm ../wgs/cohort.filtered.hf.tmp.gt.vcf.gz
rm ../wgs/cohort.filtered.hf.gt.vcf.gz
rm ../wgs/cohort.filtered.hf.DP10.GQ20.gt.vcf.gz
rm ../wgs/cohort.filtered.hf.DP10.GQ20.allele.gt.vcf.gz
rm ../wgs/cohort.filtered.hf.DP10.GQ20.allele.repmask.gt.vcf
rm ../wgs/cohort.filtered.hf.DP10.GQ20.allele.repmask.trf.gt.vcf.gz
rm ../wgs/cohort.filtered.hf.DP10.GQ20.allele.repmask.trf.sc.tmp.gt.vcf.gz
rm ../wgs/cohort.filtered.hf.DP10.GQ20.allele.repmask.trf.sc.gt.vcf.gz
rm ../wgs/cohort.filtered.raw.gt.vcf.gz.tbi
rm ../wgs/cohort.filtered.raw.gt.snp.vcf.gz.tbi
rm ../wgs/cohort.filtered.raw.gt.indel.vcf.gz.tbi
rm ../wgs/cohort.filtered.raw.gt.mixed.vcf.gz.tbi
rm ../wgs/cohort.filtered.hf.gt.snp.vcf.gz.tbi
rm ../wgs/cohort.filtered.hf.gt.indel.vcf.gz.tbi
rm ../wgs/cohort.filtered.hf.gt.mixed.vcf.gz.tbi
rm ../wgs/cohort.filtered.hf.tmp.gt.vcf.gz.tbi
rm ../wgs/cohort.filtered.hf.gt.vcf.gz.tbi
rm ../wgs/cohort.filtered.hf.DP10.GQ20.allele.repmask.trf.gt.vcf.gz.tbi
rm ../wgs/cohort.filtered.hf.DP10.GQ20.allele.repmask.trf.sc.tmp.gt.vcf.gz.tbi
rm ../wgs/cohort.filtered.hf.DP10.GQ20.allele.repmask.trf.sc.gt.vcf.gz.tbi
