/home/jianhuig/projects/def-leisun/jianhuig/#!/bin/bash
module load plink
cd /home/jianhuig/projects/def-leisun/jianhuig/
###############################################################
############# Obtain gentic and phenotype data#################
###############################################################
# dowmload ukbb file handlers
wget  -nd  biobank.ndph.ox.ac.uk/ukb/util/ukbmd5
wget  -nd  biobank.ndph.ox.ac.uk/ukb/util/ukbconv
wget  -nd  biobank.ndph.ox.ac.uk/ukb/util/ukbunpack
wget  -nd  biobank.ndph.ox.ac.uk/ukb/ukb/utilx/encoding.dat
wget  -nd  biobank.ndph.ox.ac.uk/ukb/util/gfetch
chmod 755 ukbmd5 ukbunpack ukbconv gfetch

# Validate and unpack main dataset
./ukbmd5 ukb47570.enc
./ukbunpack ukb47570.enc k64875r47570.key

# Subset Age, Sex, Ethinicty + Phenotype
./ukbconv ukb47570.enc_ukb r -iphenotype_list.txt -oPhenotype

# Obatin genetic data
for i in {1..22}
do
	./gfetch 22418 -c$i -ak64875r47570.key
	./gfetch 22418 -c$i -m -ak64875r47570.key
done

# Obtain SNP-index file
wget  -nd  biobank.ctsu.ox.ac.uk/ukb/ukb/auxdata/ukb_snp_bim.tar
tar -xvf ukb_snp_bim.tar

# Rename .bed and .fam file to match the name of .bim
for i in {1..22}
do
	mv ukb22418_c${i}_b0_v2.bed ukb_snp_chr${i}_v2.bed
	mv ukb22418_c${i}_b0_v2_s488175.fam ukb_snp_chr${i}_v2.fam
done


###############################################################
#################### Basic QC control##########################
###############################################################
for i in {1..22}
do
	plink --bfile ukb_snp_chr${i}_v2 \
	--mind 0.1 \ # Exclude inviduals with > 10% missing genotypes. 
	--maf 0.01 \ # Include SNPs with maf >=0.01
	--geno 0.1 \  # 90% genotype rate
	--hwe 0.00001 \ # Exclude markers fail hwe 0.00001
	--indep-pairwise 50 5 0.9 \ # LD pruning. Window size = 50 SNPs, pairwise LD > 0.9 is removed.
	--make-bed \ 
	--out chr${i}.cleaned 
done

# Merge all per chomosome genetic data to a single file
for i in {1..22}
do
echo chr$i.cleaned >> mergelist.txt
done
plink --merge-list mergelist.txt --make-bed --out allchromosome

# PCA needs to be scheduled to control memory
# sbatch pca.sh

# Remove related individuals based on kingship
# Needs to be paralleled to control memory
sbatch kingship.sh
# Combined paralled matrix into one
find . -type f -name 'split.king.bin.*' -exec cat {} + >> split.king.bin
# kingship cutoff
plink2 --bfile allchromosome --king-cutoff split 0.0625 --out final