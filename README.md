# cross-river
Scripts necessary to run analysis for manakin cross-river gene flow manuscript in progress

First, run the raw_to_VCF code to combine reads from sequencing runs, align to the vitellinus reference, and call genotypes with Stacks. Next, run the filtering script to produce a filtered, thinned VCF file. The VCF is provided here, so the first steps can be skipped.

Next run the R scripts in any order to make a PCA, admixture plot, run an ANOVA, and calculate IBD stats.

Other files included here: 

filtered4_thin5k.recode.vcf is the filtered VCF file used in most analyses.
plink.raw is the same data in raw plink format.
f4t5k_males.recode.vcf is the same data omitting females.
sampleinfo_filtered2.txt is additional sample information from Table S1.
pops_gps includes location info for IBD analyses.
