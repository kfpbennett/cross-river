# cross-river
Scripts necessary to run analysis for manakin cross-river gene flow manuscript in progress

First, run the raw_to_VCF code to combine reads from sequencing runs, align to the vitellinus reference, and call genotypes with Stacks. Next, run the filtering script to produce a filtered, thinned VCF file. The VCF is provided here, so the first steps can be skipped.

Next run the R scripts in any order to make a PCA, admixture plot, run an ANOVA, and calculate IBD stats.
