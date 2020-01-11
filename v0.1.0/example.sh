
###########################################
# Data
###########################################

# 1) CPgeno.pair, CPgeno.tag, CPgeno.tped, CPgeno.unorder.pair, CPgeno.vcf, CPpheno.Fix.sort, CPpheno.sort, CPpheno.unorder
     # this data set is from the software 'sommer' and the genotype was transformed to vcf, tped format. For simplicity non integer
	 # value was represent by '1/1' or 'G G'.
	
	 # CPpheno.unorder is disordered phenotype
	 
	 # for CPpheno.Fix.sort four covariates were simulated, two discrete and two consecutive and the effects were added to dependent variate 'yield'.

# 2) mouse_hs1940.gp.vcf.gz, mouse_hs1940.pheno.sort, mouse_hs1940_split.tar.gz 
	 # this data set is from the software 'GEMMA' and the genotype was transformed to vcf format.
	 
	 # for mouse_hs1940_split.tar.gz genotype was split by chromosome
	 
# 3) RepeatAbel.geno.tag, RepeatAbel.geno.tped, RepeatAbel.pair, RepeatAbel.pheno.sort
	 # this data set is from the software 'RepeatAbel' and the genotype was transformed to tped format.


###########################################
# GWAS and GS
###########################################

# 1) calculate genome relationship
	 # For genotype NA is allowed but not suggest. It's better to perform genotype imputation before building relationship
	 # in VCF format NA is './.'; in TPED format NA is '0 0'
			###########################################
			# additive
			###########################################
			
			# in tped format
			MGP kin --in geno.tped --oKinBin add.kin.bin --oKinTxt add.kin.txt --oGtBin add.genotype.bin --oPos add.pos --gType genotype --kinType add --addMethod 1 --tped
			
			# in vcf format, using genotype
			MGP kin --in vcf.out.gz --oKinBin add.kin.bin --oKinTxt add.kin.txt --oGtBin add.genotype.bin --oPos add.pos --gType genotype --kinType add --addMethod 1
			
			# in vcf format, using dosage 

			MGP kin --in vcf.out.gz --oKinBin add.kin.bin --oKinTxt add.kin.txt --oGtBin add.genotype.bin --oPos add.pos --gType dosage --kinType add --addMethod 1
			
			#*********************
			#genotype file is partial 
			#*********************
			
			# genotype
			# <1> intermediate kinship
			MGP kin --in chr1.vcf.out.gz --outBin chr1.add.kin.bin --oGtBin chr1.add.genotype.bin --oPos chr1.add.pos --gType genotype --kinType add --addMethod 1 --part
			MGP kin --in chr2.vcf.out.gz --outBin chr2.add.kin.bin --oGtBin chr2.add.genotype.bin --oPos chr2.add.pos --gType genotype --kinType add --addMethod 1 --part
			
			# <2.1> merge intermediate kinship
			MGP merge --ikBin chr1.add.kin.bin --ikBin chr2.add.kin.bin --okBin all.add.kin.bin -okTxt all.add.kin.txt
			
			# <2.2> merge binary genotype
			MGP merge --inHeader tag.vcf --igBin chr1.add.genotype.bin --igBin chr2.add.genotype.bin --ogBin all.add.kin.bin
			
			# dosage
			#  similar to genotype but parameter '--gType' is set to  dosage
			
			
			############################################
			# dominance
			############################################
			
			# similar to description above , parameter '--kinType' is set to dom
			
# 2) sort or trim the phenotype
			
			# sort phenotype
			
			MGP pheno --rawPheno CPpheno.unorder --id CPgeno.unorder.pair --vcf CPgeno.tag --formatPheno CPpheno.sort.out
			
			# remove closely related individuals and sort
			MGP pheno --vcf CPgeno.tag --rawPheno CPpheno.unorder --id CPgeno.unorder.pair --formatPheno CPpheno.sort.out --kinBin add.kin.bin --trimKin trim.add.kin.bin --inGbin add.genotype.bin --outGbin trim.genotype.bin

			# add PCA to phenotype file
			
			MGP pheno --rawPheno CPpheno.unorder --id CPgeno.unorder.pair --vcf CPgeno.tag --formatPheno CPpheno.sort.out --kinBin add.kin.bin --PCA

# 3) Genome wide association study (GWAS)

			# ordinary linear model
			MGP mix --GWAS --phenotype phenotype.sort --pCol 3 -covD 2 --genotype add.genotype.bin --out pvalue.out
			
			
			# efficient mix model (variance components are re-estimated for each SNP) 
			MGP mix --GWAS --phenotype phenotype.sort --pCol 3 -covD 2 --genotype add.genotype.bin --out pvalue.out --kinAdd add.kin.bin --varMethod EPAI
			
			
			# variance component method (variance components are estimated one time)
			MGP mix --GWAS --phenotype phenotype.sort --pCol 3 -covD 2 --genotype add.genotype.bin --out pvalue.out --kinAdd add.kin.bin --varMethod PAI
			
			
			# repeat measures
			MGP mix --GWAS --phenotype RepeatAbel.pheno.sort --genotype add.genotype.bin --kinAdd add.kin.bin --out pvalue.out --pCol 2 --covD 3 --covC 4 --randV 1 --varMethod PAI
			
			# *********************
			# multiple randome effects are supported but not tested sufficiently
			
			
# 4) Genome prediction or selection (GS)

			# 1 random effect
			MGP mix --GS --phenotype CPpheno.sort --kinAdd add.kin.bin --out ebv.out --pCol 6 --varMethod PAI

			# 2 random effects
			MGP mix --GS --phenotype CPpheno.sort --kinAdd add.kin.bin --randV 2 --out $out --pCol 6 --varMethod PAI
			
			
			# 3 random effects
			MGP mix --GS --phenotype CPpheno.sort --kinAdd add.kin.bin --randV 2,3 --out $out --pCol 6 --varMethod PAI
			
			# repeatability model
			MGP mix --GS --phenotype RepeatAbel.pheno.sort --genotype add.genotype.bin --kinAdd add.kin.bin --out ebv.out --pCol 2 --covD 3 --covC 4 --randV 1 --varMethod PAI
			
			# random effect(dominance + additive + others)
			MGP mix --GS --phenotype CPpheno.Fix.sort --kinAdd add.kin.bin --kinDom dom.kin.bin --out ebv.out --covD 9,10 --covC 11,12 --randV 3 --pCol 6 --varMethod PAI


			
			
			
			
			
			
			


