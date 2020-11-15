# MGP: Mixed model based Genome wide association studies and genome Prediction

## Description

There are already a lot of softwares that can be used to perform genome wide association studies (GWAS) and genome prediction (GP) but some of them are not computationally efficient or difficulty to use. The objective of the project is to develop a versatile, efficient and easy-to-use tool for GWAS and GP. Genotype files in VCF format, commonly used in resequencing or TPED format are read in fragment. The efficient mixed-model implemented in GEMMA was improved. The algorithm, average information restricted maximum likelihood (AI-REML) used to estimate variance component was implemented in a more efficient way. Univariate and repeated-measured GWAS and genome prediction of additive and dominant effects can be performed easily.

**NOTE:** It's recommended to use version >= v0.1.1. They are more precise for analysis involving low rank matrix (detail in ChangeLog).

## Version

The program was written by C++ and compiled using gcc 6.4.0 on a 64-bit Linux version. MGP_float_v0.1.3 (based on Eigen) and MGP_float_MKL_v0.1.3 (compiled based on Intel Math Kernel Library by defining EIGEN_USE_MKL_ALL, which is faster than Eigen for large matrix multiplication) are float versions and genotypes are stored as floats. MGP_char_v0.1.0 (with some bugs, will be updated in spare time) is the character version and genotypes are stored as characters, which can save storage space and improve the speed of reading genotype files in GWAS. The float version can use the genotype and dosage information. But the character version can’t use the dosage information which is a float. Although the current version may be not computationally efficient for a trait that was measured repeatedly it is still fast for small sample size and the more efficient algorithm implemented in the new version is being tested.

Starting from version 0.1.2 association analysis can be parallelized by two dimensions, by chromosome (chromosome segment) and by SNPs in one chromosome (chromosome segment). The number of threads used can be controlled by setting the value of the environment variable OMP_NUM_THREADS (export OMP_NUM_THREADS=N). For example, exploiting 2 cores, you can set OMP_NUM_THREADS=2 and then run the program MGP.

## Install

\# download MGP_float_v0.1.3    
chmod +x MGP_float_v0.1.3

\# download MGP_float_MKL_v0.1.3.zip    
unzip MGP_float_MKL_v0.1.3.zip    
chmod +x MGP_float_MKL_v0.1.3

## Command and option

**kin**　　　　 MGP kin [options]  

Either complete or chunked genotype files can be used to build additive and dominant relationship. If chunked genotype files are used ‘merge’ command is needed to merge the binary intermediate relationship files. If there are missing genotypes they will be sampled from marginal distribution, randomly. If missing genotypes were imputed by more accurate methods e.g. beagle dosage information can be used. SNPs with minimal allele frequency close to zero will be removed.  

    --in         FILE  input file (genotype file).  
    --part             genotype is partial.  
    --outBin     FILE  output intermediate binary kinship (--part)  
    --oKinBin    FILE  output binary kinship file (no --part)  
    --oKinTxt    FILE  output kinship file in text (no --part)  
    --oGtBin     FILE  output binary genotype file (necessary for GWAS but not for GS)  
    --oPos       FILE  output position file with two columns, chr and pos.  
    --oFreq      FILE  output allele frequency product (p*q), in current version this parameter is unuseful.   
    --gType      CHAR  genotype or dosage  
    --kinType    CHAR  model to calculate genome relationship,add(additive),dom(dominance), by default: add  
    --addMethod  INT   method to calculate additive genome relationship (1 or 2). Method 1- formula according to Vanraden, method 2- similar to GCTA  
    --tped             plink format tped file  
    --rSnp       INT   number of snps kept in memory,by default: 1000 

**pheno**　　　MGP pheno [options]  

For phenotype missing values must be represented by ‘NA’ and the first column must be id or ID. This command will sort the phenotype and keep the orders of IDs in phenotype consistent with that in genotype. It’s not required that number of individuals in phenotype is same to that in phenotype. The program will extract the matched individuals.  

    --vcf          FILE  individuals with genotype(vcf file which can just contain the header lines).  　
    --rawPheno     FILE  input phenotype file.  
    --id           FILE  input id pair file with two columns, id in vcf and id in phenotype.  
    --formatPheno  FILE  output phenotype file.  
    --kinBin       FILE  input binary kinship.  
    --trimKin      FILE  output binary kinship after removing closely related individuals.  　
    --rValue       FLOAT threshold value used to remove related individuals,range (0,1],by default: 0.9.  　
    --inGbin       FILE  input binary genotype.  
    --outGbin      FILE  output binary genotype after remove related individuals.  
    --PCA                if this parameter was set principal components would be added to final phenotype file  
    --nPCA         INT   number of principal components that would be added to final phenotype file,by default: 5  


**merge**　　　　MGP merge [options]  

    --igBin      FILE   input partitioned binary genotype.  
    --ikBin      FILE   input partitioned binary kinship.  
    --inHeader   FILE   vcf file or vcf header file  
    --ogBin      FILE   output merged binary genotype.  
    --okBin      FILE   output merged binary kinship.  
    --okTxt      FILE   output merged kinship in text.


**mix**　　　　　MGP mix [options]  

Perform GWAS and genome prediction/selection (GS). For univariate GWAS three algorithm can be used: EPAI, PAI and EMMA(X). EPAI is efficient mixed model, variance parameters were re-estimated for each SNP. PAI is a variance component method and variance parameters were estimated only one time by computationally efficient AI-REML. EMMA(X) is also a variance component method and variance parameters were estimated by EMMA algorithm and Newton down-hill method.  

    --phenotype  FILE  phenotype file  
    --genotype   FILE  genotype file(binary),perform GWAS or weighted GS  　
    --iFreq      FILE  allele frequency file(binary),perform weighted GWAS or weighted GS, not useful for current version  
    --kinAdd     FILE  additive relationship file(binary)  
    --kinDom     FILE  dominant relationship file(binary)  　
    --out        FILE  result file(p-value or breeding value)  
    --GWAS             genome wide association study  
    --GS               genome selection  
    --varMethod  STR   algorithm to estimate variance parameters(EPAI|PAI|EMMA)  
    --wError     FILE  weight of residual variance  
    --pCol       INT   column of dependent variable  
    --covD       INT,  columns of covariate(discrete),if include multivariate they are delimited by comma, e.g. INT1,INT2,INT3  　
    --covC       INT,  columns of covariate(consecutive),if include multivariate they are delimited by comma  
    --randV      INT,  columns of random effects other than additive or dominant,if include multivariate they are delimited by comma  
    --knownPar   FILE  variance parameters file if they are known  
    --initPar    FILE  initial values of variance parameters for iterations  
    --varIter    INT   the max number of iterations to estimate variance,by default 20  
    --tol        FLOAT convergence criteria,by default 1e-4  
    --rSnp       INT   number of SNPs kept in memory at a time(GWAS),by default 1000  
    --logML            output log maximum likelihood of the model when PAI algorithm is used
    --logREML          output log restricted maximum likelihood of the model when PAI algorithm is used
    --NoSe             not output standard error of fix and random effects



## Formats  

### Input format 

#### Genotype

VCF format or TPED format in Plink.  

#### Phenotype

Include all covariates (discrete or consecutive) and dependent variates. NA is allowed.  

### Output

#### GWAS output

Three columns: effects of SNPs, standard error, p value. Orders of SNPs are same to the position file  

#### GP output

\#Class  Effect  EffectSe        Observe Predict Residual  

#### Stdout

The program will output information of running which include the information of model and  variance component. For example:  

Yield = Mean + Rand_add + Rand_Row + Rand_Col + Rand_error  

\*********** Variance **********  
Var: 803.311 849.447 197.807 2995.05  
SE: 307.661 397.552 126.464 14.1251  
\********************************   

The order of ‘Var’ is consistent with the random effects in model.

## Example

Test data and script to run the program can be found in folder v0.1.0. New parameters of MGP were not included in the script but it is easy to appended.

## Dependency

Eigen  
boost  
zlib  
MKL(Intel Math Kernel Library, for program MGP_float_MKL_v0.1.3)

## Contact

If you find bugs please email me.  
Email:  miaozepu@genomics.cn (stopped)  
        &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;runzejiaoku@sina.com  
Institute: Applied agriculture of BGI  


## Licence

Free for research. For commercial use please contact P_agro_pmo@genomics.cn  

