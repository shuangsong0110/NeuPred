# NeuPred

**NeuPred** is an R implementation of a unified Bayesian framework for polygenic risk scores construction.

Song, S., Hou, L., & Liu, J. A flexible Bayesian regression approach for accurate polygenic risk prediction. *Submitted*, 2021.


## Table of contents
* [Getting Started](#getting-started)
* [Using NeuPred or NeuPred-I](#using-neupred--neupred-i)
* [Output](#output)
* [A toy example](#a-toy-example)

## Getting Started
NeuPred is an R package which can be installed using the command:
```r
devtools::install_github('shuangsong0110/NeuPred')
```

### Download the LD reference panel:
- Download the 1000 Genome Project reference panel:

`wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_plinkfiles.tgz`

`tar -xvzf 1000G_Phase3_plinkfiles.tgz`

- Users could also specify their own LD reference files with plink bfile format (.bim, .fam, .bed).
- Note that the LD reference plink bfile data should either be a single file (merge 22 chromosomes with plink), or by each chromosome with the names ended by the number of chromosomes, e.g., 1000G.EUR.QC.1.bim (,.fam, .bed), ..., 1000G.EUR.QC.22.bim (,.fam, .bed). 

### Other dependencies
- Install [plink](http://zzz.bwh.harvard.edu/plink/) (could be ignored if is already installed)
- We have attached ldetect files for European, African, and Asian ancestries in our software. (Berisa T, Pickrell J K. Approximately independent linkage disequilibrium blocks in human populations[J]. Bioinformatics, 2016, 32(2): 283.)

## Using NeuPred / NeuPred-I
### NeuPred (GWAS summary statistics based PRS approach)
```r
library(NeuPred)
NeuPred.run(summs=SUMMARY_STATISTICS (required),
	    LDpath=PATH_TO_LD_REFERENCE,
	    LDpath_chr=PATH_TO_LD_REFERENCE_FOR_EACH_CHROMOSOME (either LDpath or LDpath_chr is required),
	    n=GWAS_SAMPLE_SIZE (required),
	    external.ld=EXTERNAL_LD_REFERENCE (T / F , required),
	    ethnic=ETHNIC ('EUR'/'AFR'/'ASN', optional, default='EUR'),
	    plinkpath=PATH_TO_PLINK (required),
	    path=OUTPUT_DIR (required),
            [prior=PRIOR_TO_BE_USED ('auto' / 'L' / 'C' / 'HS', optional, default='auto'),
            testpath=PATH_TO_TEST_DATA (optional),
  	    MCMC=MCMC_ITERATIONS (optional),
   	    BURN=BURNING_TIMES (optional),
    	    parallel=T / F (optional, default=T),
     	    cores=NUMBER_OF_RUNNING_CORES (optional),
     	    chr=CHROM (optionall, default=1:22),
	    plot=ROC_PLOT (T / F, optional),
	    tmp=TEMP_FILES (T / F, optional)])                    
```
- SUMMARY_STATISTICS (required):
Prepare the summary statistics in the following format (including the header line):
```
   CHR        SNP      A1    A2      BETA        P
    1      rs4040617    G     A      0.013     8.42e-01
    1      rs4075116    C     T     -0.308     5.18e-01
    1      rs9442385    T     G      0.001     9.87e-01
    ...
```
Or:
```
   CHR        SNP      A1    A2         OR        P
    1      rs4040617    G     A      1.013     8.42e-01
    1      rs4075116    C     T      0.970     5.18e-01
    1      rs9442385    T     G      1.001     9.87e-01
    ...
```
where SNP is the rs ID, A1 is the effect allele, A2 is the alternative allele, BETA/OR is the effect/odds ratio of the A1 allele, P is the p-value of the effect. In fact, BETA/OR is used to determine the direction of an association. 

- PATH_TO_LD_REFERENCE (required): Either `LDpath` or `LDpath_chr` is required.  If the data are prepared as a single file, the `LDpath` should include the file name (without .bim/.fam/.bed) and the `LDpath_chr` is not required; if the data are split into chromosomes, the `LDpath_chr` should include the file name but not the exact number of chromosome (e.g., LDpath_chr='path/1000G.EUR.QC.').

- GWAS_SAMPLE_SIZE (required): Sample size of the GWAS.

- EXTERNAL_LD_REFERENCE (required): T / F, which indicated whether the LD reference is from an external LD panel (T) or from in-sample LD (F). 

- ETHNIC (optional): 'EUR' / 'AFR' / 'ASN', EUR: European ancestry; AFR: African ancestry; ASN: Asian ancestry.

- PATH_TO_PLINK (required): Full path to the plink software.

- OUTPUT_DIR (required): Full path for output files.

- PRIOR_TO_BE_USED (optional): chooose from 'auto' / 'L' / 'C' / 'HS'. If not specified, when `external.ld=T`, the prior would be set to 'L'; when `external.ld=F`, meaning that the LD is accurately estimated, the prior would be set to 'auto', and the algorithm will use a CV strategy to automatically select the best-performing prior for each chromosome using only training data (details are provided in the paper). 

- PATH_TO_TEST_DATA (optional): Full path to the test data (plink bfile format). If the test data is provided, the algorithm will calculate the overlapped SNPs between test data and training data, and then present the predictive r^2 and AUC (for binary traits), and give a ROC plot. If the test data is not provided, the algorithm will derive the posterior effect sizes for all SNPs in training dataset.

- MCMC_ITERATIONS (optional): The number of MCMC iteration times, default is 1e4.

- BURNING_TIMES (optional): The number of MCMC burning times, default is 2000.

- PARALLEL (optional): T / F. To decide whether to run the algorithm in parallel to speed up, default=T. Note that `parallel=T` will improve the computational effiency, but pay attention to the limits of memory.

- NUMBER_OF_CORES (optional): Number of cores for running MCMC when `parallel=T`, the optional is 5.

- CHROM (optional): The chromosome(s) on which the model is fitted, e.g., `chr=c(1, 3, 5)`. The default is `chr=1:22`.

- PLOT (optional): T / F. If `plot=T` and the outcome is binary, the software will provide an ROC plot.

- TEMP_FILES (optional): T / F. If `tmp=T`, the temp files will be kept, including the LD blocks, test files, etc. If `tmp=F`, the temp files will be deleted after the results have been generated.


### NeuPred-I (Individual-level data based PRS method)
The individual-level data version, NeuPred-I can be run using the following command:
```r
library(NeuPred)
NeuPred.I.run(trainpath=PATH_TO_TRAINING_DATA (required),
              plinkpath=PATH_TO_PLINK (required),
	      path=OUTPUT_DIR (required),
	      [prior=PRIOR_TO_BE_USED ('auto' / 'L' / 'C' / 'HS', optional),
              testpath=PATH_TO_TEST_DATA (optional),
	      MCMC=MCMC_ITERATIONS (optional),
   	      BURN=BURNING_TIMES (optional),
    	      parallel=T / F (optional),
     	      cores=NUMBER_OF_RUNNING_CORES (optional),
     	      chr=CHROM (optional),
	      plot=ROC_PLOT (T / F, optional),
	      tmp=TEMP_FILES (T / F, optional) ])       
             
```


- PATH_TO_TRAINING_DATA (required): Individual-level plink bfile as the training dataset.

- PATH_TO_PLINK (required): Full path to the plink software.
- OUTPUT_DIR (required): Full path for output files.

- PRIOR_TO_BE_USED (optional): chooose from 'auto' / 'L' / 'C' / 'HS'. If not specified, when `external.ld=T`, the prior would be set to 'L'; when `external.ld=F`, meaning that the LD is accurately estimated, the prior would be set to 'auto', and the algorithm will use a CV strategy to automatically select the best-performing prior for each chromosome using only training data (details are provided in the paper). 

- PATH_TO_TEST_DATA (optional): Full path to the test data (plink bfile format). If the test data is provided, the algorithm will calculate the overlapped SNPs between test data and training data, and then present the predictive r^2 and AUC (for binary traits), and give a ROC plot. If the test data is not provided, the algorithm will derive the posterior effect sizes for all SNPs in training dataset.

- MCMC_ITERATIONS (optional): The number of MCMC iteration times, default is 1e4.

- BURNING_TIMES (optional): The number of MCMC burning times, default is 2000.

- PARALLEL (optional): T / F. To decide whether to run the algorithm in parallel to speed up, default=T.

- NUMBER_OF_CORES (optional): Number of cores for running MCMC when `parallel=T`, the optional is 5.

- CHROM (optional): The chromosome(s) on which the model is fitted, e.g., `chr=c(1, 3, 5)`. The default is `chr=1:22`.

- PLOT (optional): T / F. If `plot=T` and the outcome is binary, the software will provide an ROC plot.

- TEMP_FILES (optional): T / F. If `tmp=T`, the temp files will be kept, including the LD blocks, test files, etc. If `tmp=F`, the temp files will be deleted after the results have been generated.


## Output

An R list containing:

$res: predictive r2 and the AUC (for qualitative traits)

$S: derived polygenic risk scores

$select.prior: the selected prior

The effect sizes for each variant are saved as files in folder `./posterior`


## A toy example
Clone this repository using the following git command:

$ git clone https://github.com/shuangsong0110/NeuPred.git

$ cd ./NeuPred

Run with R:
```r
library(NeuPred)
plinkpath <- 'PATH_TO_PLINK_SOFTWARE'
path <- './example/'
dir.create(path, recursive=T)
LDpath <- './test_dat/LD'
summs <- fread('./test_dat/train.txt')
testpath <- './test_dat/test'
n <- 3000
res <- NeuPred.run(summs=summs, LDpath=LDpath, n=n, prior='L', external.ld = F, ethnic='EUR',
                    plinkpath=plinkpath, path=path, testpath=testpath,
                    chr=22, parallel=T, cores=5, MCMC=1e4, BURN=2e3, plot=T)
```


## Maintainer

Please contact Shuang Song (song-s19@mails.tsinghua.edu.cn) if there are any problems or questions.



