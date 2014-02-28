
Some (hopefully useful) code for GWAS. Predominantly designed to work with the OXSTATGEN suite.

Requirements: python, numpy, scipy


rrm.py 

Calculates the realised related matrix on some PLINK binary data.  See the methods section of this paper for details:

http://www.nature.com/ng/journal/v42/n7/full/ng.608.html

	Usage: python rrm.py [-h] [-snps snps] [-out out] [--ibs out] plinkfile


lmem.py

Performs Genome Wide Association testing using the method described in (for example):

http://www.nature.com/nmeth/journal/v8/n10/full/nmeth.1681.html

	usage: python lmem.py [-h] [-genfile genfile] [-header nrows]
               [-covariates covariates] [-weights weights]
               [-linebuffer linebuffer] [--uncorrected uncorrected]
               [-nprocess nprocess] [-jobnum jobnum]
               phenotype kinship output


haps2vcf.py 

Converts SHAPEIT2 haplotypes to a phased VCF file.

	 usage: python haps2vcf.py [-h] [-output output.vcf.gz] [--flip] input chromo

	 positional arguments:
	 	    input                 shapeit2 output
	   	    chromo                chromosome

	 optional arguments:
  	 	  -h, --help            show this help message and exit
  		  -output output.vcf.gz output file
  		  --flip                flip alleles
