#!/usr/bin/python -O

import sys,time,os,glob,pickle,cPickle, csv,gzip,numpy as np,resource,argparse,string
import numpy as np
import iolib

parser = argparse.ArgumentParser(description='calculates the RRM on plink binary file set')
parser.add_argument('plinkfile', metavar='plinkfile', type=str, help='a binary plink file set')
parser.add_argument('-snps', metavar='snps', type=str,default='', help='list of snp ids to use for rrm calculation')
parser.add_argument('-out', metavar='out', type=str,default='', help='outfile')
parser.add_argument('--ibs', metavar='out', type=bool,default=False, help='use kinship coefficient rather than realised relationship')

args = parser.parse_args()

if args.out=='': outfile = csv.writer(iolib.wopen(args.plinkfile+".rrm.gz"),delimiter="\t")
else: outfile = csv.writer(iolib.wopen(args.out+".rrm.gz"),delimiter="\t")

if args.snps!='': 
    snps = iolib.scan(args.snps)
    print "Calculating RRM from subset of",len(snps),"SNPs"
else: snps = None


infile = iolib.plinkReader(args.plinkfile,snps=snps)

n = infile.nsample

rrm = np.zeros((n,n),np.float)
rrm_diag = np.zeros(n,np.float)

print "Calculating RRM..."
start_time = time.time()

nmono=0
for chrom,pos,rsid,ref,alt,genotypes in infile:
    afreq = float(2*(genotypes==2).sum() + (genotypes==1).sum())/(2*float(n))
    if afreq>0.0 and afreq<1.0:
        rrm += np.outer((genotypes-2*afreq),(genotypes - 2*afreq)/(2*afreq*(1-afreq)))
        rrm_diag += (genotypes**2 - (1+2*afreq)*genotypes + 2*afreq**2)/(2*afreq*(1-afreq))
    else:
        nmono += 1

rrm = rrm/(infile.nsnp-nmono)
rrm_diag = 1 + rrm_diag/(infile.nsnp-nmono)

rrm[np.arange(n),np.arange(n)] = rrm_diag

print "RRM calculation took ",time.time() - start_time," seconds\nWriting ",args.plinkfile+".rrm.gz"

for row in rrm:    outfile.writerow(row)
