import sys,time, csv,numpy as np,argparse,string,os,math,logging
import numpy as np
import iolib
import lmem_func
import multiprocessing as mp

parser = argparse.ArgumentParser(description='Linear mixed effect models for quantitative traits')
parser.add_argument('-genfile', metavar='genfile',default="", type=str, help='a .gen file')
#parser.add_argument('-dosage', metavar='dosage',default="", type=str, help='a .gen file')
parser.add_argument('-header', metavar='nrows',default=0, type=int, help='# rows in header')

parser.add_argument('phenotype', metavar='phenotype', type=str, help='a phenotype file')
parser.add_argument('kinship', metavar='kinship', type=str, help='n x n kinship matrix')
parser.add_argument('output', metavar='output', type=str, help='output')

parser.add_argument('-covariates', metavar='covariates',default=None, type=str, help='covariates (factors not handled, make your own damn design matrix)')
parser.add_argument('-weights', metavar='weights',default=None, type=str, help='weights')

parser.add_argument('-linebuffer', metavar='linebuffer',default=5000, type=int, help='line buffer (default=5000)')
# parser.add_argument('-start', metavar='start',default=None, type=int, help='')
# parser.add_argument('-stop', metavar='stop',default=None, type=int, help='')

parser.add_argument('--uncorrected', metavar='uncorrected',default=False, type=bool, help='apply linear regression without correction for relatedness')
parser.add_argument('-nprocess',metavar='nprocess',type=int,default=1)
parser.add_argument('-jobnum',metavar='jobnum',type=int,default=1)

args = parser.parse_args()

def lmem(args,start,stop,output):
    print "Processing SNPs",start,"-",stop,"..."
    print "Reading phenotypes and RRM..."

    Y = np.genfromtxt(iolib.ropen(args.phenotype),missing_values="NA",dtype=np.float,skip_header=args.header)
    if len(Y.shape)==1:
        tmp = np.empty((len(Y),1),np.float)
        tmp[:,0] = Y
        Y = tmp
    nsample,npheno = Y.shape

    genfile = iolib.genReader(args.genfile,nsample,args.linebuffer,start,stop)
    outf = csv.writer(iolib.wopen(output),delimiter="\t")
    K = np.genfromtxt(iolib.ropen(args.kinship),dtype=np.float)
    print "K.shape =",K.shape
    print "Y.shape =",Y.shape
    print nsample,"samples and",npheno,"phenotypes"
    missing = []
    null_model = []

    print "Calculating rotation matrices..."
    for mpheno in range(npheno):
        missing.append(np.logical_not(np.isnan(Y[:,mpheno])))
        ii = missing[mpheno]

        null_model.append(lmem_func.get_delta(K[ii][:,ii],Y[ii,mpheno]))

    #    Ut,denom,beta_null,sigma_g_null,loglik_null,delta_null = get_delta(K[ii][:,ii],Y[ii,mpheno]))
    print "done"

    print "Fitting LMMs..." 
    time0 = time.time()
    processed = 0
    for rsid1,rsid2,pos,G in genfile:
        output = []
        nsnp = G.shape[1]

        for mpheno in range(npheno):    
            ii = missing[mpheno]
            output.append(lmem_func.fitlmm(null_model[mpheno]['offset'],null_model[mpheno]['Ut'],null_model[mpheno]['denom'],G[ii],null_model[mpheno]['uty']))

        for i in range(nsnp):
            outf.writerow([rsid1[i],pos[i],rsid2[i]]
                          +sum([[null_model[mpheno]['beta_null'],null_model[mpheno]['sigma_g_null'],null_model[mpheno]['loglik_null']] 
                                +output[mpheno][0][i].tolist()+[output[mpheno][1][i],output[mpheno][2][i],output[mpheno][3][i]] for mpheno in range(npheno)],[]))

        processed+=nsnp
        print processed,"loci processed"
    print "Took",time.time() - time0,"seconds"


if __name__ == '__main__':
    # logger = mp.log_to_stderr()
    # logger.setLevel(logging.INFO)
    assert args.header>=0
    iolib.checkfile(args.kinship)
    iolib.checkfile(args.genfile)

    if args.output[-3:]==".gz": args.output = args.output[:-3]

    outf = csv.writer(iolib.wopen(args.output+".gz"),delimiter="\t")

    nl = iolib.nlines(iolib.ropen(args.genfile))
    print nl,"SNPs"
    neach = int(math.ceil(nl/args.nprocess))
    chunks = range(0,nl,neach) + [nl]
    print chunks

    pool = mp.Pool(processes=args.nprocess)

    Y = np.genfromtxt(iolib.ropen(args.phenotype),missing_values="NA",dtype=np.float,skip_header=args.header)

    if len(Y.shape)==1:
        nsample = Y.shape[0]
        npheno = 1
    else:
        nsample,npheno = Y.shape

    outf.writerow(["chrom","pos","rsid"]+sum([["pheno"+str(mpheno+1)+"_"+val for val in ["beta_null","sigma_null","loglik_null","beta_0","beta_1","beta_1_se","sigma_g","loglik"]] for mpheno in range(npheno)],[]))
    del outf

    for i in range(args.nprocess):
        pool.apply_async(lmem,(args,chunks[i],chunks[i+1],args.output+str(i)+".gz"))
        print "Submitted",chunks[i],"-",chunks[i+1],"writing to",args.output+str(i)
                         
    pool.close()
    pool.join()

    for i in range(args.nprocess):
        os.system("cat " + args.output+str(i)+".gz"  + " >> " + args.output+".gz")
        os.system("rm " + args.output+str(i)+".gz")

