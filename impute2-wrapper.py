 #Jared O'Connell's impute2 wrapper script.
#jared@well.ox.ac.uk
import argparse,sys,time,multiprocessing as mp,subprocess,os,string,io,gzip
import math

def ropen(filename,buf_in_MB=100):
    if not os.path.isfile(filename):
        sys.exit(filename+" does not exist!")

    else:
        if buf_in_MB != None:
            if filename[-3:]==".gz":    
                return  io.BufferedReader(gzip.GzipFile(fileobj=open(filename,"rb")),buf_in_MB*(2**20))
            else: return open(filename,"rb",buf_in_MB*(2**20))
        else:
           if filename[-3:]==".gz": return gzip.GzipFile(fileobj=cStringIO.StringIO(open(filename).read()))
           else: return cStringIO.StringIO(open(filename).read())


#WRAPPER FUNCTION FOR IMPUTE2
def impute(known_haps_g,h,l,m,Ne,k_hap,buf,start,stop,outfile):
    args = string.join(["impute2","-known_haps_g "+known_haps_g,
                     "-h " + h,
                     "-l " + l,
                     "-m " + m,
                     "-Ne "+str(Ne),
                     "-k_hap " +str(k_hap),
                     "-o " + outfile,
                     "-o_gz",
                     "-buffer "+str(buf),
                     "-int " + str(start)+" "+str(stop)])
    start_time =  time.time()
    subprocess.call(args,shell=True,stdout=open("/dev/null","wb"))
    print "Chunk",start,"-",stop,"finished in",time.time() - start_time,"seconds"


def add_chromosome_to_pool(chromosome,args,pool):
    o = args.o
    if o == "": o = args.known_haps_g[:-5]
    if os.path.isfile(o+"-imputed.gen.gz"):
        print o+"-imputed.gen.gz exists!\nExiting..."
        quit()

    if args.chromosome==None:
        o = string.join(o.split("chr1"),"chr"+chromosome)
        known_haps_g =  string.join(args.known_haps_g.split("chr1"),"chr"+chromosome)
        refhap =  string.join(args.refhap.split("chr1"),"chr"+chromosome)
        refleg =  string.join(args.refleg.split("chr1"),"chr"+chromosome)
        geneticmap = string.join(args.geneticmap.split("chr1"),"chr"+chromosome)
    else:
        known_haps_g =  args.known_haps_g
        refhap =  args.refhap
        refleg =  args.refleg
        geneticmap = args.geneticmap
        
        
    print "Outputing",o

    #GETS THE FIRST AND LAST POSITION OF INPUT HAPS
    if not os.path.isfile(known_haps_g):
        print known_haps_g,"does not exist!Exiting..."
        quit()

    infile = ropen(refleg)
    infile.next()
    start = int(infile.next().split()[1])
    for row in infile: pass
    stop = int(row.split()[1])
    npos = stop-start

    #CALCULATE CHUNKS
    if args.chunksize==0: neach = int(math.ceil(float(npos)/float(args.nprocess)))
    else: neach = int(math.ceil(float(npos)/(float(npos)/float(args.chunksize))))
    neach = min(5000000,max(1000000,neach)) #makes sure window size is between 1Mb and 5Mb
    print neach
    chunks = range(start,stop,neach)+[stop]
    print "chunks = ",chunks
    nchunk = len(chunks)-1
    print "Imputing position",chromosome,":",start,"-",stop,"in",nchunk,"chunks of roughly",neach
    # print "Chunks =",chunks

    workdir = o+".working"
    try: os.mkdir(workdir)
    except: pass

    outfiles = [workdir+"/"+o+"_"+str(chunks[idx])+"_"+str(chunks[idx+1]) for idx in range(nchunk)]
    print outfiles
    for i in range(len(chunks)-1):
        pool.apply_async(impute,(known_haps_g,refhap,refleg,geneticmap,args.Ne,args.k_hap,args.buffer,chunks[i],chunks[i+1],outfiles[i]))
    
    return outfiles

def concatenate_chromosome(chromosome,args):
    o = string.join(o.split("chr1"),"chr"+chromosome)
    known_haps_g =  string.join(args.known_haps_g.split("chr1"),"chr"+chromosome)
    refhap =  string.join(args.refhap.split("chr1"),"chr"+chromosome)
    refleg =  string.join(args.refleg.split("chr1"),"chr"+chromosome)
    geneticmap = string.join(args.geneticmap.split("chr1"),"chr"+chromosome)
    

buf = 2**20 * 100 # I/O buffer size used for all files

#REQUIRED
parser = argparse.ArgumentParser(description='Runs impute2 on chunks and then merges the output (assumes impute2 is in your path)')
parser.add_argument('-known_haps_g',metavar='known_haps_g',type=str,required=True)
parser.add_argument('-refhap',metavar='h',type=str,required=True)
parser.add_argument('-refleg',metavar='l',type=str,required=True)
parser.add_argument('-geneticmap',metavar='m',type=str,required=True)

#OPTIONAL
parser.add_argument('-o',metavar='o',type=str,default="")
parser.add_argument('-Ne',metavar='Ne',type=int,default=20000)
parser.add_argument('-k_hap',metavar='k_hap',type=int,default=1000)
parser.add_argument('-buffer',metavar='buffer',type=int,default=1000)
parser.add_argument('-chunksize',metavar='chunksize',type=int,default=5000000)
parser.add_argument('-nprocess',metavar='nprocess',type=int,default=1)
parser.add_argument('-chromosome',metavar='chromosome',type=int,default=None)


args = parser.parse_args()

if __name__ == '__main__':

    for f in [args.known_haps_g,args.refhap,args.refleg,args.geneticmap,args.o]:
        if args.chromosome==None:
            if 'chr1' not in f:
                sys.exit("Expecting"+f+" for chromosome 1")
            for c in range(2,23):
                if 'chr'+str(c) in f:
                    sys.exit("Expecting"+f+" for chromosome 1")


    prefix = args.known_haps_g.split("chr")[0]


    #SPAWNS IMPUTE PROCESSES FOR EACH CHUNK
    start_time =  time.time()
    pool = mp.Pool(processes=args.nprocess)
    outfiles = {}

    if args.chromosome==None: chromosomes = range(1,23)
    else: chromosomes = [args.chromosome]
    
    print "Imputing chromosomes",chromosomes
    print "prefix =",prefix

    for chromosome in chromosomes:
        outfiles[chromosome] = add_chromosome_to_pool(str(chromosome),args,pool)

    pool.close()
    pool.join()
    print "Imputing took ",time.time() - start_time," seconds" 


    #CONCATENATES THE OUTPUT
    for i in chromosomes:
        print "Concatenating output for chromsome",i
        chromosome = str(i)
        outfile = open(string.join([prefix,"chr",chromosome],"")+"-imputed.gen.gz","wb",buf)
        infofile =   open(string.join([prefix,"chr",chromosome],"")+"_info","wb",buf)
        infofile.write("snp_id rs_id position exp_freq_a1 info certainty type info_type0 concord_type0 r2_type0\n")

        for o in outfiles[i]:
            if os.path.isfile(o+".gz"):
                tempfile = open(o+".gz","rb",buf)
                while True:
                    data = tempfile.read(buf)
                    if data:   outfile.write(data)
                    else:  break
                tempfile = open(o+"_info","rb",buf)
                tempfile.readline()
                infofile.write(tempfile.read())
                os.remove(o+".gz")
                os.remove(o+"_info")
            else:
                print "WARNING:",o,"does not exist (maybe there were no loci in this region)"
                pass

    print "Took ",time.time() - start_time," seconds"
    
    print "Done."
