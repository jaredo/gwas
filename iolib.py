#a whole bunch of i/o functions
#jared@well.ox.ac.uk

import csv, sys,copy,gzip,numpy as np,os,math,time,cStringIO,io,itertools
import numpy as np
import subprocess,string
import string,random#,pysam
csv.field_size_limit(1000000000)

def read_genetic_map(fname):
    return np.genfromtxt(ropen(fname),dtype=np.dtype([('pos',np.int),('combined',np.float),('cm',np.float)]),skiprows=1)            


def getRho(fname,pos,rate=1.):
    gm = read_genetic_map(fname)
    ii = np.unique(gm['cm'],True)[1]
    cm = np.interp(pos,gm['pos'][ii],gm['cm'][ii])
    r = (cm[1:]-cm[:-1])/100
    rho = 1 - np.exp(-rate*r)
    return rho

def read_germline(fname):

    ibd_dt = np.dtype([('fid1', np.str_, 200),('id1', np.str_, 200),('fid2', np.str_, 200),('id2', np.str_, 200),('chromosome', np.str_,2),
                       ('pos1', np.uint),('pos2', np.uint),('rsid1', np.str_,100),('rsid2', np.str_,100),('bits',np.uint),('length',float),
                       ('unit',np.str_,2), ('mismatch',np.uint),('hom1',np.uint),('hom2',np.uint)])

    return np.genfromtxt(ropen(fname),ibd_dt,skip_header=1)

def read_shapeit_haps(fname,transpose=False):

    if fname[-8:] == ".haps.gz": fname = fname[:-8] 
    elif fname[-5:] == ".haps": fname = fname[:-5] 

    sample_file = ropen(fname+".sample")
    sample_file.next()
    sample_file.next()
    tmp = [val.split() for val in sample_file]
    sample = np.empty(len(tmp),[('id1',np.str_,300),('id2',np.str_,300)])
    sample['id1'] = [val[0] for val in tmp]
    sample['id2'] = [val[1] for val in tmp]
    
    if os.path.isfile(fname+".haps.gz"):
        infile = ropen(fname+".haps.gz")
    else:
        infile = ropen(fname+".haps")

    nsnp=nlines(infile)
    nhap=len(sample)*2

    if transpose:
        H = np.empty((nhap,nsnp),np.uint8)
    else:
        H = np.empty((nsnp,nhap),np.uint8)

    
    pos = []
    ref = []
    alt = []
    for i,row in enumerate(infile):
        data = row.split(" ",5)
        ref.append(data[3])
        alt.append(data[4])
        pos.append(int(data[2]))
        if transpose: H[:,i] = np.fromstring(data[5],dtype=np.uint8,sep=" ")
        else: H[i] = np.fromstring(data[5],dtype=np.uint8,sep=" ")
        
    return sample,pos,np.array(ref),np.array(alt),H

def get_arrayids(args):
    if args.broad:
        row = ropen(args.arrayfilename,1).next()
        return np.array(row.split()[1:])
    elif args.vcfin:
        return vcfReader(args.arrayfilename).getids()
    elif args.illuminus:
        row = ropen(args.arrayfilename,1).next()
        row = row.split()[3:]
        n = len(row)/2
        return np.array([row[idx][:-1] for idx in range(0,n*2,2)])
    else:
        header = ropen(args.arrayfilename,1).next().split()
        return np.array([header[idx][:-2] for idx in range(5,len(header),2)] )


def get_probeids(args):
    arrfile = ropen(args.arrayfilename,500)
    row  = arrfile.next()
    if ' ' in row: delim=' '
    elif '\t' in row: delim = '\t'
    else: 
        print row
        print sys.exit(fname+' does not appear to be delimited by tabs or space\nAre you sure your array file is in Illuminus format?')

    probeids = []
    if args.illuminus:
        def f(row): return row[:100].split()[0]
        probeiterator = itertools.imap(f,arrfile)
        for i in probeiterator: probeids.append(i)
#        for row in arrfile: probeids.append(row[:200].split(delim,1)[0])
    elif args.broad:
        def f(row): return row[:100].split()[0][:-2]
        probeiterator = itertools.imap(f,arrfile)
        for count,i in enumerate(probeiterator): 
            probeids.append(i)
            probeiterator.next()
            if count %10000==0: print ".",
        # for row in arrfile: 
        #     probeids.append(row[:200].split(delim,1)[0][:-2])
        #     arrfile.next()
    else:
        sys.exit('Not in illuminus or broad format')
    print
    return probeids


def get_positions(args):
    if args.chromosome!=None:    args.chromosome = str(args.chromosome)
    if args.arrayonly and args.seqfilename=='':
        print sys.exit("--arrayonly is not a valid argument when no sequence vcf is provided")

    if not (args.vcfin or args.broad or args.illuminus):
        positions = get_pos_from_chiamo(args.arrayfilename,args.chromosome)
    elif args.vcfin:
        positions = get_pos_from_vcf(args.arrayfilename,args.chromosome)
    elif args.snpinfo == None:
        print "snpinfo file is required for broad or illuminus input. The format is:"
        print "probid\tchromosome\tposition\tstrand\tallele1\tallele2"
        print "For example: "
        print "rs1000121       20      16223957        -       A       G"
        sys.exit()
    else:
        start_time0 = time.time()
        snpinfo = get_snpinfo(args.snpinfo)
        # print "1. Took ",time.time() - start_time0," seconds"
        # start_time0 = time.time()
        probeids = get_probeids(args)
        # print "2. Took ",time.time() - start_time0," seconds"
        positions = [[],[],[]]
        for linenum,i in enumerate(probeids):
            if i in snpinfo:
                chrom = snpinfo[i][0]
                if args.chromosome==None or args.chromosome==chrom:
                    positions[0].append(chrom)
                    positions[1].append(snpinfo[i][1])
                    positions[2].append(linenum+1)
            else:
                print "WARNING:",i,"not found in SNP information"
                                        # start_time0 = time.time()
        positions = [np.array(val) for val in positions]
        print "Initial pass of data took ",time.time() - start_time0," seconds"
        
    return positions


def get_sequence_likelihoods(filename,chromosome,positions,ids,snpinfo,logfile):
    # print "Reading sequence data from",filename
    start_time = time.time()
    nsample = len(ids)        
    ret = {}
    denominator = np.empty(len(ids),np.float)
    # for chrom,pos in zip(chromosomes,positions):
    stored = 0
    seqfile = vcfReader(filename,positions=positions)
    seqids = seqfile.getids()
    ii=match(ids,seqids)
    chromosomes = set(chromosome)
    for chrom,pos,ref,alt,info,values in seqfile:
        if stored == len(positions): break
        if chrom in chromosomes:
            if pos in snpinfo:
                stored+=1
                # if stored % 1000 == 0: print ".",
                snp = snpinfo[pos][0]
                if len(alt)==1:
                    if (ref==snp[3] and alt[0]==snp[4]) or (ref==snp[4] and alt[0]==snp[3]):
                        processrow=True
                    else:
                        processrow = False
                        logfile.write("WARNING: sequence alleles do not match array annotation at "+str(chrom)+":"+str(pos)+" (calling loci without sequence data)\n")
                else:
                    processrow = False
                    logfile.write("WARNING: multiple alternate alleles in sequence data for "+str(chrom)+" "+str(pos)+" (calling loci without sequence data)\n")
            else:
                processrow=False
            if processrow:
                l = np.empty((nsample,3),np.float,order="F")
                l[:] = np.nan
                if 'GL' in info:
                    idx = info.index('GL')
                    for i,val in enumerate(values): 
                        try: l[ii[i]] = np.power(10.,np.fromstring(values[i].split(":")[idx],np.float,3,","))
                        except: pass
                elif 'PL' in info:
                    idx = info.index('PL')
                    for i,val in enumerate(values): 
                        try: l[ii[i]] = np.power(10.,np.fromstring(values[i].split(":")[idx],np.float,3,",")/-10.)
                        except: pass
                else:
                    logfile.write('WARNING: '+str(chrom)+" "+str(pos)+' does not have a valid genotype likelihood in the sequence data (GL or PL)\n')
                denominator[:] = l.sum(1)
                l2 = l.T.view()
                l2 /= denominator
                flip = (ref==snp[4]) != (snp[8]=='ALT:REF')
                if flip:
                    l = l[:,::-1]
                    tmp=ref
                    ref=alt
                    alt=ref
                ret[pos] = {'chrom':chrom,'pos':pos,'id':snp[2],'ref':ref,'alt':alt,'lik':l,'flip':flip}

    # print "Took ",time.time() - start_time," seconds"    
    # print "\nStored likelihoods for",stored,"SNPs"
    return ret

def array_setup(args,positions=None):
    if args.arrayonly and args.seqfilename=='':
        sys.exit("--arrayonly is not a valid argument when no sequence vcf is provided")
    if positions != None:
        use_positions = True
    elif args.positionfile != "":
        positions = [int(val) for val in open(args.positionfile).read().split() if val !='']
        print "Calling",len(positions),"positions listed in",args.positionfile
        use_positions = True
    else:    
        print "Calling all positions in",args.arrayfilename
        use_positions = False
        positions = None

    print 'Reading array data from',args.arrayfilename
    if not(args.broad or args.illuminus or args.vcfin):
        arrids,snpinfo,array_chromosomes,array_positions,signal =  get_chiamo_array_data(args.arrayfilename,positions)        
    elif args.vcfin:
        arrids,pop,snpinfo,array_chromosomes,array_positions,signal = get_array_data(args.arrayfilename,positions,args.pop,args.chromosome)        
    elif args.snpinfo == None:
        print "snpinfo file is required for broad or illuminus input. The format is:"
        print "probid\tchromosome\tposition\tstrand\tallele1\tallele2"
        print "For example: "
        print "rs1000121       20      16223957        -       A       G"
        sys.exit()
    else:
        if args.broad:
            arrids,snpinfo,array_chromosomes,array_positions,signal = get_broad_array_data(args.arrayfilename,args.snpinfo,positions,args.pop,args.chromosome)
        elif args.illuminus:
            arrids,snpinfo,array_chromosomes,array_positions,signal = get_illuminus_array_data(args.arrayfilename,args.snpinfo,positions,args.pop,args.chromosome)

    if args.pop!="":
        # print "Reading populations from",args.pop
        poptable,pop = setup_populations(ropen(args.pop).read().split())
    else:
        # print "No population information provide I am assuming all samples are from the same population"
        pop = np.array(["ALL" for idx in range(len(arrids))])

    if not use_positions:
        positions = list(set(array_positions))
    else:
        positions = list(set(positions).intersection(set(array_positions)))

    positions.sort()

    return positions,arrids,pop,snpinfo,array_chromosomes,array_positions,signal

def scan(filename): return ropen(filename).read().split()

def randomstring(N=10):
    return ''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(N))

def match(a,b): return np.concatenate([np.where(a==i)[0] for i in b if i in a and i in b])
    
def table(x):
    y = np.unique(x)
    return dict([(val,(x==val).sum()) for val in y])

def print_table(t):
    print string.join(t.keys(),"\t")
    print string.join([str(val) for val in t.values()],"\t")

def checkfile(filename):
    if not os.path.isfile(filename):
        sys.exit(filename+" does not exist!")


class gzipReader:

    def __init__(self,filename,buf=500):

        checkfile(filename)

        if filename[-3:] != ".gz":
            sys.exit("input needs to be gzipped!")

        self.fileobj=open(filename,"rb",buffering=buf*(2**20))
        self.gfile = gzip.GzipFile(fileobj=self.fileobj)
        self.infile = io.BufferedReader(self.gfile)

    def __iter__(self):
        return self

    def next(self):
        return self.infile.tell(),self.infile.next()            

    def goto(self,filepos):
        self.infile.seek(filepos)



def ropen(filename,buf_in_MB=100):
    if not os.path.isfile(filename):
        sys.exit(filename+" does not exist!")

    else:
        if buf_in_MB != None:
            if filename[-3:]==".gz":    
                # p = subprocess.Popen(["zcat", filename],stdout = subprocess.PIPE)
                # return cStringIO.StringIO(p.communicate()[0])
                #                return io.BufferedReader(p.communicate()[0],buf_in_MB*(2**20))
                # return            gzip.GzipFile(fileobj=open(filename,"rb",buffering=buf_in_MB*(2**20)))
                # return  io.BufferedReader(gzip.open(filename),buf_in_MB*(2**20))
                return  io.BufferedReader(gzip.GzipFile(fileobj=open(filename,"rb")),buf_in_MB*(2**20))# GzipFile(fileobj=open(filename,"rb",buffering=buf_in_MB*(2**20))))
            else: return open(filename,"rb",buf_in_MB*(2**20))
        else:
           if filename[-3:]==".gz": return gzip.GzipFile(fileobj=cStringIO.StringIO(open(filename).read()))
           else: return cStringIO.StringIO(open(filename).read())

def wopen(filename,force=False,buf_in_MB=100):
    if not force and os.path.isfile(filename):
        sys.exit(filename+" exists!")
    else:
        if buf_in_MB != None:
            if filename[-3:]==".gz":
                # return gzip.GzipFile(fileobj=open(filename,"wb",buffering=buf_in_MB*(2**20)))
                return  io.BufferedWriter(gzip.open(filename,"wb"),buf_in_MB*(2**20))
                                          # GzipFile(fileobj=open(filename,"wb",buffering=buf_in_MB*(2**20))))

            else: return open(filename,"wb",buf_in_MB*(2**20))
        else:
            if filename[-3:]==".gz": return gzip.GzipFile(fileobj=open(filename,"wb"))
            else: return open(filename,"wb")


def checkfiledoesnotexist(filename):
    if os.path.isfile(filename):
        sys.exit(filename+" exists!") 

def genotype_likelihood(info,vals):
    def f(x):
        try: return x.split(':')[idx].split(',')
        except: return [0.,0.,0]

    if 'GL' in info:
        idx = info.index('GL')
        return  np.power(10, np.float64([f(val) for val in vals],order='F'))
    elif 'PL' in info:
        idx = info.index('PL')
        return  np.power(10,-np.float64([f(val) for val in vals],order='F')/10.)
    else:
        sys.exit("Could not find a likelihood (GL or PL) in the supplied sequence file.\nExiting...")

def get_refinfo(filename,chrom=""):
    print "Building refinfo from",filename,"..."
    infile = ropen(filename)
    refinfo={}
    row = infile.next()
    while row[0]=="#": row=infile.next()
    row=    infile.next()
    for row in infile:
        info = row[:1000].split()#"\t",5)[:5]
        refinfo[info[2]]=info
    return refinfo


def get_raf(filename,positions):
    infile = ropen(filename)
    print "Reading allele frequencing information from",filename
    raf = {}
    for row in infile:
        x = row.split()
        pos = int(x[1])
        if pos in positions:
            raf[pos] = 1.-float(x[6])
    return raf

def get_snpinfo(filename):
    print "Reading SNP information from",filename

    baselookup = {'C':'G','G':'C','A':'T','T':'A'}
    infile = ropen(filename)

    snpinfo = {}
    i=0
    for row in infile:
        x = row.split()
        try: #CHROM POS REF ALT
            if x[4]=="+": snpinfo[x[0]] = [x[1],int(x[2]),x[5][0],x[5][1]]
            else: snpinfo[x[0]] = [x[1],int(x[2]),baselookup[x[5][0]],baselookup[x[5][1]]]
            i+=1
        except:
            pass

    print "Stored probe information for",i,"SNPs"
    return snpinfo


def setup_populations(pops,verbose=False):
    popcount = {}
    for p in pops:
        if p not in popcount:  popcount[p]=1
        else:  popcount[p]+=1

    if verbose:
        for p in popcount:print p,"\t",
        print "Total"
        for p in popcount:   print popcount[p],"\t",
        print len(pops)

    pops = np.array(pops)
    pop_table=dict([(popname,i) for i,popname in enumerate(np.unique(pops))])
    pop_factor = np.array([pop_table[s] for s in pops])
    return pop_table,pop_factor

def file_len(fname):
    if fname[-3:]==".gz":
        f = gzip.open(fname,"rb")
    else:
        f = open(fname,"rb")
    for i, l in enumerate(f):
        pass
    return i + 1
        
def read_table(filename,dtypes=None,keep=None,header=False):
    f =  (ropen(filename)).read().split("\n")
    nvar = len(f[0].split())
    if dtypes !=None:
        if len(dtypes)!=nvar:
            sys.exit("Error:len(dtype)!=nvar")

    x = [[] for idx in range(nvar)]
    if header: 
        h=f[0].split()
        f=f[1:]
    for row in f:        
        for idx,val in enumerate(row.split()):
            x[idx].append(val)
    if keep==None: keep=range(len(x))
    if header:
        return dict([(k,np.array(x[idx])) for k,idx in zip(h,keep)])
    else:
        return [np.array(x[idx]) for idx in keep]
    

def read_sequence_data(filename):
    if filename[len(filename)-3:]==".gz":
        f = csv.reader(gzip.open(filename,"rb"),delimiter="\t")
    else:
        f = csv.reader(open(filename,"rb"),delimiter="\t")
    f.next()
    row = f.next()
    ids = row[9:]
    n = len(ids)
    d = {}
    ref = {}
    for row in f:
        chrom = row[0]
        if not chrom in d:
            d[chrom] = {}
            ref[chrom]={}
        pos = int(row[1])
        ref[chrom][pos] = [row[3]]+row[4].split(",")
        liks = [val.split(":")[0] for val in row[9:]]
        d[chrom][pos] = np.array([[pow(10,-float(x2)/10.0) for x2 in x1.split(",")] for x1 in liks])   
    return ids,ref,d
  
                

def read_illumina_data(bead,snpinfo,infile):
    print "Reading SNP information...."
    total_n_snps=sum([len(snpinfo[k]) for k in snpinfo.keys()])
    f = csv.reader(gzip.open(infile,"rb"),delimiter="\t")
    for i in range(6): f.next()#skip headers
    row = f.next()
    n = int(row[1])#get sample size
    for i in range(3): f.next()#skip headers
    print n,"samples..."
    nsnps=sum([sum([snpinfo[k1][k2][0]==bead for k2 in snpinfo[k1].keys()]) for k1 in snpinfo.keys()])
    print "Read SNP data for bead",bead,"(",nsnps," SNPs)"
    data = np.ndarray(shape=[n,nsnps,2],dtype=np.float32)
    snp_index=0
    sample_index=-1
    current_sample = "I am not a sample"
    for i,row in enumerate(f):
        new_sample = row[1]
        if not new_sample == current_sample:
            print sample_index,
            sample_index+=1
            snp_index=0
            current_sample = new_sample
        chrom = row[22]
        pos = int(row[23])
        if snpinfo[chrom][pos][0]==bead:
            data[sample_index][snp_index][0] = float(row[7])
            data[sample_index][snp_index][1] = float(row[8])
            snp_index+=1
            
    return data

def qqnormalise(data):
    print "Normalising..."
    nsnps=data.shape[0]
    n=data.shape[1]
    for j in range(2):
        indices = data[:,:,j].argsort(0)
        xbar=data[indices,range(n),j].mean(1)
        for i in range(n):
            data[indices[i],i,j] = xbar
    return data


def get_ids(infile):
    f = csv.reader(gzip.open(infile,"r",pow(2,20)),delimiter="\t")
    ids = f.next()[1:]
    return ids

def read_snp_table(infile,poslist): #positions of snps IN THE FILE
    f = csv.reader(gzip.open(infile,"r",pow(2,20)),delimiter="\t")
    ids = f.next()[1:]
    d = np.ndarray((len(poslist),len(ids),2),"int16")

    skips = [poslist[0]+1] + [poslist[i+1]-poslist[i] for i in range(len(poslist)-1)]
    for i,skip in enumerate(skips):
        for j in range((skip-1)*2): f.next()
        row = f.next()
        d[i,:,0] = [int(x) for x in row[1:]]
        row = f.next()
        d[i,:,1] = [int(x) for x in row[1:]]
    return ids, d

class dataFrame:
    def __init__(self, *args, **kwargs):
        if 'filename' in kwargs and 'data' in kwargs:
            sys.exit( "Only one of filename and data should be specified!")

        if 'filename' not in kwargs and 'data' not in kwargs:
            sys.exit("One of filename and data must be specified!")
            
        if 'filename' in kwargs:
            self.data = read_table(kwargs['filename'])
        if 'data' in kwargs:
            n = set([len(v) for k,v in kwargs['data'].iteritems()])
            if len(n)!=1: #check all arrays have same length
                print 'dataFrame: arrays are of different length'
                print n
                sys.exit()
            else:
                self.names = kwargs['data'].keys()
                self.n = len(kwargs['data'][self.names[0]])
                self.data = kwargs['data']

    def __getitem__(self, index):
        if type(index)==str:
            return self.data[index]
        elif type(index)==int:
            return [(v[index]) for k,v in self.data.iteritems()]
        elif type(index)==list:
            return dataFrame(dict([(k,self.data[k]) for k in index]))
        else:
            return dataFrame(data=dict([(k,v[index]) for k,v in self.data.iteritems()]))

    def __setitem__(self, index,val):
        if type(index)==str:
            self.data[index]=val
        else:
            sys.exit('dataFrame: I only take dict keys')

class vcfReader:
    def __init__(self,filename,positions=None,buf=100):
        checkfile(filename)
        if filename[-3:] == ".gz":
            if buf<=0: 
                self.seqfile = gzip.GzipFile(fileobj=cStringIO.StringIO(open(filename).read()))
            else: 
                self.seqfile = ropen(filename,buf)
        else:
            self.seqfile = ropen(filename,buf)
        row = self.seqfile.next()
        while "#CHROM" not in row: 
            row = self.seqfile.next()
        self.seqids = np.array(row.split()[9:])
        self.nsample = len(self.seqids)
        if positions !=None:
            self.positions=set(positions)
        else: self.positions = None
        self.linenum=0

    def getids(self):
        return self.seqids

    def __iter__(self):
        return self

    def next(self):
        row = self.seqfile.next()            
        pos = int(row[:20].split("\t",2)[1])
        if self.positions!=None:
            while pos not in self.positions:
                row = self.seqfile.next()            
                pos = int(row[:20].split("\t",2)[1])
        # print pos
        data = row.split("\t",9)
        chrom = data[0]
        pos = int(data[1])
        ref = data[3]
        rsid = data[2]
        alt = data[4].split(",")[0]
        return chrom,pos,rsid,ref,alt,data[5],data[6],data[7],data[8].split(":"),data[9].split()# np.array([genotype_likelihood(val,pl_index,gq_index) for val in row[9:]])
    
    def getHaps(self): #returns H: NSNP X NSAMPLE matrix of haps
        def tohap(s):
            try:
                return (int(s[0]),int(s[2]))
            except:
                return (3,3)
        positions=[]
        chroms=[]
        refs=[]
        alts=[]
        infos=[]
        rsids=[]
        assert self.linenum==0
        tmpH = []
        for chrom,pos,rsid,ref,alt,info,vals in self:
            gtidx = [idx for idx,val in enumerate(info) if val=="GT"][0]
            chroms.append(chrom)
            positions.append(pos)
            rsids.append(rsid)
            refs.append(ref)
            alts.append(alt)
            infos.append(info)
            for v in vals:
                assert v[1]=="|"

            for v in vals: 
                v.replace(".","3")
            tmpH.append([tohap(val.split(":")[gtidx]) for val in vals])

        self.nsnp = len(tmpH)
        H = np.empty((self.nsnp,self.nsample*2),np.uint8)
        for i in range(self.nsnp):
            for j in range(self.nsample):
                H[i,j*2] = tmpH[i][j][0]
                H[i,j*2+1] = tmpH[i][j][1]

        return chroms,np.array(positions,np.int),rsids,refs,alts,infos,H

class genReader:
    def __init__(self,filename,nsample,buf=5000,start=None,stop=None):
        checkfile(filename)
        if filename[-3:] == ".gz":
            self.genfile = ropen(filename,100)

        self.genotypes = np.empty((buf,nsample*3),np.float32,order='F')
        self.dosage = np.empty((nsample,buf),order='F')
        self.i1 = np.arange(1,nsample*3,3)
        self.i2 = np.arange(2,nsample*3,3)
        self.buf = buf
        self.linenum = 0
        self.stop = stop

        if start!=None:
            for row in range(start):
                self.genfile.next()
                self.linenum += 1
         
        print "Reading genotypes from",filename,"line number",start,"to",stop
        
    def __iter__(self):
        return self

    def next(self):
        nsnp = 0
        pos = []
        rsid1=[]
        rsid2=[]
        for i in range(self.buf):
            try:
                if self.linenum==self.stop: 
                    raise StopIteration

                row = self.genfile.next().split()
                rsid1.append(row[0])
                rsid2.append(row[1])
                pos.append(row[2])
                self.genotypes[i] = np.array(row[5:],np.float32)
                nsnp +=1 
                self.linenum += 1
            except StopIteration:
                if nsnp>0:
                    break
                else:
                    raise StopIteration
        
        self.dosage[:] = (self.genotypes[:,self.i1] + 2*self.genotypes[:,self.i2]).T

        return rsid1,rsid2,pos,self.dosage[:,:nsnp]

    def getRegion(start,stop):
        nsnp = 0
        pos = []
        rsid1=[]
        rsid2=[]
        g = []
        pos = 0
        while end > pos:
            row = self.genfile.next()
            pos = int(row[:100].split()[2])
            if start <= pos and pos <= end:
                data = row.split()
                rsid1.append(data[0])
                rsid2.append(data[1])
                pos.append(data[2])
                g.append(data[5:])
        g = np.array(g,np.float32)
        dosage = (g[:,self.i1] + 2*g[:,self.i2]).T

        return rsid1,rsid2,pos,dosage

class chunkReader:
    def __init__(self,filename,bufsize):
        self.f = ropen(filename)
        self.bufsize = bufsize
        self.buf = ""
        print "Reading lines from",filename,"with a",bufsize,"byte buffer"

    def __iter__(self):
        return self

    def readLines(self,nlines):
        if self.buf.count("\n") < nlines:
            self.buf = self.buf.join(self.f.read(self.bufsize))
            return self.readLines(nlines)
        else:
            cr_index = [idx for idx,val in enumerate(self.bufsize) if val=="\n"]
            buf_idx = cr_index[nlines-1]+1
            ret = self.buf[:buf_idx]
            self.buf = self.buf[buf_idx:]
            return ret

def zapsmall(x,small=1e-6): x[x<small] = 0.0
        
class vcfWriter:
    def __init__(self,filename,ids=None,significant_figures=3,header=None):
        if filename[-7:] != ".vcf.gz": filename+=".vcf.gz"
        self.out = wopen(filename)
        self.sf = significant_figures
        self.npsf = "S"+str(self.sf*2)
        self.gpformat = "{0:."+str(self.sf)+"f}"
        self.out.write("##fileformat=VCFv4.0\n")
        self.out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')            
        if header != None:
            for s in header:
                self.out.write("##"+s+"\n")                
        if ids!=None:        
            try:
                self.out.write(string.join(["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"]+ids.tolist(),'\t')+'\n')
            except AttributeError:
                self.out.write(string.join(["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"]+ids,'\t')+'\n')              
        self.gtags = np.array(['0/0','0/1','1/1','./.'])

    def format(self,x):
        if type(x)==np.string_: return x
        else: return ",".join(x)

    def writerow(self,chrom,pos,snpid,ref,alt,data,qual=".",info=".",filterstring='.'):

        if 'GP' in data and 'GT' not in data: #EXTRACTS HARD GENOTYPES FROM PROBS
            if data['GP'].shape[1]!=4 and data['GP'].shape[1]!=3:
                sys.exit( "Error: GP has invalid number of columns")

            data['GT'] = self.gtags[data['GP'].argmax(1)]
        elif 'GT' in data:
            if type(data['GT'][0])==str:
                data['GT'] = np.array(data['GT'])
            else:
                data['GT'] = self.gtags[data['GT']]
        else:
            sys.exit( 'no genotypes provided to vcfWriter:writerow')


        if len(data.keys()) > 1:#OTHER DATA AVAILABLE - GP, GL ETC.
            data2 = {}
            for k in data.keys():
                if k!='GT':
                    data2[k] = np.array(data[k].round(self.sf),self.npsf)
            FORMAT="GT:"+string.join(data2.keys(),":")
            vals = [":".join(val3) for val3 in zip(data['GT'],[string.join([self.format(val2) for val2 in val],":") for val in zip(*data2.values())])]
        else:#WE JUST HAVE GENOTYPE DATA
            FORMAT='GT'
            vals = data['GT'].tolist()

        self.out.write(string.join([str(chrom),str(pos),snpid,ref,alt,qual,filterstring,info,FORMAT]+vals,'\t')+'\n')

    def close(self):
        self.outfile.close()
        
        
#snpnfo is 
#1       723918  SNP1-713781     G       A       0       PASS    ;NIT=11;AFAIL=7e-05;SFAIL=0.00018;HWE=0.364315400333;RETRY=0    GT:GP idx
#from CHROM POS ID REF ALT STRAND
#      0    1   2   3  4   5
def get_illuminus_array_data(fname,snpinfo_filename,positions=None,pop=None,chromosome=None,seek=0):
#    if chromosomes==None:    chromosomes = set([str(val) for val in range(1,23)])
    baselookup = {'C':'G','G':'C','A':'T','T':'A'}
    arrfile = ropen(fname,500)
    snplook = get_snpinfo(snpinfo_filename)
    row  = arrfile.next()
    if ' ' in row: delim=' '
    elif '\t' in row: delim = '\t'
    else: 
        print row
        sys.exit( fname + ' does not appear to be delimited by tabs or space\nAre you sure your array file is in Illuminus format?')

    row = row.split()[3:]
    n = len(row)/2
    arrids = np.array([row[idx][:-1] for idx in range(0,n*2,2)])

    signal = []
    snpinfo = {}
    i=0
    array_positions=[]
    array_chroms=[]
    skipped = {}
    if positions != None:
        number_to_store = len(positions)
        positions = set(positions)
    else: number_to_store = -1
    if seek>0: arrfile.seek(seek)
    for row in arrfile:
        if i==number_to_store: break
        row = row.split(delim,3)
        probeid = row[0]
        if probeid in snplook:
            snp = snplook[probeid]
            if snp[0]==chromosome:
                pos = snp[1]
                if positions==None or pos in positions:
                    signal.append(np.uint16(np.fromstring(row[3],np.float,sep=" ").reshape(n,2)*2000))
                    snp = snp[:2]+[probeid]+snp[2:]+["." for idx in range(4)]+[i]
                    array_positions.append(pos)
                    array_chroms.append(snp[0])                   
                    if pos in snpinfo: snpinfo[pos].append(snp)
                    else: snpinfo[pos] = [snp]
                    i+=1
                    if i % 1000 ==0: print ".",
            else:
                if snp[1] not in autosomes:
                    if snp[1] in skipped: skipped[snp[1]] += 1
                    else: skipped[snp[1]] = 1

#    if len(skipped)>0:
#        print "WARNING: removed non-autosomal SNPs:"
#        for k in skipped:  print skipped[k],'from',k
    if i==0: 
        sys.exit( "No SNPs found!\nLooks like a problem with"+fname+"\nExiting...")

#    print "\n",i,"SNPs stored"
    return arrids,snpinfo,array_chroms,array_positions,np.array(signal)

def get_broad_array_data(fname,snpinfo_filename,positions=None,pop=None,chromosome=None,seek=0):
#    if chromosomes==None: chromosomes = [str(val) for val in range(1,23)]
    baselookup = {'C':'G','G':'C','A':'T','T':'A'}
    arrfile = ropen(fname,500)
    snplook = get_snpinfo(snpinfo_filename)
    
    row  = arrfile.next()
    if ' ' in row: delim=' '
    elif '\t' in row: delim = '\t'
    else: 
        print row
        sys.exit( fname +' does not appear to be delimited by tabs or space')

    
    arrids = np.array(row.split()[1:])
    print len(arrids),"samples in",fname
    signal = []
    snpinfo = {}
    i=0
    array_positions=[]
    array_chroms=[]    
    skipped = {}
    if positions != None:
        number_to_store = len(positions)
        positions = set(positions)
    else: number_to_store = -1
    if seek>0: arrfile.seek(seek)
    for row in arrfile:
        if i==number_to_store: break
        probeid = (row.split(delim,1)[0])[:-2]
        if probeid in snplook:
            snp = snplook[probeid]
            if snp[0]==chromosome:
                pos = snp[1]
                if positions==None or pos in positions:
                    x1 = row.split()[1:]
                    x2 = arrfile.next().split()[1:]
                    signal.append(np.vstack((np.array(x1,np.uint16),np.array(x2,np.uint16))).T)
                    snp = snp[:2]+[probeid]+snp[2:]+["." for idx in range(4)]+[i]
                    array_positions.append(pos)
                    array_chroms.append(snp[0])
                    if pos in snpinfo: snpinfo[pos].append(snp)
                    else: snpinfo[pos] = [snp]
                    i+=1
                    if i % 1000 ==0: print ".",
            else:
                if snp[0] in skipped: skipped[snp[0]] += 1
                else: skipped[snp[0]] = 1
                # if snp[1] not in autosomes: print "WARNING:",snp,"does not come from an autosome. Skipping it."

#    if len(skipped)>0:
#        print "WARNING: removed non-autosomal SNPs:"
#        for k in skipped:  print skipped[k],'from',k
    if i==0: 
        sys.exit( "Exiting" )

#    print "\n",i,"SNPs stored"
    return arrids,snpinfo,array_chroms,array_positions,np.array(signal)

# reads in array data from chiamo(ish) format - preferred method
def get_chiamo_array_data(fname,positions=None,chromosome=None,seek=0):    

    file1 = ropen(fname,500)
    header = file1.next().split()
    arrids = [header[idx][:-2] for idx in range(5,len(header),2)]
    nsample = len(arrids)
    dt = np.dtype([('chrom',np.str_,2),('pos',np.uint),('snpid',np.str_,20),('ref',np.str_,1),('alt',np.str_,1),('intensity',np.float32,(nsample,2))])
    snpinfo = {}
    intensity = []
    array_positions = []
    chromosomes = []
    i=0
    for linenum,row in enumerate(file1):
        try:
            snp = row[:100].split()[:5]
            pos = int(snp[1])
            if (chromosome==None or snp[0]==chromosome) and (positions==None or pos in positions):
                array_positions.append(pos)
                chromosomes.append(snp[0])
                x = 2000.*np.array(row.split()[5:],np.float).reshape((nsample,2))
                intensity.append(x)
                if pos in snpinfo: snpinfo[pos].append(snp+['-','-','-','-',i])
                else: snpinfo[pos] = [snp+['-','-','-','-',i]]
                i+=1
        except: sys.exit("Error reading "+fname+".  Line "+str(linenum))

    return np.array(arrids),snpinfo,np.array(chromosomes),np.array(array_positions),np.array(intensity)

    
def get_array_data(fname,positions=[],population='',chromosome=None,seek=0):    
    start_time = time.time()

    arrfile = ropen(fname,500)
    
    if positions != None: use_positions=True
    else: use_positions=False
    row = arrfile.next()#.strip("\r\n").split("\t")#arrfile[rownum].strip("\r").split("\t")
    rownum=0
    pops = []
    try:
        while row[:6] != "#CHROM":
            if "##population" in row:
                pops = (row.split("=")[1]).strip().split(",")
            row = arrfile.next()#.strip().split("\t")#arrfile[rownum].strip("\r").split("\t")
            rownum+=1
        arrids = np.array(row.strip().split()[9:])
    except:
        sys.exit( fname + " does not appear to be a VCF\nMaybe you want to try the --broad or --illuminus arguments for different input formats?")


    if len(pops)==0:
        popcount = {'ALL':len(arrids)}
        pops = ["ALL" for idx in range(len(arrids))]
        "Warning: no populations listed in array file.  Assuming all samples are from same population."
    else:
        pop_table,pops =  setup_populations(pops)

    if len(pops)>0 and len(pops) != len(arrids):
        print "Error: len(pops) != len(ids) in ",fname
        sys.exit( "Exiting..." )

    
    if seek>0:    arrfile.seek(seek)#SEEKS TO POSITION PROVIDED BY MASTER PROCESS
    array_positions=[]
    array_chromosomes=[]
    if use_positions:       signal = np.ndarray((np.int(np.ceil(len(positions)*1.1)),len(arrids),2),np.uint16)
    else:       signal = np.ndarray((nlines(arrfile),len(arrids),2),np.uint16)
    snpinfo = {}
#    import pdb;pdb.set_trace()
    i=0
    for linenum,row in enumerate(arrfile):#[(rownum+1):]:    
        if linenum % 10000 == 0: print ".",
        if row!= '':
            tmp  = row[:200].split('\t',2)# row.strip("\r\n").split("\t")
            pos = int(tmp[1])
            if not use_positions or pos in positions:
                if chromosome==None or tmp[0]==chromosome:
                    tmp = row.strip().split()
                    array_chromosomes.append(tmp[0])
                    array_positions.append(pos)
                    if pos in snpinfo: snpinfo[pos].append(tmp[:9]+[i])
                    else: snpinfo[pos] = [tmp[:9]+[i]]
                    signal[i] = np.array([np.fromstring(val,np.uint16,sep=":") for val in tmp[9:]])
                    i+=1        
        else:
            break
    if population!= "":
            print "Only calling individuals from",population,"population"
            ii=pops==pop_table[population]
            signal = signal[:,ii]
            arrids = arrids[ii]
            pops = pops[ii]
    print "Took ",time.time() - start_time," seconds"    
    return arrids,pops,snpinfo,array_chromosomes,array_positions,signal[:i]


def nlines(file):
    return_to = file.tell()
    linenum = 0
    for row in file: linenum+=1
    file.seek(return_to)
    return linenum

def get_pos_from_vcf(filename,chromosome=None): #dumps all (chrom,positions) from a vcf nlines
    if chromosome != None: chromosome = str(chromosome)
    inf = ropen(filename,500)
    pos = []
    chrom = []
    index = []
    row = inf.next()
    while row[:6] != '#CHROM': 
        row = inf.next()
    inf.next()
    idx = inf.tell()
    for row in inf:
        tmp = row[:100].split("\t",2)
        if tmp[0]==chromosome or chromosome==None:
            chrom.append(tmp[0])
            pos.append(tmp[1])
            index.append(idx)
            idx = inf.tell()
    return  [np.array(chrom),np.array(pos,np.int),np.array(index,np.int)]

def get_pos_from_chiamo(filename,chromosome=None): #dumps all (chrom,positions) from a vcf nlines
    if chromosome != None: chromosome = str(chromosome)
    inf = ropen(filename,500)
    pos = []
    chrom = []
    index = []
    row = inf.next()
#    idx = inf.tell()
    idx=0
    for linenum,row in enumerate(inf):
        tmp = row[:100].split()
        if len(tmp)<=5: sys.exit(filename+" does not appear to be in the correct format (check line "+str(linenum)+")")
        if tmp[0]==chromosome or chromosome==None:
            chrom.append(tmp[0])
            try: pos.append(int(tmp[1]))
            except ValueError: sys.exit(filename+" does not appear to be in the correct format (check line "+str(linenum)+")")
            index.append(idx)
#            idx = inf.tell()
    return  [np.array(chrom),np.array(pos,np.int),np.array(index,np.int)]


def vcfdump(filename,var=None,positions=None,chromosome=None):
    file1 = ropen(filename)
    if positions==None:
        if filename[-3:]==".gz": 
            nl = 0
            for row in file1: nl+=1
            file1.seek(0)
        else: 
            nl = nlines(file1)                

    def getprob(val,idx): 
        try:
            return np.fromstring(val.split(":")[idx],np.float,sep=",")
        except: 
            return 0.
    
    def getcall(val,idx):
        try:
            tmp = val.split(":")[idx]
            return int(tmp[0])+int(tmp[2])
        except: 
            return 3

    rownum=1
    row = file1.next()
    while not row.startswith("#CHROM"):
        row = file1.next()
        rownum+=1

    ids1 = np.array(row[:-1].split()[9:])
    n = len(ids1)
    if positions == None: nsnp = nl-rownum
    else: nsnp = len(positions)
    try:    poslook = set(positions.tolist())
    except: poslook = positions
    print filename,"has",nsnp,"SNPs and",n,"individuals"
    if var != None: probs1 = np.zeros((nsnp,n),np.float)
    calls1 = 3*np.ones((nsnp,n),np.uint8)
    pos1 = -1*np.ones(nsnp,np.uint32)
    ref1 = np.empty(nsnp,"|S1")
    alt1 = np.empty(nsnp,"|S1")    
    chrom1 = []
    snpnames = []
    info = []
    i = 0

    if var=="GP":
        tmpprobs = np.empty((n,4),np.float)
    if var=="GC":
        tmpprobs = np.empty(n,np.float)
    #    print positions,chromosomes
    for rownum,row in enumerate(file1):
#        print i,rownum,
        processrow=False
        tmp = row[:200].split()            

        if tmp[1] == "-":
            print 'WARNING: blank rows found'
            processrow=False
        elif positions != None:            
            chrom = tmp[0]
            snpid = tmp[2]
            pos = int(tmp[1])
            if pos in poslook and chrom==chromosome:
                processrow = True
                vals = row[:-1].split()
        else:            
            processrow = True
            vals = row[:-1].split()
            pos = int(vals[1])

        # if pos in poslook:        
        #     print pos, pos in poslook

        if processrow:
            form =  vals[8].split(":")
            if "GT" not in form:
                processrow=False
                import pdb;pdb.set_trace()

        if processrow:
 #           print "stored"
            ref1[i] = vals[3]
            alt1[i] = vals[4]
            snpnames.append(vals[2])
            info.append(vals[7])
            pos1[i] = pos
            chrom1.append(vals[0])
            gtidx = form.index("GT")

            if var==None:
                for j,value in enumerate(vals[9:]):
                    spvalue = value.split(":")
                    try:
                        calls1[i,j] = int(spvalue[gtidx][0])+int(spvalue[gtidx][2])
                    except:
                        pass
                
            elif var=="GP": 
                varidx = form.index(var)
                for j,value in enumerate(vals[9:]):
                    spvalue = value.split(":")
                    tmpprobs[j] = np.fromstring(spvalue[varidx],np.float,sep=",")
                    try:   calls1[i,j] = int(spvalue[gtidx][0])+int(spvalue[gtidx][2])
                    except:  pass

                probs1[i] = tmpprobs.max(1)

            elif var=="GC": 
                varidx = form.index(var)
                for j,value in enumerate(vals[9:]):
                    spvalue = value.split(":")
                    try:  
                        calls1[i,j] = int(spvalue[gtidx][0])+int(spvalue[gtidx][2])
                        probs1[i,j] = float(spvalue[varidx])
                    except:   
                        pass
            i+=1

    print "vcfdump:",i,'rows stored'
    if var != None: 
        return ids1,{'chromosome':np.array(chrom1),
                     'position':pos1[:i],
                     'info':np.array(info),
                     'rsid':np.array(snpnames),
                     'ref':ref1[:i],
                     'alt':alt1[:i],
                     'genotype':calls1[:i],
                     'probability':probs1[:i]}
    
    else:           
        return ids1,{'chromosome':np.array(chrom1),
                     'position':pos1[:i],
                     'info':np.array(info),
                     'rsid':np.array(snpnames),
                     'ref':ref1[:i],
                     'alt':alt1[:i],
                     'genotype':calls1[:i]}


def decompressPlinkGenotype(x):
    g = np.zeros(4,np.uint8)
    raw = [(x%(2**val[1]))/(2**val[0]) for val in zip(range(7,-1,-1),range(8,0,-1))]
    raw.reverse()
    for i in range(0,7,2):
        if raw[i]==1 and raw[i+1]==0: g[i/2]=3
        else: g[i/2]=raw[i]+raw[i+1]
    return g

class plinkReader:
    def __init__(self,filename,snps=None,buf=False):
        self.rawdata = np.fromfile(filename+".bed",np.uint8)
        fam = read_table(filename+".fam")
        self.ids = fam[1]
        self.sex = fam[4]
        self.snpinfo = read_table(filename+".bim")
        self.nsample = len(self.ids)
        self.nlines = len(self.snpinfo[1])
        self.rowlength = int(math.ceil(self.nsample/4.))
        self.counter = 3 #where to read bytes from
        self.linenum = 0 
        self.lookup=dict([(val,decompressPlinkGenotype(val)) for val in range(256)])
        print filename,"contains",self.nlines,"SNPs and",self.nsample,"samples"
        if snps==None: 
            self.snps = None
            self.nsnp = len(self.snpinfo[0])
        else: 
            self.snps = set(np.intersect1d(snps,self.snpinfo[1]))
            if len(self.snps) < len(snps):
                print "WARNING: Only",len(self.snps),"of",len(snps),"listed in -snps were present in plink file"
            self.nsnp = len(self.snps)

    def getids(self):
        return self.ids

    def __iter__(self):
        return self
    
    def getSnp(self,snpid):
        location = np.where(snpid==self.snpinfo[1])[0][0]*self.rowlength + 3
        ret = np.empty(self.rowlength*4,np.uint8)
        for i in range(self.rowlength): ret[(i*4):(4*i+4)]=self.lookup[self.rawdata[location+i]]
        return ret[:self.nsample]

    def next(self):
        if self.snps!=None:
            while self.snpinfo[1][self.linenum] not in self.snps:
                self.counter+=self.rowlength
                self.linenum+=1

        if self.linenum==self.nlines: raise StopIteration
            

        if self.linenum<self.nlines:
            if self.linenum % 10000==0:                  
                print "SNP",self.linenum,"of",self.nlines,"processed"
            ret = np.empty(self.rowlength*4,np.uint8)
            for i in range(self.rowlength): ret[(i*4):(4*i+4)]=self.lookup[self.rawdata[self.counter+i]]
            self.counter+=self.rowlength
            self.linenum+=1
            return self.snpinfo[0][self.linenum-1],self.snpinfo[3][self.linenum-1],self.snpinfo[1][self.linenum-1],self.snpinfo[4][self.linenum-1],self.snpinfo[5][self.linenum-1],ret[:self.nsample]
        else: 
            raise StopIteration  

    def getGenotypes(self):
        ret = np.empty((self.nsnp,self.nsample),np.uint8)
        for i in range(self.nsnp):
            ret[i] = self.next()[5]
        return ret

    def getBytes(self):
        return self.rawdata

class SNPError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

def plainTextReader(filename,buf=100):
    def f(x): return x.split()
    return itertools.imap(f,ropen(filename,buf))



