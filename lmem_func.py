import numpy as np,scipy.optimize,sys
import iolib


def gen2dosage(fname,nsample):
    gen_dt = np.dtype([('rsid1',np.str_,50),('rsid2',np.str_,50),('pos',np.uint),('ref',np.str_,1),('alt',np.str_,1)]+
                      [("g"+str(idx),np.float) for idx in range(3*nsample)])
                       
    tmp = np.genfromtxt(iolib.ropen(fname),dtype=gen_dt)
    dosage = np.empty((tmp.shape[0],nsample),np.float32)
    for i in range(nsample):
        i1 = "g"+str(i*3+1)
        i2 = "g"+str(i*3+2)
        dosage[:,i] = tmp[i1] + 2*tmp[i2]
    return tmp['rsid1'],tmp['rsid2'],tmp['pos'],dosage


def get_delta(K,y):
    assert K.shape[0]==K.shape[1] 
    assert K.shape[0]==len(y)
    nsample = K.shape[0]
    S,U = np.linalg.eig(K)
    ii = S.argsort()[::-1]
    U = -U[:,ii]
    S = S[ii]
    Ut = U.T
    offset = np.dot(Ut,np.ones(nsample,np.float).reshape(nsample,1))
    uty = np.dot(Ut,y.reshape(nsample,1,))

    def loglik_delta(delta):
        denom = (S+delta).reshape((nsample,1))
        beta = (offset*uty/denom).sum() / (np.power(offset,2)/denom).sum()
        yhat = offset*beta
        sigma_g_null = (np.power(uty-yhat,2) / denom).mean()
        return 0.5 * (nsample * np.log(2*np.pi) + np.log(denom).sum() + nsample + nsample * np.log(sigma_g_null))

    solved = scipy.optimize.minimize_scalar(loglik_delta,bounds=[0.001,1000],method="bounded",tol=1e-12)
    delta_null = solved['x']
    loglik_null = -solved['fun']
    denom = (S+delta_null).reshape(nsample,1)    
    beta_null = (offset*uty/denom).sum() / (np.power(offset,2)/denom).sum()
    yhat = offset*beta_null
    sigma_g_null = (np.power(uty-yhat,2) / denom).mean()
    loglik_null = -0.5 * (nsample * np.log(2*np.pi) + np.log(denom).sum() + nsample + nsample * np.log(sigma_g_null))
    sys.stderr.flush()

    return {'offset':offset,'Ut':Ut,'denom':denom,'uty':uty,'beta_null':beta_null,'sigma_g_null':sigma_g_null,'loglik_null':loglik_null,'delta_null':delta_null}
    
def fitlmm(offset,Ut,denom,dosage,uty): 


    if Ut.shape[0]!=dosage.shape[0]:
        print "Ut.shape[0]!=dosage.shape[0]"
        sys.exit()

    nsample = dosage.shape[0]
    nsnp = dosage.shape[1]

    utx =  np.empty((nsample,2),np.float)
    utx[:,0] = offset[:,0]
#    uty /= denom
    constant1 = nsample * np.log(2*np.pi) + np.log(denom).sum() + nsample#constant to be added to the LL at each iteration
    beta_alt = np.empty((nsnp,2),np.float)
    sigma_g = np.empty(nsnp)
    beta_1_se = np.empty(nsnp)
    loglik = np.empty(nsnp)
#    print "Delta =",delta_null
    UTX = np.dot(Ut,dosage)
    residuals = np.empty(nsample,np.float)
    yhat = np.empty(nsample,np.float)
    Sxx = np.empty((2,2),np.float)
    Sxy = np.empty(2,np.float)
    bot = np.empty((2,2),np.float)
    for i in range(nsnp):
        utx[:,1] = UTX[:,i]
        Sxx[:] = np.dot(utx.T,utx)
        # Sxy[:] = (utx*uty/denom).sum(0)
        a=Sxx[0,0];b=Sxx[0,1];c=Sxx[1,0];d=Sxx[1,1]
        #bot[:] = np.array([[d,-b],[-c,a]])/(a*d-b*c)
        # beta = np.dot(bot,Sxy)
        sqrtdenom = np.sqrt(denom)
        tmp = np.linalg.lstsq(utx/sqrtdenom,uty/sqrtdenom)
        beta=tmp[0]
        rss=tmp[1]
        yhat[:] = np.dot(utx,beta)[:,0]
        residuals[:] = np.power(uty[:,0]-yhat,2)
        rss = (residuals/denom[:,0]).sum()
        # print beta,rss
        # print tmp[0],tmp[1]
        beta_alt[i] = beta.reshape((1,2))
        sigma_g[i] = rss / nsample
        loglik[i] = -0.5 * (constant1 + nsample * np.log(sigma_g[i]))
        beta_1_se[i] = np.sqrt((residuals.sum()/(nsample-2))* (a/(a*d-b*c)))

    return beta_alt,beta_1_se,sigma_g,loglik

def fitlm(dosage,y):
    nsample = dosage.shape[0]
    nsnp = dosage.shape[1]

    beta_alt = np.empty((nsnp,2),np.float)
    dosage_centred = dosage - dosage.mean(0)
    Sxx = np.power(dosage_centred,2).sum(0)
    y_mean = y.mean()
    y_centred = y - y_mean
    beta_alt[:,1] = np.dot(y_centred,dosage_centred) / Sxx
    beta_alt[:,0] = (y - (dosage * beta_alt[:,1]).T).mean(1)
    yhat = (beta_alt[:,0] + (dosage * beta_alt[:,1]))
    residuals =  np.power(y-yhat.T,2)
    rss = residuals.sum(1)
    sigma =  rss / (nsample - 2)
    loglik = -0.5 * (nsample * np.log(2*np.pi) + nsample * np.log(sigma) + rss/sigma )
    beta_1_se = np.sqrt(sigma/Sxx)

    return beta_alt,beta_1_se,sigma,loglik
