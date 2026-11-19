"""Reference-coded genomic features and training-only imputation."""
import numpy as np
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

class GenomicEncoder:
    def fit(self, x, paternal=None):
        if paternal is None:
            self.mode='diploid'
            self.p=np.asarray(x,dtype=float).mean(axis=0)/2
            self.keep=(self.p>.01)&(self.p<.99)
        else:
            self.mode='interpool'
            self.p=np.asarray(x,dtype=float).mean(axis=0)/2
            self.q=np.asarray(paternal,dtype=float).mean(axis=0)/2
            self.keep=(self.p>.01)&(self.p<.99)&(self.q>.01)&(self.q<.99)
        if not self.keep.any():raise ValueError('No segregating training markers')
        return self

    def transform(self,x,paternal=None,additive_only=False):
        p=self.p[self.keep];a=np.asarray(x,dtype=float)[:,self.keep]
        if self.mode=='interpool':
            q=self.q[self.keep];b=np.asarray(paternal,dtype=float)[:,self.keep]
            z1=a/2-p;z2=b/2-q
            one=z1/np.sqrt(np.sum(p*(1-p)))
            two=z2/np.sqrt(np.sum(q*(1-q)))
            if additive_only:return np.c_[one,two]
            w=-2*z1*z2
            d=w/np.sqrt(np.sum(4*p*(1-p)*q*(1-q)))
            return np.c_[one,two,d]
        z=(a-2*p)/np.sqrt(np.sum(2*p*(1-p)))
        if additive_only:return z
        w=np.where(a==0,-2*p*p,np.where(a==1,2*p*(1-p),-2*(1-p)**2))
        w=w/np.sqrt(np.sum((2*p*(1-p))**2))
        return np.c_[z,w]

def fit_ridge(x,y,alpha=1.0):
    model=make_pipeline(SimpleImputer(strategy='median',add_indicator=True),StandardScaler(with_std=False),Ridge(alpha=alpha))
    return model.fit(x,y)
