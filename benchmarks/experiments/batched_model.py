#!/usr/bin/env python3
"""Experimental batched cubic/PCA model evaluator: NumPy validation or CuPy GPU.
Not a FERRE optimizer replacement. Supports prepared 4D/5D float32 PCA slices,
default reverse-axis interpolation order, no continuum/LSF/resampling changes.
"""
import argparse
import json
from pathlib import Path
import re
import statistics
import time
import numpy as np


class BatchedModel:
    def __init__(self, header, xp=np):
        self.xp = xp
        text = Path(header).read_text()
        blocks = list(re.finditer(r'&SYNTH.*?\n\s*/',text,re.S|re.I))
        if len(blocks)!=4:
            raise ValueError('Expected the prepared DR17 PCA header')
        first=blocks[0][0]
        def numbers(key,dtype=float):
            m=re.search(r'\b'+key+r'\s*=([^\n]+)',first)
            return np.array(m[1].split(),dtype=dtype)
        self.shape=tuple(numbers('N_P',int))
        self.low=numbers('LLIMITS',np.float32)
        self.step=numbers('STEPS',np.float32)
        self.sections=tuple(numbers('NPCA',int))
        self.ncoef=int(numbers('NPIX',int)[0])
        if len(self.shape) not in (4,5) or self.ncoef!=900 or len(self.sections)!=12:
            raise ValueError('Only the prepared DR17 4D/5D slices are supported')
        if re.search(r'\b(CONSTANT|TRANSPOSED)\s*=',text[:blocks[-1].end()],re.I):
            raise ValueError('Nondefault CONSTANT/TRANSPOSED are unsupported')
        self.nvar=self.ncoef//len(self.sections)
        data=np.fromstring(text[blocks[-1].end():],sep=' ',dtype=np.float32)
        data=data.reshape(self.nvar+2,sum(self.sections))
        self.mean=xp.asarray(data[0])
        self.basis=xp.asarray(data[2:].T.copy())
        binary=Path(header).with_suffix('.unf')
        grid=np.fromfile(binary,dtype='<f4')
        if grid.size!=np.prod(self.shape)*self.ncoef:
            raise ValueError('Grid size disagrees with header')
        self.grid=xp.asarray(grid.reshape(-1,self.ncoef))
        self.strides=np.array([np.prod(self.shape[i+1:],dtype=np.int64)
                               for i in range(len(self.shape))],dtype=np.int64)
        self.digits=xp.asarray(np.array(list(np.ndindex(*(4,)*len(self.shape))),dtype=np.int64))
        self.low_device=xp.asarray(self.low)
        self.range_device=xp.asarray(self.step*(np.array(self.shape,dtype=np.float32)-1))
        self.shape_device=xp.asarray(np.array(self.shape,dtype=np.int64))
        self.strides_device=xp.asarray(self.strides)

    def evaluate(self, physical):
        xp=self.xp
        physical=xp.asarray(physical,dtype=xp.float32)
        p=(physical-self.low_device)/self.range_device
        # Caller must keep all points inside [0,1); validation is outside timing.
        t=p*(self.shape_device.astype(xp.float32)-1)+xp.float32(1)
        cell=xp.floor(t).astype(xp.int64)
        offset=xp.where(t<2,1,xp.where(t>=self.shape_device-1,-1,0))
        coords=cell[:,None,:]-2+offset[:,None,:]+self.digits[None,:,:]
        indices=xp.sum(coords*self.strides_device,axis=2)
        work=self.grid[indices]
        half=xp.float32(.5);third=xp.float32(1/3);three=xp.float32(3)
        for axis in range(len(self.shape)-1,-1,-1):
            d=(t[:,axis]-cell[:,axis].astype(xp.float32))[:,None,None]
            o=xp.float32(1)-d
            y1,y2,y3,y4=(work[:,j::4,:] for j in range(4))
            c0=y2+third*(half*(y3-y1))
            c1=y3-third*(half*(y4-y2))
            cubic=y2*(o**3)+y3*(d**3)+c0*(three*d*o**2)+c1*(three*d**2*o)
            lower=y1*(o**2)+y2*(d**2)+(y2-half*(half*(y3-y1)))*(xp.float32(2)*d*o)
            upper=y3*(o**2)+y4*(d**2)+(y3+half*(half*(y4-y2)))*(xp.float32(2)*d*o)
            flag=offset[:,axis,None,None]
            work=xp.where(flag==0,cubic,xp.where(flag==1,lower,upper))
        coef=work[:,0,:]
        out=xp.empty((len(physical),sum(self.sections)),dtype=xp.float32)
        start=0
        for section,count in enumerate(self.sections):
            out[:,start:start+count]=coef[:,section*self.nvar:(section+1)*self.nvar] @ self.basis[start:start+count].T
            start+=count
        return out+self.mean


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--header',type=Path,required=True)
    p.add_argument('--points',type=Path,required=True,help='Numeric physical parameter matrix (no IDs)')
    p.add_argument('--reference',type=Path,required=True,help='FERRE NOV=0 OFFILE for those points')
    p.add_argument('--backend',choices=['numpy','cupy'],default='numpy')
    p.add_argument('--batch-size',type=int,default=8)
    p.add_argument('--repeats',type=int,default=3)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    if a.batch_size<1 or a.repeats<1:p.error('Positive batch size/repeats required')
    xp=np
    if a.backend=='cupy':
        try:
            import cupy as xp
            if xp.cuda.runtime.getDeviceCount()<1:raise RuntimeError('No CUDA GPU')
        except (ImportError,RuntimeError) as e:
            p.error('CuPy with a working CUDA GPU is required: '+str(e))
    def sync():
        if a.backend=='cupy':xp.cuda.get_current_stream().synchronize()
    begin=time.perf_counter();model=BatchedModel(a.header,xp);sync()
    setup=time.perf_counter()-begin
    points=np.atleast_2d(np.loadtxt(a.points,dtype=np.float32))
    normal=(points-model.low)/(model.step*(np.array(model.shape,dtype=np.float32)-1))
    if not np.isfinite(normal).all() or np.any(normal<0) or np.any(normal>=1):
        p.error('All parameter vectors must be inside grid bounds, excluding upper endpoints')
    ref=np.atleast_2d(np.loadtxt(a.reference))
    if ref.shape!=(len(points),sum(model.sections)):p.error('Reference shape differs')
    # Prepare parameter batches once: persistent grid and inputs on device.
    batches=[xp.asarray(points[i:i+a.batch_size]) for i in range(0,len(points),a.batch_size)]
    sync()
    output=[]
    begin=time.perf_counter()
    for batch in batches:output.append(model.evaluate(batch))
    sync();first=time.perf_counter()-begin
    got=np.concatenate([xp.asnumpy(x) if a.backend=='cupy' else x for x in output])
    passed=bool(np.isfinite(got).all() and np.allclose(got,ref,atol=1e-5,rtol=1e-5))
    times=[]
    for repeat in range(a.repeats):
        sync();begin=time.perf_counter()
        for batch in batches:result=model.evaluate(batch)
        sync();times.append(time.perf_counter()-begin)
    # Include per-batch host input/output transfers separately, excluding grid setup.
    sync();begin=time.perf_counter()
    for i in range(0,len(points),a.batch_size):
        result=model.evaluate(xp.asarray(points[i:i+a.batch_size]))
        if a.backend=='cupy':xp.asnumpy(result)
    sync();transfer_inclusive=time.perf_counter()-begin
    report={'backend':a.backend,'points':len(points),'batch_size':a.batch_size,
            'grid_setup_seconds':setup,'first_pass_seconds':first,'warm_seconds':times,
            'warm_median_seconds':statistics.median(times),
            'transfer_inclusive_seconds':transfer_inclusive,
            'max_absolute_flux_difference':float(np.max(np.abs(got-ref))),
            'agreement_passed':passed,
            'scope':'Batched cubic interpolation and PCA reconstruction only; no optimizer/continuum/errors',
            'gpu_executed':a.backend=='cupy'}
    if a.backend=='cupy':
        report['device']=str(xp.cuda.runtime.getDeviceProperties(xp.cuda.Device().id)['name'])
        report['cupy_version']=xp.__version__
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
    raise SystemExit(0 if passed else 1)


if __name__=='__main__':main()
