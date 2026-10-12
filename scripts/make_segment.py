"""Extract a time segment of the raw AP recording plus its matching motion trace,
so the MATLAB interpolation can be run on it exactly as it would run on the full file.

    python make_segment.py <t0_s> <t1_s> <outdir_name>

The displacement REFERENCE is taken from the whole recording, not the segment, so a
segment run is equivalent to the corresponding batches of a full-length run.
"""
import numpy as np, sys, shutil, os, json
from pathlib import Path
from scipy.signal import resample_poly

import recinfo
_R=recinfo.get()
AP=_R['ap']
NCH=_R['nch']; FS=30000; NFRAMES=_R['nframes']; MOTION_FS=250

def main(t0,t1,name,motion=None):
    S=AP/name; S.mkdir(exist_ok=True)
    F0,F1=int(t0*FS),min(int(t1*FS),NFRAMES)   # clamp to the recording
    t1=F1/FS
    link_from=os.environ.get('NPX_LINK_FROM')
    if link_from:
        # Reuse an already-extracted segment covering the SAME time range: a hard link, no
        # 20+ GB copy. Safe because MATLAB only ever opens the source read-only.
        srcf=AP/link_from/'continuous.dat'
        want=(F1-F0)*NCH*2
        if srcf.stat().st_size!=want:
            raise SystemExit('link source %s has %d bytes, segment needs %d' % (srcf, srcf.stat().st_size, want))
        dst=S/'continuous.dat'
        if dst.exists(): dst.unlink()
        os.link(srcf, dst)
        print('hard-linked continuous.dat from', link_from)
    else:
        src=np.memmap(AP/'continuous.dat',dtype='int16',mode='r',shape=(NFRAMES,NCH))
        out=np.memmap(S/'continuous.dat',dtype='int16',mode='w+',shape=(F1-F0,NCH))
        for i in range(0,F1-F0,500_000):
            j=min(F1-F0,i+500_000); out[i:j]=src[F0+i:F0+j]
        out.flush(); del out,src
    shutil.copy(_R['chanmap'],S/'chanMap.mat')
    if motion is None:
        dsrc, csrc = AP/'nonrigid_displacement.csv', AP/'dredge_window_centers.csv'
    else:
        M=_R['motion']
        dsrc, csrc = M/(motion+'_displacement.csv'), M/(motion+'_window_centers.csv')
    shutil.copy(csrc, S/'dredge_window_centers.csv')
    D=np.loadtxt(dsrc,delimiter=',')
    if D.ndim==1: D=D[None,:]
    ref=np.median(resample_poly(D.T,120,128,axis=0),axis=0)
    np.savetxt(S/'ref_displacement.csv',ref,delimiter=',',fmt='%.6f')
    # the motion estimate may not start at file t=0 (msit004 starts at 805 s)
    mt0=0.0
    tf=dsrc.parent/dsrc.name.replace('_displacement.csv','_t0.txt')
    if tf.exists(): mt0=float(tf.read_text().strip())
    c0,c1=int(round((t0-mt0)*MOTION_FS)),int(round((t1-mt0)*MOTION_FS))+4
    if c0<0 or c1>D.shape[1]+4:
        raise SystemExit('segment %.0f-%.0f s outside the motion estimate (%.0f-%.0f s)'
                         %(t0,t1,mt0,mt0+D.shape[1]/MOTION_FS))
    sl=D[:,c0:min(c1,D.shape[1])]
    if sl.shape[1]==0: raise SystemExit('empty motion slice')
    np.savetxt(S/'nonrigid_displacement.csv',sl,delimiter=',',fmt='%.6f')
    (S/'segment_t0.txt').write_text('%.6f' % t0)
    # everything matlab/run_interp.m needs, so no MATLAB file has to be edited per run
    (S/'interp_params.json').write_text(json.dumps(dict(
        motion_fs=MOTION_FS,      # sample rate of nonrigid_displacement.csv
        batch_samples=128,        # interpolation batch length, in 30 kHz samples
        nchan_raw=NCH,            # channels in the FILE (sync included) - NOT the neural count
        m_tol_um=0.5,             # rebuild the kriging operator when the shift moves this far
        out_name='interp'+name,
    ), indent=1))
    print(f'{name}: t={t0:.0f}-{t1:.0f} s, {F1-F0} frames -> {S}')

if __name__=='__main__':
    main(float(sys.argv[1]),float(sys.argv[2]),sys.argv[3],
         sys.argv[4] if len(sys.argv)>4 else None)
