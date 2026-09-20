#!/usr/bin/env python3
"""
build_precip_onestep.py -- one-step precipitation reconstruction (drip -> P/T -> P) for Fig 6.

Replaces legacy/precip_recon/build_onestep_dense.py. Reads everything from the workbook and
external/: the age-propagated drip posterior (Fig5_driprate_AP; mu already anchored to the
2004-2023 annual-baseflow mean, so NO rescaling is applied), the calibration table
(external/calibration_onestep.csv: annual baseflow, P, T -> P/T_h), and Wang et al. (2018)
RAN15-MAAT (external/T_recon_Wang_et_al.xlsx). Monte Carlo (N = 10,000 per horizon):
lognormal drip from med/pc25/pc75; normal a_h, b_h (regression SEs) and residual; normal T (SD 2.6).
Writes sheets Fig6_P_recon_dense (every 5th drip-grid year) and Fig6_P_recon_LR (73 Wang horizons).
"""
import os, numpy as np, pandas as pd, openpyxl
from scipy import stats
from scipy.stats import norm
HERE=os.path.dirname(os.path.abspath(__file__)); EXT=os.path.join(HERE,'external'); WB=os.path.join(HERE,'HS4_SourceData.xlsx')
N=10000; rng=np.random.default_rng(0)
cal=pd.read_csv(os.path.join(EXT,'calibration_onestep.csv'))
a_h,b_h,r,pv,a_se=stats.linregress(cal.baseflow,cal.PT_h); resid=(cal.PT_h-(a_h*cal.baseflow+b_h)).std(ddof=2); b_se=a_se*np.sqrt((cal.baseflow**2).mean())
print(f'calibration: a_h={a_h:.5f}+/-{a_se:.5f} b_h={b_h:+.5f}+/-{b_se:.5f} resid={resid:.4f} R2={r**2:.3f} n={len(cal)}')
raw=pd.read_excel(WB,sheet_name='Fig5_driprate_AP',header=None); h=[i for i in range(10) if str(raw.iloc[i,0]).strip()=='age_calBP'][0]
ap=pd.read_excel(WB,sheet_name='Fig5_driprate_AP',header=h).apply(pd.to_numeric,errors='coerce').dropna()
T=pd.read_excel(os.path.join(EXT,'T_recon_Wang_et_al.xlsx'))[['Age (yr, BP)','RAN15-MAAT(°C)']].dropna(); T.columns=['age','T']
def mc(rows):
    out=[]
    for _,row in rows.iterrows():
        sg=(np.log(row.DR_pc75)-np.log(row.DR_pc25))/(2*norm.ppf(0.75)); drip=np.exp(rng.normal(np.log(row.DR_med),sg,N))
        PT=np.maximum(rng.normal(a_h,a_se,N)*drip+rng.normal(b_h,b_se,N)+rng.normal(0,resid,N),0)
        P=PT*rng.normal(row['T'],2.6,N)*365.25
        out.append({'age':row.age_calBP,'P_med':np.median(P),'P_pc25':np.percentile(P,25),'P_pc75':np.percentile(P,75)})
    return pd.DataFrame(out).sort_values('age')
d=ap.sort_values('age_calBP').iloc[::5].copy(); d['T']=np.interp(d.age_calBP,T.age,T['T']); d=d[(d.age_calBP>=T.age.min())&(d.age_calBP<=T.age.max())&(d.DR_pc25<d.DR_pc75)]
dense=mc(d)
lr_rows=[]
for _,t in T.iterrows():
    i=(ap.age_calBP-t.age).abs().idxmin(); lr_rows.append(dict(age_calBP=t.age,DR_med=ap.loc[i,'DR_med'],DR_pc25=ap.loc[i,'DR_pc25'],DR_pc75=ap.loc[i,'DR_pc75'],T=t['T']))
lr=mc(pd.DataFrame(lr_rows))
early=dense[(dense.age>8000)&(dense.age<9500)]; late=dense[(dense.age>500)&(dense.age<2000)]; pk=early.loc[early.P_med.idxmax()]
print(f'dense ({len(dense)} horizons): peak {pk.P_med:.0f} ({pk.P_pc25:.0f}-{pk.P_pc75:.0f}) at {pk.age:.0f} BP; late {late.P_med.median():.0f}; decline {100*(pk.P_med-late.P_med.median())/pk.P_med:.0f}%')
wb=openpyxl.load_workbook(WB)
for name,df,note in [('Fig6_P_recon_dense',dense,'dense, every 5th drip-grid year'),('Fig6_P_recon_LR',lr,'73 Wang-T horizons')]:
    ws=wb[name]; ws['A2']=f'one-step P reconstruction (drip->P/T->P), sigma = pi/sqrt(6); {note}. Regenerated 2026-09-21 by build_precip_onestep.py from Fig5_driprate_AP (568-point input, 20 second-laboratory samples excluded; mu anchored to 16.66 drips/min; no rescaling), calibration_onestep.csv (a_h {a_h:.5f}, b_h {b_h:+.5f}, resid {resid:.4f}) and Wang et al. 2018 T; MC N=10,000.'
    for rr in range(7,ws.max_row+1):
        for c in range(1,5): ws.cell(row=rr,column=c,value=None)
    for i,row in enumerate(df[['age','P_med','P_pc25','P_pc75']].itertuples(index=False),7):
        for j,v in enumerate(row,1): ws.cell(row=i,column=j,value=float(v))
wb.save(WB); print('workbook sheets written')
