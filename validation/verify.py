"""Independent numerical verification of the documented QP (not MATLAB)."""
from pathlib import Path
import numpy as np
import scipy
from scipy.optimize import minimize, LinearConstraint, nnls
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import csv, sys, platform
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'/'independent'; OUT.mkdir(parents=True,exist_ok=True)
N,M,ds=121,120,.5
s=np.arange(N)*ds; stops=np.array([0,40,80,120]); weights=[0,1,10,100]; w=.08
B=np.diff(np.eye(N),axis=0)/(2*ds); D=np.diff(np.eye(M),axis=0)
Aeq=np.vstack([np.hstack([B,-np.eye(M)]),np.hstack([np.eye(N)[stops],np.zeros((4,M))])])
lb=np.r_[np.zeros(N),np.full(M,-2.)]; ub=np.r_[np.full(N,36.),np.full(M,1.5)]
# Eliminate exactly the dynamics and mandatory stops; alpha = B beta.
free=np.setdiff1d(np.arange(N),stops); P=np.eye(N)[:,free]; T=np.vstack([P,B@P]); G=np.vstack([T,-T]); h=np.r_[ub,-lb]
# Omit zero rows from eliminated stop constraints.
keep=np.linalg.norm(G,axis=1)>0; G=G[keep]; h=h[keep]
rows=[]; profiles=[]; active=[]
def save_csv(name,header,rows):
    with (OUT/name).open('w',newline='') as f:
        writer=csv.writer(f);writer.writerow(header);writer.writerows(rows)
for ws in weights:
    H=np.zeros((N+M,N+M)); H[:N,:N]=2*w*np.eye(N); H[N:,N:]=2*ws*D.T@D
    f=np.r_[np.full(N,-2*w*25),np.zeros(M)]; Q=T.T@H@T; c=T.T@f
    # Scale objective for SLSQP only, preserving the minimizer.
    scale=max(1,ws)
    sol=minimize(lambda z:(.5*z@Q@z+c@z)/scale,np.zeros(len(free)),jac=lambda z:(Q@z+c)/scale,constraints=[LinearConstraint(G,-np.inf,h)],method='SLSQP',options={'ftol':1e-12,'maxiter':3000})
    # Polish the candidate active set by solving the equality-constrained QP.
    ids=np.flatnonzero(h-G@sol.x<1e-5); A=G[ids]
    K=np.block([[Q,A.T],[A,np.zeros((len(ids),len(ids)))]])
    polished=np.linalg.lstsq(K,np.r_[-c,h[ids]],rcond=None)[0][:len(free)]
    x=T@polished
    # Independent full-space KKT certificate: solve for nonnegative active
    # bound multipliers in the nullspace, then recover free equality multipliers.
    low=np.flatnonzero(x-lb<1e-7); up=np.flatnonzero(ub-x<1e-7)
    E=np.eye(N+M); C=np.column_stack([-E[:,low],E[:,up]])
    grad=H@x+f; mu,_=nnls(T.T@C,-T.T@grad,maxiter=10000)
    nu=np.linalg.lstsq(Aeq.T,-grad-C@mu,rcond=None)[0]
    lamL=np.zeros(N+M);lamU=np.zeros(N+M);lamL[low]=mu[:len(low)];lamU[up]=mu[len(low):]
    eq=np.max(np.abs(Aeq@x)); bound=max(0,np.max(lb-x),np.max(x-ub)); stat=np.max(np.abs(grad+Aeq.T@nu-lamL+lamU)); dual=max(0,-min(lamL.min(),lamU.min())); comp=max(np.max(np.abs(lamL*(x-lb))),np.max(np.abs(lamU*(ub-x))))
    assert max(eq,bound,stat,dual,comp)<1e-7, (ws,eq,bound,stat,dual,comp,sol.message)
    beta=x[:N]; a=x[N:]; v=np.sqrt(np.maximum(beta,0)); tr=np.sum((beta-25)**2); sm=np.sum(np.diff(a)**2); acc=np.abs(a-1.5)<1e-4; brake=np.abs(a+2)<1e-4
    rows.append([ws,tr,w*tr,sm,ws*sm,w*tr+ws*sm,acc.sum(),brake.sum(),acc.sum()+brake.sum(),eq,bound,stat,dual,comp,sol.success,sol.message])
    profiles.append((v,a)); active.append(acc.astype(int)-brake.astype(int))
    save_csv(f'ws_{ws}_nodes.csv',['node','s_m','beta','speed_mps'],zip(range(N),s,beta,v));save_csv(f'ws_{ws}_segments.csv',['segment','s_start_m','s_end_m','acceleration_mps2','active_upper','active_lower'],zip(range(M),s[:-1],s[1:],a,acc.astype(int),brake.astype(int)))
    save_csv(f'ws_{ws}_multipliers.csv',['variable_index','lambda_lower','lambda_upper'],zip(range(N+M),lamL,lamU));save_csv(f'ws_{ws}_equality_multipliers.csv',['row_index','nu'],enumerate(nu))
    print(ws, rows[-1],flush=True)
header=['ws','raw_tracking','weighted_tracking','raw_smoothness','weighted_smoothness','J','active_acceleration','active_braking','active_total','equality_inf','bound_violation','stationarity_inf','dual_violation','complementarity_inf','slsqp_success','slsqp_message'];save_csv('summary.csv',header,rows)
assert Aeq.shape==(124,241) and np.linalg.matrix_rank(Aeq)==124 and D.shape==(119,120)
(OUT/'environment.txt').write_text(f'Python {sys.version}\nNumPy {np.__version__}\nSciPy {scipy.__version__}\nPlatform {platform.platform()}\nSolver: scipy SLSQP + equality QP polish + full-space KKT certificate\nAeq rank: {np.linalg.matrix_rank(Aeq)}\n')
fig,axs=plt.subplots(2,1,figsize=(9,6),sharex=True,layout='constrained')
for ws,(v,a) in zip(weights,profiles):
    axs[0].plot(s,v,label=f'ws = {ws}');axs[1].stairs(a,s,baseline=None,label=f'ws = {ws}')
axs[0].axhline(5,ls=':',color='gray');axs[1].axhline(1.5,ls=':',color='gray');axs[1].axhline(-2,ls=':',color='gray')
for ax in axs:
    for st in s[stops]:ax.axvline(st,ls=':',color='gray',alpha=.5)
    ax.grid(alpha=.2)
axs[0].set(ylabel='Speed (m/s)',title='Independent reconstruction: speed and acceleration');axs[1].set(xlabel='Path distance (m)',ylabel='Acceleration (m/s²)');axs[0].legend(ncol=4);fig.savefig(OUT/'profiles.png',dpi=180);plt.close(fig)
from matplotlib.colors import ListedColormap,BoundaryNorm
from matplotlib.patches import Patch
fig,axs=plt.subplots(1,2,figsize=(10,3.8),layout='constrained',gridspec_kw={'width_ratios':[2,1]}); colors=['#d65b50','#e6ebf0','#168e7c'];axs[0].imshow(active,extent=[0,60,3.5,-.5],aspect='auto',cmap=ListedColormap(colors),norm=BoundaryNorm([-1.5,-.5,.5,1.5],3),interpolation='nearest');axs[0].set(yticks=range(4),yticklabels=weights,ylabel='Smoothing weight',xlabel='Path distance (m)',title='Active bounds by segment');axs[0].legend(handles=[Patch(color=colors[2],label='Acceleration'),Patch(color=colors[0],label='Braking'),Patch(color=colors[1],label='Inactive')],loc='upper center',bbox_to_anchor=(.5,-.2),ncol=3,fontsize=8)
i=np.arange(4);axs[1].bar(i-.18,[r[6] for r in rows],.36,label='Acceleration',color=colors[2]);axs[1].bar(i+.18,[r[7] for r in rows],.36,label='Braking',color=colors[0]);axs[1].set(xticks=i,xticklabels=weights,xlabel='Smoothing weight',ylabel='Active segments',title='Active counts');fig.savefig(OUT/'active_set.png',dpi=180);plt.close(fig)
fig,axs=plt.subplots(1,2,figsize=(9,3.8),layout='constrained');wt=[r[2] for r in rows];axs[0].bar(i,wt,label='Weighted tracking');axs[0].bar(i,[r[4] for r in rows],bottom=wt,label='Weighted smoothing');axs[0].legend();axs[0].set(ylabel='Objective contribution',title='Weighted objective');axs[1].bar(i,[r[3] for r in rows]);axs[1].set(ylabel='Raw smoothness',title='Acceleration differences');
for ax in axs: ax.set(xticks=i,xticklabels=weights,xlabel='Smoothing weight')
fig.savefig(OUT/'objective.png',dpi=180);plt.close(fig)
# Preserve the original table verbatim and compare only explicitly reported endpoints.
with (ROOT/'reference'/'original_report_table3.csv').open() as f:
    original=list(csv.DictReader(f))
comparison=[]
for entry in original:
    r=next(r for r in rows if r[0]==int(entry['ws']))
    old_t=float(entry['reported_raw_tracking_as_labeled']);old_s=float(entry['reported_raw_smoothness']);old_a=int(entry['reported_active_total'])
    comparison.append([r[0],old_t,r[1],r[2],r[2]-old_t,old_s,r[3],r[3]-old_s,old_a,r[8],r[8]-old_a])
save_csv('comparison_to_report.csv',['ws','reported_raw_tracking_as_labeled','reconstructed_raw_tracking','reconstructed_weighted_tracking','weighted_minus_reported','reported_smoothness','reconstructed_smoothness','smoothness_delta','reported_active','reconstructed_active','active_delta'],comparison)
