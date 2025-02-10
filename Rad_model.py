import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
import matplotlib.gridspec as gridspec
from pylab import *
from sklearn.feature_selection import f_regression, mutual_info_regression
from scipy import stats
import seaborn as sns
import xarray as xr
import pdb
from Casicus import exp_rrtmg, daily_states, real_rrtmg, battery_exp
import sys

def plot_vert_rrtmg(revs,orgis,xlabel):
    linestyle = ['--','-','-.']
    titles = ['(a)','(b)','(c)']
    fig = plt.figure(figsize=(12,8))
    gs = gridspec.GridSpec(1,3, left=0.05, right=0.975, hspace=0.25, wspace=0.1, top=0.9, bottom=0.1)
    for i,rev in enumerate(revs):
        ax=subplot(gs[i])
        plt.plot(rev, pres_long, color = 'darkorchid', linestyle = linestyle[i], label = 'Reversals')
        plt.plot(orgis[i], pres_long, color = 'blue', linestyle = linestyle[i], label = 'Reversals')
        plt.ylim(1000,100)
        plt.xlabel(xlabel)
        plt.grid(linestyle=':')
        plt.title(titles[i], fontweight='bold', loc = 'left')
        plt.axvline(0,linestyle = '--', color = 'k', linewidth=0.75)
        if i == 0:
            plt.ylabel('Pressure')
        else:
            ax.yaxis.set_major_formatter(NullFormatter())
        if i == 1:
            plt.legend(frameon=False, ncol = 2, loc = 'upper center', bbox_to_anchor=(0.5,1.075))
            
def plot_sct_RRTMG(revs,orgs):
    xaxis = np.arange(0,3)
    titles = ['(a)','(b)','(c)','(d)']
    ylabel = ['LW Top Net Total (W m$^{-2}$)', 'LW Sfc Upw Total (W m$^{-2}$)', 
              'LWclr Top Net Total (W m$^{-2}$)', 'LWclr Sfc Upw Total (W m$^{-2}$)']

    fig = plt.figure(figsize=(8,8))
    gs = gridspec.GridSpec(2,2, left=0.05, right=0.975, hspace=0.25, wspace=0.4, top=0.9, bottom=0.1)

    for i,rev in enumerate(revs):
        ax=subplot(gs[i])
        plt.scatter(xaxis,rev, color = 'darkorchid', label = 'Reversals', marker = '*')
        plt.scatter(xaxis,orgs[i], color = 'blue', label = 'Organised', marker = '*')
        plt.xticks(xaxis, ['2N-9N','3N-10N','3S-4N'])
        plt.title(titles[i], fontweight = 'bold', loc = 'left', fontsize = 12)
        if i == 0:
            plt.legend(frameon=False, ncol = 2, loc = 'upper center', bbox_to_anchor=(1.2,1.2))
        plt.grid(linestyle = ':', linewidth = 0.5)
        plt.ylabel(ylabel[i])

### See Nomenclature above!!!
path = '/home/tompkins-archive/acasallas/RRTMG_data/'
# '2-9' and '135-145'
mt2,mlo2,mhi2,mne2, st2,slo2,shi2,sne2, ctt2,ctlo2,cthi2,ctne2 = daily_states('2-9','135-145',0.5, True)
# '3-10' and '147-157'
mt3,mlo3,mhi3,mne3, st3,slo3,shi3,sne3, ctt3,ctlo3,cthi3,ctne3 = daily_states('3-10','147-157',0.5, True)
# '3S-4N' and '156-166'
mt3s,mlo3s,mhi3s,mne3s, st3s,slo3s,shi3s,sne3s, ctt3s,ctlo3s,cthi3s,ctne3s = daily_states('3S-4N','156-166', 0.5, True)

pres_long = [100,  125,  150,  175,  200,  225,  250,  300,  350,  400, 450,
             500,  550,  600,  650,  700,  750,  775,  800,  825,
             850,  875,  900,  925,  950,  975, 1000]
#### Air temperature, specifict humidity and cloud cover
temp2 = mt2.columns[178:205]
sphum = mt2.columns[205:232]
cover = mt2.columns[147:174]

### Battery of experiments
print('-----------------------------------------------------------------------------')
print('------------- Starting to calculate the battery of experiments! -------------')
print('-----------------------------------------------------------------------------')

rads = ['OLR','OLRclr','SFC','SFCclr']

jump = 1

if sys.argv[2] == 'exp':
    varis = [sys.argv[4]]
    sta = sys.argv[5]
    end = sys.argv[6]
else:
    sta = 0
    end = 175

days = np.arange(int(sta),int(end)+0.1,jump)
### System arguments, 1:area; 2:real or exp?; 3:rev/org, r_in_o/o_in_r; 4:SST/Hum/Temp/Clouds 
############# 2-9
if sys.argv[1] == '2-9':
    if sys.argv[2] == 'real':
        if sys.argv[3] == 'rev':
            exp_inp = {'Real': [np.array(mlo2['SST']), np.array(mlo2[sphum]),
                                np.array(mlo2[temp2]), np.array(mlo2[cover]),
                                np.array(mlo2['ciwp']), np.array(mlo2['clwp'])]
                      }
            df = real_rrtmg(days,exp_inp,'2-9','Real_rev',pres_long)
        elif sys.argv[3] == 'org':
            exp_inp = {'Real': [np.array(mhi2['SST']), np.array(mhi2[sphum]),
                                np.array(mhi2[temp2]), np.array(mhi2[cover]),
                                np.array(mhi2['ciwp']), np.array(mhi2['clwp'])]
                      }
            df = real_rrtmg(days,exp_inp,'2-9','Real_org',pres_long)
    elif sys.argv[2] == 'exp':
        if sys.argv[3] == 'r_in_o':
            exp_inp_r_in_o = {'SST': [np.array(mlo2['SST']), np.array(mhi2[sphum]),
                                      np.array(mhi2[temp2]), np.array(mhi2[cover]),
                                      np.array(mhi2['ciwp']), np.array(mhi2['clwp'])],
                              'Hum': [np.array(mhi2['SST']),np.array(mlo2[sphum]),
                                      np.array(mhi2[temp2]), np.array(mhi2[cover]),
                                      np.array(mhi2['ciwp']), np.array(mhi2['clwp'])],
                              'Temp': [np.array(mhi2['SST']),np.array(mhi2[sphum]),
                                       np.array(mlo2[temp2]), np.array(mhi2[cover]),
                                       np.array(mhi2['ciwp']), np.array(mhi2['clwp'])],
                              'Clouds': [np.array(mhi2['SST']),np.array(mhi2[sphum]),
                                         np.array(mhi2[temp2]), np.array(mlo2[cover]),
                                         np.array(mlo2['ciwp']), np.array(mlo2['clwp'])]
                             }
            df = battery_exp(varis,rads,days,exp_inp_r_in_o,'2-9','r_in_o',pres_long,sta,end)    
        elif sys.argv[3] == 'o_in_r':
            exp_inp_o_in_r = {'SST': [np.array(mhi2['SST']), np.array(mlo2[sphum]),
                                      np.array(mlo2[temp2]), np.array(mlo2[cover]),
                                      np.array(mlo2['ciwp']), np.array(mlo2['clwp'])],
                              'Hum': [np.array(mlo2['SST']),np.array(mhi2[sphum]),
                                      np.array(mlo2[temp2]), np.array(mlo2[cover]),
                                      np.array(mlo2['ciwp']), np.array(mlo2['clwp'])],
                              'Temp': [np.array(mlo2['SST']),np.array(mlo2[sphum]),
                                       np.array(mhi2[temp2]), np.array(mlo2[cover]),
                                       np.array(mlo2['ciwp']), np.array(mlo2['clwp'])],
                              'Clouds': [np.array(mlo2['SST']),np.array(mlo2[sphum]),
                                         np.array(mlo2[temp2]), np.array(mhi2[cover]),
                                         np.array(mhi2['ciwp']), np.array(mhi2['clwp'])]
                               }
            df = battery_exp(varis,rads,days,exp_inp_o_in_r,'2-9','o_in_r',pres_long,sta,end)
############# 3-10  
elif sys.argv[1] == '3-10':
   if sys.argv[2] == 'real':
        if sys.argv[3] == 'rev':
            exp_inp = {'Real': [np.array(mlo3['SST']), np.array(mlo3[sphum]),
                                np.array(mlo3[temp2]), np.array(mlo3[cover]),
                                np.array(mlo3['ciwp']), np.array(mlo3['clwp'])]
                      }
            df = real_rrtmg(days,exp_inp,'3-10','Real_rev',pres_long)
        elif sys.argv[3] == 'org':
            exp_inp = {'Real': [np.array(mhi3['SST']), np.array(mhi3[sphum]),
                                np.array(mhi3[temp2]), np.array(mhi3[cover]),
                                np.array(mhi3['ciwp']), np.array(mhi3['clwp'])]
                      }
            df = real_rrtmg(days,exp_inp,'3-10','Real_org',pres_long)
   elif sys.argv[2] == 'exp':
       if sys.argv[3] == 'r_in_o':
           exp_inp_r_in_o = {'SST': [np.array(mlo3['SST']), np.array(mhi3[sphum]),
                                     np.array(mhi3[temp2]), np.array(mhi3[cover]),
                                     np.array(mhi3['ciwp']), np.array(mhi3['clwp'])],
                             'Hum': [np.array(mhi3['SST']),np.array(mlo3[sphum]),
                                     np.array(mhi3[temp2]), np.array(mhi3[cover]),
                                     np.array(mhi3['ciwp']), np.array(mhi3['clwp'])],
                             'Temp': [np.array(mhi3['SST']),np.array(mhi3[sphum]),
                                      np.array(mlo3[temp2]), np.array(mhi3[cover]),
                                      np.array(mhi3['ciwp']), np.array(mhi3['clwp'])],
                             'Clouds': [np.array(mhi3['SST']),np.array(mhi3[sphum]),
                                        np.array(mhi3[temp2]), np.array(mlo3[cover]),
                                        np.array(mlo3['ciwp']), np.array(mlo3['clwp'])]
                            }
           df = battery_exp(varis,rads,days,exp_inp_r_in_o,'3-10','r_in_o',pres_long,sta,end)
       elif sys.argv[3] == 'o_in_r':
           exp_inp_o_in_r = {'SST': [np.array(mhi3['SST']), np.array(mlo3[sphum]),
                                     np.array(mlo3[temp2]), np.array(mlo3[cover]),
                                     np.array(mlo3['ciwp']), np.array(mlo3['clwp'])],
                             'Hum': [np.array(mlo3['SST']),np.array(mhi3[sphum]),
                                     np.array(mlo3[temp2]), np.array(mlo3[cover]),
                                     np.array(mlo3['ciwp']), np.array(mlo3['clwp'])],
                             'Temp': [np.array(mlo3['SST']),np.array(mlo3[sphum]),
                                      np.array(mhi3[temp2]), np.array(mlo3[cover]),
                                      np.array(mlo3['ciwp']), np.array(mlo3['clwp'])],
                             'Clouds': [np.array(mlo3['SST']),np.array(mlo3[sphum]),
                                        np.array(mlo3[temp2]), np.array(mhi3[cover]),
                                        np.array(mhi3['ciwp']), np.array(mhi3['clwp'])]
                              }
           df = battery_exp(varis,rads,days,exp_inp_o_in_r,'3-10','o_in_r',pres_long,sta,end)
############# 3S-4N
elif sys.argv[1] == '3S-4N':
    if sys.argv[2] == 'real':
        if sys.argv[3] == 'rev':
            exp_inp = {'Real': [np.array(mlo3s['SST']), np.array(mlo3s[sphum]),
                                np.array(mlo3s[temp2]), np.array(mlo3s[cover]),
                                np.array(mlo3s['ciwp']), np.array(mlo3s['clwp'])]
                      }
            df = real_rrtmg(days,exp_inp,'3S-4N','Real_rev',pres_long)
        elif sys.argv[3] == 'org':
            exp_inp = {'Real': [np.array(mhi3s['SST']), np.array(mhi3s[sphum]),
                                np.array(mhi3s[temp2]), np.array(mhi3s[cover]),
                                np.array(mhi3s['ciwp']), np.array(mhi3s['clwp'])]
                      }
            df = real_rrtmg(days,exp_inp,'3S-4N','Real_org',pres_long)
    elif sys.argv[2] == 'exp':
        if sys.argv[3] == 'r_in_o':
            exp_inp_r_in_o = {'SST': [np.array(mlo3s['SST']), np.array(mhi3s[sphum]),
                                      np.array(mhi3s[temp2]), np.array(mhi3s[cover]),
                                      np.array(mhi3s['ciwp']), np.array(mhi3s['clwp'])],
                              'Hum': [np.array(mhi3s['SST']),np.array(mlo3s[sphum]),
                                      np.array(mhi3s[temp2]), np.array(mhi3s[cover]),
                                      np.array(mhi3s['ciwp']), np.array(mhi3s['clwp'])],
                              'Temp': [np.array(mhi3s['SST']),np.array(mhi3s[sphum]),
                                       np.array(mlo3s[temp2]), np.array(mhi3s[cover]),
                                       np.array(mhi3s['ciwp']), np.array(mhi3s['clwp'])],
                              'Clouds': [np.array(mhi3s['SST']),np.array(mhi3s[sphum]),
                                         np.array(mhi3s[temp2]), np.array(mlo3s[cover]),
                                         np.array(mlo3s['ciwp']), np.array(mlo3s['clwp'])]
                             }
            df = battery_exp(varis,rads,days,exp_inp_r_in_o,'3S-4N','r_in_o',pres_long,sta,end)
        elif sys.argv[3] == 'o_in_r':
            exp_inp_o_in_r = {'SST': [np.array(mhi3s['SST']), np.array(mlo3s[sphum]),
                                      np.array(mlo3s[temp2]), np.array(mlo3s[cover]),
                                      np.array(mlo3s['ciwp']), np.array(mlo3s['clwp'])],
                              'Hum': [np.array(mlo3s['SST']),np.array(mhi3s[sphum]),
                                      np.array(mlo3s[temp2]), np.array(mlo3s[cover]),
                                      np.array(mlo3s['ciwp']), np.array(mlo3s['clwp'])],
                              'Temp': [np.array(mlo3s['SST']),np.array(mlo3s[sphum]),
                                       np.array(mhi3s[temp2]), np.array(mlo3s[cover]),
                                       np.array(mlo3s['ciwp']), np.array(mlo3s['clwp'])],
                              'Clouds': [np.array(mlo3s['SST']),np.array(mlo3s[sphum]),
                                         np.array(mlo3s[temp2]), np.array(mhi3s[cover]),
                                         np.array(mhi3s['ciwp']), np.array(mhi3s['clwp'])]
                               }
            df = battery_exp(varis,rads,days,exp_inp_o_in_r,'3S-4N','o_in_r',pres_long,sta,end)

