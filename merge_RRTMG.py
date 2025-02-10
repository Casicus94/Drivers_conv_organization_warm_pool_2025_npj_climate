import pandas as pd
import numpy as np
import pdb 
import os

areas = ['2-9','3-10','3S-4N']
tipos = ['r_in_o','o_in_r']
varis = ['SST','Hum','Temp','Clouds']
stas = ['0'] #['0','50','100','150']
end = ['90'] #['50','100','150','175']

path = '/home/tompkins-archive/acasallas/RRTMG_data/'

### First lets append all the files per variable
for area in areas:
    for tipo in tipos:
        data = {}
        for var in varis:
            df = pd.DataFrame()
            for i,sta in enumerate(stas):
                filename = path+'RRTMG_exp_'+tipo+'_'+area+'_'+var+'_'+sta+'-'+end[i]+'.csv'
                df_tmp = pd.read_csv(filename); df_tmp.drop('Unnamed: 0',axis=1, inplace=True)
                df = df.append(df_tmp, ignore_index=True)
                os.system('rm '+filename)
            ### Second save them in a dictionary
            data[area+'_'+var+'_'+tipo+'_all'] = df
        ### Third concat all the variables in 1 dataframe
        area_tot = pd.concat([data[area+'_SST_'+tipo+'_all'], data[area+'_Hum_'+tipo+'_all'],data[area+'_Temp_'+tipo+'_all'],data[area+'_Clouds_'+tipo+'_all']],axis = 1)
        ### Fourth put the dataframe in a file
        area_tot.to_csv(path+area+'_'+tipo+'_all.csv', index=False)
        
