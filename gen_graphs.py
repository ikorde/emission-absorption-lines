import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


df_all = pd.read_pickle('df_all.pkl')


# dz1 = df_all.loc[df_all['Z'] <= 1]
# dz1_LYA = len(dz1[dz1['LYA','SNR'].notnull()])
# dz1_OII = len(dz1[dz1['OII','SNR'].notnull()])
# dz1_HeII = len(dz1[dz1['HeII','SNR'].notnull()])
# dz1_HBeta = len(dz1[dz1['HBeta','SNR'].notnull()])
# dz1_NII = len(dz1[dz1['NII','SNR'].notnull()])
# dz1_HAlpha = len(dz1[dz1['HAlpha','SNR'].notnull()])
# dz1_CIV = len(dz1[dz1['CIV','SNR'].notnull()])
# dz1_OIII = len(dz1[dz1['OIII','SNR'].notnull()])



# dz2 = df_all.loc[(df_all['Z'] > 1) & (df_all['Z'] <=2) ]
# dz2_LYA = len(dz2[dz2['LYA','SNR'].notnull()])
# dz2_OII = len(dz2[dz2['OII','SNR'].notnull()])
# dz2_HeII = len(dz2[dz2['HeII','SNR'].notnull()])
# dz2_HBeta = len(dz2[dz2['HBeta','SNR'].notnull()])
# dz2_NII = len(dz2[dz2['NII','SNR'].notnull()])
# dz2_HAlpha = len(dz2[dz2['HAlpha','SNR'].notnull()])
# dz2_CIV = len(dz2[dz2['CIV','SNR'].notnull()])
# dz2_OIII = len(dz2[dz2['OIII','SNR'].notnull()])

# dz3 = df_all.loc[(df_all['Z'] > 2) & (df_all['Z'] <=3) ]
# dz3_LYA = len(dz3[dz3['LYA','SNR'].notnull()])
# dz3_OII = len(dz3[dz3['OII','SNR'].notnull()])
# dz3_HeII = len(dz3[dz3['HeII','SNR'].notnull()])
# dz3_HBeta = len(dz3[dz3['HBeta','SNR'].notnull()])
# dz3_NII = len(dz3[dz3['NII','SNR'].notnull()])
# dz3_HAlpha = len(dz3[dz3['HAlpha','SNR'].notnull()])
# dz3_CIV = len(dz3[dz3['CIV','SNR'].notnull()])
# dz3_OIII = len(dz3[dz3['OIII','SNR'].notnull()])

# dz4 = df_all.loc[(df_all['Z'] > 3) & (df_all['Z'] <=4) ]
# dz4_LYA = len(dz4[dz4['LYA','SNR'].notnull()])
# dz4_OII = len(dz4[dz4['OII','SNR'].notnull()])
# dz4_HeII = len(dz4[dz4['HeII','SNR'].notnull()])
# dz4_HBeta = len(dz4[dz4['HBeta','SNR'].notnull()])
# dz4_NII = len(dz4[dz4['NII','SNR'].notnull()])
# dz4_HAlpha = len(dz4[dz4['HAlpha','SNR'].notnull()])
# dz4_CIV = len(dz4[dz4['CIV','SNR'].notnull()])
# dz4_OIII = len(dz4[dz4['OIII','SNR'].notnull()])

# dz5 = df_all.loc[(df_all['Z'] > 4) ]
# dz5_LYA = len(dz5[dz5['LYA','SNR'].notnull()])
# dz5_OII = len(dz5[dz5['OII','SNR'].notnull()])
# dz5_HeII = len(dz5[dz5['HeII','SNR'].notnull()])
# dz5_HBeta = len(dz5[dz5['HBeta','SNR'].notnull()])
# dz5_NII = len(dz5[dz5['NII','SNR'].notnull()])
# dz5_HAlpha = len(dz5[dz5['HAlpha','SNR'].notnull()])
# dz5_CIV = len(dz5[dz5['CIV','SNR'].notnull()])
# dz5_OIII = len(dz5[dz5['OIII','SNR'].notnull()])

# n_groups = 5
# d1 = (dz1_LYA, dz2_LYA, dz3_LYA, dz4_LYA,dz5_LYA)
# d2 = (dz1_HeII,dz2_HeII,dz3_HeII,dz4_HeII,dz5_HeII)
# d3 = (dz1_OII,dz2_OII,dz3_OII,dz4_OII,dz5_OII)
# d4 = (dz1_HBeta, dz2_HBeta, dz3_HBeta,dz4_HBeta,dz5_HBeta)
# d5 = (dz1_NII,dz2_NII,dz3_NII,dz4_NII,dz5_NII)
# d6 = (dz1_HAlpha, dz2_HAlpha,dz3_HAlpha,dz4_HAlpha,dz5_HAlpha)
# d7 = (dz1_CIV,dz2_CIV,dz3_CIV,dz4_CIV,dz5_CIV)
# d8 = (dz1_OIII,dz2_OIII,dz3_OIII,dz4_OIII,dz5_OIII)


# # create plot
# fig, ax = plt.subplots()
# index = np.arange(n_groups)
# bar_width = 0.1
# opacity = 0.8


# rects1 = plt.bar(index, d4, bar_width,alpha=opacity,label='HBeta')
# rects2 = plt.bar(index+ bar_width, d8, bar_width,alpha=opacity,label='OIII')
# rects3 = plt.bar(index+ bar_width*2, d3, bar_width,alpha=opacity,label='OII')
# rects4 = plt.bar(index+ bar_width*3, d5, bar_width,alpha=opacity,label='NII')
# rects5 = plt.bar(index+ bar_width*4, d6, bar_width,alpha=opacity,label='HAlpha')
# rects6 = plt.bar(index+ bar_width*5, d2, bar_width,alpha=opacity,label='HeII')
# rects7 = plt.bar(index+ bar_width*6, d7, bar_width,alpha=opacity,label='CIV')
# rects8 = plt.bar(index+ bar_width*7, d1, bar_width,alpha=opacity,label='LYA')


# plt.xlabel('Redshift')
# plt.ylabel('Frequency')
# plt.title('Lines at redshifts')
# plt.xticks(index + bar_width, ('0-1', '1-2', '2-3', '3-4', '4+'))
# plt.legend()

# plt.tight_layout()
# plt.show()


dz1 = df_all.loc[df_all['Z'] <= 1]
dz1_LYA = dz1[dz1['LYA','EWcenter'].notnull()]
dz1_OII = dz1[dz1['OII','EWcenter'].notnull()]
dz1_HeII = dz1[dz1['HeII','EWcenter'].notnull()]
dz1_HBeta = dz1[dz1['HBeta','EWcenter'].notnull()]
dz1_NII = dz1[dz1['NII','EWcenter'].notnull()]
dz1_HAlpha = dz1[dz1['HAlpha','EWcenter'].notnull()]
dz1_CIV = dz1[dz1['CIV','EWcenter'].notnull()]
dz1_OIII = dz1[dz1['OIII','EWcenter'].notnull()]
breakpoint()

dz2 = df_all.loc[(df_all['Z'] > 1) & (df_all['Z'] <=2) ]
dz2_LYA = dz2[dz2['LYA','EWcenter'].notnull()]
dz2_OII = dz2[dz2['OII','EWcenter'].notnull()]
dz2_HeII = dz2[dz2['HeII','EWcenter'].notnull()]
dz2_HBeta = dz2[dz2['HBeta','EWcenter'].notnull()]
dz2_NII = dz2[dz2['NII','EWcenter'].notnull()]
dz2_HAlpha = dz2[dz2['HAlpha','EWcenter'].notnull()]
dz2_CIV = dz2[dz2['CIV','EWcenter'].notnull()]
dz2_OIII = dz2[dz2['OIII','EWcenter'].notnull()]

dz3 = df_all.loc[(df_all['Z'] > 2) & (df_all['Z'] <=3) ]
dz3_LYA = dz3[dz3['LYA','EWcenter'].notnull()]
dz3_OII = dz3[dz3['OII','EWcenter'].notnull()]
dz3_HeII = dz3[dz3['HeII','EWcenter'].notnull()]
dz3_HBeta = dz3[dz3['HBeta','EWcenter'].notnull()]
dz3_NII = dz3[dz3['NII','EWcenter'].notnull()]
dz3_HAlpha = dz3[dz3['HAlpha','EWcenter'].notnull()]
dz3_CIV = dz3[dz3['CIV','EWcenter'].notnull()]
dz3_OIII = dz3[dz3['OIII','EWcenter'].notnull()]

dz4 = df_all.loc[(df_all['Z'] > 3) & (df_all['Z'] <=4) ]
dz4_LYA = dz4[dz4['LYA','EWcenter'].notnull()]
dz4_OII = dz4[dz4['OII','EWcenter'].notnull()]
dz4_HeII = dz4[dz4['HeII','EWcenter'].notnull()]
dz4_HBeta = dz4[dz4['HBeta','EWcenter'].notnull()]
dz4_NII = dz4[dz4['NII','EWcenter'].notnull()]
dz4_HAlpha = dz4[dz4['HAlpha','EWcenter'].notnull()]
dz4_CIV = dz4[dz4['CIV','EWcenter'].notnull()]
dz4_OIII = dz4[dz4['OIII','EWcenter'].notnull()]

dz5 = df_all.loc[(df_all['Z'] > 4) ]
dz5_LYA = dz5[dz5['LYA','EWcenter'].notnull()]
dz5_OII = dz5[dz5['OII','EWcenter'].notnull()]
dz5_HeII = dz5[dz5['HeII','EWcenter'].notnull()]
dz5_HBeta = dz5[dz5['HBeta','EWcenter'].notnull()]
dz5_NII = dz5[dz5['NII','EWcenter'].notnull()]
dz5_HAlpha = dz5[dz5['HAlpha','EWcenter'].notnull()]
dz5_CIV = dz5[dz5['CIV','EWcenter'].notnull()]
dz5_OIII = dz5[dz5['OIII','EWcenter'].notnull()]

n_groups = 5
d1 = (dz1_LYA, dz2_LYA, dz3_LYA, dz4_LYA,dz5_LYA)
d2 = (dz1_HeII,dz2_HeII,dz3_HeII,dz4_HeII,dz5_HeII)
d3 = (dz1_OII,dz2_OII,dz3_OII,dz4_OII,dz5_OII)
d4 = (dz1_HBeta, dz2_HBeta, dz3_HBeta,dz4_HBeta,dz5_HBeta)
d5 = (dz1_NII,dz2_NII,dz3_NII,dz4_NII,dz5_NII)
d6 = (dz1_HAlpha, dz2_HAlpha,dz3_HAlpha,dz4_HAlpha,dz5_HAlpha)
d7 = (dz1_CIV,dz2_CIV,dz3_CIV,dz4_CIV,dz5_CIV)
d8 = (dz1_OIII,dz2_OIII,dz3_OIII,dz4_OIII,dz5_OIII)


# create plot
fig, ax = plt.subplots()
index = np.arange(n_groups)
bar_width = 0.1
opacity = 0.8


rects1 = plt.bar(index, d4, bar_width,alpha=opacity,label='HBeta')
rects2 = plt.bar(index+ bar_width, d8, bar_width,alpha=opacity,label='OIII')
rects3 = plt.bar(index+ bar_width*2, d3, bar_width,alpha=opacity,label='OII')
rects4 = plt.bar(index+ bar_width*3, d5, bar_width,alpha=opacity,label='NII')
rects5 = plt.bar(index+ bar_width*4, d6, bar_width,alpha=opacity,label='HAlpha')
rects6 = plt.bar(index+ bar_width*5, d2, bar_width,alpha=opacity,label='HeII')
rects7 = plt.bar(index+ bar_width*6, d7, bar_width,alpha=opacity,label='CIV')
rects8 = plt.bar(index+ bar_width*7, d1, bar_width,alpha=opacity,label='LYA')


plt.xlabel('Redshift')
plt.ylabel('Frequency')
plt.title('Lines at redshifts')
plt.xticks(index + bar_width, ('0-1', '1-2', '2-3', '3-4', '4+'))
plt.legend()

plt.tight_layout()
plt.show()