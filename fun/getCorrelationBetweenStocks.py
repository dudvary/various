# Converted from juypter notebook script
# Calculate correlation between stocks. Here the following methodology is used:
# Correlate daily change (X[i]-X[i-1])/X[i-1] (this seems to be the most reasonable, and is what the internet suggests)
# Historic stock performance downloaded as wkn_[wkn]_historic.csv with semicolon delimiter from https://www.ariva.de/
# All .csv in one folder 

# Other options, not used/advised
# Option 2: Correlate daily change X[i]/X[i-1]
# Option 3: Correlate time series
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import date, timedelta
import glob, os

foldername = ''
prefix = ''

# Iterate over all csv files and store table in dictionary
sharesDict = {}
filelist = glob.glob(foldername + '/*.csv')
col_list = ['Datum','Schlusskurs'] # No difference between Schlusskurs or Hoechstkurs 
for file in filelist:
    ID = file.split("_")
    df = pd.read_csv(file,sep=';',index_col='Datum',usecols=col_list, decimal=",", parse_dates=True)
    df = df.rename(columns={"Schlusskurs": "X"})
    
    change_tmp = np.empty((len(df.X),))
    change_tmp[0] = 0
    for i in range(1,len(df.X)):
        change_tmp[i] = (df.X[i]-df.X[i-1])/df.X[i-1]*100

    df["change"] = change_tmp
    
    sharesDict.update({ID[2]:df})

# COMPUTE CORRELATION (by day)
def computeCorrelation(df1,df2):
    
    # Find date range of both data sets
    startDate = min([min(df1.index), min(df2.index)])
    endDate = max([max(df1.index), max(df2.index)])
    dates = pd.date_range(start=startDate, end=endDate)
    #print ('from %s to %s' % (startDate.strftime("%Y-%m-%d"),endDate.strftime("%Y-%m-%d")))

    # Error handling
    if len(dates)<1:
        print('WARNING! Only ',len(dates),' found!')
        R = np.empty((3,3,))
        R[:] = np.nan
        return R,np.nan
    
    x1 = []
    x2 = []
    
    # Iterate over all dates and check for overlap
    for d in dates:
        b1 = (d in df1.index)
        b2 = (d in df2.index)
        if b1 & b2:
            x1.append(df1.change[d])
            x2.append(df2.change[d])

    x1 = np.array(x1)
    x2 = np.array(x2)
    n = len(x2)

    # Error handling
    # Check if at least 2 overlapping values (necessary for correlation)
    if n<2:
        print('WARNING! Cannot compute correlation. Less than ',n, ' overlapping values!')
        R = np.empty((3,3,))
        R[:] = np.nan
        return R,n
    
    # Calculate correlation coefficient
    R = np.corrcoef(x1,x2)
    return R, n

# Compute change in percentage
def computeChange(df):
    change = (df.X[0]-df.X[-1])/df.X[-1]*100
    return change

R = np.empty((len(sharesDict),len(sharesDict),))
change = np.empty((len(sharesDict),))
x = 0
labels = []
Rlist = []
IDpairs = []

for key1 in sharesDict:
    labels.append(key1)
    change[x] = computeChange(sharesDict[key1])
    print(' %s (n=%d) %.2f%%' % (key1,sharesDict[key1].size,change[x]))
    y = 0
    for key2 in sharesDict:
        if x==y:
            R[x,y] = 1.0
        elif x<y:
            rtmp, n = computeCorrelation(sharesDict[key1],sharesDict[key2])
            #print('%s-%s: R=%.2f n=%d' % (key1,key2,rtmp[1,0],n))
            R[x,y] = rtmp[1,0]
            R[y,x] = rtmp[1,0]
            Rlist.append(rtmp[1,0])
            IDpairs.append((x,y))
        y = y+1
    x = x+1

mask = np.triu(np.ones_like(R, dtype=bool))
sns.set_theme()
f, ax = plt.subplots(figsize=R.shape)
heatmap = sns.heatmap(R, annot=True, linewidths=.5, ax=ax, cmap="vlag", 
                      xticklabels=labels, yticklabels=labels,center=0,
                     cbar_kws={'label': 'correlation R','shrink':0.5}, 
                     mask=mask, square=True)
heatmap.set_xticklabels(heatmap.get_xticklabels(), rotation=45);
heatmap.set_yticklabels(heatmap.get_yticklabels(), rotation=0);

plt.savefig(foldername + '_CorrMatrix' + prefix + '.png',bbox_inches='tight')


plt.barh(y=range(0,len(change)),width=change,tick_label=labels);
plt.xlabel('change(%)');
plt.gca().invert_yaxis()
plt.savefig(foldername + '_Change' + prefix + '.png',bbox_inches='tight')

for key in sharesDict:
    x = sharesDict[key].X
    x = x/x[-1]
    plt.plot(x,label=key);
plt.legend(bbox_to_anchor=(1.05,1),loc='upper left',borderaxespad=0);
plt.savefig(foldername + '_Chart' + prefix + '.png',bbox_inches='tight')

# Write csv with results
f = open(foldername + '_table' + prefix + '.csv','w')

str = 'WKN,'
for key1 in sharesDict:
    str = str + key1 + ','
f.write(str + ',Change(%),\n')

x = 0
for key1 in sharesDict:
    y = 0
    str = key1 + ','
    for key2 in sharesDict:
        str = str + ('%.2f,' % R[x,y])   
        y = y+1
    str = str + (',%.3f,' % change[x])
    f.write(str + '\n')
    x = x+1
f.close()

# Write csv with sorted Correlations
idxSort = sorted(range(len(Rlist)), key=lambda k: Rlist[k],reverse=True)
f = open(foldername + '_sortedCorrelations' + prefix + '.csv','w')

f.write('WKN1,WKN2,R,Change1(%),Change2(%),\n')
for i in idxSort:
    tmp = IDpairs[i]
    label1 = labels[tmp[0]]
    label2 = labels[tmp[1]]
    c1 = change[tmp[0]]
    c2 = change[tmp[1]]
    str = ('%s,%s,%.2f,%.2f,%.2f,' % (label1,label2,Rlist[i],c1,c2))
    f.write(str + '\n')
f.close()