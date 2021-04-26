# LOAD all packages
# see https://gist.github.com/josef-pkt/c932904296270d75366a24ee92a4eb2f
# https://www.statsmodels.org/stable/generated/statsmodels.discrete.count_model.ZeroInflatedNegativeBinomialP.html
# PYTHON 3 script
import numpy as np
import statsmodels.api as sm
import pandas as pd
import statsmodels.discrete._diagnostics_count as dia
from pandas import DataFrame

# Define fitting Function
def fitZINB(preCellType,postCellType):
    
    # Get data from file
    filename = "data_dense_model\%s_%s.csv" % (preCellType,postCellType)
    df = pd.read_csv(filename,header=None,names=["data"])
    
    # Prepare data for fitting
    X = df.data
    nobs = len(X)
    exog = np.ones(nobs)
    freq = np.bincount(X) / nobs
    binValue = list(range(0,len(freq)))
    
    # Fit Data
    mod_ZINB = sm.ZeroInflatedNegativeBinomialP(X, exog)
    res_ZINB = mod_ZINB.fit(disp=False)
    
    # Get fitting results
    probs_zinb = res_ZINB.predict(which='prob')
    probsm_zinb = probs_zinb.mean(0)
    
    # Export freq and probsm_zinb
    values = {'x': freq,
                'xFit': probsm_zinb}
    outputDF = DataFrame(values, columns= ['x', 'xFit'])
    outputfilename = "fit_dense_model\%s_%s_ZINB.csv" % (preCellType,postCellType)
    export_csv = outputDF.to_csv (outputfilename,index=None,header=True)
    
    # Export fit results
    X = res_ZINB.summary().as_csv()
    outputfilenameFit = "fit_dense_model\%s_%s_ZINB_FitResults.csv" % (preCellType,postCellType)
    text_file = open(outputfilenameFit, "w")
    n = text_file.write(X)
    text_file.close()

# Define fitting Function
def fitZIP(preCellType,postCellType):
    
    # Get data from file
    filename = "data_dense_model\%s_%s.csv" % (preCellType,postCellType)
    df = pd.read_csv(filename,header=None,names=["data"])
    
    # Prepare data for fitting
    X = df.data
    nobs = len(X)
    exog = np.ones(nobs)
    freq = np.bincount(X) / nobs
    binValue = list(range(0,len(freq)))
    
    # Fit Data
    mod_ZIP = sm.ZeroInflatedPoisson(X, exog)
    res_ZIP = mod_ZIP.fit(disp=False)
    
    # Get fitting results
    probs_zip = res_ZIP.predict(which='prob')
    probsm_zip = probs_zip.mean(0)
    
    # Export freq and probsm_zinb
    values = {'x': freq,
                'xFit': probsm_zip}
    outputDF = DataFrame(values, columns= ['x', 'xFit'])
    outputfilename = "fit_dense_model\%s_%s_ZIP.csv" % (preCellType,postCellType)
    export_csv = outputDF.to_csv (outputfilename,index=None,header=True)
    
    # Export fit results
    X = res_ZIP.summary().as_csv()
    outputfilenameFit = "fit_dense_model\%s_%s_ZIP_FitResults.csv" % (preCellType,postCellType)
    text_file = open(outputfilenameFit, "w")
    n = text_file.write(X)
    text_file.close()

# Iterate over all cell type combinations
preCellTypes = ['L2PY','L3PY','L4PY','L4sp','L4ss','L5IT','L5PT','L6ACC','L6BCC','L6CT','INH','VPM']
postCellTypes = ['L2PY','L3PY','L4PY','L4sp','L4ss','L5IT','L5PT','L6ACC','L6BCC','L6CT','INH']

for preCellType in preCellTypes:
    for postCellType in postCellTypes:
        print(preCellType + ' ' + postCellType + ' ZINB fit')
        fitZINB(preCellType,postCellType)
        print(preCellType + ' ' + postCellType + ' ZIP fit')
        fitZIP(preCellType,postCellType)

