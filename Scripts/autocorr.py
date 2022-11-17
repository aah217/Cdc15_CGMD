import numpy as np
import pandas as pd
import time
from scipy import optimize
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-f","--file", help="file",type=str)
parser.add_argument("-c","--col", help="column for autocorrelation analysis",type=int)
parser.add_argument("-o","--outfile", help="output file",type=str)
args = parser.parse_args()

startTime = time.time()

file = args.file
outfile = args.outfile
col = args.col - 1

def test_func(x, a):
    return np.exp(-x / a)
print("Reading file...",flush=True)
my_csv = pd.read_csv(file, sep="\s+", header=None)
s = my_csv[col]
x_list = []
y_list = []
print("Autocorrelation...",flush=True)
mod_check = round(float(s.size)/100.0) #gives 1% progress increments
for i in range(0,s.size-2):
    if i % mod_check == 0.0: print("...",round(100.0*float(i)/float(s.size)),"%",end='\r',flush=True)
    x_list.append(i)
    y_list.append(s.autocorr(i))
x_data = np.array(x_list)    
y_data = np.array(y_list)
print("Curve fit...",flush=True)
params, params_covariance = optimize.curve_fit(test_func, x_data, y_data)
print(params[0],flush=True)
stop = np.ceil(2*params[0]*np.log(2))
print("split on 2*tau_0.5 = ",stop,flush=True)

print("Writing out file...",flush=True)
with open(file,'r') as f, open(outfile,'w') as o:        
    for i,line in enumerate(f):
        if i+1 >= stop:
            o.write(line)
executionTime = (time.time() - startTime)
print('Execution time in seconds: ' + str(executionTime),flush=True)     