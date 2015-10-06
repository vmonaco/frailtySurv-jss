'''
Created on Jul 24, 2015

@author: vinnie
'''

import os,sys
import pandas as pd
from IPython import embed

RAW_COLS = [
'subject_id',
'laser_type',
'treated_eye',
'age_at_onset',
'diabetes_type',
'treated_risk',
'treated_status',
'treated_time',
'untreated_risk',
'untreated_status',
'untreated_time',
]

TREATED_COLS = ['treated_risk', 'treated_status', 'treated_time']
UNTREATED_COLS = ['untreated_risk', 'untreated_status', 'untreated_time']

def load_data(fname='http://www.mayo.edu/research/documents/diabetesdat/DOC-10027862'):
    
    df = pd.read_fwf(fname, header=None, widths=[5,5,5,5,5,5,5,10,5,5,10])
    df.columns = RAW_COLS
    return df

def preprocess(df):
    
    treated_eye = df['treated_eye'].squeeze()
    
    right_eye = df[TREATED_COLS] if treated_eye == 1 else df[UNTREATED_COLS]
    left_eye = df[TREATED_COLS] if treated_eye == 2 else df[UNTREATED_COLS]
    
    right_eye.columns = ['risk','status','time']
    left_eye.columns = ['risk','status','time']
    
    right_eye['eye'] = 1
    left_eye['eye'] = 2
    
    right_eye['treated'] = 1 if treated_eye == 1 else 0
    left_eye['treated'] = 1 if treated_eye == 2 else 0
    
    newdf = pd.concat([right_eye, left_eye], axis=0)
    newdf['age_at_onset'] = df['age_at_onset']
    newdf['laser_type'] = df['laser_type']
    newdf['diabetes_type'] = df['diabetes_type']
    newdf = newdf[['eye','time','status','treated','age_at_onset','laser_type','diabetes_type']]
    return newdf

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage: python preprocess_drs.py [output.csv]')
        sys.exit(1)
        
    df = load_data()
    df = df.groupby('subject_id').apply(preprocess).reset_index(level=1, drop=True)
    df.to_csv(sys.argv[1])