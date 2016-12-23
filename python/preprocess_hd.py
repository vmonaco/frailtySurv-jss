'''
Created on Jul 14, 2015

@author: vinnie
'''

import os, sys
import pandas as pd
from glob import glob
from IPython import embed

# DATA_DIR = os.getenv('HOME') + '/workspace/data/hard-drive-failure/data/'

DATE_COLS = ['date','day','month','dayofyear','dayofweek','weekofmonth']

SMART_COLS = {
'smart_1_raw':'rer',
'smart_5_raw':'rsc',
'smart_9_raw':'time',
'smart_194_raw':'temp',
'smart_197_raw':'psc',

# 'smart_3_raw':'spin_up_time', 
# 'smart_4_raw':'start_stop_count',
# 'smart_7_raw':'seek_error_rate', 
# 'smart_10_raw':'spin_retry_count',
# 'smart_12_raw':'power_cycle_count',
# 'smart_192_raw':'unsafe_shutdown_count',
# 'smart_193_raw':'load_cycle_count',
# 'smart_198_raw':'uncorrectable_sector_count',
# 'smart_199_raw':'crc_error_count',
}

SMART_COLS_KEYS = list(SMART_COLS.keys())

def parse_date(df):
    date = pd.DatetimeIndex(df['date'], format='%Y-%M-%D')
    df['day'] = date.day
    df['month'] = date.month
    df['dayofyear'] = date.dayofyear
    df['dayofweek'] = date.dayofweek
    df['weekofmonth'] = date.day//7
    return df

def load_data(fnames):
    
    # Should sort by capture date, which is the filename
    fname = sorted(fnames)
    
    # Load the most recent. Assumes the fnames are sorted by date
    df = pd.read_csv(fnames[-1], index_col=[1])
    
    # Make sure serials are unique
    df = df.groupby(df.index).last()
    
    # Load additional data: new drives
    for fname in fnames[-2::-1]:
        tmp = pd.read_csv(fname, index_col=[1])
        
        # if tmp contains a new drive, add it
        df = df.append(tmp[tmp.index.get_level_values('serial_number').isin(df.index.get_level_values('serial_number').values)==False])
        
        # Get the most recent smart stats before a hd failed
        # if tmp contains an existing drive that failed, and has SMART data, add it
        idx = (df['failure']==1)&df.index.get_level_values('serial_number').isin(tmp.index.get_level_values('serial_number'))&df[SMART_COLS_KEYS].isnull().any(axis=1)
        serials = df[idx].index.get_level_values('serial_number').values
        for serial in serials:
            df.loc[serial, SMART_COLS_KEYS] = tmp.loc[serial, SMART_COLS_KEYS]
    
    df.index.name = 'serial'
    
    # HD age in days
    df['time'] = df['smart_9_raw']/24
    df['status'] = df['failure'].astype(int)
    df['capacity'] = df['capacity_bytes']
    
    df = df.dropna(subset=SMART_COLS.keys()).dropna(axis=1)
    
    # Creates new columns
#     df = parse_date(df)
    
    # Select interesting columns and set the idx
    df = df[['model', 'time','status','capacity'] + [c for c in df.columns if 'raw' in c]]
    df.columns = [SMART_COLS[colname] if colname in SMART_COLS.keys() else colname for colname in df.columns]
    
    df['rer'] = (df['rer'] > 0).astype(int)
    df['rsc'] = (df['rsc'] > 0).astype(int)
    df['psc'] = (df['psc'] > 0).astype(int)
    
    df = df.sort('model')
    
    return df

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('Usage: python preprocess_hd.py [datadir] [output.csv]')
        sys.exit(1)
        
    datadir = sys.argv[1]
    output = sys.argv[2]
    fnames = sorted(glob(datadir + '*.csv'))
    df = load_data(fnames)
    df.to_csv(output)