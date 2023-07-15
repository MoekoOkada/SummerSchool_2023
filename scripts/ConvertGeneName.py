# Script name: ConvertGeneName.py
# Project: SummerSchool
# Date: 14-Jul-2023
# usage: MoekoOkada

from os import sep
import pandas as pd
import io
import sys
args = sys.argv

# read files
rbh = pd.read_table(args[1], sep = '\t', header = 0)
rbh['Ahal'] = rbh['Ahal'].str.replace('scaffold\_\d*\.', '')
rbh['TAIR10_cds'] = rbh['TAIR10_cds'].str.replace('\.\d+', '')
print(rbh.head())

metal = pd.read_table(args[2], sep='\t', header = 0)

metal = metal.rename(columns={'GeneID': 'TAIR10_cds'})
print(metal.head())

metal2 = pd.merge(metal, rbh, on = 'TAIR10_cds', how = 'left')
print(metal2.head())

metal2.to_csv(args[3], header = True, index=None, sep="\t")
