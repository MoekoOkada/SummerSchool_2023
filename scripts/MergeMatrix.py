# Script name: MergeMatrix.py
# Project: SummerSchool
# Date: 11-Jul-2023
# usage: python MergeMatrix.py rbh hal lyr output


from os import sep
import pandas as pd
import io
import sys
args = sys.argv

# read files
rbh = pd.read_table(args[1], sep = '\t', header = 0)
rbh['Ahal'] = rbh['Ahal'].str.replace('g0000', 'g')
rbh['Ahal'] = rbh['Ahal'].str.replace('g000', 'g')
rbh['Ahal'] = rbh['Ahal'].str.replace('g00', 'g')
rbh['Ahal'] = rbh['Ahal'].str.replace('g0', 'g')

print(rbh.head())

rbh['Alyr'] = rbh['Alyr'].str.replace('g0000', 'g')
rbh['Alyr'] = rbh['Alyr'].str.replace('g000', 'g')
rbh['Alyr'] = rbh['Alyr'].str.replace('g00', 'g')
rbh['Alyr'] = rbh['Alyr'].str.replace('g0', 'g')
print(rbh.head())

hal = pd.read_table(args[2], sep='\t', header = 0)
lyr = pd.read_table(args[3], sep='\t', header = 0)

hal = hal.rename(columns={'gene_id': 'Ahal'})
hal['Ahal'] = hal['Ahal'].str.replace('\.t1', '')
print(hal.head())

lyr = lyr.rename(columns={'gene_id': 'Alyr'})
lyr['Alyr'] = lyr['Alyr'].str.replace('\.t1', '')
print(lyr.head())

mat1 = pd.merge(rbh, hal, on = 'Ahal', how = 'left')
print(mat1.head())

mat2 = pd.merge(mat1, lyr, on = 'Alyr', how = 'left')
print(mat2.head())

mat2.to_csv(args[4], header = True, index=None, sep="\t")
