#! /usr/bin/env python3
import os
import pandas as pd
import numpy as np
from functools import reduce
from itertools import combinations
import argparse

parser = argparse.ArgumentParser(description='mergebarcode')

parser.add_argument(
    '-r', '--reference', default='reference.txt'
)

parser.add_argument(
    '--controllabel', nargs='*', required=True
)

parser.add_argument(
    '--controlinput', nargs='*', required=True
)

parser.add_argument(
    '--treatlabel', nargs='*', required=True
)

parser.add_argument(
    '--treatinput', nargs='*', required=True
)

parser.add_argument(
    '--output', default='result.csv'
)

args = parser.parse_args()

reference = pd.read_table(args.reference, header=0)

labels = args.controllabel + args.treatlabel

files = args.controlinput + args.treatinput

rawdata = dict()
bestdict = dict()

# reading data
for l, f in zip(labels, files):
    data = pd.read_table(f, header=None, sep=' ')
    data.columns = ['guide', 'bar1', 'bar2', 'count'] # make title
    data.bar2 = data.bar2.str.strip()
    rawdata[l] = data

    # join data with library by bar1 and bar2 separately
    # join data with library by bar1
    # colnames: guide, bar1, bar2, count, gene, barcode
    data = pd.merge(
        data,
        reference[['gene', 'guide', 'barcode']],
        left_on=['guide', 'bar1'],
        right_on=['guide', 'barcode'],
        how='left'
    )

    # join data with library by bar2
    # colnames: guide, bar1, bar2, count, gene_x, barcode_x, gene_y, barcode_y
    data = pd.merge(
        data,
        reference[['gene', 'guide', 'barcode']],
        left_on=['guide', 'bar2'],
        right_on=['guide', 'barcode'],
        how='left'
    )

    # merge gene_x, gene_y, and barcode_x, barcode_y
    # make all the gene_y which the gene_x not null to empty string
    data.loc[data['gene_x'].notnull(), 'gene_y'] = ''
    # fill na in gene_x to empty string
    data['gene_x'].fillna('', inplace=True)
    # fill rest na in gene_y to empty string
    data['gene_y'].fillna('', inplace=True)
    # concatinate the strings in gene_x and gene_y to form gene
    data['gene'] = data['gene_x'].str.cat(data['gene_y'])
    # make all the barcode_y which the relative barcode_x is not null to empty string
    data.loc[data['barcode_x'].notnull(), 'barcode_y'] = ''
    # fill na in barcode_x to empty string
    data['barcode_x'].fillna('', inplace=True)
    # fill rest na in barcode_y to empty string
    data['barcode_y'].fillna('', inplace=True)
    # concatinate the strings in barcode_x and barcode_y to form gene
    data['barcode'] = data['barcode_x'].str.cat(data['barcode_y'])

    # remove gene_x gene_y barcode_x barcode_y bar1 bar2 columns
    del data['gene_x']
    del data['gene_y']
    del data['barcode_x']
    del data['barcode_y']
    del data['bar1']
    del data['bar2']

    bestdata = data.loc[data['gene'] != ''].groupby(
        ['gene', 'guide', 'barcode']
    ).sum().reset_index()
    bestdata.columns = ['gene', 'guide', 'barcode', l]
    bestdict[l] = bestdata

    # make good data
    del data['gene']
    del data['barcode']


# merge all the columns into one matrix
def mapmerge(x, y):
    return pd.merge(
        x, y,
        on=['gene', 'guide', 'barcode'],
        how='outer'
    )

clab = args.controllabel

# print(x)
bestresult = reduce(mapmerge, [bestdict[i] for i in labels])
bestresult.sort_values(['gene', 'guide'], inplace=True)
# controlcol = [i for i in labels if i in clab]
# bestresult = bestresult[['gene', 'guide', 'barcode'] + labels].loc[
#     reduce(np.logical_and, [bestresult[i] > 10 for i in controlcol])
# ]
bestresult.fillna(0).to_csv(
    args.output,
    index=False
)

print('finished')

####################
