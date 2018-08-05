#! /bin/env python3
# ------------------
# Library
# ------------------

import pandas as pd
import logging

from .decorator import helpstring
from .decorator import AppendHelp

# ------------------
# Function
# ------------------

# checking functions

def iscsv(filepath):
    return filepath.split('.')[-1] == 'csv'

# ------------------


def istsv(filepath):
    return filepath.split('.')[-1] == 'tsv'

# ------------------


def istxt(filepath):
    return filepath.split('.')[-1] == 'txt'

# ------------------


_helpdoc = dict()

_helpdoc['readdata'] = helpstring(
    describe='',
    parameterdicts={
        'filepath': 'string, indicate the inputdata file, should be csv or tsv (txt)',
        'genelab': 'string, the column name of gene.',
        'guidelab': 'string, the column name of guide.',
        'barcodelab': 'string, the column name of barcode.',
        'controlids': 'list, list contains the column names of control data.',
        'treatids': 'list, list contains the column names of treatment data.',
    },
    returns='pd.DataFrame, the DataFrame containing the input data.',
    examplecodelists=[
        "inputdata = analysis.readdata(",
        "    inputpath,",
        "    genelab='gene',",
        "    guidelab='guide',",
        "    barcodelab='barcode',",
        "    controlids='control',",
        "    treatids='treat',",
        "    barcode=True",
        ")"
    ]
)

@AppendHelp(_helpdoc['readdata'], join='')
def readdata(
        filepath,
        genelab,
        guidelab,
        barcodelab,
        controlids,
        treatids,
        hasbarcode=True):
    '''
    Reading input data.
    The input data should be saved as csv or tsv (txt),
    each columns should have a header and separate by "," or tab.
    The columns should contain:
        <gene>, <guide>, <barcode>,
        <control1>, [<control2>, ...],
        <treat1>, [<treat2>, ...]
    '''
    # check whether the input file as csv or tsv
    logging.info('Reading data: {0:s}.'.format(filepath))

    if (iscsv(filepath)):
        inputdata = pd.read_csv(filepath, header=0)
    elif (istsv(filepath) or istxt(filepath)):
        inputdata = pd.read_table(filepath, header=0, sep='\t')
    else:
        logging.error('Input file should be csv, tsv or txt, column should be separate by , or tab.')
        raise ValueError('Wrong input file type.')

    logging.info(
        'Data with {0:d} Controls, {1:d} Treatments,'.format(
            len(controlids), len(treatids)
        ) + ' {0:d} guide RNAs'.format(inputdata.shape[0])
    )

    # if there is no barcode ignore the column of barcode
    if hasbarcode:
        colnm1 = [
            genelab, guidelab, barcodelab
        ] + controlids + treatids
    else:
        colnm1 = [
            genelab, guidelab, guidelab
        ] + controlids + treatids
    colnm2 = [
        'gene', 'guide', 'barcode'
    ] + controlids + treatids
    colnm3 = [
        'gene', 'guide', 'gid', 'barcode', 'bid'
    ] + controlids + treatids

    # select columns
    data = inputdata[colnm1].copy()
    # rename column names
    data.columns = colnm2

    # make guide level id
    data['gid'] = data[['gene', 'guide']].apply(
        lambda x: '.'.join(x), axis=1
    )
    # make barcode level id if barcode exists
    if hasbarcode:
        data['bid'] = inputdata[['gene', 'guide', 'barcode']].apply(
            lambda x: '.'.join(x), axis=1
        )
    else:
        data['bid'] = data['gid']

    # zero count to 1
    for x in controlids + treatids:
        data.loc[data[x] == 0, x] = 1

    return data[colnm3]

# ------------------
