#! /bin/env python3
# ------------------
# Library
# ------------------

import pandas as pd
import numpy as np
import logging

from .decorator import helpstring
from .decorator import AppendHelp

# ------------------
# Function
# ------------------

_helpdoc = dict()

_helpdoc['df_geomean'] = helpstring(
    describe='',
    parameterdicts={
        'dat': 'pd.DataFrame, data to process.',
        'label': 'list, column names of data to calculate mean.'
    },
    returns='pd.Series, containing the mean of each row of columns indicated by label.',
    examplecodelists=[
        "datgm = df_geomean(dataframe, ['one', 'two'])"
    ]
)

@AppendHelp(_helpdoc['df_geomean'], join='')
def df_geomean(dat, label):
    '''
    Function used to calculate geometric mean.
    The geometric mean of each row is calculated.
    '''
    if (len(label) == 0):
        logging.error('Length of label should be at least 1.')
        raise ValueError('Length of label should be at least 1.')
    elif (len(label) == 1):
        return dat[label]
    else:
        datlog = np.log(dat[label] + 1)
        datlogsum = datlog.sum(axis=1)
        gm = np.exp(datlogsum / len(label)) - 1
        return gm

# ------------------

_helpdoc['df_median_ratio_normfactor'] = helpstring(
    describe='',
    parameterdicts={
        'dat': 'pd.DataFrame, data to process.',
        'label': 'list, column names of data to calculate normalize.'
    },
    returns='pd.DataFrame, the corresponding normalized data.',
    examplecodelists=[
        "datnorm = df_normalization(dataframe, ['one', 'two'])"
    ]
)

@AppendHelp(_helpdoc['df_median_ratio_normfactor'], join='')
def df_median_ratio_normfactor(dat, label):
    '''
    Normalize input data indicated by the label.
    Now only median ratio normalization is available.
    '''
    datgm = np.exp(np.log(dat[label] + 1.0).sum(axis=1) / len(label))
    datgm[datgm <= 0] = 1
    meanfactor = dat[label].div(datgm, axis=0)
    normfactor = 1 / meanfactor.median(axis=0)
    return normfactor

# ------------------

def df_total_count_normfactor(dat, label):
    colsum = dat[label].sum(axis=0)
    colsummean = colsum.sum() / len(label)
    normfactor = colsummean / colsum
    return normfactor

# ------------------

def df_normalization(dat, label, method):
    normfactor = df_total_count_normfactor(dat, label)
    if method == 'none':
        normfactor = np.array([1]*len(normfactor))
    elif method == 'median':
        medianfactor = df_median_ratio_normfactor(dat, label)
        if (medianfactor == 0).any():
            logging.warning('Median factor is zero, using total count normalization')
        elif ((dat[label] == 0).sum(axis=0) / dat.shape[0] > 0.45).any():
            logging.warning('Too many zeros in counts, using total count normalization')
        else:
            normfactor = medianfactor
    result = dat[label].mul(normfactor, axis=1)
    return result

# ------------------

_helpdoc['df_leastsquare'] = helpstring(
    describe='',
    parameterdicts={
        'dat': 'pd.DataFrame, data to process.',
        'xlabel': 'string, column name of x in dat.',
        'ylabel': 'string, column name of y in dat.',
        'weight': 'string, column name of weight in dat, default is None means no weight considered.'
    },
    returns='tuple, (b, a) of the y = a + bx function.',
    examplecodelists=[
        "k, b = df_leastsquare(data['x'], data['y'], data['w'])"
    ]
)

@AppendHelp(_helpdoc['df_leastsquare'], join='')
def df_leastsquare(dat, xlabel, ylabel, weightlabel=None):
    '''
    Least squares Fitting
    coefficients from y = a + bx
    reference: http://mathworld.wolfram.com/LeastSquaresFitting.html
    For weighted least square: http://goo.gl/pGpTZ6
    '''
    n = dat.shape[0]
    if weightlabel is None:
        sy = dat[ylabel].sum()
        sx = dat[xlabel].sum()
        sx2 = dat[xlabel].mul(dat[xlabel], axis=0).sum()
        sxy = dat[xlabel].mul(dat[ylabel], axis=0).sum()
        a = (sy*sx2-sx*sxy) / (n*sx2-sx*sx)
        b = (n*sxy-sx*sy) / (n*sx2-sx*sx)
        return (b,a)
    else:
        nw = dat[weightlabel].sum()
        sy = dat[ylabel].mul(dat[weightlabel], axis=0).sum()
        sx = dat[xlabel].mul(dat[weightlabel], axis=0).sum()
        sx2 = dat[xlabel].mul(
            dat[xlabel], axis=0
        ).mul(
            dat[weightlabel], axis=0
        ).sum()
        sxy = dat[xlabel].mul(
            dat[ylabel], axis=0
        ).mul(
            dat[weightlabel], axis=0
        ).sum()
        a = (sy * sx2 - sx * sxy) / (nw * sx2 - sx * sx)
        b = (nw * sxy - sx * sy) / (nw * sx2 - sx * sx)
        return (b,a)

# ------------------

_helpdoc['df_modelmeanvar'] = helpstring(
    describe='',
    parameterdicts={
        'dat': 'pd.DataFrame, data to process.',
        'label': 'list, column names of data to calculate the model.'
    },
    returns='tuple, (k, b) of the ln(Var - Mean) = ln(b) + k * ln(Mean).',
    examplecodelists=[
        "k, b = df_modelmeanvar(data, ['R1', 'R2'])"
    ]
)

@AppendHelp(_helpdoc['df_modelmeanvar'], join='')
def df_modelmeanvar(dat, label):
    '''
    Calculate the parameters in  Negative Binomial model of sequence tag data.
    Model: Var = Mean + alpha * Mean ^ beta
    with: log2(Var - Mean) = log2(b) + k * log2(Mean)
          Var = Mean + 2 ^ b * Mean ^ k
          alpha = 2 ^ b, beta = k
    '''
    datgm = df_geomean(dat, label)
    datvar = dat[label].var(axis=1)
    goodidx = datgm < datvar
    lgm = np.log2(datgm[goodidx] + 1)
    lvar = np.log2(datvar[goodidx] - datgm[goodidx] + 1)
    data = pd.DataFrame(
        {
            'lgm': lgm,
            'lvar': lvar,
            'w': datgm[goodidx]
        }
    )
    k, b = df_leastsquare(data, 'lgm', 'lvar', 'w')
    k = max(k, 1)
    b = max(b, 0)
    return (k, b)

# ------------------

_helpdoc['df_estvar'] = helpstring(
    describe='',
    parameterdicts={
        'dat': 'pd.DataFrame, data to process.',
        'meanlabel': 'string, column name of mean data in dat.',
        'k': 'numeric, k in Var = Mean + 2 ^ b * Mean ^ k.',
        'b': 'numeric, b in Var = Mean + 2 ^ b * Mean ^ k.'
    },
    returns='pd.Series, adjusted variance of the mean inputed.',
    examplecodelists=[
        "estvar = df_estvar(data, ['mean'], 2.5, 1)"
    ]
)

@AppendHelp(_helpdoc['df_estvar'], join='')
def df_estvar(dat, meanlabel, k, b):
    '''
    Calculate the modeled variance by mean and k, b parameters in Negative Binomial model of sequence tag data.
    Model: Var = Mean + alpha * Mean ^ beta
    with: log2(Var - Mean) = log2(b) + k * log2(Mean)
          Var = Mean + 2 ^ b * Mean ^ k
          alpha = 2 ^ b, beta = k
    '''
    var = dat[meanlabel].map(
        lambda x: x ** k
    ) * 2 ** b + dat[meanlabel]
    return var

# ------------------
# EOF
# ------------------
