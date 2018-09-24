#! /bin/env python3
# ------------------
# Library
# ------------------

import pandas as pd
import numpy as np
import logging
import os
from scipy.stats import norm

from .decorator import helpstring
from .decorator import AppendHelp
from .dfcalculate import df_geomean
from .dfcalculate import df_normalization
from .dfcalculate import df_leastsquare
from .dfcalculate import df_modelmeanvar
from .dfcalculate import df_estvar

from .sysrun import robustrank

# ------------------
# Function
# ------------------

_helpdoc = dict()

_helpdoc['analysis'] = helpstring(
    describe='',
    parameterdicts={
        'inputdata': 'pd.DataFrame, data to process.',
        'outprefix': 'string, output data name prefix.',
        'controlids': 'list, column names of control data',
        'treatids': 'list, column names of treatment data',
        'hasbarcode': 'bool, whether the screening using barcode',
        'normthreshold': 'numeric, threshold used in scoring, the normalized data less than the score will be punished.',
        'test': 'string, test method, "norm" for normal test.',
        'rrapath': 'string, path of RobustRankAggregation program.'
    },
    returns='No specific returns.',
    examplecodelists=[
        "analysis(",
        "    inputdata,",
        "    outprefix=outprefix,",
        "    controlids=['control1', 'control2'],",
        "    treatids=['treat1', 'treat2'],",
        "    hasbarcode=True,",
        "    normthreshold=10,",
        "    test='norm',",
        "    rrapath='RRA'",
        ")"
    ]
)

@AppendHelp(_helpdoc['analysis'], join='')
def analysis(inputdata,
             outprefix,
             controlids,
             treatids,
             hasbarcode=True,
             normthreshold=10,
             test='norm',
             rrapath='RRA'):
    '''
    Pipeline function in testing of the CRISPR/Cas9 screening data.
    '''
    # setting logging level

    # output file names
    files = {
        'sgrnaout': outprefix + '.sgrna.txt',
        'plowout': outprefix + '.plow.txt',
        'phighout': outprefix + '.phigh.txt',
        'genelow': outprefix + '.gene.low.txt',
        'genehigh': outprefix + '.gene.high.txt'
    }

    # make the column names
    infocolnm = ['gene', 'guide', 'gid', 'barcode', 'bid']

    # make the labels for all consider counts columns in calculation
    conlabels = controlids + treatids
    if (len(controlids) > 1):
        conlabels = controlids

    # normalization
    logging.info('Normalizing data.')

    datanorm = pd.concat(
        [
            inputdata[infocolnm],
            df_normalization(inputdata, controlids + treatids, 'median')
        ],
        axis=1
    )

    datanorm.columns = infocolnm + controlids + treatids

    data = datanorm

    # linear regression
    logging.info('Estimating parameters of mean and var of controls.')

    k, b = df_modelmeanvar(
        data,
        conlabels
    )

    logging.info(
        'Estimated: Var = Mean + {0:.2f} * Mean ^ {1:.2f}.'.format(
            2 ** b, k
        )
    )
    # calculate mean and var
    logging.info('Calculating means and variance of data.')

    data['controlmean'] = df_geomean(data, controlids)
    data['treatmean'] = df_geomean(data, treatids)
    data['controlvar'] = data[conlabels].var(axis=1)
    data['estvar'] = df_estvar(data, 'controlmean', k, b)

    # log fold change
    logging.info('Calculating guide log2 fold change.')

    data['lfc'] = np.log2(
        data['treatmean'] + 1
    ) - np.log2(
        data['controlmean'] + 1
    )

    # lfc direction
    data['lfc_bin'] = data['lfc'].map(
        lambda x: 1 if x > -0.1 else 0
    ) + data['lfc'].map(
        lambda x: -1 if x < 0.1 else 0
    )

    # large norm data
    data['large'] = (
        data[controlids + treatids] > normthreshold
    ).sum(axis=1) > len(controlids)

    # sgRNA in barcode with same direction
    data['direction'] = data['lfc_bin'].mul(
        data['large'], axis=0
    )

    # calculate adjusted variance
    if hasbarcode:
        logging.info('Adjusting variance of data in guide level.')
        # aggregate variance
        agg_guide_var = data[
            ['guide', 'controlvar']
        ].groupby('guide').mean()

        agg_guide_samedirection = data[['guide', 'direction']].groupby(
            'guide'
        ).apply(
            lambda x: max(x['direction']) * min(x['direction'])
        )

        agg_guide_guidecount = data['guide'].value_counts()

        agg_guide = pd.DataFrame(
            {
                'controlvar': agg_guide_var['controlvar'],
                'guidecount': agg_guide_guidecount[agg_guide_var.index],
                'samedirection': agg_guide_samedirection[agg_guide_var.index] != -1
            }
        )

        # adjvar
        adjust_var = agg_guide.loc[
            data['guide'], 'controlvar'
        ].mul(
            (1 - agg_guide.loc[data['guide'], 'samedirection']),
            axis=0
        )
        adjust_var.index = data['guide'].index
        data['adjvar'] = data['estvar'] + adjust_var
    else:
        data['adjvar'] = data['estvar']


    # normalize treatment mean value
    logging.info('Normalizing treatment values.')
    theta = (
        data['treatmean'] - data['controlmean']
    ).div(np.sqrt(data['adjvar']), axis=0)

    data = data.assign(
        treat_z=theta
    )

    if test == 'norm':
        data['p.low'] = data['treat_z'].map(norm.cdf)
        data['p.high'] = data['treat_z'].map(norm.sf)

    data.to_csv(files['sgrnaout'], index=False, sep='\t')

    # prepare for Robust Rank Aggregation
    pcolnm = ['sgrna', 'symbol', 'pool', 'plow', 'prob', 'chosen']
    # lower direction
    logging.info('Robust Rank Aggregation of lower direction data.')
    plowout = pd.DataFrame(
        {
            'sgrna': data['bid'],
            'symbol': data['gene'],
            'pool': ['list'] * data['gene'].size,
            'plow': data['treat_z'],
            'prob': [1] * data['gene'].size,
            'chosen': [1] * data['gene'].size
        }
    )

    plowout[pcolnm].sort_values(
        'plow'
    ).to_csv(files['plowout'], index=False, sep='\t')

    percentilelow = (data['p.low'] < 0.25).sum() / data['p.low'].size
    robustrank(
        rrapath,
        infile=files['plowout'],
        outfile=files['genelow'],
        percentile=percentilelow
    )

    # higher direction
    logging.info('Robust Rank Aggregation of higher direction data.')
    phighout = pd.DataFrame(
        {
            'sgrna': data['bid'],
            'symbol': data['gene'],
            'pool': ['list'] * data['gene'].size,
            'plow': -data['treat_z'],
            'prob': [1] * data['gene'].size,
            'chosen': [1] * data['gene'].size
        }
    )

    phighout[pcolnm].sort_values(
        'plow'
    ).to_csv(files['phighout'], index=False, sep='\t')

    percentilehigh = (data['p.high'] < 0.25).sum() / data['p.high'].size

    robustrank(
        rrapath,
        infile=files['phighout'],
        outfile=files['genehigh'],
        percentile=percentilehigh
    )

# ------------------
####################
