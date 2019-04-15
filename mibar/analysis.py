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
from .dfcalculate import array_fdr
from .programio import read_rra
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
             gene_test_threshold=0.25,
             test='norm',
             tworra=False,
             rrapath='RRA'):
    '''
    Pipeline function in testing of the CRISPR/Cas9 screening data.
    '''
    # setting logging level

    # output file names
    files = {
        'barcodeout': outprefix + '.barcode.txt',
        'sgrnaout': outprefix + '.sgrna.txt',
        'geneout': outprefix + '.gene.txt',
        'plowout': outprefix + '.plow.txt',
        'phighout': outprefix + '.phigh.txt',
        'sgrnalow': outprefix + '.sgrna.low.txt',
        'sgrnahigh': outprefix + '.sgrna.high.txt',
        'genelow': outprefix + '.gene.low.txt',
        'genehigh': outprefix + '.gene.high.txt'
    }
    files['firstlevel'] = files['sgrnaout']
    if hasbarcode or tworra:
        files['firstlevel'] = files['barcodeout']
    files['rra_low_in'] = files['plowout']
    files['rra_low_out'] = files['genelow']
    files['rra_high_in'] = files['phighout']
    files['rra_high_out'] = files['genehigh']
    if tworra:
        files['rra_low_in'] = files['plowout']
        files['rra_low_out'] = files['sgrnalow']
        files['rra_high_in'] = files['phighout']
        files['rra_high_out'] = files['sgrnahigh']
        files['rra2_low_in'] = files['sgrnalow']
        files['rra2_low_out'] = files['genelow']
        files['rra2_high_in'] = files['sgrnahigh']
        files['rra2_high_out'] = files['genehigh']

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
        data['treatmean'] + 1.0
    ) - np.log2(
        data['controlmean'] + 1.0
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

    if tworra:
        hasbarcode = False

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
        treat_zscore=theta
    )

    if test == 'norm':
        data['p.low'] = data['treat_zscore'].map(norm.cdf)
        data['p.high'] = data['treat_zscore'].map(norm.sf)

    data['p.twoside'] = data[['p.low', 'p.high']].apply(
        lambda x: 2 * x['p.low'] if x['p.low'] < x['p.high'] else 2 * x['p.high'],
        axis=1
    )
    data['fdr'] = array_fdr(data['p.twoside'])
    data.to_csv(files['firstlevel'], index=False, sep='\t')

    # fold change
    foldchange = data.groupby(['gene'])['lfc'].mean().reset_index()
    data['symbol'] = data['gene']
    if tworra:
        foldchange = data.groupby(['gene', 'guide'])['lfc'].mean().reset_index()
        data['symbol'] = data['gid']

    # prepare for Robust Rank Aggregation
    pcolnm = ['sgrna', 'symbol', 'pool', 'p', 'prob', 'chosen']
    # lower direction
    logging.info('Robust Rank Aggregation of lower direction data.')
    plowout = pd.DataFrame(
        {
            'sgrna': data['bid'],
            'symbol': data['symbol'],
            'pool': ['list'] * data['gene'].size,
            'p': data['treat_zscore'],
            'prob': [1] * data['gene'].size,
            'chosen': [1] * data['gene'].size
        }
    )

    plowout[pcolnm].sort_values(
        'p'
    ).to_csv(files['rra_low_in'], index=False, sep='\t')

    percentilelow = (
        data['p.low'] < gene_test_threshold
    ).sum() / data['p.low'].size
    robustrank(
        rrapath,
        infile=files['rra_low_in'],
        outfile=files['rra_low_out'],
        percentile=percentilelow
    )
    rralow = read_rra(files['rra_low_out'])
    # columns: group_id, items_in_group, beta, p, FDR, goodsgrna

    # higher direction
    logging.info('Robust Rank Aggregation of higher direction data.')
    phighout = pd.DataFrame(
        {
            'sgrna': data['bid'],
            'symbol': data['symbol'],
            'pool': ['list'] * data['gene'].size,
            'p': data['treat_zscore'] * -1,
            'prob': [1] * data['gene'].size,
            'chosen': [1] * data['gene'].size
        }
    )

    phighout[pcolnm].sort_values(
        'p'
    ).to_csv(files['rra_high_in'], index=False, sep='\t')

    percentilehigh = (
        data['p.high'] < gene_test_threshold
    ).sum() / data['p.high'].size

    robustrank(
        rrapath,
        infile=files['rra_high_in'],
        outfile=files['rra_high_out'],
        percentile=percentilehigh
    )
    rrahigh = read_rra(files['rra_high_out'])
    # columns: group_id, items_in_group, beta, p, FDR, goodsgrna
    mresult = pd.merge(
        rralow, rrahigh, how='inner',
        on=['group_id'], suffixes=['.low', '.high']
    )

    if tworra:
        # low
        plowout2 = pd.DataFrame(
            {
                'sgrna': rralow['group_id'],
                'symbol': rralow['group_id'].str.split('.').map(lambda x: x[0]),
                'pool': ['list'] * rralow['group_id'].size,
                'p': rralow['beta'],
                'prob': [1] * rralow['group_id'].size,
                'chosen': [1] * rralow['group_id'].size,
            }
        )
        plowout2[pcolnm].sort_values(
            'p'
        ).to_csv(files['rra2_low_in'], index=False, sep='\t')

        rra2percentilelow = (
            rralow['FDR'] < gene_test_threshold
        ).sum() / rralow['FDR'].size
        robustrank(
            rrapath,
            infile=files['rra2_low_in'],
            outfile=files['rra2_low_out'],
            percentile=rra2percentilelow
        )
        rralow2 = read_rra(files['rra2_low_out'])
        # high
        phighout2 = pd.DataFrame(
            {
                'sgrna': rrahigh['group_id'],
                'symbol': rrahigh['group_id'].str.split('.').map(lambda x: x[0]),
                'pool': ['list'] * rrahigh['group_id'].size,
                'p': rrahigh['beta'],
                'prob': [1] * rrahigh['group_id'].size,
                'chosen': [1] * rrahigh['group_id'].size,
            }
        )
        plowout2[pcolnm].sort_values(
            'p'
        ).to_csv(files['rra2_high_in'], index=False, sep='\t')
        rra2percentilehigh = (
            rrahigh['FDR'] < gene_test_threshold
        ).sum() / rrahigh['FDR'].size
        robustrank(
            rrapath,
            infile=files['rra2_high_in'],
            outfile=files['rra2_high_out'],
            percentile=rra2percentilehigh
        )
        rrahigh2 = read_rra(files['rra2_high_out'])
        mresult = pd.merge(
            rralow2, rrahigh2, how='inner',
            on=['group_id'], suffixes=['.low', '.high']
        )
    return mresult



# ------------------
####################
