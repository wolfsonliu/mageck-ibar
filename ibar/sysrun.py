# ------------------
# Library
# ------------------

import logging
import os

from .decorator import helpstring
from .decorator import AppendHelp

# ------------------
# Function
# ------------------

_helpdoc = dict()

_helpdoc['robustrank'] = helpstring(
    describe='',
    parameterdicts={
        'rrapath': '',
        'infile': 'string, the file path of input data. Format: <item id> <group id> <list id> <value> [<probability>] [<chosen>]',
        'outfile': 'string, the file path of output data. Format: <group id> <number of items in the group> <lo-value> <false discovery rate>',
        'percentile': 'numeric, RRA only consider the items with percentile smaller than this parameter. Default=0.1'
    },
    returns='output txt files with the RRA results',
    examplecodelists=[
        "robustrank(",
        "    'RRA',",
        "    infile=out.plow.txt,",
        "    outfile=out.gene.low.txt,",
        "    percentile=0.1",
        ")"
    ]
)

@AppendHelp(_helpdoc['robustrank'], join='')
def robustrank(rrapath, infile, outfile, percentile):
    '''
    Wrapper function of RRA, which was writen by Wei Li.

    Input file:
        Format: <item id> <group id> <list id> <value> [<probability>] [<chosen>]
    Output file:
        Format: <group id> <number of items in the group> <lo-value> <false discovery rate>
    Percentile:
        Maximum percentile.
        RRA only consider the items with percentile smaller than this parameter. Default=0.1.
    '''
    logging.info(
        'RRA start: maximum percentile is {0:.6f}'.format(percentile)
    )
    cmd = ' '.join(
        [
            rrapath, '-i', infile, '-o', outfile,
            '-p', str(percentile),
            '--skip-gene NA --skip-gene na'
        ]
    )
    os.system(cmd)
    logging.info('RRA finished.')

# ------------------
# EOF
# ------------------
