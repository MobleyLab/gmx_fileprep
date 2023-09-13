# utils.py
#  * this module contains utility functions that may be helpful
#    for any module

import sys


def printerr(*args, **kwargs):
    ''' print arguments to stderr
        use the same as builtin print()
    '''
    print(*args, file=sys.stderr, **kwargs)