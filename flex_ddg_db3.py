#!/usr/bin/env python3

'''Small helpers shared by analyze_flex_ddG.py and extract_structures.py.

Standard library only, so that extract_structures.py keeps running without numpy or pandas.
'''

import re
import sqlite3

# Rosetta records the full option list it was invoked with in the "protocols" table of every
# features database it writes, so the stride a run actually used can be read back out of the
# output instead of being remembered and re-typed by hand.
_stride_re = re.compile( r'backrub_trajectory_stride[= ]+(\d+)' )

def trajectory_stride_from_db3( db3_file ):
    '''The backrub_trajectory_stride the run in db3_file was launched with.

    Returns None if it cannot be determined, which happens for databases written by Rosetta
    versions that predate the protocols table, or if the stride was left at the XML default
    instead of being passed on the command line. Callers should fall back to a default and say
    so, because the stride does not affect any energy -- it only labels the checkpoints, so a
    wrong value silently mislabels every row (and misnames every extracted PDB) rather than
    causing a visible failure.
    '''
    try:
        conn = sqlite3.connect( 'file:%s?mode=ro' % db3_file, uri = True )
    except sqlite3.Error:
        return None

    try:
        try:
            rows = conn.execute( 'SELECT specified_options, command_line FROM protocols' ).fetchall()
        except sqlite3.Error:
            return None

        strides = set()
        for row in rows:
            for field in row:
                if field:
                    strides.update( int(m) for m in _stride_re.findall( field ) )
    finally:
        conn.close()

    # More than one distinct value means the database mixes runs; we cannot label it correctly.
    if len( strides ) == 1:
        return strides.pop()
    return None
