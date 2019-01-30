#!/usr/bin/env python3

import os
import sys
import subprocess
import re
import shutil

use_multiprocessing = True
if use_multiprocessing:
    import multiprocessing

from Reporter import Reporter

struct_db3_file = 'struct.db3'

# Important - to correctly name extracted structures by the stride, this trajectory_stride must be used
trajectory_stride = 5

def recursive_find_struct_dbs( input_dir ):
    return_list = []

    for path in [os.path.join(input_dir, x) for x in os.listdir( input_dir )]:
        if os.path.isdir( path ):
            return_list.extend( recursive_find_struct_dbs( path ) )
        elif os.path.isfile( path ) and os.path.basename( path ) == struct_db3_file:
            return_list.append( path )

    return return_list

def extract_structures( struct_db, rename_function = None ):
    args = [
        os.path.expanduser( '~/rosetta/source/bin/score_jd2' ),
        '-inout:dbms:database_name', struct_db3_file,
        '-in:use_database',
        '-out:pdb',
    ]

    working_directory = os.path.dirname( struct_db )
    rosetta_outfile_path = os.path.join(working_directory, 'structure_output.txt' )
    rosetta_outfile = open( rosetta_outfile_path, 'w')
    if not use_multiprocessing:
        print( ' '.join( args ) )
    rosetta_process = subprocess.Popen(
        ' '.join(args),
        stdout=rosetta_outfile, stderr=subprocess.STDOUT, close_fds = True, cwd = working_directory, shell = True,
    )
    return_code = rosetta_process.wait()
    rosetta_outfile.close()

    if return_code == 0:
        os.remove( rosetta_outfile_path )

    if rename_function != None:
        output_pdbs = []
        for path in [ os.path.join( working_directory, x ) for x in os.listdir( working_directory ) ]:
            m = re.match( '(\d+)_0001.pdb', os.path.basename(path) )
            if m:
                dest_path = os.path.join( working_directory, rename_function( int(m.group(1)) ) )
                shutil.move( path, dest_path )

    return return_code

def flex_ddG_rename(struct_id):
    steps = [
        'backrub',
        'wt',
        'mut',
    ]

    return '%s_%05d.pdb' % ( steps[ (struct_id-1) % len(steps) ], (((struct_id-1) // len(steps)) + 1) * trajectory_stride )

def main( input_dir ):
    struct_dbs = recursive_find_struct_dbs( input_dir )
    print( 'Found {:d} structure database files to extract'.format( len(struct_dbs) ) )

    if use_multiprocessing:
        pool = multiprocessing.Pool()
    r = Reporter('extracting structure database files', entries = '.db3 files')
    r.set_total_count( len(struct_dbs) )

    for struct_db in struct_dbs:
        if use_multiprocessing:
            pool.apply_async(
                extract_structures,
                args = (struct_db,),
                kwds = {'rename_function' : flex_ddG_rename},
                callback = r.increment_report_callback
            )
        else:
            r.increment_report_callback(
                extract_structures( struct_db, rename_function = flex_ddG_rename )
            )

    if use_multiprocessing:
        pool.close()
        pool.join()
    r.done()

if __name__ == '__main__':
    for x in sys.argv[1:]:
        if os.path.isdir(x):
            main( x )
        else:
            print( 'ERROR: %s is not a valid directory' % x )
