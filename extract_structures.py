#!/usr/bin/env python3

import os
import sys
import argparse
import functools
import subprocess
import re
import shutil
import datetime
import math
import collections
import threading

import flex_ddg_db3

use_multiprocessing = False
if use_multiprocessing:
    import multiprocessing

# The Reporter class is useful for printing output for tasks which will take a long time
# Really, you should just use tqdm now, but I used this before I knew about tqdm and it removes a dependency

# Time in seconds function
# Converts datetime timedelta object to number of seconds
def ts(td):
    return (td.microseconds + (td.seconds + td.days * 24 * 3600) * 1e6) / 1e6

def mean(l):
    # Not using numpy mean to avoid dependency
    return float( sum(l) ) / float( len(l) )

class Reporter:
    def __init__( self, task, entries = 'files', print_output = True, eol_char = '\r' ):
        self._lock = threading.Lock()
        self.print_output = print_output
        self.start = datetime.datetime.now()
        self.entries = entries
        self.lastreport = self.start
        self.task = task
        self.report_interval = datetime.timedelta( seconds = 1 ) # Interval to print progress
        self.n = 0
        self.completion_time = None
        if self.print_output:
            print('\nStarting ' + task)
        self.total_count = None # Total tasks to be processed
        self.maximum_output_string_length = 0
        self.rolling_est_total_time = collections.deque( maxlen = 50 )
        self.kv_callback_results = {}
        self.list_results = []
        self.eol_char = eol_char

    def set_total_count(self, x):
        self.total_count = x
        self.rolling_est_total_time = collections.deque( maxlen = max(1, int( .05 * x )) )

    def decrement_total_count(self):
        if self.total_count:
            self.total_count -= 1

    def report(self, n):
        with self._lock:
            self.n = n
            time_now = datetime.datetime.now()
            if self.print_output and self.lastreport < (time_now - self.report_interval):
                self.lastreport = time_now
                if self.total_count:
                    percent_done = float(self.n) / float(self.total_count)
                    est_total_time_seconds = ts(time_now - self.start) * (1.0 / percent_done)
                    self.rolling_est_total_time.append( est_total_time_seconds )
                    est_total_time = datetime.timedelta( seconds = mean(self.rolling_est_total_time) )
                    time_remaining = est_total_time - (time_now - self.start)
                    eta = time_now + time_remaining
                    time_remaining_str = 'ETA: %s Est. time remaining: ' % eta.strftime("%Y-%m-%d %H:%M:%S")

                    time_remaining_str += str( datetime.timedelta( seconds = int(ts(time_remaining)) ) )

                    output_string = "  Processed: %d %s (%.1f%%) %s" % (n, self.entries, percent_done*100.0, time_remaining_str)
                else:
                    output_string = "  Processed: %d %s" % (n, self.entries)

                output_string += self.eol_char

                if len(output_string) > self.maximum_output_string_length:
                    self.maximum_output_string_length = len(output_string)
                elif len(output_string) < self.maximum_output_string_length:
                    output_string = output_string.ljust(self.maximum_output_string_length)
                sys.stdout.write( output_string )
                sys.stdout.flush()

    def increment_report(self):
        self.report(self.n + 1)

    def increment_report_callback(self, cb_value):
        self.increment_report()

    def increment_report_keyval_callback(self, kv_pair):
        key, value = kv_pair
        self.kv_callback_results[key] = value
        self.increment_report()

    def increment_report_list_callback(self, new_list_items):
        self.list_results.extend(new_list_items)
        self.increment_report()

    def decrement_report(self):
        self.report(self.n - 1)

    def add_to_report(self, x):
        self.report(self.n + x)

    def done(self):
        self.completion_time = datetime.datetime.now()
        if self.print_output:
            print('Done %s, processed %d %s, took %s\n' % (self.task, self.n, self.entries, self.completion_time-self.start))

    def elapsed_time(self):
        if self.completion_time:
            return self.completion_time - self.start
        else:
            return datetime.datetime.now() - self.start


struct_db3_file = 'struct.db3'

# Extraction uses the score_jd2 binary, not rosetta_scripts. It is built alongside
# rosetta_scripts by the standard Rosetta build, but is a separate executable.
score_jd2_path = os.path.expanduser( '~/rosetta/source/bin/score_jd2' )

# Only a fallback. Extracted structures are named by how many backrub steps produced them, so
# the stride each run used is read back out of its own struct.db3. This value is used only when
# the database does not record it, and a warning is printed.
default_trajectory_stride = 5

def resolve_trajectory_stride( struct_db, stride_override = None ):
    '''Stride to name this database's extracted PDBs with, preferring what the run recorded.'''
    if stride_override is not None:
        return stride_override

    stride = flex_ddg_db3.trajectory_stride_from_db3( struct_db )
    if stride is not None:
        return stride

    print( 'WARNING: %s does not record backrub_trajectory_stride; assuming %d.' % (
        struct_db, default_trajectory_stride ) )
    print( '         If the run used a different stride, pass --stride, or the extracted PDBs' )
    print( '         will be named with the wrong backrub step counts.' )
    return default_trajectory_stride

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
        os.path.abspath( score_jd2_path ),
        '-inout:dbms:database_name', struct_db3_file,
        '-in:use_database',
        '-out:pdb',
    ]

    working_directory = os.path.dirname( struct_db )
    rosetta_outfile_path = os.path.join(working_directory, 'structure_output.txt' )
    if not use_multiprocessing:
        print(rosetta_outfile_path)
    rosetta_outfile = open( rosetta_outfile_path, 'w')
    if not use_multiprocessing:
        print( ' '.join( args ) )
    # No shell: joining the arguments into a string breaks as soon as a path contains a space.
    rosetta_process = subprocess.Popen(
        args,
        stdout=rosetta_outfile, stderr=subprocess.STDOUT, close_fds = True, cwd = working_directory,
    )
    return_code = rosetta_process.wait()
    rosetta_outfile.close()

    if return_code == 0:
        os.remove( rosetta_outfile_path )
    else:
        print( 'ERROR: score_jd2 failed on %s (exit %d) -- see %s' % (
            struct_db, return_code, rosetta_outfile_path ) )
        return return_code

    if rename_function != None:
        for path in [ os.path.join( working_directory, x ) for x in os.listdir( working_directory ) ]:
            m = re.match( r'(\d+)_0001\.pdb$', os.path.basename(path) )
            if m:
                dest_path = os.path.join( working_directory, rename_function( int(m.group(1)) ) )
                shutil.move( path, dest_path )

    return return_code

def flex_ddG_rename(struct_id, trajectory_stride):
    steps = [
        'backrub',
        'wt',
        'mut',
    ]

    return '%s_%05d.pdb' % ( steps[ (struct_id-1) % len(steps) ], (((struct_id-1) // len(steps)) + 1) * trajectory_stride )

def main( input_dir, stride_override = None ):
    struct_dbs = recursive_find_struct_dbs( input_dir )
    print( 'Found {:d} structure database files to extract'.format( len(struct_dbs) ) )

    if use_multiprocessing:
        pool = multiprocessing.Pool()
    r = Reporter('extracting structure database files', entries = '.db3 files')
    r.set_total_count( len(struct_dbs) )

    for struct_db in struct_dbs:
        # Each database is named using the stride its own run was launched with.
        # functools.partial rather than a lambda, so that this stays picklable for the
        # multiprocessing path below.
        stride = resolve_trajectory_stride( struct_db, stride_override )
        rename_function = functools.partial( flex_ddG_rename, trajectory_stride = stride )
        if use_multiprocessing:
            pool.apply_async(
                extract_structures,
                args = (struct_db,),
                kwds = {'rename_function' : rename_function},
                callback = r.increment_report_callback
            )
        else:
            r.increment_report_callback(
                extract_structures( struct_db, rename_function = rename_function )
            )

    if use_multiprocessing:
        pool.close()
        pool.join()
    r.done()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = 'Extract PDBs from the struct.db3 files under a flex ddG output folder.' )
    parser.add_argument( 'output_folders', nargs = '+', help = 'flex ddG output folder(s)' )
    parser.add_argument( '--stride', type = int, default = None,
                         help = 'override backrub_trajectory_stride instead of reading it from'
                                ' each struct.db3. Affects extracted PDB names only.' )
    parsed_args = parser.parse_args()

    if not os.path.isfile( score_jd2_path ):
        print( 'ERROR: "score_jd2_path" variable must be set to the location of the "score_jd2" binary executable' )
        print( 'This file might look something like: "score_jd2.linuxgccrelease"' )
        print( 'Note that this is a different executable from the "rosetta_scripts" binary used to run flex ddG' )
        raise Exception( 'score_jd2 missing' )

    for x in parsed_args.output_folders:
        if os.path.isdir(x):
            main( x, parsed_args.stride )
        else:
            print( 'ERROR: %s is not a valid directory' % x )
