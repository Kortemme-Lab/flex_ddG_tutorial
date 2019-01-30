from __future__ import print_function

import datetime
import math
import sys
import collections
import threading

# The Reporter class is useful for printing output for tasks which will take a long time
# It can even predict a finish time!

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
            return time.time() - self.start
