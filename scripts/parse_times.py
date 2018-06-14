#!/usr/bin/env python

"""
Script to generate a timeline and stats per partition and worker from
controller log file:
# It uses input-statement-id to mark the beginning of an iteration.
# Total iteration time for ith iteration = T for the next stage with
input-statement-id minus T for the current stage with input-statement-id.
# Total compute time for ith iteration = sum of Ws starting after current
stage id for given input-statement-id minus current W.
# Output is 4 tables.
  1. partition_time: compute time for each partition for each iteration
  2. worker_time: average compute time per core for each worker for each
    iteration
  3. iteration_time: iteration time for each iteration
  4. marker: partition to worker mapping for each iteration
"""

import argparse
import pickle
import numpy as np
import os

def ParseArguments():
    parser = argparse.ArgumentParser(description="Process controller log file")
    parser.add_argument('--cores', type=int, required=True,
                        help="Cores per worker")
    parser.add_argument('--offset', type=int, required=True,
                        help="Offset for main workers")
    parser.add_argument('--input', type=str, required=True,
                        help="Input dir")
    parser.add_argument('--output', type=str, required=False, default='output_',
                        help="Output dir")
    parser.add_argument('--marker', type=str, required=False,
                        default='ReceiveGlobalDt',
                        help="Statement that marks beginning of an iteration")
    parser.add_argument('--count-also', action='store_true',
                        help="Also parse local fluid count from worker 0")
    parser.add_argument('--last-iter', type=int, required=False, default='-1',
                        help="Number of iterations")
    parser.add_argument('--scale', type=float, required=False, default='1.0',
                        help="Normalizing factor")
    return parser.parse_args()

"""
Stores and adds information to tables.
"""
class IterationData(object):
    def __init__(self, cores, num_partitions, num_workers, num_iterations):
        self.num_cores      = cores
        self.num_partitions = num_partitions
        self.num_workers    = num_workers
        self.num_iterations = num_iterations
        # Compute time for each partition
        self.partition_time = np.zeros(
            shape=(num_partitions, num_iterations), dtype=float)
        # Average per core compute time for each worker
        self.worker_time = np.zeros(
            shape=(num_workers, num_iterations), dtype=float)
        # Iteration time for each partition
        self.iteration_time = np.zeros(shape=(num_iterations), dtype=float)
        # Partition to worker map
        self.partition_to_worker = np.zeros(
            shape=(num_partitions, num_iterations), dtype=int) 

    """
    Save data.
    """
    def Save(self, out, scale):
        sum_time = self.partition_time.sum(axis=0)
        result = sum_time[0:-1].sum()
        print("Cumulative average time = %f" % (result))
        np.save(os.path.join(out, 'partition_time'), self.partition_time/scale)
        np.save(os.path.join(out, 'worker_time'), self.worker_time/scale)
        np.save(os.path.join(out, 'iteration_time'), self.iteration_time/scale)
        np.save(os.path.join(out, 'partition_to_worker'), self.partition_to_worker)
        parameters = {
            "num_cores"      : self.num_cores,
            "num_partitions" : self.num_partitions,
            "num_workers"    : self.num_workers,
            "num_iterations" : self.num_iterations
        }
        pickle.dump(parameters, open(os.path.join(out, 'parameters.txt'), 'wt'))


    """
    Record iteration time.
    """
    def RecordIterationTime(self, itr, time):
        self.iteration_time[itr] = time

    """
    Record compute time logged for a partition for an iteration.
    """
    def RecordComputeTime(self, partition_id, itr, compute_time):
        self.partition_time[partition_id, itr] = compute_time

    """
    Record worker over which an iteration ran.
    """
    def RecordWorker(self, partition_id, itr, worker_id):
        self.partition_to_worker[partition_id, itr] = worker_id

    """
    Compute average time per partition for partitions that ran on a worker,
    for every iteration.
    """
    def ComputeWorkerComputeTime(self):
        for w in range(self.num_workers):
            partitions = (self.partition_to_worker == w)
            partition_time = self.partition_time.copy()
            partition_time[~partitions] = 0
            self.worker_time[w,:] = partition_time.sum(axis=0)/self.num_cores

"""
Parses and stores iteration data.
"""
class LogParser(object):
    def __init__(self, offset, cores, controller_log_file_name, stats_file_name,
                 num_iterations):
        # Controller log file that dumps all of log, to determine the statement
        # ID for ReceiveGlobalDt.
        self.controller_log_file_name = controller_log_file_name
        # Statement IDs for the statement that marks start of an iteration.
        self.iteration_marker = set()
        # Stats file that controller dumps at the end of application run.
        self.stats_file_name = stats_file_name
        # Iteration data
        self.iteration_data = None
        # Variable group to partitions.
        self.variable_group_to_partitions = dict()
        self.variable_group_main = None  # main variable group
        self.num_partitions_main = None  # number of partitions fo main variable
        # Total number of workers, and cores per worker.
        self.offset = offset
        self.num_cores = cores
        self.num_workers = None
        # Number of iterations.
        self.num_iterations = num_iterations

    """
    Go over log file and get Statement IDs for given statement, to mark the
    beginning of an iteration.
    """
    def InitializeIterationMarkers(self, statement_name):
        with open(self.controller_log_file_name, 'r') as controller_log_file:
            for line in controller_log_file:
                tokens = line.split()
                idx = -1
                for i,t in enumerate(tokens):
                    if statement_name in t:
                        idx = i
                if idx != -1:
                    if idx < 2:
                        continue
                    subtokens = tokens[idx-2].strip(',').split('#')
                    if len(subtokens) !=2 or subtokens[0] != 'statement':
                        continue
                    self.iteration_marker.add(int(subtokens[1]))
        print("Markers", self.iteration_marker)

    """
    Determine the main variable group, number of partitions, number of workers,
    and number of iterations. Initialize iteration data.
    """
    def InitializeIterationData(self):
        worker_set = set()
        main_group = None
        total_iterations = 0
        # Make a pass over stats file, recording information about partitions,
        # iterations and workers.
        with open(self.stats_file_name, 'r') as stats_file:
            line_num = 0
            current_variable_group = None
            for line in stats_file:
                tokens = line.split()
                category = tokens[0]
                if category == 'P':
                    current_variable_group = int(tokens[2])
                    partition_id = int(tokens[3])
                    partitions = self.variable_group_to_partitions.setdefault(
                        current_variable_group, set())
                    partitions.add(partition_id)
                    worker_id = int(tokens[5])
                    worker_set.add(worker_id)
                elif category == 'T':
                    if current_variable_group == None:
                        continue
                    statement_id = int(tokens[2])
                    if statement_id in self.iteration_marker:
                        total_iterations += 1
                        if main_group is None:
                            main_group = current_variable_group
                            print("Main group is %d" % main_group)
                        else:
                            assert(main_group == current_variable_group)
                line_num += 1
                if line_num%(10**7) == 0:
                    print("Initializing, processed %d lines" % line_num)
        # Count total number of workers.
        self.num_workers = len(worker_set) - self.offset
        # Main variable group is the variable group with non-zero iterations.
        # Verify that this is the one with maximum number of partitions.
        max_num_partitions = max(
          len(partitions) for partitions in
          self.variable_group_to_partitions.values())
        assert(len(self.variable_group_to_partitions[main_group]) ==
               max_num_partitions)
        self.variable_group_main = main_group
        self.num_partitions_main = max_num_partitions
        # Total number of iterations equals number of iterations on main
        # variable group.
        print("Total iterations X %d partitions = %d" % (max_num_partitions, total_iterations))
        #assert(total_iterations % max_num_partitions == 0)
        # Number of iterations equals number of times marker appears for
        # an iteration, minus -- we're counting iteration time as time between
        # markers.
        if self.num_iterations == -1:
            self.num_iterations = total_iterations/max_num_partitions - 1
        print("Will use %d iterations" % (self.num_iterations))
        # Finally initialize iteration data.
        self.iteration_data = IterationData(
          self.num_cores, self.num_partitions_main,
          self.num_workers, self.num_iterations)
        print("Number of partitions", self.num_partitions_main)
        print("Number of workers", self.num_workers)
        print("Number of iterations", self.num_iterations)
        print("Main variable group", self.variable_group_main)

    """
    Now go over the log, and for every partition over main variable group,
    record total time and compute time for every iteration, and also, the
    worker over which the iteration ran (began).
    """
    def ParseIterationTimes(self):
        # Store time stamp for each iteration, to compute iteration time.
        iteration_timestamps = dict()  # Stage ID to list of time stamps
        # Store total busy time for each partition, iteration.
        compute_times = list()
        # Store worker for each partition, iteration.
        iteration_workers = list()
        with open(self.stats_file_name, 'r') as stats_file:
            recording = False
            current_partition = -1
            current_worker = -1
            for p in range(self.num_partitions_main):
                compute_times.append(dict())
                iteration_workers.append(dict())
            # Go over log file and record times for each iteration.
            line_num = 0
            for line in stats_file:
                tokens = line.split()
                category = tokens[0]
                if category == 'P':
                    # Record current partition information.
                    variable_group = int(tokens[2])
                    partition_id = int(tokens[3])
                    worker_id = int(tokens[5]) - self.offset
                    if variable_group == self.variable_group_main:
                        recording = True
                        current_partition = partition_id
                        current_worker = worker_id
                    else:
                        recording = False
                elif category == 'T' and recording:
                    # Record iteration time stamp.
                    stage_id = int(tokens[1])
                    statement_id = int(tokens[2])
                    time = float(tokens[3])
                    if statement_id not in self.iteration_marker:
                        continue
                    timestamps = iteration_timestamps.setdefault(
                      stage_id, list())
                    timestamps.append(time)
                    # Record worker.
                    iteration_workers[current_partition][stage_id] = current_worker
                elif category == 'W' and recording:
                    # Record busy time.
                    stage_id = int(tokens[1])
                    time = float(tokens[3])
                    compute_times[current_partition][stage_id] = time
                line_num += 1
                if line_num%(10**7) == 0:
                    print("Processed %d lines" % line_num)
        stage_ids = iteration_timestamps.keys()
        stage_ids.sort()
        num_stages = len(stage_ids)
        # Average the timestamps for each iteration.
        #for t in iteration_timestamps.values():
        #    assert(len(t) == self.num_partitions_main)
        for i in range(self.num_iterations+1):
            print(i, len(iteration_timestamps[stage_ids[i]]))
            assert(len(iteration_timestamps[stage_ids[i]]) == self.num_partitions_main)
        iteration_timestamps = {
            sid:sum(t)/self.num_partitions_main
            for sid, t in iteration_timestamps.iteritems() }
        # Generate a sorted list of stage ids that mark beginning of an
        # iteration.
        print("Number of stages recorded = %d" % (num_stages))
        # Record iteration times.
        for i in range(min(num_stages-1, self.num_iterations)):
            sid1 = stage_ids[i]
            sid2 = stage_ids[i+1]
            self.iteration_data.RecordIterationTime(
              i, iteration_timestamps[sid2] - iteration_timestamps[sid1])
        # Record worker and busy time for each partition and iteration.
        for p in range(self.num_partitions_main):
            all_stages = compute_times[p].keys()
            all_stages.sort()
            for i in range(min(num_stages-1, self.num_iterations)):
                sid1 = stage_ids[i]
                sid2 = stage_ids[i+1]
                self.iteration_data.RecordWorker(
                  p, i, iteration_workers[p][sid1])
                time = 0
                # Find the stage just greater than sid1
                for s, sid in enumerate(all_stages):
                    if sid > sid1:
                        break
                # Add busy times till stage = sid2
                time += compute_times[p][sid]
                while sid < sid2:
                    s += 1
                    sid = all_stages[s]
                    time += compute_times[p][sid]
                if i < self.num_iterations:
                    self.iteration_data.RecordComputeTime(
                      p, i, time)
                    assert(sid == sid2)
        # Finally computer worker compute times.
        self.iteration_data.ComputeWorkerComputeTime()

    """
    Save data.
    """
    def Save(self, prefix, scale):
        self.iteration_data.Save(prefix, scale)

    """
    Run parser.
    """
    def Run(self, statement_name):
        self.InitializeIterationMarkers(statement_name)
        self.InitializeIterationData()
        self.ParseIterationTimes()

if __name__ == '__main__':
    args = ParseArguments()
    stats_file = os.path.join(args.input, 'controller.txt')
    log_file = os.path.join(args.input, 'canary_controller.INFO')
    log_parser = LogParser(args.offset, args.cores, log_file, stats_file, args.last_iter)
    log_parser.Run(args.marker)
    log_parser.Save(args.output, args.scale)
