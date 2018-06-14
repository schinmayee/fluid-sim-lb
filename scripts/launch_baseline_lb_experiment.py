#!/usr/bin/env python

import argparse
import commands
import os
import subprocess
import time

import config

def parse_arguments():
    parser = argparse.ArgumentParser(description="Process controller log file")
    parser.add_argument('--ssh-key', type=str, required=True,
                        help="SSH Key to use")
    parser.add_argument('--project', type=str, default="",
                        help="Project ID")
    parser.add_argument('--config', type=str, required=True,
                        help="Experiment configuration file")
    return parser.parse_args()

def launch_experiment():
    launch_args = parse_arguments()

    project = ""
    if launch_args.project:
        project = "--project %s" % launch_args.project

    args = config.ParseConfiguration(launch_args.config)

    # Get instance names and select a controller instance and worker instances.
    response = commands.getstatusoutput(
      "gcloud compute instance-groups list-instances %s --zone %s %s" %
      (args.group, args.zone, project))
    assert(response[0] == 0)
    print(response[1])
    instances = list()
    instances_string = response[1].split('\n')[1:]
    for line in instances_string:
        name = line.split()[0]
        instances.append(name)
    num_instances = len(instances)
    controller = instances[0]
    workers = instances[1:num_instances]

    # Get ip address for controller and worker instances.
    response = commands.getstatusoutput("gcloud compute instances list %s" %
                                        project)
    assert(response[0] == 0)
    print(response[1])
    controller_internal_ip = None
    controller_external_ip = None
    workers_internal_ip = list()
    workers_external_ip = list()
    instances_string = response[1].split('\n')
    header = instances_string[0].split()
    assert('INTERNAL_IP' in header)
    assert('EXTERNAL_IP' in header)
    i_internal = header.index('INTERNAL_IP')  # hack, preemptible is empty
    i_external = header.index('EXTERNAL_IP')  # hack, preemptible is empty
    header_len = len(header)
    i_internal = i_internal - header_len
    i_external = i_external - header_len
    print(i_internal, i_external)
    for line in instances_string[1:]:
        tokens = line.split()
        if tokens[0] == controller:
            controller_internal_ip = tokens[i_internal]
            controller_external_ip = tokens[i_external]
        if tokens[0] in workers:
            workers_internal_ip.append(tokens[i_internal])
            workers_external_ip.append(tokens[i_external])
    assert(len(workers) == len(workers_internal_ip))
    assert(controller_internal_ip is not None)

    # Launch controller.
    print("Launching controller on %s" % controller)
    print(commands.getstatusoutput(
      "gcloud compute scp --zone %s %s ./launch_controller.sh %s:~/" %
      (args.zone, project, controller)))
    subprocess.call(
      "ssh -o StrictHostKeyChecking=no -i %s -f %s ./launch_controller.sh" %
      (launch_args.ssh_key, controller_external_ip), shell=True)
    time.sleep(2)

    # Launch workers.
    for i, w in enumerate(workers):
        print("Launching worker on %s" % w)
        print(commands.getstatusoutput(
          "gcloud compute scp --zone %s %s ./launch_worker.sh %s:~/" %
          (args.zone, project, w)))
        print(commands.getstatusoutput(
          "gcloud compute scp --zone %s %s %s %s:~/init_file" %
          (args.zone, project, args.init_file, w)))
        subprocess.call(
          "ssh -o StrictHostKeyChecking=no -i %s -f %s ./launch_worker.sh %s %d" %
          (launch_args.ssh_key, workers_external_ip[i],
           controller_internal_ip, args.worker_threads), shell=True)
    time.sleep(2)

    # Launch application.
    print("Launching application on %s" % controller)
    print(commands.getstatusoutput(
      "gcloud compute scp --zone %s %s ./launch_baseline_lb_application.sh %s:~/" %
      (args.zone, project, controller)))
    print(commands.getstatusoutput(
      "gcloud compute scp --zone %s %s %s %s:~/init_file" %
      (args.zone, project, args.init_file, controller)))
    params  = "PLACEMENT=%s " % args.placement
    params += "NUM_WORKERS=%d " % len(workers)
    params += "INIT=%d " % args.init
    params += "FRAMES=%d " % args.frames
    params += "NGRID=%d " % args.ngrid
    params += "X=%d " % args.x
    params += "Y=%d " % args.y
    params += "Z=%d " % args.z
    params += "P=%d " % args.p
    params += "AX=%d " % args.ax
    params += "AY=%d " % args.ay
    params += "AZ=%d " % args.az
    params += "F=%f " % args.current_load_factor
    params += "HORIZON=%d " % args.horizon
    params += "FLIP=%f " % args.flip
    params += "CFL=%d " % args.cfl
    params += "RATE=%d " % args.frame_rate
    params += "SOLVER_ITERS=%d " % args.solver_iterations
    if args.use_geometric:
        params += "GEOMETRIC=true "
    print("Using launch parameters:")
    print(params)
    subprocess.call(
      "ssh -o StrictHostKeyChecking=no -i %s -f %s '%s ./launch_baseline_lb_application.sh'" %
      (launch_args.ssh_key, controller_external_ip, params), shell=True)

if __name__ == '__main__':
    launch_experiment()
