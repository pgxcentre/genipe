
# This file is part of gwip.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


import os
import sys
import shlex
import logging
import traceback
from os.path import isfile
from multiprocessing import Pool
from subprocess import Popen, PIPE
from tempfile import NamedTemporaryFile

from ..db import *
from ..error import ProgramError


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = "Copyright 2014, Beaulieu-Saucier Pharmacogenomics Centre"
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


__all__ = ["launch_tasks", ]


def launch_tasks(to_process, nb_threads, check_rc=True, hpc=False,
                 hpc_options=None, out_dir=None):
    """Executes commands."""
    # Do we need to check the return code?
    to_run = []
    for i in range(len(to_process)):
        assert "name" in to_process[i]
        assert "task_id" in to_process[i]
        assert "task_db" in to_process[i]
        assert "o_files" in to_process[i]

        task_name = to_process[i]["name"]
        task_id = to_process[i]["task_id"]
        db_name = to_process[i]["task_db"]
        o_files = to_process[i]["o_files"]

        # Checking if we need to run this task
        if check_task_completion(task_id, db_name):
            if _check_output_files(o_files):
                run_time = get_task_runtime(task_id, db_name)
                logging.info("Task '{}': already performed in {:,d} "
                             "seconds".format(task_name, run_time))
                continue

        # The name of the task
        task_id = to_process[i]["task_id"]

        # Some options to add
        to_process[i]["check_retcode"] = check_rc
        to_process[i]["out_dir"] = out_dir

        # Setting the task options
        if hpc:
            assert hpc_options is not None
            walltime = None
            nodes = None
            if task_id in hpc_options:
                if "walltime" in hpc_options[task_id]:
                    walltime = hpc_options[task_id]["walltime"]
                if "nodes" in hpc_options[task_id]:
                    nodes = hpc_options[task_id]["nodes"]
            to_process[i]["walltime"] = walltime
            to_process[i]["nodes"] = nodes

        # Adding to list to run
        to_run.append(to_process[i])

    # The execution command
    execute_func = _execute_command
    if hpc:
        execute_func = _execute_command_drmaa

    # Launching the command
    if nb_threads > 1:
        # Running all the processes
        pool = Pool(processes=nb_threads)
        results = None

        try:
            results = pool.map(execute_func, to_run)

        except Exception as e:
            pool.terminate()
            traceback.print_exc(file=sys.stdout)
            raise

        finally:
            # Closing the pool
            pool.close()

        # Checking the results
        problems = []
        for result in results:
            if not result[0]:
                problems.append(result[1])
                logging.error("Task '{}': did not finish...".format(result[1]))
            else:
                logging.info("Task '{}': {} in {:,d} "
                             "seconds".format(result[1], result[2], result[3]))
        if len(problems) > 0:
            raise ProgramError("the following task did not work: " +
                               repr(problems))

    else:
        for data in to_run:
            logging.info("Executing {}".format(data["name"]))
            result = execute_func(data)
            if result[0]:
                logging.info("Task '{}': {} in {:,d} "
                             "seconds".format(result[1], result[2], result[3]))
            else:
                raise ProgramError("problem executing {}".format(data["name"]))


def _check_output_files(o_files):
    """Check that the files exist."""
    for filename in o_files:
        if filename.endswith(".impute2"):
            # IMPUTE2 files might be gzipped
            if not (isfile(filename) or isfile(filename + ".gz")):
                return False

        elif not isfile(filename):
            return False

    return True


def _execute_command(command_info):
    """Executes a single command."""
    # Some assertions
    assert "task_id" in command_info
    assert "name" in command_info
    assert "command" in command_info
    assert "check_retcode" in command_info
    assert "task_db" in command_info
    assert "o_files" in command_info

    # Getting the command's information
    name = command_info["name"]
    command = command_info["command"]
    check_rc = command_info["check_retcode"]
    task_id = command_info["task_id"]
    db_name = command_info["task_db"]

    logging.debug("Checking status for '{}'".format(task_id))
    # Checking if the command was completed
    if check_task_completion(task_id, db_name):
        if _check_output_files(command_info["o_files"]):
            logging.debug("'{}' completed".format(task_id))
            runtime = get_task_runtime(task_id, db_name)
            return True, name, "already performed", runtime
        else:
            logging.debug("'{}' problem with output files".format(task_id))
    logging.debug("'{}' to run".format(task_id))

    # Creating a new entry in the database
    create_task_entry(task_id, db_name)

    # Launching the command
    proc = Popen(command, stdout=PIPE, stderr=PIPE)
    logging.debug("'{}' finished".format(task_id))

    # Waiting for the process to terminate
    outs, errs = proc.communicate()
    rc = proc.returncode
    if check_rc and rc != 0:
        # There was a problem...
        logging.debug("'{}' exit status problem".format(task_id))
        return False, name, "problem", None

    # The task was performed correctly, so we update to completed
    mark_task_completed(task_id, db_name)

    # Everything when well
    logging.debug("'{}' everything was fine".format(task_id))
    return True, name, "performed", get_task_runtime(task_id, db_name)


def _execute_command_drmaa(command_info):
    """Executes a command using DRMAA (usually on a HPC)."""
    import drmaa

    # Some assertions
    assert "out_dir" in command_info
    assert "command" in command_info
    assert "walltime" in command_info
    assert "nodes" in command_info
    assert "task_id" in command_info
    assert "task_db" in command_info
    assert "name" in command_info
    assert "check_retcode" in command_info
    assert "o_files" in command_info

    # Getting the command's information
    name = command_info["name"]
    command = command_info["command"]
    task_id = command_info["task_id"]
    db_name = command_info["task_db"]
    out_dir = command_info["out_dir"]
    check_rc = command_info["check_retcode"]

    # Checking if the command was completed
    logging.debug("Checking status for '{}'".format(task_id))
    if check_task_completion(task_id, db_name):
        if _check_output_files(command_info["o_files"]):
            logging.debug("'{}' completed".format(task_id))
            runtime = get_task_runtime(task_id, db_name)
            return True, name, "already performed", runtime
        else:
            logging.debug("'{}' problem with output files".format(task_id))
    else:
        logging.debug("'{}' to run because not completed".format(task_id))
    logging.debug("'{}' to run".format(task_id))

    # Creating the script
    tmp_file = NamedTemporaryFile(mode="w", suffix="_execute.sh", delete=False,
                                  dir=out_dir)
    print("#!/usr/bin/env bash", file=tmp_file)
    print(command[0], end=" ", file=tmp_file)
    for chunk in command[1:]:
        print(shlex.quote(chunk), end=" ", file=tmp_file)
    print("", file=tmp_file)
    tmp_file.close()

    # Making the script executable
    os.chmod(tmp_file.name, 0o755)

    # Creating the job template
    s = drmaa.Session()
    s.initialize()

    # Creating the job template
    job = s.createJobTemplate()
    job.remoteCommand = tmp_file.name
    job.jobName = "_{}".format(task_id)
    job.workingDirectory = os.getcwd()
    job.jobEnvironment = os.environ
    if command_info["walltime"] is not None:
        job.hardWallclockTimeLimit = command_info["walltime"]
    if command_info["nodes"] is not None:
        job.nativeSpecification = command_info["nodes"]

    # Creating a new entry in the database
    create_task_entry(task_id, db_name)

    # Running the job
    job_id = s.runJob(job)

    # Waiting for the job
    ret_val = s.wait(job_id, drmaa.Session.TIMEOUT_WAIT_FOREVER)
    logging.debug("'{}' finished".format(task_id))

    # Deleting the job
    s.deleteJobTemplate(job)

    # Closing the connection
    s.exit()

    # Removing the temporary file
    os.remove(tmp_file.name)

    # Checking the task's return values
    if ret_val.hasCoreDump or ret_val.wasAborted or ret_val.hasSignal:
        logging.debug("'{}' problems ({}, {}, {})".format(
            task_id,
            ret_val.hasCoreDump,
            ret_val.wasAborted,
            ret_val.hasSignal,
        ))
        return False, name, "problem", None
    if check_rc and ret_val.exitStatus != 0:
        logging.debug("'{}' exit status problem".format(task_id))
        return False, name, "problem", None

    # Getting the time
    launch_time = float(ret_val.resourceUsage["submission_time"])
    start_time = float(ret_val.resourceUsage["start_time"])
    end_time = float(ret_val.resourceUsage["end_time"])

    # The task was performed correctly, so we update to completed
    mark_drmaa_task_completed(task_id, launch_time, start_time, end_time,
                              db_name)

    # Everything when well
    logging.debug("'{}' everything was fine".format(task_id))
    return True, name, "performed", get_task_runtime(task_id, db_name)
