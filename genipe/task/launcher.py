
# This file is part of genipe.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


import os
import re
import sys
import time
import shlex
import logging
import traceback
from os.path import isfile
from datetime import datetime
from multiprocessing import Pool
from subprocess import Popen, PIPE
from tempfile import NamedTemporaryFile
from multiprocessing.pool import ThreadPool

from ..db import utils as db
from ..error import GenipeError


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = "Copyright 2014, Beaulieu-Saucier Pharmacogenomics Centre"
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


__all__ = ["launch_tasks", ]


def launch_tasks(to_process, nb_threads, check_rc=True, hpc=False,
                 hpc_options=None, out_dir=None, preamble=""):
    """Executes commands.

    Args:
        to_process (list): a list of tasks to process
        nb_threads (int): the number of processes that is required
        check_rc (bool): whether or not to check the return code of the task
        hpc (bool): whether or not to execute the tasks on a cluster (DRMAA)
        hpc_options (dict): the DRMAA options
        out_dir (str): the output directory
        preamble (str): the script preamble (for DRMAA)

    """
    # Do we need a DRMAA session?
    drmaa_session = None
    if hpc:
        import drmaa
        drmaa_session = drmaa.Session()
        drmaa_session.initialize()

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
        if db.check_task_completion(task_id, db_name):
            if _check_output_files(o_files, task_id):
                run_time = db.get_task_runtime(task_id, db_name)
                logging.info("Task '{}': already performed in {:,d} "
                             "seconds".format(task_name, run_time))
                continue

            else:
                # The DB said the task was completed, but there is a missing
                # output files. Setting this task completion to '0'
                db.mark_task_incomplete(task_id, db_name)

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
            to_process[i]["preamble"] = preamble
            to_process[i]["drmaa_session"] = drmaa_session

        # Adding to list to run
        to_run.append(to_process[i])

    # The execution function and the execution type
    execute_func = _execute_command
    GenipePool = Pool
    if hpc:
        execute_func = _execute_command_drmaa
        GenipePool = ThreadPool

    # Launching the command
    if nb_threads > 1:
        # Running all the processes
        pool = GenipePool(processes=nb_threads)
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

            if drmaa_session is not None:
                drmaa_session.exit()

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
            raise GenipeError("the following task did not work: " +
                              repr(problems))

    else:
        try:
            for data in to_run:
                logging.info("Executing {}".format(data["name"]))
                result = execute_func(data)
                if result[0]:
                    logging.info("Task '{}': {} in {:,d} " "seconds".format(
                        result[1],
                        result[2],
                        result[3],
                    ))
                else:
                    raise GenipeError(
                        "problem executing {}".format(data["name"])
                    )

        finally:
            if drmaa_session is not None:
                drmaa_session.exit()


def _check_output_files(o_files, task):
    """Check that the files exist.

    Args:
        o_files (list): the list of files got check
        task (str): the name of the task

    Returns:
        bool: ``True`` if all files exist, ``False`` otherwise

    If the file to check is an impute2 file, and that this file is missing, we
    check for further statistics using the :py:func:`_check_impute2_file`.

    Note
    ----
        If the file name ends with ``.impute2`` and the file doesn't exist, we
        look for the compressed file (``.impute2.gz``) instead.

    """
    for filename in o_files:
        if filename.endswith(".impute2"):
            # IMPUTE2 files might be gzipped
            if not (isfile(filename) or isfile(filename + ".gz")):
                if not _check_impute2_file(filename, task):
                    return False

        elif filename.endswith(".snp.strand"):
            # SHAPEIT alignment file might not exits
            if not isfile(filename):
                if not _check_shapeit_align_file(filename, task):
                    return False

        elif not isfile(filename):
            return False

    return True


def _check_shapeit_failed_rc(fn, task=None):
    """Checks the log to explain a failure return code.

    Args:
        fn (str): the name of the file to check
        task (str): the name of the task

    Returns:
        bool: ``True`` if everything is norma, ``False`` otherwise

    This function looks for a known error message in the log file. If the
    message ``ERROR: Reference and Main panels are not well aligned:`` appears
    in the log file, then it's normal that the job failed.

    """
    # The name of the log file
    log_fn = fn.replace(".snp.strand", "") + ".log"
    if not os.path.isfile(log_fn):
        # The log file doesn't exist...
        return False

    # Reading the content of the file
    log = None
    with open(log_fn, "r") as i_file:
        log = i_file.read()

    # Checking that the SNPs were read in the legend file (the strand error
    # message is just after this notice).
    match = re.search(
        r"\sERROR: Reference and Main panels are not well aligned:\n",
        log,
    )

    if match is None:
        return False
    return True


def _check_shapeit_align_file(fn, task=None):
    """Checks the log to explain the absence of an .snp.strand file.

    Args:
        fn (str): the name of the file to check
        task (str): the name of the task

    Returns:
        bool: ``True`` if everything is normal, ``False`` otherwise.

    This function looks for known message in the log file. If the SNPs were
    read from the legend file and the haplotypes were read from the hap file,
    then there were no SNPs flip issue.

    """
    # The name of the log file
    log_fn = fn.replace(".snp.strand", "") + ".log"
    if not os.path.isfile(log_fn):
        # The log file doesn't exist...
        return False

    # Reading the content of the file
    log = None
    with open(log_fn, "r") as i_file:
        log = i_file.read()

    # Checking that the SNPs were read in the legend file (the strand error
    # message is just after this notice).
    match = re.search(r"\sReading SNPs in \[.+\]\n", log)
    if match is None:
        return False

    # Checking if the step after was run (i.e. meaning there were no strand
    # issue)
    match = re.search(r"\sReading reference haplotypes in \[.+\]\n", log)
    if match is None:
        return False

    # We are here, so it is normal that no *.snp.strand file exists
    if task:
        logging.info("{}: there are no flip issue".format(task))
    return True


def _check_impute2_file(fn, task=None):
    """Checks the summary to explain the absence of an .impute2 file.

    Args:
        fn (str): the name of the file to check
        task (str): the name of the task

    Returns:
        bool: ``True`` if everything is normal, ``False`` otherwise.

    This function looks for known message in the summary file. Three possible
    ways that an impute2 file is missing:

    1. there are no SNPs in the imputation interval;
    2. there are no type 2 SNPs after applying the settings;
    3. there are no SNPs for output.

    """
    # The name of the summary file
    summary_fn = fn + "_summary"
    if not os.path.isfile(summary_fn):
        # The summary file doesn't exist...
        return False

    # Reading the file content
    summary = None
    with open(summary_fn, "r") as i_file:
        summary = i_file.read()

    # Checking if there are no SNPs in the imputation interval?
    match = re.search(
        r"\sThere are no SNPs in the imputation interval, so there is "
        "nothing for IMPUTE2 to analyze; the program will quit now.",
        summary,
    )
    if match is not None:
        if task:
            logging.warning("{}: there are no SNPs in the imputation "
                            "interval".format(task))
        return True

    # Checking if there are not type 2 SNPs
    match = re.search(
        r"\sERROR: There are no type 2 SNPs after applying the command-line "
        "settings for this run, which makes it impossible to perform "
        "imputation.",
        summary,
    )
    if match is not None:
        if task:
            logging.warning("{}: there are no type 2 SNPs for this "
                            "run".format(task))
        return True

    # Checking if there are no output SNPs
    match = re.search(
        r"\sYour current command-line settings imply that there will not be "
        "any SNPs in the output file, so IMPUTE2 will not perform any "
        "analysis or print output files.",
        summary,
    )
    if match is not None:
        if task:
            logging.warning("{}: no SNPs in the output file".format(task))
        return True

    # If attained, there is a problem
    return False


def _execute_command(command_info):
    """Executes a single command.

    Args:
        command_info (dict): information about the command

    Returns:
        tuple: a tuple containing 4 entries: whether the task completed (bool),
               the name of the task (str), the status of the run (str) and the
               execution time in seconds (int)

    """
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
    if db.check_task_completion(task_id, db_name):
        if _check_output_files(command_info["o_files"], task_id):
            logging.debug("'{}' completed".format(task_id))
            runtime = db.get_task_runtime(task_id, db_name)
            return True, name, "already performed", runtime
        else:
            logging.debug("'{}' problem with output files".format(task_id))
            db.mark_task_incomplete(task_id, db_name)
    logging.debug("'{}' to run".format(task_id))

    # Creating a new entry in the database
    db.create_task_entry(task_id, db_name)

    # Launching the command
    proc = Popen(command, stdout=PIPE, stderr=PIPE)
    logging.debug("'{}' finished".format(task_id))

    # Waiting for the process to terminate
    outs, errs = proc.communicate()
    rc = proc.returncode
    if check_rc and rc != 0:
        if task_id.startswith("impute2"):
            # Task is IMPUTE2, and it might be normal according to message in
            # the summary file
            impute2_fn = None
            for fn in command_info["o_files"]:
                if fn.endswith(".impute2"):
                    impute2_fn = fn
                    break
            if not _check_impute2_file(impute2_fn):
                logging.debug("'{}' exit status problem".format(task_id))
                return False, name, "problem", None

        elif task_id.startswith("shapeit_check"):
            shapeit_fn = None
            for fn in command_info["o_files"]:
                if fn.endswith(".alignments.snp.strand"):
                    shapeit_fn = fn
                    break
            if not _check_shapeit_failed_rc(shapeit_fn):
                logging.debug("'{}' exit status problem".format(task_id))
                return False, name, "problem", None

        else:
            # There was a problem...
            logging.debug("'{}' exit status problem".format(task_id))
            return False, name, "problem", None

    # Checking all the required files were generated
    if not _check_output_files(command_info["o_files"], task_id):
        logging.debug("'{}' exit status problem".format(task_id))
        return False, name, "problem", None

    # The task was performed correctly, so we update to completed
    db.mark_task_completed(task_id, db_name)

    # Everything when well
    logging.debug("'{}' everything was fine".format(task_id))
    return True, name, "performed", db.get_task_runtime(task_id, db_name)


def _execute_command_drmaa(command_info):
    """Executes a command using DRMAA (usually on a HPC).

    Args:
        command_info (dict): information about the command

    Returns:
        tuple: a tuple containing 4 entries: whether the task completed (bool),
               the name of the task (str), the status of the run (str) and the
               execution time in seconds (int)

    Note
    ----
        The preamble (if required) is inserted between the shebang line and the
        actual command.

    """
    # Some import
    from drmaa import Session, JobControlAction

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
    assert "preamble" in command_info
    assert "drmaa_session" in command_info

    # Getting the command's information
    name = command_info["name"]
    command = command_info["command"]
    task_id = command_info["task_id"]
    db_name = command_info["task_db"]
    out_dir = command_info["out_dir"]
    check_rc = command_info["check_retcode"]
    preamble = command_info["preamble"]
    drmaa_session = command_info["drmaa_session"]

    # Checking if the command was completed
    logging.debug("Checking status for '{}'".format(task_id))
    if db.check_task_completion(task_id, db_name):
        if _check_output_files(command_info["o_files"], task_id):
            logging.debug("'{}' completed".format(task_id))
            runtime = db.get_task_runtime(task_id, db_name)
            return True, name, "already performed", runtime
        else:
            logging.debug("'{}' problem with output files".format(task_id))
            db.mark_task_incomplete(task_id, db_name)
    else:
        logging.debug("'{}' to run because not completed".format(task_id))
    logging.debug("'{}' to run".format(task_id))

    # Creating the script
    tmp_file = NamedTemporaryFile(mode="w", suffix="_execute.sh", delete=False,
                                  dir=out_dir)

    # Writing the shebang
    print("#!/usr/bin/env bash", file=tmp_file)

    # Writing the preamble
    print(preamble, file=tmp_file)

    # Writing the command
    print(command[0], end=" ", file=tmp_file)
    for chunk in command[1:]:
        print(shlex.quote(chunk), end=" ", file=tmp_file)
    print("", file=tmp_file)

    # Closing the temporary file
    tmp_file.close()

    # Making the script executable
    os.chmod(tmp_file.name, 0o755)

    # Creating the job template
    job = drmaa_session.createJobTemplate()
    job.remoteCommand = tmp_file.name
    job.jobName = "_{}".format(task_id)
    job.workingDirectory = os.getcwd()
    if command_info["walltime"] is not None:
        job.hardWallclockTimeLimit = command_info["walltime"]
    if command_info["nodes"] is not None:
        job.nativeSpecification = command_info["nodes"]

    # Creating a new entry in the database
    db.create_task_entry(task_id, db_name)

    # Running the job
    job_id = drmaa_session.runJob(job)

    # Waiting for the job
    ret_val = None
    try:
        ret_val = drmaa_session.wait(job_id, Session.TIMEOUT_WAIT_FOREVER)

    except KeyboardInterrupt:
        drmaa_session.control(job_id, JobControlAction.TERMINATE)
        logging.warning("{}: terminated".format(task_id))
        raise

    finally:
        drmaa_session.deleteJobTemplate(job)

    # The job is done
    logging.debug("'{}' finished".format(task_id))

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
        if task_id.startswith("impute2"):
            # Task is IMPUTE2, and it might be normal according to message in
            # the summary file
            impute2_fn = None
            for fn in command_info["o_files"]:
                if fn.endswith(".impute2"):
                    impute2_fn = fn
                    break
            if not _check_impute2_file(impute2_fn):
                logging.debug("'{}' exit status problem".format(task_id))
                return False, name, "problem", None

        elif task_id.startswith("shapeit_check"):
            shapeit_fn = None
            for fn in command_info["o_files"]:
                if fn.endswith(".alignments.snp.strand"):
                    shapeit_fn = fn
                    break
            if not _check_shapeit_failed_rc(shapeit_fn):
                logging.debug("'{}' exit status problem".format(task_id))
                return False, name, "problem", None

        else:
            # There was a problem...
            logging.debug("'{}' exit status problem".format(task_id))
            return False, name, "problem", None

    # Checking all the required files were generated
    if not _check_output_files(command_info["o_files"], task_id):
        logging.debug("'{}' exit status problem".format(task_id))
        return False, name, "problem", None

    # Getting the launch time (should always been present)
    launch_time = float(ret_val.resourceUsage["submission_time"])

    # Getting the start time. If 'start_time' is missing from the dictionary
    # (e.g. using Slurm), we use the launch time as the start time
    start_time = float(ret_val.resourceUsage.get("start_time", launch_time))

    # Getting the end time. If 'end_time' is missing from the dictionary
    # (e.g. using Slurm), we use the current time as the end time
    end_time = float(ret_val.resourceUsage.get(
        "end_time", time.mktime(datetime.now().timetuple()),
    ))

    # The task was performed correctly, so we update to completed
    db.mark_drmaa_task_completed(task_id, launch_time, start_time, end_time,
                                 db_name)

    # Everything when well
    logging.debug("'{}' everything was fine".format(task_id))
    return True, name, "performed", db.get_task_runtime(task_id, db_name)
