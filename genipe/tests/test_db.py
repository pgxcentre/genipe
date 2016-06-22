
# This file is part of genipe.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


import time
import logging
import sqlite3
import unittest
from datetime import datetime
from tempfile import TemporaryDirectory

from ..db import utils as db_utils
from ..db.utils import _create_db_connection


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = "Copyright 2014, Beaulieu-Saucier Pharmacogenomics Centre"
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


__all__ = ["TestDB"]


class TestDB(unittest.TestCase):

    def setUp(self):
        """Setup the tests."""
        # Creating the temporary directory
        self.output_dir = TemporaryDirectory(prefix="genipe_test_")

        # We need to create an empty database
        self.db_name = db_utils.create_task_db(self.output_dir.name)

        # We're going to add thee entries
        self.creation_times = []
        self.task_names = []
        for i in range(4):
            task_name = "dummy_task_{}".format(i + 1)
            self.creation_times.append(datetime.now())
            db_utils.create_task_entry(task_name, self.db_name)
            self.task_names.append(task_name)

    def tearDown(self):
        """Finishes the test."""
        # Deleting the output directory
        self.output_dir.cleanup()

    def test_create_task_db(self):
        """Tests the 'create_task_db' function."""
        # The DB should already be created, so we connect to it
        conn = sqlite3.connect(
            self.db_name,
            timeout=360,
            detect_types=sqlite3.PARSE_DECLTYPES | sqlite3.PARSE_COLNAMES,
        )
        c = conn.cursor()

        # Getting all the tables
        c.execute("SELECT name FROM sqlite_master WHERE type='table'")
        tables = c.fetchall()

        # There should only be one table (and one column)
        self.assertEqual(1, len(tables))
        self.assertEqual(1, len(tables[0]))

        # The only value should have name=genipe_task
        self.assertEqual("genipe_task", tables[0][0])

        # Checking the columns (name, type, notnull, default, primary)
        expected_columns = {
            "name": ("name", "TEXT", 0, None, 1),
            "launch": ("launch", "TIMESTAMP", 0, None, 0),
            "start": ("start", "TIMESTAMP", 0, None, 0),
            "end": ("end", "TIMESTAMP", 0, None, 0),
            "completed": ("completed", "INT", 0, None, 0),
        }
        c.execute("PRAGMA table_info(genipe_task)")
        col_names = set()
        for result in c.fetchall():
            # Getting the name of the column
            col_name = result[1]
            col_names.add(col_name)

            # Checking the content of the columns
            self.assertEqual(expected_columns[col_name], result[1:])

        # Checking we have all columns
        col_diff = (col_names & set(expected_columns.keys())) - col_names
        if len(col_diff) != 0:  # pragma: no cover
            self.fail("not all DB columns are present")

    def test_create_db_connection(self):
        """Tests the '_create_db_connection' function."""
        # Creating the connection
        conn, c = _create_db_connection(self.db_name)

        # Checking that the table 'genipe_task' exists
        c.execute("SELECT name FROM sqlite_master WHERE type='table'")
        self.assertEqual("genipe_task", c.fetchone()[0])

        # Closing the connection
        conn.close()

    def test_check_task_completion(self):
        """Tests the 'check_task_completion' function."""
        # Marking the first and fourth tasks as completed
        start = self.creation_times[3].timestamp()
        now = datetime.now().timestamp()
        db_utils.mark_task_completed(self.task_names[0], self.db_name)
        db_utils.mark_drmaa_task_completed(self.task_names[3], start, start,
                                           now, self.db_name)

        # Marking the second one as incomplete
        db_utils.mark_task_incomplete(self.task_names[1], self.db_name)

        # Checking the values
        self.assertTrue(db_utils.check_task_completion(self.task_names[0],
                                                       self.db_name))
        self.assertFalse(db_utils.check_task_completion(self.task_names[1],
                                                        self.db_name))
        self.assertFalse(db_utils.check_task_completion(self.task_names[2],
                                                        self.db_name))
        self.assertTrue(db_utils.check_task_completion(self.task_names[3],
                                                       self.db_name))

        # Setting a completely random value for the third task
        conn, c = _create_db_connection(self.db_name)
        c.execute("UPDATE genipe_task SET completed='foo' WHERE name=?",
                  (self.task_names[3], ))
        conn.commit()
        conn.close()
        self.assertFalse(db_utils.check_task_completion(self.task_names[3],
                                                        self.db_name))

        # The logging capability might be disable...
        disable_lvl = logging.Logger.manager.disable
        logging.disable(logging.NOTSET)

        # Checking the status of a missing task
        with self.assertLogs(level="DEBUG") as cm:
            db_utils.check_task_completion("DUMMY_NAME", self.db_name)
        log_m = ("DEBUG:root:'DUMMY_NAME' no entry")
        self.assertEqual(1, len(cm.output))
        self.assertEqual(log_m, cm.output[0])

        # Setting back to disabled level (if required)
        logging.disable(disable_lvl)

    def test_create_task_entry(self):
        """Tests the 'create_task_entry' function."""
        # The task that will be modified
        modified_task = self.task_names[0]

        # Three tasks should have been created
        conn, c = _create_db_connection(self.db_name)

        # Fetching all the tasks
        c.execute("SELECT name FROM genipe_task")
        results = [r[0] for r in c.fetchall()]
        self.assertEqual(self.task_names, results)

        # Checking that the times are the same
        c.execute(
            "SELECT name, launch, start, end, completed FROM genipe_task"
        )
        results = c.fetchall()
        iteration = zip(self.task_names, self.creation_times, results)
        for task_name, expected_time, result in iteration:
            # The observed values
            o_name, o_launch, o_start, o_end, o_completed = result

            # The name should be the same
            self.assertEqual(task_name, o_name)

            # The launch and start times should be the same (max 1 second diff)
            self.assertEqual(o_launch, o_start)
            t_delta = abs((o_launch - expected_time).total_seconds())
            self.assertTrue(t_delta >= 0 and t_delta <= 1)

            # End and completed should be none
            self.assertTrue(o_end is None)
            self.assertTrue(o_completed is None)

        # We are going to relaunch a task after 3 seconds
        time.sleep(3)
        now = datetime.now()
        db_utils.create_task_entry(modified_task, self.db_name)
        c.execute(
            "SELECT name, launch, start, end, completed FROM genipe_task"
        )
        results = c.fetchall()
        iteration = zip(self.task_names, self.creation_times, results)
        for task_name, expected_time, result in iteration:
            # The expected values
            o_name, o_launch, o_start, o_end, o_completed = result

            # The name should be the same
            self.assertEqual(task_name, o_name)

            # The expected time is different for the first task
            if task_name == modified_task:
                expected_time = now

            # The launch and start times should be the same (max 1 second diff)
            self.assertEqual(o_launch, o_start)
            t_delta = abs((o_launch - expected_time).total_seconds())
            self.assertTrue(t_delta >= 0 and t_delta <= 1)

            # End and completed should be none (except for first task)
            self.assertTrue(o_end is None)

            if task_name == modified_task:
                self.assertEqual(0, o_completed)
            else:
                self.assertTrue(o_completed is None)

        # Closing the connection
        conn.close()

    def test_mark_task_completed(self):
        """Tests the 'mark_task_completed' function."""
        # The task that will be modified
        modified_task = self.task_names[0]

        # We are going to mark the first task as completed after 3 seconds
        time.sleep(3)
        completion_time = datetime.now()
        db_utils.mark_task_completed(modified_task, self.db_name)

        # Creating the connection
        conn, c = _create_db_connection(self.db_name)

        # Checking that the times are the same
        c.execute(
            "SELECT name, launch, start, end, completed FROM genipe_task"
        )
        results = c.fetchall()
        iteration = zip(self.task_names, self.creation_times, results)
        for task_name, expected_time, result in iteration:
            # The observed values
            o_name, o_launch, o_start, o_end, o_completed = result

            # The name should be the same
            self.assertEqual(task_name, o_name)

            # The launch and start times should be the same (max 1 second diff)
            self.assertEqual(o_launch, o_start)
            t_delta = abs((o_launch - expected_time).total_seconds())
            self.assertTrue(t_delta >= 0 and t_delta <= 1)

            # End and completed should be none unless it's the first task
            if task_name != modified_task:
                self.assertTrue(o_end is None)
                self.assertTrue(o_completed is None)
            else:
                # The task should be completed
                self.assertEqual(1, o_completed)

                # Time difference between completion times (max 1 second diff)
                t_delta = abs((o_end - completion_time).total_seconds())
                self.assertTrue(t_delta >= 0 and t_delta <= 1)

                # Time difference between start and end should be 3 (max 1
                # second diff)
                t_delta = abs((o_end - o_start).total_seconds())
                self.assertTrue(t_delta >= 3 and t_delta <= 4)

        conn.close()

    def test_mark_task_incomplete(self):
        """Tests the 'mark_task_incomplete' function."""
        # Marking an incomplete task shouldn't change its values (except
        # completed)
        modified_task = self.task_names[0]
        db_utils.mark_task_incomplete(modified_task, self.db_name)

        # Creating the connection
        conn, c = _create_db_connection(self.db_name)

        # Checking that the times are the same
        c.execute(
            "SELECT name, launch, start, end, completed FROM genipe_task"
        )
        results = c.fetchall()
        iteration = zip(self.task_names, self.creation_times, results)
        for task_name, expected_time, result in iteration:
            # The observed values
            o_name, o_launch, o_start, o_end, o_completed = result

            # The name should be the same
            self.assertEqual(task_name, o_name)

            # The launch and start times should be the same (max 1 second diff)
            self.assertEqual(o_launch, o_start)
            t_delta = abs((o_launch - expected_time).total_seconds())
            self.assertTrue(t_delta >= 0 and t_delta <= 1)

            # End and completed should be none unless it's the first task
            self.assertTrue(o_end is None)
            if task_name != modified_task:
                self.assertTrue(o_completed is None)
            else:
                # The task should be completed
                self.assertEqual(0, o_completed)

        # Now, marking a task as completed, waiting for 3 seconds and mark it
        # as incomplete
        modified_task = self.task_names[1]
        completion_time = datetime.now()
        db_utils.mark_task_completed(modified_task, self.db_name)
        time.sleep(3)
        db_utils.mark_task_incomplete(modified_task, self.db_name)

        # Checking that the times are the same
        c.execute(
            "SELECT name, launch, start, end, completed FROM genipe_task"
        )
        results = c.fetchall()
        iteration = zip(self.task_names, self.creation_times, results)
        for task_name, expected_time, result in iteration:
            # The observed values
            o_name, o_launch, o_start, o_end, o_completed = result

            # The name should be the same
            self.assertEqual(task_name, o_name)

            # The launch and start times should be the same (max 1 second diff)
            self.assertEqual(o_launch, o_start)
            t_delta = abs((o_launch - expected_time).total_seconds())
            self.assertTrue(t_delta >= 0 and t_delta <= 1)

            # End and completed should be none unless it's the first task
            if task_name != modified_task:
                self.assertTrue(o_end is None)
                self.assertTrue((o_completed is None) or (o_completed != 1))
            else:
                # The task should be completed
                self.assertEqual(0, o_completed)

                # Time difference between completion times (max 1 second diff)
                t_delta = abs((o_end - completion_time).total_seconds())
                self.assertTrue(t_delta >= 0 and t_delta <= 1)

        # Closing the connection
        conn.close()

    def test_mark_drmaa_task_completed(self):
        """Tests the 'mark_drmaa_task_completed' function."""
        # The task that will be modified
        modified_task = self.task_names[0]

        # Waiting 1 second and "launch" task
        time.sleep(1)
        launch_time = datetime.now()

        # Waiting 1 second and "start" task
        time.sleep(1)
        start_time = datetime.now()

        # Waiting 3 seconds and "ending" task
        time.sleep(3)
        end_time = datetime.now()
        db_utils.mark_drmaa_task_completed(
            modified_task, launch_time.timestamp(), start_time.timestamp(),
            end_time.timestamp(), self.db_name,
        )

        # Creating the connection
        conn, c = _create_db_connection(self.db_name)

        # Checking that the times are the same
        c.execute(
            "SELECT name, launch, start, end, completed FROM genipe_task"
        )
        results = c.fetchall()
        iteration = zip(self.task_names, self.creation_times, results)
        for task_name, expected_time, result in iteration:
            # The observed values
            o_name, o_launch, o_start, o_end, o_completed = result

            # The name should be the same
            self.assertEqual(task_name, o_name)

            # The launch and start times should be the same (max 1 second diff)
            # for the other tasks
            if o_name != modified_task:
                self.assertEqual(o_launch, o_start)
                t_delta = abs((o_launch - expected_time).total_seconds())
                self.assertTrue(t_delta >= 0 and t_delta <= 1)
            else:
                t_delta = abs(launch_time.timestamp() - o_launch.timestamp())
                self.assertTrue(t_delta >= 0 and t_delta <= 1)

                t_delta = abs(start_time.timestamp() - o_start.timestamp())
                self.assertTrue(t_delta >= 0 and t_delta <= 1)

                t_delta = abs(end_time.timestamp() - o_end.timestamp())
                self.assertTrue(t_delta >= 0 and t_delta <= 1)

            # End and completed should be none unless it's the first task
            if task_name != modified_task:
                self.assertTrue(o_end is None)
                self.assertTrue(o_completed is None)
            else:
                # The task should be completed
                self.assertEqual(1, o_completed)

                # Time difference between completion times (max 1 second diff)
                t_delta = abs((end_time - o_end).total_seconds())
                self.assertTrue(t_delta >= 0 and t_delta <= 1)

                # Time difference between launch and start should be 1 (max 1
                # second diff)
                t_delta = abs((o_start - o_launch).total_seconds())
                self.assertTrue(t_delta >= 1 and t_delta <= 2)

                # Time difference between start and end should be 3 (max 1
                # second diff)
                t_delta = abs((o_end - o_start).total_seconds())
                self.assertTrue(t_delta >= 3 and t_delta <= 4)

        conn.close()

    def test_get_task_runtime(self):
        """Tests the 'task_runtime' function."""
        # Those two tasks will be modified
        modified_task_1 = self.task_names[0]
        modified_task_2 = self.task_names[1]

        # Waiting 1 second and "launch" task
        time.sleep(1)
        launch_time = datetime.now()

        # Waiting 1 second and "start" task
        time.sleep(1)
        start_time = datetime.now()

        # Waiting 3 seconds and "ending" two task
        time.sleep(3)
        end_time = datetime.now()
        db_utils.mark_task_completed(modified_task_1, self.db_name)
        db_utils.mark_drmaa_task_completed(
            modified_task_2, launch_time.timestamp(), start_time.timestamp(),
            end_time.timestamp(), self.db_name,
        )

        # Getting the first task time
        task_time_1 = db_utils.get_task_runtime(modified_task_1, self.db_name)
        task_time_2 = db_utils.get_task_runtime(modified_task_2, self.db_name)

        # Comparing the time
        self.assertEqual(5, task_time_1)
        self.assertEqual(3, task_time_2)

    def test_get_all_runtimes(self):
        """Tests the 'get_all_runtimes' function."""
        # The time that the task started
        start = self.creation_times[-1].timestamp()

        # Waiting one second for each task
        end_times = []
        for task_name in self.task_names:
            time.sleep(1)
            if task_name != self.task_names[-1]:
                db_utils.mark_task_completed(task_name, self.db_name)

            else:
                now = datetime.now().timestamp()
                db_utils.mark_drmaa_task_completed(task_name, start, start,
                                                   now, self.db_name)
            end_times.append(time.time())

        # The expected time
        expected_time = {
            task_name: int(round(elapsed - start, 0))
            for task_name, elapsed in zip(self.task_names, end_times)
        }

        # Getting the time for all tasks
        observed_time = db_utils.get_all_runtimes(self.db_name)

        # Comparing the results
        self.assertEqual(set(expected_time.keys()), set(observed_time.keys()))

        for task_name in expected_time.keys():
            t_delta = abs(expected_time[task_name] - observed_time[task_name])
            self.assertTrue(t_delta >= 0 and t_delta <= 1)

        # Setting one of the task's end time to None
        conn, c = _create_db_connection(self.db_name)
        c.execute("UPDATE genipe_task SET end=NULL WHERE name=?",
                  (self.task_names[0], ))
        conn.commit()
        conn.close()
        with self.assertLogs(level="WARNING") as cm:
            db_utils.get_all_runtimes(self.db_name)
        log_m = "WARNING:root:{}: no execution time for task"
        self.assertEqual(1, len(cm.output))
        self.assertEqual(log_m.format(self.task_names[0]), cm.output[0])

        # Setting one of the task's start time to None
        conn, c = _create_db_connection(self.db_name)
        c.execute("UPDATE genipe_task SET start=NULL WHERE name=?",
                  (self.task_names[0], ))
        conn.commit()
        conn.close()
        with self.assertLogs(level="WARNING") as cm:
            db_utils.get_all_runtimes(self.db_name)
        log_m = "WARNING:root:{}: no execution time for task"
        self.assertEqual(1, len(cm.output))
        self.assertEqual(log_m.format(self.task_names[0]), cm.output[0])
