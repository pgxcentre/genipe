
import os
import logging
import sqlite3
from datetime import datetime


__all__ = ["create_task_db", "check_task_completion", "create_task_entry",
           "mark_task_completed", "get_task_runtime",
           "mark_drmaa_task_completed"]


def create_task_db(out_dir):
    """Creates a task DB."""
    # The name
    db_name = os.path.join(out_dir, "tasks.db")
    logging.info("Connecting to DB '{}'".format(db_name))

    # The DB
    conn, c = _create_db_connection(db_name)

    # Creating the table if it doesn't exists
    c.execute("""CREATE TABLE IF NOT EXISTS gwip_task (
                    name TEXT PRIMARY KEY,
                    launch TIMESTAMP,
                    start TIMESTAMP,
                    end TIMESTAMP,
                    completed INT)""")

    # Committing the changes
    conn.commit()
    conn.close()

    return db_name


def _create_db_connection(db_name):
    """Creates a DB connection."""
    conn = sqlite3.connect(
        db_name,
        timeout=360,
        detect_types=sqlite3.PARSE_DECLTYPES | sqlite3.PARSE_COLNAMES,
    )
    c = conn.cursor()

    return conn, c


def check_task_completion(task_id, db_name):
    """Checks if the task exists and is completed."""
    conn, c = _create_db_connection(db_name)

    # Retrieving the task information
    c.execute("SELECT completed FROM gwip_task WHERE name=?", (task_id, ))
    r = c.fetchone()

    conn.close()

    if r is None:
        # There is not entry with this task ID
        logging.debug("'{}' no entry".format(task_id))
        return False

    if r[0] is None or r[0] != 1:
        # There is an entry, but it wasn't completed
        logging.debug("'{}' not completed ({})".format(task_id, r[0]))
        return False

    return True


def create_task_entry(task_id, db_name):
    """Creates a update a task entry."""
    conn, c = _create_db_connection(db_name)

    # Checking if the entry already exists
    c.execute("SELECT name FROM gwip_task WHERE name=?", (task_id, ))
    r = c.fetchone()

    # The time of launch
    time = datetime.now()

    if r is None:
        # This is the first time we see this task, so we create an new entry
        c.execute("INSERT INTO gwip_task (name, launch, start) "
                  "VALUES (?, ?, ?)",
                  (task_id, time, time))

    else:
        # We saw this entry, but we need to relaunch it
        c.execute("UPDATE gwip_task "
                  "SET launch=?, start=? WHERE name=?",
                  (time, time, task_id))

    conn.commit()
    conn.close()


def mark_task_completed(task_id, db_name):
    """Marks the task as completed."""
    conn, c = _create_db_connection(db_name)

    # Updating the end time
    c.execute("UPDATE gwip_task SET end=?, completed=1 WHERE name=?",
              (datetime.now(), task_id))

    conn.commit()
    conn.close()


def mark_drmaa_task_completed(task_id, launch_time, start_time, end_time,
                              db_name):
    """Marks a task run by DRMAA as completed (while updating times)."""
    conn, c = _create_db_connection(db_name)

    # The time
    launch_time = datetime.fromtimestamp(launch_time)
    start_time = datetime.fromtimestamp(start_time)
    end_time = datetime.fromtimestamp(end_time)

    # Updating
    c.execute("UPDATE gwip_task SET launch=?, start=?, end=?, completed=1 "
              "WHERE name=?", (launch_time, start_time, end_time, task_id))

    conn.commit()
    conn.close()


def get_task_runtime(task_id, db_name):
    """Gets the task run time."""
    conn, c = _create_db_connection(db_name)

    # Getting the start and end time
    c.execute("SELECT start, end FROM gwip_task WHERE name=?", (task_id, ))
    r = c.fetchone()

    conn.close()

    return int((r[1] - r[0]).total_seconds())
