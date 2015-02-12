
# This file is part of gwip.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


import logging
import configparser

from .. import chromosomes


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = "Copyright 2014, Beaulieu-Saucier Pharmacogenomics Centre"
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


__all__ = ["parse_drmaa_config", ]


def parse_drmaa_config(configfile):
    """Parses the tasks' configuration file for DRMAA."""
    # This will save the final DRMAA configuration
    final_config = {}

    # Loading the configuration file
    drmaa_config = configparser.ConfigParser()
    drmaa_config.read(configfile)

    # Going through file task by task
    # First is 'plink_exclude'
    config = _generate_default_values("plink_exclude", drmaa_config)
    final_config.update(config)


def _generate_default_values(task_name, config):
    """Generates default values for missing DRMAA configuration."""
    ## "walltime": bytes("00:15:00", encoding="ascii"),
    ## "nodes": bytes("-l nodes=1:ppn=1", encoding="ascii"),
    # The final tool configuration
    final_tool_config = {}

    # The default walltime, nodes and ppn values
    walltime = "00:15:00"
    nodes = "1"
    ppn = "1"

    # The default task configuration
    task_config = {
        "walltime": walltime,
        "nodes":    nodes,
        "ppn":      ppn,
    }

    # Is there a configuration for this task?
    if task_name in config:
        # The task configuration
        task_config = config[task_name]

    # Getting the general walltime
    walltime = task_config.get("walltime", walltime)
    nodes = task_config.get("nodes", nodes)
    ppn = task_config.get("ppn", ppn)

    # Getting the final values for each chromosomes
    for chrom in chromosomes:
        # Getting the chromosome specif walltime
        chrom_walltime = task_config.get("chr{}_walltime".format(chrom),
                                         walltime)
        chrom_walltime = bytes(chrom_walltime, encoding="ascii")

        # Getting the chromosome specif nodes/ppn
        chrom_nodes = task_config.get("chr{}_nodes".format(chrom), nodes)
        chrom_ppn = task_config.get("chr{}_ppn".format(chrom), ppn)
        chrom_nodes = bytes("-l nodes={}:ppn={}".format(chrom_nodes,
                                                        chrom_ppn),
                            encoding="ascii")

        # Saving the final config
        final_tool_config["{}_chr{}".format(task_name, chrom)] = {
            "walltime": chrom_walltime,
            "nodes":    chrom_nodes,
        }

    return final_tool_config
