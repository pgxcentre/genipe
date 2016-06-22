
# This file is part of genipe.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


import configparser

from .. import chromosomes


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = "Copyright 2014, Beaulieu-Saucier Pharmacogenomics Centre"
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


__all__ = ["parse_drmaa_config", ]


def parse_drmaa_config(configfile):
    """Parses the tasks' configuration file for DRMAA.

    Args:
        configfile (str): the name of the configuration file

    Returns:
        dict: the DRMAA configuration for each task

    """
    # This will save the final DRMAA configuration
    final_config = {}

    # Loading the configuration file
    drmaa_config = configparser.ConfigParser()
    drmaa_config.read(configfile)

    # Is SKIP is set?
    if "main" in drmaa_config:
        if drmaa_config["main"].get("skip_drmaa_config", "no") == "yes":
            return {"skip_drmaa_config": True}

    # First is 'plink_exclude'
    config = _generate_default_values("plink_exclude", drmaa_config)
    final_config.update(config)

    # Second is 'plink_missing_rate' (only one value)
    config = _generate_default_values("plink_missing_rate", drmaa_config,
                                      only_one=True)
    final_config.update(config)

    # Third is 'shapeit_check_1'
    config = _generate_default_values("shapeit_check_1", drmaa_config,
                                      template="shapeit_check_chr{chrom}_1")
    final_config.update(config)

    # Fourth is 'plink_flip'
    config = _generate_default_values("plink_flip", drmaa_config)
    final_config.update(config)

    # Fifth is 'shapeit_check_2'
    config = _generate_default_values("shapeit_check_2", drmaa_config,
                                      template="shapeit_check_chr{chrom}_2")
    final_config.update(config)

    # Sixth is 'plink_final_exclude'
    config = _generate_default_values("plink_final_exclude", drmaa_config)
    final_config.update(config)

    # Seventh is 'shapeit_phase'
    config = _generate_default_values("shapeit_phase", drmaa_config)
    final_config.update(config)

    # Eight is 'impute2_chr{}_{}_{}'
    config = _generate_default_values("impute2", drmaa_config)
    final_config.update(config)

    # Ninth  is 'merge_impute2'
    config = _generate_default_values("merge_impute2", drmaa_config)
    final_config.update(config)

    # Tenth is 'bgzip'
    config = _generate_default_values("bgzip", drmaa_config)
    final_config.update(config)

    return final_config


def _generate_default_values(task_name, config, walltime="00:15:00", nodes="1",
                             ppn="1", only_one=False, template=None):
    """Generates default values for missing DRMAA configuration.

    Args:
        task_name (str): the name of the task
        config (dict): the configuration
        walltime (str): the default execution time
        nodes (str): the default number of nodes
        ppn (str): the default number of processes
        only_one (bool): if there is only on task (and not one per chromosome)
        template (str): task name template (for each chromosome)

    Returns:
        dict: the final configuration for this task

    """
    # The final tool configuration
    final_tool_config = {}

    # The nodes configuration string
    nodes_string = "-l nodes={n}:ppn={p}"

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
    walltime = task_config.pop("walltime", walltime)
    nodes = task_config.pop("nodes", nodes)
    ppn = task_config.pop("ppn", ppn)

    # Is there only one value?
    if only_one:
        final_tool_config[task_name] = {
            "walltime": bytes(walltime, encoding="ascii"),
            "nodes":    bytes(nodes_string.format(n=nodes, p=ppn),
                              encoding="ascii"),
        }
        return final_tool_config

    # The template for chromosome specific tasks
    if template is None:
        template = task_name + "_chr{chrom}"

    # Getting the final values for each autosomes
    for chrom in chromosomes + ("25_1", "25_2"):
        # Getting the chromosome specif walltime
        chrom_walltime = task_config.pop("chr{}_walltime".format(chrom),
                                         walltime)
        chrom_walltime = bytes(chrom_walltime, encoding="ascii")

        # Getting the chromosome specif nodes/ppn
        chrom_nodes = task_config.pop("chr{}_nodes".format(chrom), nodes)
        chrom_ppn = task_config.pop("chr{}_ppn".format(chrom), ppn)
        chrom_nodes = bytes(nodes_string.format(n=chrom_nodes, p=chrom_ppn),
                            encoding="ascii")

        # Saving the final config
        final_tool_config[template.format(chrom=chrom)] = {
            "walltime": chrom_walltime,
            "nodes":    chrom_nodes,
        }

    # Saving what's remaining values
    suffixes = ["_walltime", "_nodes", "_ppn"]
    for suffix in suffixes:
        remaining = [
            i for i in task_config.keys() if i.endswith(suffix)
        ]
        for k in remaining:
            # The name of the remaining task (removing the _walltime)
            remaining_task_name = k[:-len(suffix)]

            # The walltime
            walltime_k = remaining_task_name + "_walltime"
            remaining_walltime = task_config.pop(walltime_k, walltime)
            remaining_walltime = bytes(remaining_walltime, encoding="ascii")

            # The nodes and ppn values
            nodes_k = remaining_task_name + "_nodes"
            ppn_k = remaining_task_name + "_ppn"
            remaining_nodes = task_config.pop(nodes_k, nodes)
            remaining_ppn = task_config.pop(ppn_k, ppn)
            remaining_nodes = bytes(
                nodes_string.format(n=remaining_nodes, p=remaining_ppn),
                encoding="ascii",
            )

            # Saving the final config
            final_tool_config[task_name + "_" + remaining_task_name] = {
                "walltime": remaining_walltime,
                "nodes":    remaining_nodes,
            }

    return final_tool_config
