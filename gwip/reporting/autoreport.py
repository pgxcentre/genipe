

__all__ = ["generate_report", ]


import os
from datetime import date

from .utils import *
from ..error import ProgramError


def generate_report(out_dir, run_opts, run_info):
    """Generate the report."""
    # Configuring Jinja2
    jinja2_env = config_jinja2()

    # Gathering the report data
    report_data = {
        "report_number": sanitize_tex(run_opts.report_number),
        "title":         sanitize_tex(run_opts.report_title),
        "author":        sanitize_tex(run_opts.report_author),
        "date":          sanitize_tex("{:%B %d, %Y}".format(date.today())),
    }

    # Gathering the report content
    report_content = ""
    report_content += _generate_methods(jinja2_env, run_opts, run_info)
    report_data["report_content"] = report_content

    # Getting the template
    main_template = jinja2_env.get_template("main_template.tex")

    # Writing the report
    report_filename = os.path.join(out_dir, "report.tex")
    try:
        with open(report_filename, "w") as o_file:
            print(main_template.render(**report_data), file=o_file)

    except FileNotFoundError:
        raise ProgramError("{}: cannot write file".format(report_filename))


def _generate_methods(templates, run_options, run_information):
    """Generate the method section of the report."""
    # Some assertions
    required_variables = ["shapeit_version", "impute2_version",
                          "initial_nb_markers", "initial_nb_samples",
                          "nb_duplicates", "nb_ambiguous", "nb_flip",
                          "nb_exclude", "nb_phasing_markers"]
    for required_variable in required_variables:
        assert required_variable in run_information

    # Loading the templates
    section_template = templates.get_template("section_template.tex")
    itemize_template = templates.get_template("iterate_template.tex")

    # The text of the introduction
    content = wrap_tex(sanitize_tex(
        "The following (cleaned) files provided information about the study "
        "cohort dataset for {nb_samples:,d} samples and {nb_markers:,d} "
        "markers (including {nb_special:,d} markers located on sexual or "
        "mitochondrial chromosomes):".format(
            nb_samples=run_information["initial_nb_samples"],
            nb_markers=run_information["initial_nb_markers"],
            nb_special=run_information["nb_special_markers"],
        )
    ))

    # The input files
    input_files = [
        "{}.{}".format(run_options.bfile, ext) for ext in ["bed", "bim", "fam"]
    ]
    input_files = [
        format_tex(sanitize_tex(text), "texttt") for text in input_files
    ]

    # Creating the iteration
    content += "\n" + itemize_template.render(iteration_type="itemize",
                                              iteration_list=input_files)

    # The values for the different steps
    content += "\n" + wrap_tex(sanitize_tex(
        "IMPUTE2's pre-phasing approach can work with phased haplotypes from "
        "SHAPEIT, a highly accurate phasing algorithm that can handle "
        "mixtures of unrelated samples, duos or trios. The usage of SHAPEIT "
        "is highly recommended to infer haplotypes underlying the study "
        "genotypes. The phased haplotypes are then passed to IMPUTE2 for "
        "imputation. Although pre-phasing allows for very fast imputation, it "
        "leads to a small loss in accuracy since the estimation uncertainty "
        "in the study haplotypes is ignored. SHAPEIT version {shapeit} and "
        "IMPUTE2 version {impute2} were used for this analysis.".format(
            shapeit=run_information["shapeit_version"],
            impute2=run_information["impute2_version"],
        )
    ))

    content += r"\\" + "\n\n" + wrap_tex(sanitize_tex(
        "To speed up the pre-phasing and imputation steps, the dataset was "
        "split by chromosome. The following quality steps were then performed "
        "on each chromosome:"
    ))

    # The text for the different steps
    steps = []

    # The ambiguous and duplicated markers that were removed
    steps.append(wrap_tex(sanitize_tex(
        "Ambiguous markers with alleles "
    ) + format_tex("A", "texttt") + "/" + format_tex("T", "texttt") + " and "
        + format_tex("C", "texttt") + "/" + format_tex("G", "texttt")
        + (
        ", duplicated markers (same position), and markers located on special "
        "chromosomes (sexual or mitochondrial chromosomes) were excluded from "
        "the imputation. "
    ) + format_tex(
        sanitize_tex(
            "In total, {ambiguous:,d} ambiguous, {duplicated:,d} duplicated "
            "and {special:,d} special markers were excluded.".format(
                ambiguous=run_information["nb_ambiguous"],
                duplicated=run_information["nb_duplicates"],
                special=run_information["nb_special_markers"],
            )
        ),
        "textbf"
    )))

    # The number of markers that were flipped
    steps.append(wrap_tex(sanitize_tex(
        "Markers' strand was checked using the SHAPEIT algorithm and "
        "IMPUTE2's reference files. "
    ) + format_tex(
        sanitize_tex(
            "In total, {nb_markers:,d} markers had an incorrect strand and "
            "were flipped using Plink.".format(
                nb_markers=run_information["nb_flip"],
            )
        ),
        "textbf",
    )))

    # The number of excluded markers because of strand problem
    steps.append(wrap_tex(sanitize_tex(
        "The strand of each marker was checked again using SHAPEIT against "
        "IMPUTE2's reference files. "
    ) + format_tex(
        sanitize_tex(
            "In total, {nb_markers:,d} markers were found to still be on the "
            "wrong strand, and were hence excluded from the final dataset "
            "using Plink.".format(
                nb_markers=run_information["nb_exclude"],
            )
        ),
        "textbf",
    )))

    content += "\n" + itemize_template.render(iteration_type="enumerate",
                                              iteration_list=steps)

    content += "\n" + wrap_tex(format_tex(
        sanitize_tex("In total, {:,d} were used for phasing using "
                     "SHAPEIT.".format(run_information["nb_phasing_markers"])),
        "textbf",
    ) + sanitize_tex(
        " IMPUTE2 was then used to impute markers genome-wide using its "
        "reference file."
    ))

    # Returning the section
    return section_template.render(section_name="Methods",
                                   section_content=content)
