

__all__ = ["generate_report", ]


import os
import shutil
from datetime import date

from pkg_resources import resource_filename

from .utils import *
from .. import __version__
from ..error import ProgramError


def generate_report(out_dir, run_opts, run_info):
    """Generate the report."""
    # Configuring Jinja2
    jinja2_env = config_jinja2()

    # Gathering the report data
    report_data = {
        "report_number":   sanitize_tex(run_opts.report_number),
        "title":           sanitize_tex(run_opts.report_title),
        "author":          sanitize_tex(run_opts.report_author),
        "date":            sanitize_tex("{:%B %d, %Y}".format(date.today())),
        "package_name":    sanitize_tex(__name__.split(".")[0]),
        "package_version": sanitize_tex(__version__),
    }

    # Gathering the report content
    report_content = ""
    report_content += _generate_background(jinja2_env, run_opts, run_info)
    report_content += _generate_methods(jinja2_env, run_opts, run_info)
    report_content += _generate_results(jinja2_env, run_opts, run_info)
    report_content += _generate_conclusions(jinja2_env, run_opts, run_info)
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

    # Copying the bibliography file
    bib_file = resource_filename(__name__, "templates/biblio/references.bib")
    shutil.copy(bib_file, out_dir)

    # Copying the bibliography style
    bib_style = resource_filename(__name__, "templates/biblio/references.bst")
    shutil.copy(bib_style, out_dir)


def _generate_background(templates, run_options, run_information):
    """Generates the background section of the report."""
    # Loading the template
    section_template = templates.get_template("section_template.tex")
    background = templates.get_template("parts/background.tex")

    # Returning the section
    return section_template.render(
        section_name="Background",
        section_type="section",
        section_label="sec:background",
        section_content=background.render(**run_information),
    )


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
    methods = templates.get_template("parts/methods.tex")

    # The input files
    data_files = [
        "{}.{}".format(run_options.bfile, ext) for ext in ["bed", "bim", "fam"]
    ]
    data_files = [
        format_tex(sanitize_tex(text), "texttt") for text in data_files
    ]

    # Creating the iteration
    data_files = itemize_template.render(iteration_type="itemize",
                                         iteration_list=data_files)

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
            "In total, {ambiguous} ambiguous, {duplicated} duplicated and "
            "{special} special markers were excluded.".format(
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

    steps = itemize_template.render(iteration_type="enumerate",
                                    iteration_list=steps)

    # Returning the section
    return section_template.render(
        section_name="Methods",
        section_type="section",
        section_label="sec:methods",
        section_content=methods.render(
            data_files=data_files,
            steps_data=steps,
            **run_information
        ),
    )


def _generate_results(templates, run_options, run_information):
    """Generates the results section of the report."""
    # Some assertions
    required_variables = ["cross_validation_final_nb_genotypes",
                          "cross_validation_nb_genotypes_chrom",
                          "cross_validation_table_1",
                          "cross_validation_table_2",
                          "cross_validation_table_1_chrom",
                          "cross_validation_table_2_chrom", "prob_threshold",
                          "nb_imputed", "average_comp_rate", "rate_threshold",
                          "nb_good_sites", "average_comp_rate_cleaned",
                          "mean_missing", "nb_samples", "nb_genotyped",
                          "nb_genotyped_not_complete",
                          "pct_genotyped_not_complete", "nb_geno_now_complete",
                          "pct_geno_now_complete", "nb_site_now_complete",
                          "pct_good_sites", "nb_missing_geno"]
    for required_variable in required_variables:
        assert required_variable in run_information

    # Loading the templates
    section_template = templates.get_template("section_template.tex")
    tabular_template = templates.get_template("tabular_template.tex")
    float_template = templates.get_template("float_template.tex")
    cross_validation = templates.get_template("parts/cross_validation.tex")
    completion_rate = templates.get_template("parts/completion_rate.tex")

    # The header of the two kind of tables
    header_table_1 = [
        format_tex(sanitize_tex("Interval"), "textbf"),
        format_tex(sanitize_tex("Nb Geno"), "textbf"),
        format_tex(sanitize_tex("Concordance (%)"), "textbf"),
    ]
    header_table_2 = [
        format_tex(sanitize_tex("Interval"), "textbf"),
        format_tex(sanitize_tex("Called (%)"), "textbf"),
        format_tex(sanitize_tex("Concordance (%)"), "textbf"),
    ]

    # Creating the tables
    tables = ""

    # Adding the table for each of the chromosomes
    for chrom in range(1, 23):
        # Getting the table 1
        table_1 = run_information["cross_validation_table_1_chrom"][chrom]
        for i in range(len(table_1)):
            table_1[i][0] = tex_inline_math(table_1[i][0])
        table_1 = create_tabular(
            template=tabular_template,
            header=header_table_1,
            col_align=["c", "r", "r"],
            data=table_1,
        )

        # Getting the table 2
        table_2 = run_information["cross_validation_table_2_chrom"][chrom]
        for i in range(len(table_2)):
            table_2[i][0] = tex_inline_math(
                table_2[i][0].replace(">=", r"\geq "),
            )
        table_2 = create_tabular(
            template=tabular_template,
            header=header_table_2,
            col_align=["c", "r", "r"],
            data=table_2,
        )

        # The number of genotypes
        nb_genotypes = run_information["cross_validation_nb_genotypes_chrom"]
        nb_genotypes = nb_genotypes[chrom]

        # Adding the float
        tables += create_float(
            template=float_template,
            float_type="table",
            caption=wrap_tex(sanitize_tex(
                "IMPUTE2's internal cross-validation for chromosome {}. "
                "Tables show the percentage of concordance between genotyped "
                "calls and imputed calls for {:,d} "
                "genotypes.".format(chrom, nb_genotypes)
            )),
            label="tab:cross_validation_chr_{}".format(chrom),
            placement="H",
            content=table_1 + r"\hfill" + table_2,
        )

    # Adding the table for all the autosomes (Table 1)
    table_1 = run_information["cross_validation_table_1"]
    for i in range(len(table_1)):
        table_1[i][0] = tex_inline_math(table_1[i][0])
    table_1 = create_tabular(template=tabular_template, header=header_table_1,
                             col_align=["c", "r", "r"], data=table_1)

    # Adding the table for all the autosomes (Table 2)
    table_2 = run_information["cross_validation_table_2"]
    for i in range(len(table_2)):
        table_2[i][0] = tex_inline_math(table_2[i][0].replace(">=", r"\geq "))
    table_2 = create_tabular(template=tabular_template, header=header_table_2,
                             col_align=["c", "r", "r"], data=table_2)

    # The number of genotypes
    nb_genotypes = run_information["cross_validation_final_nb_genotypes"]

    # Adding the float
    tables += "\n\n" + create_float(
        template=float_template,
        float_type="table",
        caption=wrap_tex(sanitize_tex(
            "IMPUTE2's internal cross-validation across the genome. Tables "
            "show the percentage of concordance between genotyped calls and "
            "imputed calls for {:,d} genotypes.".format(nb_genotypes)
        )),
        label="tab:cross_validation",
        placement="H",
        content=table_1 + r"\hfill" + table_2,
    )

    # Creating the cross-validation subsection
    cross_validation_content = section_template.render(
        section_name="Cross-validation",
        section_type="subsection",
        section_content=cross_validation.render(tables=tables),
        section_label="subsec:cross_validation",
    )

    completion_rate_content = section_template.render(
        section_name="Completion rate",
        section_type="subsection",
        section_content=completion_rate.render(**run_information),
        section_label="subsec:completion_rate",
    )

    # The final content
    content = cross_validation_content + completion_rate_content

    return section_template.render(section_name="Results",
                                   section_type="section",
                                   section_content=content,
                                   section_label="sec:results")


def _generate_conclusions(templates, run_options, run_information):
    """Generates the background section of the report."""
    # Loading the template
    section_template = templates.get_template("section_template.tex")
    conclusions = templates.get_template("parts/conclusions.tex")

    # Returning the section
    return section_template.render(
        section_name="Conclusions",
        section_type="section",
        section_label="sec:conclusions",
        section_content=conclusions.render(**run_information),
    )
