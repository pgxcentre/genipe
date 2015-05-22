
# This file is part of genipe.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


import re
from textwrap import wrap

import jinja2


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = "Copyright 2014, Beaulieu-Saucier Pharmacogenomics Centre"
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


__all__ = ["config_jinja2", "sanitize_tex", "format_tex", "wrap_tex",
           "create_tabular", "create_float", "tex_inline_math", "format_time",
           "colorize_time"]


_char_mod = {
    "~": r"$\sim$",
}
_escaped_char = ["$", "%", "_", "}", "{", "&", "#"]
_valid_tex_formats = {"texttt", "emph", "textbf", "textit"}


def config_jinja2():
    """Configure the jinja2 environment for LaTeX.

    Note
    ----
        The configuration used is for LaTeX documents. Hence, a block command
        is done using ``\BLOCK{}`` and variables using ``\VAR{}`` in the Jinja2
        template.

    """
    return jinja2.Environment(
        block_start_string='\BLOCK{',
        block_end_string='}',
        variable_start_string='\VAR{',
        variable_end_string='}',
        comment_start_string='\#{',
        comment_end_string='}',
        line_statement_prefix='%-',
        line_comment_prefix='%#',
        trim_blocks=True,
        autoescape=False,
        loader=jinja2.PackageLoader(__name__, "templates"),
    )


def sanitize_tex(original_text):
    """Sanitize TeX text.

    Args:
        original_text (str): the text to sanitize for LaTeX

    Text is sanitized by following these steps:

    1. Replaces ``\\\\`` by ``\\textbackslash``
    2. Escapes certain characters (such as ``$``, ``%``, ``_``, ``}``, ``{``,
       ``&`` and ``#``) by adding a backslash (*e.g.* from ``&`` to ``\\&``).
    3. Replaces special characters such as ``~`` by the LaTeX equivalent
       (*e.g.* from ``~`` to ``$\\sim$``).

    """
    # The backslashes
    sanitized_tex = original_text.replace("\\", r"\textbackslash ")

    # Escaping
    sanitized_tex = re.sub(r"([{}])".format("".join(_escaped_char)),
                           r"\\\g<1>", sanitized_tex)

    # Replacing
    for character, mod in _char_mod.items():
        sanitized_tex = sanitized_tex.replace(character, mod)

    return sanitized_tex


def wrap_tex(original_text):
    """Wraps the text.

    Args:
        original_text (str): the text to wrap

    Returns:
        str: a string where the original text was wrapped

    Wraps the text so that lines are no longer than 80 characters. Uses the
    :py:func:`str.join` function on the results of the :py:func:`wrap`
    function, so that a single string is returned.

    """
    return "\n".join(wrap(original_text))


def format_tex(text, tex_format):
    """Change the TeX text format.

    Args:
        text (str): the text for which the format needs to be specified
        tex_format (str): the format of the text to return

    Returns:
        str: the formatted text

    This will change the format by adding the LaTeX format command (*e.g.* from
    ``text`` to ``\\texttt{text}``).

    Note
    ----
        Only the following format are available:

        * ``texttt``
        * ``emph``
        * ``textbf``
        * ``textit``

    """
    assert tex_format in _valid_tex_formats, "invalid format"
    assert _is_sanitized(text), "text not sanitized"

    return r"\%s{%s}" % (tex_format, text)


def tex_inline_math(content):
    """Creates an inline mathematical formula in TeX.

    Args:
        content (str): the content of the mathematical formula

    Returns:
        str: the formatted mathematical formula

    The function only adds ``$`` symbols before and after the content (*e.g.*
    from ``\\pi`` to ``$\\pi$``).

    """
    return "${}$".format(content)


def _is_sanitized(text):
    """Check if text is sanitized.

    Args:
        text (str): the text to check

    Returns:
        bool: ``True`` if the text is sanitized, ``False`` otherwise

    """
    # Checking the escaped characters
    sanitized = re.search(r"[^\\][{}]".format("".join(_escaped_char)), text)
    sanitized = sanitized is None

    # Checking the characters to replace
    for character in _char_mod.keys():
        sanitized = sanitized and (character not in text)

    return sanitized


def create_tabular(template, header, data, header_multicol=None,
                   col_align=None):
    """Creates a TeX tabular.

    Args:
        template (jinja2.Template): the tabular template
        header (list): the header of the tabular
        data (list): the tabular data
        header_multicol (list): the number of columns for the header
        col_align (list): the column alignement

    Returns:
        str: a string representation of a LaTeX tabular

    """
    if header_multicol is None:
        header_multicol = [1 for i in header]

    # Getting the number of columns
    nb_col = sum(header_multicol)

    if col_align is None:
        col_align = ["c"] * nb_col

    # Checking that the number of columns holds
    assert len(header) == len(header_multicol), "len(header) != len(multicol)"
    assert len(col_align) == nb_col, "len(align) != number of columns"

    # Generating the tabular data
    tabular_data = {
        "col_alignments": "".join(col_align),
        "header_data":    zip(header, header_multicol),
        "tabular_data":   data,
    }

    # Rendering
    return template.render(**tabular_data)


def create_float(template, float_type, caption, label, content, placement="H"):
    """Creates a TeX float.

    Args:
        template (jinja2.Template): the float template
        float_type (str): the type of float (``figure`` or  ``table``)
        caption (str): the caption of the float
        label (str): the label of the float
        content (str): the content of the float
        placement (str): the float placement (*e.g.* ``H``)

    Returns:
        str: a string representation of a LaTeX float

    """
    # Some assertions
    assert float_type in ["figure", "table"], "invalid float type"
    for character in placement:
        assert character in "htbp!H", "invalid placement"
    if "H" in placement:
        assert placement == "H", "placement 'H' should be alone"

    # Generating the float data
    float_data = {
        "float_type":      float_type,
        "float_placement": placement,
        "float_caption":   caption,
        "float_label":     label,
        "float_content":   content,
    }

    # Rendering
    return template.render(**float_data)


def format_time(total_seconds, written_time=False):
    """Format time (either "HH:MM:SS" or "H hours, M minutes and S seconds".

    Args:
        total_seconds (int): the total number of seconds
        written_time (bool): whether to write time in written language

    Returns:
        str: a string representation of the total time

    If ``written_time`` is ``True``, time will be displayed as "H hours, M
    minutes and S seconds". Otherwise, the time will be represented as
    HH:MM:SS.

    """
    # The format for the time
    time_fmt = "{hours:02d}:{minutes:02d}:{seconds:02d}"

    # Getting the number of seconds, minutes and hours
    minutes, seconds = divmod(total_seconds, 60)
    hours, minutes = divmod(minutes, 60)

    # Formatting
    if not written_time:
        return time_fmt.format(seconds=seconds, minutes=minutes, hours=hours)

    # The written time
    written_time = []

    # The hours
    if hours > 0:
        written_time.append(
            "{} hour{}".format(hours, "s" if hours > 1 else "")
        )

    # The minutes
    if minutes > 0:
        written_time.append(
            "{} minute{}".format(minutes, "s" if minutes > 1 else "")
        )

    # The seconds
    if seconds > 0:
        written_time.append(
            "{} second{}".format(seconds, "s" if seconds > 1 else "")
        )

    # Formatting the written time
    if len(written_time) == 0:
        return "no time"

    if len(written_time) == 1:
        return written_time[0]

    return ", ".join(written_time[:-1]) + " and " + written_time[-1]


def colorize_time(total_seconds):
    """Colorize the time.

    Args:
        total_seconds (int): the total number of seconds

    Returns:
        str: a colorized LaTeX string representation of time

    The time is displayed as ``HH:MM:SS``, but insignificant zeros are
    grayed-out.

    """
    # Formatting the time
    formatted_time = format_time(total_seconds)

    # Setting the color
    colored_time = formatted_time
    to_color = re.match("([0:]+)", formatted_time)
    if to_color is not None:
        colored_time = r"{\color{light_gray}"
        colored_time += formatted_time[:to_color.end()]
        colored_time += "}" + formatted_time[to_color.end():]

    # Formatting the time
    return colored_time
