
# This file is part of gwip.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = "Copyright 2014, Beaulieu-Saucier Pharmacogenomics Centre"
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


__all__ = ["config_jinja2", "sanitize_tex", "format_tex", "wrap_tex",
           "create_tabular", "create_float", "tex_inline_math", ]


import re
from textwrap import wrap

import jinja2


_char_mod = {
    "~": r"$\sim$",
}
_escaped_char = ["$", "%", "_", "}", "{", "&", "#"]
_valid_tex_formats = {"texttt", "emph", "textbf", "textit"}


def config_jinja2():
    """Configure the jinja2 environment."""
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
    """Sanitize TeX text."""
    # Escaping
    sanitized_tex = re.sub(r"([{}])".format("".join(_escaped_char)),
                           r"\\\g<1>", original_text)

    # Replacing
    for char, mod in _char_mod.items():
        sanitized_tex = sanitized_tex.replace(char, mod)

    return sanitized_tex


def wrap_tex(original_text):
    """Wraps the text."""
    return "\n".join(wrap(original_text))


def format_tex(text, tex_format):
    """Change the TeX text format."""
    assert tex_format in _valid_tex_formats
    assert _is_sanitized(text)

    return r"\%s{%s}" % (tex_format, text)


def tex_inline_math(content):
    """Creates an inline mathematical formula in TeX."""
    return "${}$".format(content)


def _is_sanitized(text):
    """Check if text is sanitized."""
    return re.search(r"[^\\][{}]".format("".join(_escaped_char)), text) is None


def create_tabular(template, header, data, header_multicol=None,
                   col_align=None):
    """Creates a TeX tabular."""
    if header_multicol is None:
        header_multicol = [1 for i in header]

    # Getting the number of columns
    nb_col = sum(header_multicol)

    if col_align is None:
        col_align = ["c"] * nb_col

    # Checking that the number of columns holds
    assert len(header) == nb_col
    assert len(col_align) == nb_col

    # Generating the tabular data
    tabular_data = {
        "col_alignments": "".join(col_align),
        "header_data":    zip(header, header_multicol),
        "tabular_data":   data,
    }

    # Rendering
    return template.render(**tabular_data)


def create_float(template, float_type, caption, label, content, placement="H"):
    """Creates a TeX float."""
    # Some assertions
    assert float_type in ["figure", "table"]

    # Generating the float data
    float_data = {
        "float_type":      float_type,
        "float_placement": placement,
        "float_caption":   caption,
        "float_label":     label,
        "table_content":   content,
    }

    # Rendering
    return template.render(**float_data)
