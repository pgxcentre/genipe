

__all__ = ["config_jinja2", "sanitize_tex", "format_tex", "wrap_tex"]


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


def _is_sanitized(text):
    """Check if text is sanitized."""
    return re.search(r"[^\\][{}]".format("".join(_escaped_char)), text) is None
