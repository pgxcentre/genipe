
# This file is part of gwip.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


import unittest
from tempfile import TemporaryDirectory

from ..reporting.utils import *
from ..reporting.utils import _is_sanitized


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = "Copyright 2014, Beaulieu-Saucier Pharmacogenomics Centre"
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


class TestReportingUtils(unittest.TestCase):

    def setUp(self):
        """Setup the tests."""
        # Creating the temporary directory
        self.output_dir = TemporaryDirectory(prefix="gwip_test_")

    def tearDown(self):
        """Finishes the test."""
        # Deleting the output directory
        self.output_dir.cleanup()

    def test_sanitize_tex(self):
        """Tests the 'sanitize_tex' function."""
        # Sanitize all required characters
        string = "$%_}{&#~"
        sanitized = r"\$\%\_\}\{\&\#$\sim$"
        self.assertEqual(sanitized, sanitize_tex(string))

        # Sanitize only required characters (and not others)
        string = "ad$sf]35[fsdfs%f_s3f}546t%?%$634{&#dfs&?%%$/B~DF"
        sanitized = (r"ad\$sf]35[fsdfs\%f\_s3f\}546t\%?\%\$634\{\&\#dfs"
                     r"\&?\%\%\$/B$\sim$DF")
        self.assertEqual(sanitized, sanitize_tex(string))

        # Sanitizing a LaTeX command
        string = r"\texttt{test}\\"
        sanitized = (r"\textbackslash texttt\{test\}\textbackslash "
                     r"\textbackslash ")
        self.assertEqual(sanitized, sanitize_tex(string))

    def test_wrap_tex(self):
        """Tests the 'wrap_tex' function."""
        string = ". " * 80
        expected = ((". " * 35).rstrip() + "\n") * 2 + (". " * 10).rstrip()
        self.assertEqual(expected, wrap_tex(string))

    def test_format_tex(self):
        """Tests the 'format_tex' function."""
        # Testing all valid format
        valid_formats = ["texttt", "emph", "textbf", "textit"]
        string = "dummy string"
        for valid_format in valid_formats:
            expected = "\\" + valid_format + "{" + string + "}"
            self.assertEqual(expected, format_tex(string, valid_format))

        # Testing an invalid format
        with self.assertRaises(AssertionError) as cm:
            format_tex(string, "invalid_format")
        self.assertEqual(str(cm.exception), "invalid format")

        # Testing an un-sanitized string
        with self.assertRaises(AssertionError) as cm:
            format_tex("50%", valid_formats[0])
        self.assertEqual(str(cm.exception), "text not sanitized")

    def test_tex_inline_math(self):
        """Tests the 'tex_inline_math' function."""
        string = r"\sum_{n=1}^{100} \frac{1}{n}"
        expected = "$" + string + "$"
        self.assertEqual(expected, tex_inline_math(string))

    def test_is_sanitized(self):
        """Tests the '_is_sanitized' function."""
        # Not sanitized
        characters = "$%_}{&#~"
        padding = "pad"
        for character in characters:
            self.assertFalse(_is_sanitized(padding + character + padding))
        self.assertFalse(_is_sanitized(r"\ %"))

        # Sanitized
        string = r"ksj\$h\}f2\%o39\_875n\{9nt\&98\#34"
        self.assertTrue(_is_sanitized(string))

    def test_create_tabular(self):
        """Tests the 'create_tabular' function."""
        # Just testing the assertions for now
        with self.assertRaises(AssertionError) as cm:
            create_tabular(
                template=None,
                header=[1, 2, 3, 4, 5],
                data=None,
                header_multicol=[1, 1, 1, 1],
                col_align=None,
            )
        self.assertEqual(str(cm.exception), "len(header) != len(multicol)")

        with self.assertRaises(AssertionError) as cm:
            create_tabular(
                template=None,
                header=[1, 2, 3, 4, 5],
                data=None,
                header_multicol=[2, 1, 1, 1, 1],
                col_align=["c"] * 5,
            )
        self.assertEqual(str(cm.exception), "len(align) != number of columns")

    def test_create_float(self):
        """Tests the 'create_float' function."""
        # Just testing the assertions for now
        with self.assertRaises(AssertionError) as cm:
            create_float(
                template=None,
                float_type="unknown_type",
                caption=None,
                label=None,
                content=None,
                placement="H",
            )
        self.assertEqual(str(cm.exception), "invalid float type")

        with self.assertRaises(AssertionError) as cm:
            create_float(
                template=None,
                float_type="figure",
                caption=None,
                label=None,
                content=None,
                placement="Hb",
            )
        self.assertEqual(str(cm.exception), "placement 'H' should be alone")

        with self.assertRaises(AssertionError) as cm:
            create_float(
                template=None,
                float_type="figure",
                caption=None,
                label=None,
                content=None,
                placement="xb!",
            )
        self.assertEqual(str(cm.exception), "invalid placement")

    @unittest.skip("Test not implemented")
    def test_format_time(self):
        """Tests the 'format_time' function."""
        self.fail("Test not implemented")

    @unittest.skip("Test not implemented")
    def test_colorize_time(self):
        """Tests the 'colorize_time' function."""
        self.fail("Test not implemented")
