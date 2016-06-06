
# This file is part of genipe.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


import unittest

from ..reporting import utils as report_utils


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = "Copyright 2014, Beaulieu-Saucier Pharmacogenomics Centre"
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


__all__ = ["TestReportingUtils"]


class TestReportingUtils(unittest.TestCase):

    def test_sanitize_tex(self):
        """Tests the 'sanitize_tex' function."""
        # Sanitize all required characters
        string = "$%_}{&#~"
        sanitized = r"\$\%\_\}\{\&\#$\sim$"
        self.assertEqual(sanitized, report_utils.sanitize_tex(string))

        # Sanitize only required characters (and not others)
        string = "ad$sf]35[fsdfs%f_s3f}546t%?%$634{&#dfs&?%%$/B~DF"
        sanitized = (r"ad\$sf]35[fsdfs\%f\_s3f\}546t\%?\%\$634\{\&\#dfs"
                     r"\&?\%\%\$/B$\sim$DF")
        self.assertEqual(sanitized, report_utils.sanitize_tex(string))

        # Sanitizing a LaTeX command
        string = r"\texttt{test}\\"
        sanitized = (r"\textbackslash texttt\{test\}\textbackslash "
                     r"\textbackslash ")
        self.assertEqual(sanitized, report_utils.sanitize_tex(string))

    def test_wrap_tex(self):
        """Tests the 'wrap_tex' function."""
        string = ". " * 80
        expected = ((". " * 35).rstrip() + "\n") * 2 + (". " * 10).rstrip()
        self.assertEqual(expected, report_utils.wrap_tex(string))

    def test_format_tex(self):
        """Tests the 'format_tex' function."""
        # Testing all valid format
        valid_formats = ["texttt", "emph", "textbf", "textit"]
        string = "dummy string"
        for valid_format in valid_formats:
            expected = "\\" + valid_format + "{" + string + "}"
            self.assertEqual(
                expected, report_utils.format_tex(string, valid_format),
            )

        # Testing an invalid format
        with self.assertRaises(AssertionError) as cm:
            report_utils.format_tex(string, "invalid_format")
        self.assertEqual(str(cm.exception), "invalid format")

        # Testing an un-sanitized string
        with self.assertRaises(AssertionError) as cm:
            report_utils.format_tex("50%", valid_formats[0])
        self.assertEqual(str(cm.exception), "text not sanitized")

    def test_tex_inline_math(self):
        """Tests the 'tex_inline_math' function."""
        string = r"\sum_{n=1}^{100} \frac{1}{n}"
        expected = "$" + string + "$"
        self.assertEqual(expected, report_utils.tex_inline_math(string))

    def test_is_sanitized(self):
        """Tests the '_is_sanitized' function."""
        # Not sanitized
        characters = "$%_}{&#~"
        padding = "pad"
        for character in characters:
            self.assertFalse(report_utils._is_sanitized(
                padding + character + padding,
            ))
        self.assertFalse(report_utils._is_sanitized(r"\ %"))

        # Sanitized
        string = r"ksj\$h\}f2\%o39\_875n\{9nt\&98\#34"
        self.assertTrue(report_utils._is_sanitized(string))

    def test_create_tabular(self):
        """Tests the 'create_tabular' function."""
        # Just testing the assertions for now
        with self.assertRaises(AssertionError) as cm:
            report_utils.create_tabular(
                template=None,
                header=[1, 2, 3, 4, 5],
                data=None,
                header_multicol=[1, 1, 1, 1],
                col_align=None,
            )
        self.assertEqual(str(cm.exception), "len(header) != len(multicol)")

        with self.assertRaises(AssertionError) as cm:
            report_utils.create_tabular(
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
            report_utils.create_float(
                template=None,
                float_type="unknown_type",
                caption=None,
                label=None,
                content=None,
                placement="H",
            )
        self.assertEqual(str(cm.exception), "invalid float type")

        with self.assertRaises(AssertionError) as cm:
            report_utils.create_float(
                template=None,
                float_type="figure",
                caption=None,
                label=None,
                content=None,
                placement="Hb",
            )
        self.assertEqual(str(cm.exception), "placement 'H' should be alone")

        with self.assertRaises(AssertionError) as cm:
            report_utils.create_float(
                template=None,
                float_type="figure",
                caption=None,
                label=None,
                content=None,
                placement="xb!",
            )
        self.assertEqual(str(cm.exception), "invalid placement")

    def test_format_time(self):
        """Tests the 'format_time' function."""
        seconds = [0, 1, 30]
        minutes = [0, 1, 30]
        hours = [0, 1, 30, 101]

        written_times = ["no time", "1 hour", "30 hours", "101 hours",
                         "1 minute", "1 hour and 1 minute",
                         "30 hours and 1 minute", "101 hours and 1 minute",
                         "30 minutes", "1 hour and 30 minutes",
                         "30 hours and 30 minutes", "101 hours and 30 minutes",
                         "1 second", "1 hour and 1 second",
                         "30 hours and 1 second", "101 hours and 1 second",
                         "1 minute and 1 second",
                         "1 hour, 1 minute and 1 second",
                         "30 hours, 1 minute and 1 second",
                         "101 hours, 1 minute and 1 second",
                         "30 minutes and 1 second",
                         "1 hour, 30 minutes and 1 second",
                         "30 hours, 30 minutes and 1 second",
                         "101 hours, 30 minutes and 1 second", "30 seconds",
                         "1 hour and 30 seconds", "30 hours and 30 seconds",
                         "101 hours and 30 seconds", "1 minute and 30 seconds",
                         "1 hour, 1 minute and 30 seconds",
                         "30 hours, 1 minute and 30 seconds",
                         "101 hours, 1 minute and 30 seconds",
                         "30 minutes and 30 seconds",
                         "1 hour, 30 minutes and 30 seconds",
                         "30 hours, 30 minutes and 30 seconds",
                         "101 hours, 30 minutes and 30 seconds"]

        i = 0
        for second in seconds:
            for minute in minutes:
                for hour in hours:
                    # Computing the total number of seconds
                    total_seconds = second + (minute * 60) + (hour * 3600)

                    # Checking the time
                    self.assertEqual(
                        "{:02d}:{:02d}:{:02d}".format(hour, minute, second),
                        report_utils.format_time(total_seconds),
                    )

                    # Checking the written time
                    self.assertEqual(
                        written_times[i],
                        report_utils.format_time(total_seconds,
                                                 written_time=True),
                    )
                    i += 1

    def test_colorize_time(self):
        """Tests the 'colorize_time' function."""
        # Testing 1 second
        self.assertEqual(r"{\color{light_gray}00:00:0}1",
                         report_utils.colorize_time(1))

        # Testing 10 seconds
        self.assertEqual(r"{\color{light_gray}00:00:}10",
                         report_utils.colorize_time(10))

        # Testing 1 minute
        self.assertEqual(r"{\color{light_gray}00:0}1:00",
                         report_utils.colorize_time(60))

        # Testing 10 minutes
        self.assertEqual(r"{\color{light_gray}00:}10:00",
                         report_utils.colorize_time(600))

        # Testing 1 hour
        self.assertEqual(r"{\color{light_gray}0}1:00:00",
                         report_utils.colorize_time(3600))

        # Testing 10 hours
        self.assertEqual("10:00:00",
                         report_utils.colorize_time(36000))

        # Testing 100 hours
        self.assertEqual("100:00:00",
                         report_utils.colorize_time(360000))
