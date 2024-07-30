
# This file is part of genipe.
#
# GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = "Copyright 2014, Beaulieu-Saucier Pharmacogenomics Centre"
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


__all__ = ["GenipeError"]


class GenipeError(Exception):
    """An Exception raised in case of a problem."""
    def __init__(self, msg):
        """Construction of the GenipeError class."""
        self.message = str(msg)

    def __str__(self):
        return self.message
