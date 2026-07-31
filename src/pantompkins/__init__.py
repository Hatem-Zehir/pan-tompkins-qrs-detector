"""
pantompkins — Python implementation of the Pan-Tompkins QRS detection algorithm.

Python port by Dr. Hatem Zehir.
Original MATLAB implementation by Hooman Sedghamiz (2018).
"""

from .pan_tompkins import pan_tompkins

__version__ = "1.0.0"
__author__ = "Hatem Zehir"
__all__ = ["pan_tompkin"]