"""
ExpOpt - A Python package for RGroup replacement and CoreHopping using RDKit.

This package provides tools for:
- RGroup replacement workflow
- CoreHopping workflow  
- Fragment/Scaffold splicing
"""

__version__ = "0.1.0"

from expopt.rgroup.replacement import RGroupReplacement
from expopt.core_hopping.hopping import CoreHopping
from expopt.splicing.rgroup_splicing import RGroupSplicing
from expopt.splicing.core_splicing import CoreSplicing

__all__ = [
    "RGroupReplacement",
    "CoreHopping", 
    "RGroupSplicing",
    "CoreSplicing",
]
