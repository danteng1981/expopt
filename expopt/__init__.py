"""
ExpOpt: Molecular Fragment Replacement and Splicing for Drug Discovery

This package provides functionality for:
- RGroup replacement: Finding candidate fragments to replace target fragments
- CoreHopping: Replacing the core scaffold of molecules
- RGroup splicing: Connecting replacement fragments at attachment points
- Core splicing: Enumerating connection point combinations for cores
"""

from expopt.rgroup import RGroupReplacer
from expopt.core_hopping import CoreHopper
from expopt.splicing import RGroupSplicer, CoreSplicer

__version__ = "0.1.0"
__all__ = ["RGroupReplacer", "CoreHopper", "RGroupSplicer", "CoreSplicer"]
