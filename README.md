# expopt

Molecular fragment replacement and splicing for drug discovery and chemistry.

## Features

### 1. RGroup Replacement
Find candidate fragments to replace target fragments while maintaining similar properties and ligand site matching.

**Workflow:**
1. **Priority: Exact Matching** - First attempts to find exact matches in the database
2. **Fallback: Candidate Generation** - If no exact match, generates similar candidates from the database

```python
from expopt import RGroupReplacer
from expopt.rgroup import FragmentDatabase

# Create and populate database
db = FragmentDatabase()
db.add_fragments_from_list(["c1ccccc1", "c1ccc(C)cc1", "c1ccc(O)cc1"])

# Find replacements
replacer = RGroupReplacer(database=db)
candidates = replacer.find_replacements(
    "c1ccc(F)cc1",  # target fragment
    similarity_threshold=0.3,
    max_candidates=10
)
```

### 2. Core Hopping
Replace the core scaffold of molecules while preserving side chains and key functional groups.

```python
from expopt import CoreHopper
from expopt.core_hopping import CoreDatabase

# Create and populate database
db = CoreDatabase()
db.add_cores_from_list([
    "c1ccc2ccccc2c1",   # naphthalene
    "c1ccc2ncccc2c1",   # quinoline
])

# Find hopping candidates
hopper = CoreHopper(database=db)
candidates = hopper.find_hopping_candidates(
    "c1ccc2ccccc2c1C",  # query molecule
    similarity_threshold=0.3,
    max_candidates=10
)
```

### 3. RGroup Splicing
Connect replacement fragments at specified attachment points via single bonds.

```python
from expopt import RGroupSplicer

splicer = RGroupSplicer()

# Basic splicing with atom indices
result = splicer.splice(
    "c1ccccc1",  # scaffold
    "C",          # fragment (methyl)
    0,            # scaffold attachment index
    0             # fragment attachment index
)

# Splicing using dummy atoms
result = splicer.splice_with_dummy_atoms(
    "*c1ccccc1",  # scaffold with attachment point
    "*C"          # fragment with attachment point
)
```

### 4. Core Splicing
Enumerate all reasonable connection point combinations between candidate cores and R-groups, scoring by fingerprint similarity.

```python
from expopt import CoreSplicer

splicer = CoreSplicer()

# Find best connections by similarity to target
results = splicer.find_best_connections(
    "c1ccccc1",       # core scaffold
    ["C", "C"],       # R-groups
    "Cc1ccc(C)cc1",   # target molecule
    max_results=10
)

# Batch processing multiple cores
results = splicer.batch_splice_and_score(
    ["c1ccccc1", "c1ccncc1"],  # candidate cores
    ["C"],                      # R-groups
    "c1ccc(C)cc1",              # target molecule
    top_n=10
)
```

## Installation

```bash
pip install -r requirements.txt
```

### Requirements
- Python >= 3.8
- RDKit >= 2023.3.1
- NumPy >= 1.24.0
- Pandas >= 2.0.0

## Running Tests

```bash
pip install pytest
pytest tests/ -v
```

## API Reference

### Classes

#### `RGroupReplacer`
- `exact_match(target_smiles)` - Find exact match in database
- `generate_candidates(target_smiles, ...)` - Generate similar candidates
- `find_replacements(target_smiles, ...)` - Find replacements (exact first, then similar)

#### `CoreHopper`
- `extract_core(mol_smiles)` - Extract Murcko scaffold
- `decompose_molecule(mol_smiles)` - Decompose into core and R-groups
- `find_hopping_candidates(mol_smiles, ...)` - Find scaffold hopping candidates

#### `RGroupSplicer`
- `splice(scaffold, fragment, ...)` - Splice fragment onto scaffold
- `splice_with_dummy_atoms(scaffold, fragment)` - Splice using dummy atoms
- `splice_multiple_fragments(scaffold, fragments)` - Splice multiple fragments

#### `CoreSplicer`
- `enumerate_connections(core, rgroups, ...)` - Enumerate connection combinations
- `splice_core_with_rgroups(core, rgroups, connections)` - Splice with specified connections
- `find_best_connections(core, rgroups, target, ...)` - Find best connections by similarity
- `batch_splice_and_score(cores, rgroups, target, ...)` - Process multiple cores
