from setuptools import setup, find_packages

setup(
    name="expopt",
    version="0.1.0",
    description="Molecular fragment replacement and splicing for drug discovery",
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[
        "rdkit>=2023.3.1",
        "numpy>=1.24.0",
        "pandas>=2.0.0",
    ],
    extras_require={
        "test": ["pytest>=7.0.0"],
    },
)
