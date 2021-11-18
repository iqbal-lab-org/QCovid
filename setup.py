from setuptools import setup, find_packages

setup(
    name="qcovid",
    version="0.1.2",
    description="amplicon based QC data extraction for sars-cov-2",
    package_dir={"": "."},
    packages=find_packages(where="."),
    scripts=[
        "qcovid/bin_amplicons.py",
        "qcovid/self_qc.py",
        "qcovid/detect_primers.py",
    ],
)
