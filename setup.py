from setuptools import setup

setup(
        name='qcovid',
        version='0.1.1',
        description='amplicon based QC data extraction for sars-cov-2',
        scripts=['qcovid/bin_amplicons.py', 'qcovid/self_qc.py'],
)
