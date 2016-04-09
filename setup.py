from setuptools import setup
from mykrobe.version import __version__
__version__ = __version__[1:].split("-")[0]

setup(
    name='mykrobe',
    version=__version__,
    packages=[
        'mykrobe',
        'mykrobe.cmds',
        'mykrobe.predict',
        'mykrobe.metagenomics'],
    license='https://github.com/iqbal-lab/Mykrobe-predictor/blob/master/LICENSE',
    url='https://github.com/iqbal-lab/Mykrobe-predictor',
    description='.',
    author='Phelim Bradley, Zamin Iqbal',
    author_email='wave@phel.im, zam@well.ox.ac.uk',
    install_requires=["mykatlas"],
    entry_points={
            'console_scripts': [
                'mykrobe = mykrobe.mykrobe_predictor:main'
            ]},
    package_data={
        'mykrobe': [
            'data/predict/*/*',
            'data/predict/taxon_coverage_threshold.json',
            'data/phylo/*',
            'data/panels/tb-species-160330.fasta.gz',
            'data/panels/staph-species-160227.fasta.gz',
            'data/panels/tb-amr-walker_2015.fasta.gz',
            'data/panels/tb-amr-bradley_2015.fasta.gz',
            'data/panels/staph-amr-probe_set_v0_3_1_2.fasta.gz',
        ]})
