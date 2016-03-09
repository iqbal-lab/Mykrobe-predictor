from setuptools import setup
from pip.req import parse_requirements
import atlas.version
install_reqs = parse_requirements('requirements.txt')
reqs = [str(ir.req) for ir in install_reqs]

setup(
    name='mykrobe',
    version='0.3.0.1',
    packages=[
        'atlas',
        'atlas.cmds',
        'atlas.pheno',
        'atlas.typing',
        'atlas.typing.typer',
        'atlas.typing.models',
        'atlas.schema',
        'atlas.schema.models',
        'atlas.stats',
        'atlas.cortex',
        'atlas.metagenomics'],
    license='https://github.com/iqbal-lab/Mykrobe-predictor/blob/master/LICENSE',
    url='http://github.com/phelimb/atlas',
    description='.',
    author='Phelim Bradley, Zamin Iqbal',
    author_email='wave@phel.im, zam@well.ox.ac.uk',
    install_requires=reqs,
    entry_points={
            'console_scripts': [
                'mykrobe = atlas.main:main',
            ]},
    package_data={
        'atlas': [
            'data/predict/*/*',
            'data/predict/taxon_coverage_threshold.json',
            'data/phylo/*',
            'data/panels/tb-species-160227.fasta',
            'data/panels/tb-amr-walker_2015.fasta',
            'data/panels/tb-amr-bradley_2015.fasta',
        ]})
