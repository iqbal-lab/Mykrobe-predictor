from setuptools import setup
from pip.req import parse_requirements
import atlas.version
install_reqs = parse_requirements('requirements.txt')
reqs = [str(ir.req) for ir in install_reqs]

setup(
    name='mykrobe',
    version=str(atlas.version.__version__),
    packages=['atlas', 'atlas.cmds', 'atlas.typing', 'atlas.typing.typer'],
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
    package_data={'atlas': ['data/*']}
)
