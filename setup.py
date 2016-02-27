from setuptools import setup
from pip.req import parse_requirements

install_reqs = parse_requirements('requirements.txt')
reqs = [str(ir.req) for ir in install_reqs]

setup(
    name='atlas',
    version='0.1',
    packages=['atlas'],
    license='MIT',
    url='http://github.com/phelimb/atlas',
    description='.',
    author='Phelim Bradley',
    author_email='wave@phel.im',
    install_requires=reqs,
    entry_points={
        'console_scripts': [
            'mykrobe = atlas.main:main',
        ]},
    package_data={'atlas': ['data/*']}        
)