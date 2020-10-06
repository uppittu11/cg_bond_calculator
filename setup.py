from setuptools import setup, find_packages
import sys

try:
    import mdtraj
except ImportError:
    print('Building and running cg_bond_calculator requires mdtraj. See '
          'http://mdtraj.org/latest/installation.html for help!')
    sys.exit(1)

setup(name='cg_bond_calculator',
      version='0.1',
      description=('Using MDTraj libraries and user-defined mapping '
        'files, convert an atomistic trajectory to a coarse-grained '
        'trajectory'),
      url='https://github.com/uppittu11/cg_bond_calculator',
      author='Parashara Shamaprasad',
      author_email='p.shama@vanderbilt.edu',
      license='MIT',
      packages=find_packages(),
      package_dir={'cg_bond_calculator': 'cg_bond_calculator'},
      include_package_data=True,
      install_requires=["mdtraj", "numpy", "scipy", "unyt", "networkx"],
)
