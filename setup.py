from setuptools import setup, find_packages
from src.version import __version__

setup(
    name="ReactantGenerator",
    version=__version__,
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    entry_points={'console_scripts': [
        'reactant_generator=ReactantGenerator.main:main',
    ]},
    python_requires='=3.10',
    install_requires=[
        'numpy==1.24.3',
        'geatpy==2.7.0',
    ],
    description="ReactantGenerator",
    author="TBA",
)