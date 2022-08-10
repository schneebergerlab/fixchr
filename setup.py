from setuptools import setup, Extension
from fixchr import __version__

setup(
    name="fixchr",
    version='{}'.format(__version__),
    description='Filter and reorient genomes to get homologous chromosomes',
    author='Manish Goel',
    author_email='mnshgl0110@gmail.com',
    url='https://github.com/schneebergerlab/fixchr',
    license='MIT License',
    license_files=('LICENSE',),
    packages=["fixchr", "fixchr.scripts"],
    py_modules=["fixchr.scripts.fixchr",
                "fixchr.scripts.func",
                "fixchr.scripts.dotplot"],
    scripts=['bin/fixchr', 'bin/dotplot'],
    long_description=open('README.rst').read(),
)