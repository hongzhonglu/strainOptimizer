from setuptools import setup, find_packages
import sys

# files = ["model/data/*"]

with open('README.rst') as readme_file:
    readme = readme_file.read()

requirements = ['cobra','pytfa', 'matplotlib', 'pandas', 'tdqm']

setup_requirements = requirements + ['pytest-runner']
test_requirements = requirements + ['pytest', 'cplex']
install_requirements = requirements

setup(
    name='ETFLdesigner',
    version='0.0.0',
    python_requires='>=3.8',
    packages=find_packages(),
    include_package_data=True,
    zip_safe=False,
    install_requires=install_requirements,
    setup_requires=setup_requirements,
    tests_require=test_requirements,
    author='Haoyu Wang, Hongzhong Lu',
    author_email='wanghy@dicp.ac.cn',
    description='ETFLdesigner - Strain design algorithm for ETFL models',
    license='Apache License Version 2.0',
    keywords='strain design, ME-model',
    url='https://github.com/hongzhonglu/ETFLdesigner',
    long_description=readme,
    test_suite='tests',
)
