#! /usr/bin/env python

import os
import setuptools


def read_requirements():
    requirements = []
    with open('./requirements.txt', 'rt', encoding='utf-8') as file:
        for line in file:
            # exclude comments
            line = line[: line.find("#")] if "#" in line else line
            # clean
            line = line.strip()
            if line:
                requirements.append(line)
    return requirements


data_files_names = ["docs/README.md", "LICENSE.txt"]
data_files_locations = [
    ('.', [f]) if os.path.exists(f) else ('.', ["../" + f])
    for f in data_files_names
]
deps = read_requirements()
excalibur = os.path.join('excalibur', '__init__.py')
version = 'esp-git-rev'
with open(os.path.join(os.path.dirname(__file__), excalibur)) as f:
    t = f.read()
t = t.replace('${ESP_GIT_REV}', version)
with open(os.path.join(os.path.dirname(__file__), excalibur), 'tw') as f:
    f.write(t)
with open(data_files_names[0], "rt", encoding='utf-8') as f:
    description = f.read()
with open(data_files_names[1], "rt", encoding='utf-8') as f:
    license = f.read()
setuptools.setup(
    name='excalibur',
    version='2.0.0',
    packages=setuptools.find_packages(),
    setup_requires=deps,
    src_root=os.path.abspath(os.path.dirname(__file__)),
    install_requires=deps,
    package_data={},
    author='Gael Roudier',
    author_email='Gael.Roudier@jpl.caltech.edu',
    description='Exoplanet Calibration Bayesian Unified Retrieval Pipeline (ExCaliBUR).',
    license=license,
    keywords='baysian mcmc',
    url='https://github-fn.jpl.nasa.gov/EXCALIBUR/esp',
    zip_safe=False,
    classifiers=[
        "Programming Language :: Python :: 3",
        'Programming Language :: Python :: 3.12',
        "Operating System :: OS Independent",
        'License :: Free To Use But Restricted',
        'Development Status :: 5 - Production/Stable',
    ],
    data_files=data_files_locations,
    long_description=description,
    long_description_content_type="text/markdown",
    python_requires='>3.10',
)
