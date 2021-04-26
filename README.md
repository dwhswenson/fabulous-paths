 [![Documentation Status](https://readthedocs.org/projects/fabulous-paths/badge/?version=latest)](https://fabulous-paths.readthedocs.io/en/latest/?badge=latest)
 [![tests](https://github.com/dwhswenson/fabulous-paths/actions/workflows/tests.yml/badge.svg)](https://github.com/dwhswenson/fabulous-paths/actions/workflows/tests.yml)
[![codecov](https://codecov.io/gh/dwhswenson/fabulous-paths/branch/main/graph/badge.svg?token=Mhtza0eAID)](https://codecov.io/gh/dwhswenson/fabulous-paths)

# fabulous-paths

*Tools to integrate OpenPathSampling and FABULOUS.*

``fabulous-paths`` is a set of tools that bridge between
[OpenPathSampling](http://openpathsampling.org), especially the
[OpenPathSampling CLI](http://openpathsampling-cli.readthedocs.org), and the
analysis tools provided by
[FABULOUS](https://github.com/Ensing-Laboratory/FABULOUS/).

To learn more about FABULOUS, read the original paper here:

* [F. Hooft, A. Pérez de Alba Ortíz, B. Ensing. "Discovering Collective
  Variables of Molecular Transitions via Genetic Algorithms and Neural
  Networks." J. Chem. Theory Comput. **17**, 2294
  (2021).](https://dx.doi.org/10.1021/acs.jctc.0c00981)

## Installation

``fabulous-paths`` can be installed ... (eventually this will give ``pip`` and
``conda`` instructions). Installing ``fabulous-paths`` will automatically also
install OpenPathSampling, the OpenPathSampling CLI, and FABULOUS, if you do not
already have these installed.

## Citing

When using ``fabulous-paths``, please cite the following papers (for various
functionality included):

* **FABULOUS**: https://doi.org/10.1021/acs.jctc.0c00981
* **OpenPathSampling**: https://doi.org/10.1021/acs.jctc.8b00626
<!--* **OpenPathSampling CLI**: ???-->

## Support and Development

``fabulous-paths`` is an open source project, released under the GNU LGPL,
version 3.0 or (at your option) any later version. Development takes place in
public at http://github.com/openpathsampling/paths-fabulous.

As this package is mostly a thin adapter, problems are more likely to arise
through either OpenPathSampling or FABULOUS. Please first contact those
packages through their help forums:

* Getting help with OpenPathSampling
* Getting help with FABULOUS
