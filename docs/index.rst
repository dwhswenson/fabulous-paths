.. fabulous-paths documentation master file, created by
   sphinx-quickstart on Thu Apr  8 18:13:45 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

fabulous-paths
==============

``fabulous-paths`` connects OpenPathSampling, a package for path sampling
simulations, to FABULOUS, a machine learning approach to identify important
collective variables in transitions.

This package is just the glue between those two projects. Most of the
documentation can be found in the individual packages. The documentation
includes details on how to use this code as a library and specifics of the
OpenPathSampling CLI plugin.

To install ``fabulous-paths`` just ???. You can verify that the installation
worked by running the command ``openpathsampling fabulous --help``.

Details on how to use OpenPathSampling to perform path sampling simulations
can be found in `the OpenPathSampling documentation
<http://openpathsampling.org>`_. Details on how to use FABULOUS to extract
important CVs from path sampling data can be found in the `documentation for
FABULOUS <https://github.com/Ensing-Laboratory/FABULOUS>`_.

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   usage
   api

