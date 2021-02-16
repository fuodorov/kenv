.. kenv documentation master file, created by
   sphinx-quickstart on Tue Feb 16 11:59:36 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

KENV documentation
======================================================

KENV (Kapchinscky ENVelope) is a solver for the Kapchinsky-Vladimirsky envelope equation for electron beam with space charge.

It is well suited for electron beam posting in an accelerator path with solenoidal and quadrupole focusing.

The distinctive feature of canvases for finding the envelope in comparison with other programs (various PIC-codes) is its speed.

Contents of the documentation
-----------------------------

If you are new to KENV, we **strongly recommend** that you read the
section :doc:`overview/overview` first, so as to have a basic understanding of
what the code does.

You can then see the section :doc:`install/installation` and
:doc:`tutorial/tutorial`, to get started with using KENV. For more
information, the section :doc:`api_reference/api_reference` lists the main objects
that are accessible through KENV.

.. toctree::
  :maxdepth: 1

  overview/overview
  install/installation
  tutorial/tutorial
  api_reference/api_reference

Contributing to KENV
----------------------

KENV is open-source, and the source code is hosted `here <http://github.com/fuodorov/kenv>`_, on
Github.

We welcome contributions to the code!

Research & Attribution
----------------------

KENV was developed by `Vyacheslav Fedorov <https://fuodorov.github.io>`_, `Danila Nikiforov <mailto:nikdanila@bk.ru>`_ and `Alexey Petrenko <https://inp.nsk.su/~petrenko/>`_  at `Budker Institute of Nuclear Physics <https://inp.nsk.su>`_,

KENV's algorithms are documented in following scientific publications:

   * High-Current Electron-Beam Transport in the LIA-5 Linear Induction Accelerator (original paper):
     `D. A. Nikiforov et al., Phys. Part. Nuclei Lett. 17, 197-203 (2020) <https://link.springer.com/article/10.1134/S1547477120020156>`_

If you use KENV for your research project: that's great! We are
very pleased that the code is useful to you!

If your project even leads to a scientific publication, please consider citing at least KENV's original paper.
If your project uses the more advanced algorithms, please consider citing the respective publications in addition.
