# Kapchinsky ENVelope (KENV)
[![PyPI version](https://badge.fury.io/py/kenv.svg)](https://badge.fury.io/py/kenv)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/fuodorov/kenv/master?filepath=notebooks%2Fintroduction.ipynb)
## The solver of Kapchinsky-Vladimirsky envelope equation for electron beam with space charge.

<a href=mailto:fuodorov1998@gmail.com>V. Fedorov</a>, <a href=mailto:nikdanila@bk.ru>D. Nikiforov</a>, <a href=http://www.inp.nsk.su/~petrenko/>A. Petrenko</a>, (Novosibirsk, 2019)

## Overview

KENV is a solver code for the equation of the envelope of an electron beam with the Kapchinsky-Vladimirsky distribution for accelerator physics.

It is particularly suitable for accelerating an electron beam in direct channels with solenoidal and quadrupole focusing.

In order to use the KENV code correctly, it is **important to read the [Wiki](https://github.com/fuodorov/kenv/wiki).**

## Algorithm

The algorithm reduces to lowering the order of the Kapchinsky-Vladimirsky differential equations to the first and subsequent integration.

## Language

KENV completely written in Python.

## Installation

In order to quickly install all the required [Python](https://www.python.org/downloads/) libraries in the new environment, just download [requirements](https://github.com/fuodorov/kenv/blob/master/requirements.txt) and run the command on the command line:

```
pip install -r requirements.txt
```

## Publications

[Publication](http://www1.jinr.ru/Pepan_letters/panl_2020_2/13_nikifor.pdf) in Particles and Nuclei. 
