# Kapchinsky ENVelope (KENV)
[![PyPI version](https://badge.fury.io/py/kenv.svg)](https://badge.fury.io/py/kenv)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/fuodorov/kenv/dev?filepath=notebooks%2F00_introduction.ipynb)
## Solver of the Kapchinsky-Vladimirsky envelope equation

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

```
pip install kenv
```

## Test

Comparison of PIC-codes Astra (blue) and WARP (red), SAM (black) with KENV 
![VS](notebooks/output/demo/vs_kenv.gif)
