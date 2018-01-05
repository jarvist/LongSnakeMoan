LongSnakeMoan
=============

Codes to predict DoS of a polymer chain in vacuum with a one-dimensional tightbinding (Huckel model).

## Installation

The codes are still slightly dependent on the 'pythtb' module from Rutgers university (http://physics.rutgers.edu/pythtb/). Almost none of this library is used, but it's definitely recommended if you're interested in tight binding!

You can try and install this library automatically by pip,

```
pip install pythtb
```

Codes written by Jarvist Moore Frost 2012-2014.

Work based on this code was presented at MRS Boston 2012 (O8.08), slides are
here:
https://speakerdeck.com/jarvist/beta-phase

The code should produce live plots like this:

![Screenshot](https://raw.github.com/jarvist/LongSnakeMoan/master/2012-10-25-13h03-BetaPhasePyTB.png)

The continuation of this work (2013 onwards), including a Sturm-sequency based linear-in-time and constant-in-memory 1D Hamiltonian solver was continued in the Julia language elsewhere: 
https://github.com/jarvist/Teclo.jl
