# distancecalc

Python script to calculate distances/angles from a VASP structure file.

**Usage**

distancecalc.py [dist/angle] [options] [arguments]

**Mandatory arguments**

dist  : Distance calculation mode

angle : Angle calculation mode

**Optional arguments**

-i [input] : Designating the input file. (Default : POSCAR)

-o [output] : Designating the name of output file. If not set, results will be printed out on the screen.

-s [list of str] : Setting the atomic species to calculate distance/angle. In angle calculation, be aware of the sequence of items. (ex: dist -s Si Si, angle -s C O C)

-t [float]  : Setting the tolerance value to distinguish duplicates. (Default : 0.1 Å for distance, 0.01 degree for angle)

-d [float] : Setting the distance tolerance value to calculate angles. If given atoms are positioned near than the tolerance, angles are not calculated.
