Simoji -- Because an emoji is better than your atom name
========================================================

Likely by accident, Gromacs has no problem using emojis where it expect
strings. This means you can use emojis as atom types, atom names, residue names
or molecule types for much more insighful simulations. Simoji replaces all
names in a topology and in a structure by emoji to run such emoji molecular
dynamics simulations.

Usage: simoji.py [options]

Options:

* `-f`: input structure (.gro)
* `-p`: input topology (.top)
* `-t`: input correspondence table
* `-o`: output structure (.gro)
* `-op`: output topology (.top)
* `-ot`: output correspondence table

The emojis are chosen at random but coincide between the output structure and
the output topology so the files can be used together in a simulation. If the
`-ot` option is provided, simoji writes a correspondence table between the
original names ans the chosen emojis. This correspondence table can be provided
as an input with the `-t` option; this bypasses the random choice of the
emojis.

If the topology file is provided with the `-p` option, and if it includes ITP
files, the translated content of these ITP file is concatenated in the
resulting topology.

You should probably not use emojis in your serious simulations 😝.
