Simoji -- Because an emoji is better than your atom name
========================================================

Likely by accident, Gromacs has no problem using emojis where it expects
strings. This means you can use emojis as atom types, atom names, residue names
or molecule types for much more insightful simulations. Simoji replaces all
names in a topology and in a structure by emojis to run such emoji molecular
dynamics simulations.

Usage: `python3 simoji.py [options]`

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
original names and the chosen emojis. This correspondence table can be provided
as an input with the `-t` option; this bypasses the random choice of the
emojis.

If the topology file is provided with the `-p` option, and if it includes ITP
files, the translated content of these ITP file is concatenated in the
resulting topology.

As an example, here is a extract of the POPC topology in the Martini coarse
grained force field [1]:

```
[moleculetype]
; molname      nrexcl
   POPC          1

[atoms]
; id    type    resnr   residu  atom    cgnr    charge
1       Q0      1       POPC    NC3     1       1.0
2       Qa      1       POPC    PO4     2       -1.0
3       Na      1       POPC    GL1     3       0
4       Na      1       POPC    GL2     4       0
5       C1      1       POPC    C1A     5       0
6       C1      1       POPC    C2A     6       0
7       C3      1       POPC    D3A     7       0
8       C1      1       POPC    C4A     8       0
9       C1      1       POPC    C1B     9       0
10      C1      1       POPC    C2B     10      0
11      C1      1       POPC    C3B     11      0
12      C1      1       POPC    C4B     12      0
```

Here is the same fraction emojified:

```
[moleculetype]
; molname      nrexcl
   ⚛          1

[atoms]
; id    type    resnr   residu  atom    cgnr    charge
1       ⛽      1       ⚛       🍲      1       1.0
2       📝      1       ⚛       🌭      2       -1.0
3       🐍      1       ⚛       🍹      3       0
4       🐍      1       ⚛       🤐      4       0
5       🕷       1       ⚛       🕓      5       0
6       🕷       1       ⚛       🚊      6       0
7       🚕      1       ⚛       🚓      7       0
8       🕷       1       ⚛       🔇      8       0
9       🕷       1       ⚛       🕝      9       0
10      🕷       1       ⚛       🗨       10      0
11      🕷       1       ⚛       🎑      11      0
12      🕷       1       ⚛       👂      12      0
```

Note that simoji requires python 3. It will not work on python 2.

You should probably not use emojis in your serious simulations 😝.

    [1] The MARTINI Force Field:  Coarse Grained Model for Biomolecular Simulations
        Siewert J. Marrink, H. Jelger Risselada, Serge Yefimov, D. Peter Tieleman, and, and Alex H. de Vries
        The Journal of Physical Chemistry B 2007 111 (27), 7812-7824
        DOI: http://dx.doi.org/10.1021/jp071097f 
