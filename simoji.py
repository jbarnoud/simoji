#!/usr/bin/env python

"""
Replace names by emojis in GROMACS input files.
"""

import argparse
import collections
import configparser
import itertools
import os
import random

__author__ = 'Jonathan Barnoud'
__version__ = '0.1.0dev0'

# List of emoji originally obtained from "emoji_list":
# <https://github.com/vincentmwong/emoji_list>.
# The "emoji_list" package is released under the MIT licence.
EMOJI_RANGES = [
    (8205, 8206), (8252, 8253), (8265, 8266), (8482, 8483), (8505, 8506),
    (8596, 8602), (8617, 8619), (8986, 8988), (9000, 9001), (9167, 9168),
    (9193, 9204), (9208, 9211), (9410, 9411), (9642, 9644), (9654, 9655),
    (9664, 9665), (9723, 9727), (9728, 9733), (9742, 9743), (9745, 9746),
    (9748, 9750), (9752, 9753), (9757, 9758), (9760, 9761), (9762, 9764),
    (9766, 9767), (9770, 9771), (9774, 9776), (9784, 9787), (9792, 9793),
    (9794, 9795), (9824, 9825), (9827, 9828), (9829, 9831), (9832, 9833),
    (9851, 9852), (9855, 9856), (9874, 9880), (9881, 9882), (9883, 9885),
    (9888, 9890), (9898, 9900), (9904, 9906), (9917, 9919), (9924, 9926),
    (9928, 9929), (9935, 9936), (9937, 9938), (9939, 9941), (9961, 9963),
    (9968, 9974), (9975, 9979), (9981, 9982), (9986, 9987), (9989, 9990),
    (9992, 9998), (9999, 10000), (10002, 10003), (10004, 10005),
    (10006, 10007), (10013, 10014), (10017, 10018), (10024, 10025),
    (10035, 10037), (10052, 10053), (10055, 10056), (10060, 10061),
    (10062, 10063), (10067, 10070), (10071, 10072), (10083, 10085),
    (10133, 10136), (10145, 10146), (10160, 10161), (10175, 10176),
    (10548, 10550), (11013, 11016), (11035, 11037), (11088, 11089),
    (11093, 11094), (12336, 12337), (12349, 12350), (12951, 12952),
    (12953, 12954), (65039, 65040), (126980, 126981), (127183, 127184),
    (127344, 127346), (127358, 127360), (127374, 127375), (127377, 127387),
    (127462, 127488), (127489, 127491), (127514, 127515), (127535, 127536),
    (127538, 127547), (127568, 127570), (127744, 127778), (127780, 127892),
    (127894, 127896), (127897, 127900), (127902, 127985), (127987, 127990),
    (127991, 128254), (128255, 128318), (128329, 128335), (128336, 128360),
    (128367, 128369), (128371, 128379), (128391, 128392), (128394, 128398),
    (128400, 128401), (128405, 128407), (128420, 128422), (128424, 128425),
    (128433, 128435), (128444, 128445), (128450, 128453), (128465, 128468),
    (128476, 128479), (128481, 128482), (128483, 128484), (128488, 128489),
    (128495, 128496), (128499, 128500), (128506, 128592), (128640, 128710),
    (128715, 128723), (128736, 128742), (128745, 128746), (128747, 128749),
    (128752, 128753), (128755, 128761), (129296, 129339), (129340, 129343),
    (129344, 129350), (129351, 129357), (129360, 129388), (129408, 129432),
    (129472, 129473), (129488, 129511),
]

EMOJIS = [chr(int('{:08x}'.format(x), 16))
          for x in itertools.chain(
                  *(range(begin, end) for begin, end in EMOJI_RANGES)
          )]


class EmojiTransTable(object):
    def __init__(self):
        self._emojis = collections.deque(random.sample(EMOJIS, k=len(EMOJIS)))
        self.correspondence = collections.defaultdict(self.next_emoji)
        # 'X' as an atom type is considered a joker, it has to stay
        self.correspondence['X'] = 'X'

    def __getitem__(self, key):
        return self.correspondence[key]

    def __setitem__(self, key, value):
        self.correspondence[key] = value

    def next_emoji(self):
        return self._emojis.pop()

    def load(self, fname):
        with open(fname) as infile:
            for linenum, line in enumerate(infile):
                if not line:
                    continue
                splitted = line.split('=')
                if len(splitted) == 2:
                    name, emoji = splitted
                else:
                    raise IOError('Invalid line {} in "{}".'
                                  .format(linenum, fname))
                name = name.strip()
                emoji = emoji.strip()
                self.correspondence[name] = emoji

    def write(self, fname):
        with open(fname, 'w') as outfile:
            for name, emoji in self.correspondence.items():
                print('{} = {}'.format(name, emoji), file=outfile)


def _change_n_names(nfields):
    def wraped(self, line):
        names = line.split()[:nfields]
        for name in names:
            emoji = self.correspondence[name]
            line = line.replace(name, emoji, 1)
        return line
    return wraped


class TopReplacer(object):
    def __init__(self, lines, correspondence=None, include=True):
        if include:
            self._lines = self._recursive_top_lines(lines)
        else:
            self._lines = lines
        if correspondence is None:
            self.correspondence = EmojiTransTable()
        else:
            self.correspondence = correspondence
        self._transformers = {
            'atomtypes': self._one_name,
            'nonbond_params': self._two_names,
            'moleculetype': self._one_name,
            'atoms': self._atoms,
            'molecules': self._one_name,
            'pairtypes': self._two_names,
            'bondtypes': self._two_names,
            'constrainttypes': self._two_names,
            'angletypes': self._three_names,
            'dihedraltypes': self._four_names,
            'implicit_genborn_params': self._one_name,
            'cmaptypes': self._five_names,
        }

    def __iter__(self):
        context = None
        continuation = False
        # Maybe it is safer to have a disclaimer. We never know.
        yield '; This file has been emojified with simoji:\n'
        yield '; <https://github.com/jbarnoud/simoji>\n'
        yield '; You probably do not want to use it...\n'
        yield '\n'
        for line in self._lines:
            uncommented = self._uncomment(line).strip()
            section = self._section_name_if_any(uncommented)
            if continuation or not uncommented or uncommented[0] == '#':
                yield line
            elif section is not None:
                context = section
                yield line
            else:
                yield self._transformers.get(context, self._neutral)(line)

            if len(line) >=2 and line[-2] == '\\':
                continuation = True
            else:
                continuation = False

    _one_name = _change_n_names(1)
    _two_names = _change_n_names(2)
    _three_names = _change_n_names(3)
    _four_names = _change_n_names(4)
    _five_names = _change_n_names(5)
        
    def _neutral(self, line):
        return line

    def _atoms(self, line):
        uncommented = self._uncomment(line).strip()
        _, atomtype, _, resname, atomname, *_ = uncommented.split()
        emoji_atomtype = self.correspondence[atomtype]
        emoji_resname = self.correspondence[resname]
        emoji_atomname = self.correspondence[atomname]
        new_line = line.replace(atomtype, emoji_atomtype, 1)
        new_line = new_line.replace(resname, emoji_resname, 1)
        new_line = new_line.replace(atomname, emoji_atomname, 1)
        return new_line

    @staticmethod
    def _include_fname_if_any(line, parent):
        if line.startswith('#include'):
            parent_dir = os.path.dirname(parent)
            _, *fname = line.split()
            fname = ' '.join(fname)
            if not fname:
                raise IOError('Include with no file name')
            if not (fname[0] == fname[-1] == '"'):
                raise IOError('Missformated include')
            fname = os.path.join(parent_dir, fname[1:-1])
            if os.path.exists(fname):
                return fname
            if 'GMXDATA' in os.environ:
                fname = os.path.join(os.environ['GMXDATA'], 'top', fname)
                if os.path.exists(fname):
                    return fname
        return None

    @staticmethod
    def _section_name_if_any(line):
        if not line or line[0] != '[' or line[-1] != ']':
            return None
        return line[1:-1].strip().lower()

    @staticmethod
    def _uncomment(line, comment_char=';'):
        pos = line.find(comment_char)
        if pos > -1:
            line = line[:pos]
        return line

    def _recursive_top_lines(self, infile):
        name = infile.name
        for line in infile:
            uncommented = self._uncomment(line).strip()
            fname = self._include_fname_if_any(uncommented, name)
            if fname is not None:
                with open(fname) as infile:
                    yield from self._recursive_top_lines(infile)
            else:
                yield line


def write(fname, lines):
    with open(fname, 'w') as outfile:
        for line in lines:
            print(line, end='', file=outfile)


def gro_replacer(lines, correspondence=None):
    if correspondence is None:
        correspondence = EmojiTransTable()
    yield next(lines)  # Title
    num_str = next(lines)  # Number of atoms
    yield num_str
    # We want to stop iterating before the last line so we can get the box
    # untouched latter. To do that we cap the iteration with a range.
    for _, line in zip(range(int(num_str)), lines):
        resname = line[5:10].strip()
        atom_name = line[10:15].strip()
        res_emoji = correspondence[resname]
        atom_emoji = correspondence[atom_name]
        # Names should be 5 bytes to fit the GRO format, so we need to know
        # how many bytes the emojis are and fill the gaps with spaces.
        res_emoji = res_emoji + ' ' * (5 - len(res_emoji.encode('utf-8')))
        atom_emoji = atom_emoji + ' ' * (5 - len(atom_emoji.encode('utf-8')))
        yield line[:5] + res_emoji + atom_emoji + line[15:]
    yield next(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--version', action='version', version=__version__)
    parser.add_argument('-f', dest='structure', help='input structure (.gro)')
    parser.add_argument('-p', dest='topology', help='input topology (.top/itp)')
    parser.add_argument('-t', dest='transtable', help='input correspondance table')
    parser.add_argument('-o', dest='out_structure', default='emoji.gro',
                        help='output structure (.gro)')
    parser.add_argument('-op', dest='out_topology', default='emoji.top',
                        help='output topology (.top/itp)')
    parser.add_argument('-ot', dest='out_transtable',
                        help='output correspondence table')
    args = parser.parse_args()

    if (args.structure, args.topology) == (None, None):
        parser.error('Hey! I do not have anything to do.')

    correspondence = EmojiTransTable()
    if args.transtable is not None:
        correspondence.load(args.transtable)

    if args.topology is not None:
        with open(args.topology) as infile:
            lines = TopReplacer(infile, correspondence=correspondence)
            write(args.out_topology, lines)

    if args.structure is not None:
        with open(args.structure) as infile:
            lines = gro_replacer(infile, correspondence=correspondence)
            write(args.out_structure, lines)

    if args.out_transtable is not None:
        correspondence.write(args.out_transtable)


if __name__ == '__main__':
    main()
