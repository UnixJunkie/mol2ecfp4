#!/usr/bin/python

# get the ECFP4 bits set for each molecule in a file; using rdkit's
# python bindings

from __future__ import print_function

import sys

from rdkit import Chem
from rdkit.Chem import AllChem

# ECFP4 is used in ZINC15
# In ECFP4, 4 stands for the diameter of the atom environment
fp_diameter = 4
# but rdkit wants a radius
fp_radius = fp_diameter / 2

def RobustSmilesMolSupplier(filename):
    with open(filename) as f:
        for line in f:
            split = line.split()
            smile = split[0]
            mol_name = split[1]
            mol = Chem.MolFromSmiles(smile)
            yield mol, mol_name

def get_mol_reader(filename):
    if filename.endswith(".smi"):
        return RobustSmilesMolSupplier(filename)
    else:
        print("get_mol_reader: fatal: unsupported file format: %s" % filename,
              file=sys.stderr)
        exit(1)

def get_name(mol):
    return mol.GetProp('_Name')
    
argc = len(sys.argv)
if argc != 2:
    #                           0          1
    print("ecfp6: fatal: usage: ./ecfp6.py filename.{smi|sdf}",
          file=sys.stderr)
    exit(1)

mol_reader = get_mol_reader(sys.argv[1])
i = 0
print("#mol_name,IC50 in mol/L (0.0 means unknown),ECFP4 set bits indexes");
for mol, mol_name in mol_reader:
    try:
        fp = AllChem.GetMorganFingerprint(mol, fp_radius)
        set_bits = fp.GetNonzeroElements()
        # set_bits: a dictionary; values are all ones
        # keys are the bit index in the fp
        print("%s,0.0,[" % mol_name, end='')
        first_time = True
        for key in set_bits:
            if first_time:
                print("%d" % key, end='')
                first_time = False
            else:
                print(";%d" % key, end='')
        print("]\n", end='')
    except:
        print("ecfp6: error: molecule at index %d" % i,
              file=sys.stderr)
    i = i + 1
