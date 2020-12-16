#!/usr/bin/env python
"""
generates a matrix.csv file from one .var file for the retrained pyrazinamide RandomForest model
"""

import re
import os
import sys

from argparse import ArgumentParser, ArgumentTypeError

def file_type(fname):
    """Confirm file input is in a file"""
    if os.path.isfile(fname):
        return fname
    raise ArgumentTypeError("File not found: %s" % fname)

def main():
    """Main body of the script"""
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('mutations', help='CSV File of mutation names', type=file_type)
    parser.add_argument('strain', help='VAR File for a strain to process', type=file_type)
    args = parser.parse_args()

    with open(args.mutations, 'r') as fhl:
        variants = list(set(fhl.read().rstrip().split(",")))
    generate_matrix(variants, args.strain)

def generate_matrix(variants, filename):
    """Generate the matrix"""
    out = {}
    output2 = sys.stdout
    output2.write('strain,')
    output2.write(",".join(variants))
    output2.write("\n")
    output2.write(filename.rsplit(".", 1)[0]) # strain name

    for line in open(filename):
        items = line.rstrip().split("\t") #variant name i.e. snpname
        snpname = items[5]
        if snpname != "varname":
            out[snpname] = "1"
 
    for variant in variants:
        output2.write("," + str(int(variant in out)))
    output2.write("\n")

if __name__ == '__main__':
    main()
