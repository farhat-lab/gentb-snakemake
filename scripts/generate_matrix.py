#!/usr/bin/env python
"""
generates a matrix.csv file from one .var file
"""

import re
import os
import sys
import csv

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
    #strains = set()
    output2 = sys.stdout
    output2.write('strain,')
    output2.write(",".join(variants))
    output2.write("\n")
    output2.write(filename.rsplit(".", 1)[0]) # strain name
    output2.write(",")

    with open(filename, 'r') as fhl:
        input_csv = csv.reader(fhl, delimiter="\t")
        header = next(input_csv)
        varname_index = header.index('varname')

        out = set([])
        for record in input_csv:
            out.add(process_record(record[varname_index], "\t".join(variants)))

    for variant in variants:
        output2.write(str(int(variant in out))+',')
    output2.write("\n")

def process_record(snpname, variants):
    parts = snpname.split('_')
    if snpname in variants:
        return snpname
    elif parts[1] in ['CN', 'CS', 'CZ']:
        # coding snps pool changes that cause the same aachange
        if re.search(r'\*$', parts[4]): #for to stop changes
            parts[4] = parts[4][0:(len(parts[4])-1)] + "."
        if re.search(r'^\*', parts[4]): #for stop to changes
            parts[4] = "." + parts[4][1:]
        if re.search(r"'$", parts[5]): #for oxyR'
            parts[5] = parts[5][0:(len(parts[5])-1)] + "."
        pattern = parts[4] + "_" + parts[5]
        pattern = r"\s?([\w\.]+%s[\w\.]*)\s?" % pattern
        ita = tuple(re.finditer(pattern, variants, re.IGNORECASE))
        if len(ita) > 1:
            return ita[0].group()
            #plan to add more matching to the right mutation here in a future version
        elif len(ita) == 1:
            return ita[0].group()
    elif parts[1] in ['CF']:
        #coding frameshifts pool all that occur at the same nucleotide start
        pattern = parts[1] + '_' + parts[2] + r'_[^\s\,]+_' + parts[5]
        pattern = r"\s?([\w\.]+%s[\w\.]*)\s?" % pattern
        itb = tuple(re.finditer(pattern, variants, re.IGNORECASE))
        if len(itb) > 1:
            return itb[0].group()
            # plan to add more matching to the right mutation here in a future version
        elif len(itb) == 1:
            return itb[0].group()
    elif parts[1] == 'P':
        # promoter (to maintain compatibility with old naming used in
        # randomforest built from MIP data
        operon = parts[len(parts)-1].split('-')
        if operon[0] == "promoter":
            pattern = parts[3] + '_' + operon[0] + '_' + operon[1]
        else:
            pattern = parts[3] + '_' + parts[4] + '_' + operon[0]

        pattern = r"\s?([\w\.]+%s[\w\.]*)\s?" % pattern
        itc = tuple(re.finditer(pattern, variants, re.IGNORECASE))
        if len(itc) > 1:
            return itc[0].group()
            #plan to add more matching to the right mutation here in a future version
        elif len(itc) == 1:
            return itc[0].group()


if __name__ == '__main__':
    main()
