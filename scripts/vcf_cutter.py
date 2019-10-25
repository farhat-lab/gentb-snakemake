#!/bin/env python3
"""
A simple script to strip out non-useful references from vcf files to 'clean' them.
"""

import os
from argparse import ArgumentParser, ArgumentTypeError

def vcf_cutter(inp, output):
    """Strip out all non-mutation references from the vcf"""
    column_ref = -1
    column_alt = -1

    for line in inp:
        if line.startswith(b'#'):
            output.write(line)

            if line.startswith(b'##'):
                continue

            # extract which column REF and ALT are in
            headers = line.rstrip(b'\n').split(b'\t')
            column_ref = headers.index(b'REF')
            column_alt = headers.index(b'ALT')

            continue

        data = line.rstrip(b'\n').split(b'\t')

        # if column_ref and column_alt aren't set then this will fail
        if len(data[column_ref]) != 1 or data[column_alt] != '.':
            output.write(line)

def file_type(fname):
    """Confirm file input is in a file"""
    if os.path.isfile(fname):
        return fname
    raise ArgumentTypeError("File not found: %s" % fname)

def main():
    """Main body of the script"""
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--input', '-i', help='Filename of the input vcf file', type=file_type)
    parser.add_argument('--output', '-o', help='File of the output vcf file')
    args = parser.parse_args()
    with open(args.input, 'rb') as fhl:
        with open(args.output, 'wb') as fhl:
            vcf_cutter(args.input, fhl)

if __name__ == '__main__':
    main()
