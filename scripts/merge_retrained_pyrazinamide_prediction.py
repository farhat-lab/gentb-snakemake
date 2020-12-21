#!/usr/bin/env python

import argparse
import json

#define arguments
parser = argparse.ArgumentParser(
    description="This script processes final prediction probabilities in json files from the retrained pyrazinamide model"
)

#argument 1: final prediction file
parser.add_argument("final_json_file", help="final prediction json file")

#argument 2: prediction file from the retrained pyrazinamide model
parser.add_argument("pza_json_file", help="son file from the retrained PZA random Forest model")

#specify the output name
parser.add_argument("output", help="Output filename for the merged json")


args = parser.parse_args()


#open both json files
with open(args.final_json_file, "r") as f:
    final_json_data = json.load(f)


with open(args.pza_json_file, "r") as f:
    pza_json_data = json.load(f)


#define the probability value from the pza JSON output
pza_probability = pza_json_data[0][0][2]

#write PZA probability to final JSON file and remove values for FP and FN rate (not computed with the new PZA random Forest)
final_json_data[0][2][2] = pza_probability
final_json_data[0][2][3] = ''
final_json_data[0][2][4] = ''


#write merged output into a final merged json file in the same format
with open(args.output, 'w') as json_file:
    json_file.write(json.dumps(final_json_data, indent = 2))
