#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Generate table for efm data graphing from .gd files
"""

import argparse
import os
import re
import sys
from collections import defaultdict

parser = argparse.ArgumentParser(description="generate EFM data table from annotated .gd files")
parser.add_argument("-g", "--gd", help="folder containing annotated .gd files")
parser.add_argument("-o", "--output", default='summary_output.tsv', help="output filename")
args = parser.parse_args()

master_dict = {}

# Reference file used for breseq run uses biobrick igem part numbers, need dictionary for conversion. Dennis maintains the biobrick numbers are interesting in and of themselves to people
name_replacements = {
    "BBa_K864100": "SYFP2_CDS-BBa_K864100",  # http://parts.igem.org/Part:BBa_K864100
    "BBa_J23110": "SYFP2_Promoter-BBa_J23110",  # http://parts.igem.org/Part:BBa_J23110
    "BBa_K608006": "Cam_resistance_gene_Promoter&RBS-BBa_K608006",  # http://parts.igem.org/Part:BBa_K608006
    "BBa_B0032": "Cam_RBS-BBa_B0032"  # http://parts.igem.org/Part:BBa_B0032
}

def name_replace(string1, replacement_dictionary=name_replacements):
    # Verify search stings not used in any other solutions. This is not perfect, still possible that replacement adds new search string at front or beginning but unlikely
    for k in replacement_dictionary:
        assert not re.search(k, str(replacement_dictionary.values()).replace(replacement_dictionary[k], ""))

    for k in replacement_dictionary:
        string1 = string1.replace(k, replacement_dictionary[k])
    return string1



# read gd files in
for gd_file in os.listdir(args.gd):
    if gd_file.endswith(".gd"):

        # All plasmid sequencing files have the naming convention Day_strain_replicate.gd
        filename_parts = gd_file.split("_")
        assert len(filename_parts) == 3, "All plasmid sequencing files have the naming convention Day_strain_replicate.gd\n%s does not conform to this naming convention" % gd_file

        master_dict[gd_file] = {
            "Time": filename_parts[0],
            "Strain": filename_parts[1],
            "Replicate": filename_parts[2].replace(".gd", ""),
            "Mutations": defaultdict(lambda: defaultdict())  # default dictionary of default dictionary, 2 keys deep on the fly assignment
        }

        with open(args.gd + "/" + gd_file, "r") as f:
            for line in f:
                line=line.rstrip().split("\t")

                if re.match("^[A-Z][A-Z]$", line[0]) or re.match("#", line[0]):  # ignore evidence lines and commented lines
                    continue

                assert line[1] not in master_dict[gd_file]['Mutations'], "ID number duplicated:\n%s\n%s" % (line, master_dict[gd_file]["Mutations"])

                master_dict[gd_file]['Mutations'][line[1]]["Type"] = line[0]
                master_dict[gd_file]['Mutations'][line[1]]["Reference"] = line[3]
                master_dict[gd_file]['Mutations'][line[1]]["Position"] = line[4]

                # Mutation Size
                if line[0] == "DEL":
                    master_dict[gd_file]['Mutations'][line[1]]["Size"] = line[5]
                elif line[0] == "INS":
                    master_dict[gd_file]['Mutations'][line[1]]["Size"] = str(len(line[5]))
                else:
                    master_dict[gd_file]['Mutations'][line[1]]["Size"] = "0"  #includes SNP and MOB

                # Mobile element name
                if line[0] == "MOB":
                    master_dict[gd_file]['Mutations'][line[1]]["Mobile Element Name"] = line[5]
                else:
                    master_dict[gd_file]['Mutations'][line[1]]["Mobile Element Name"] = "NA"


                freq = "Unknown"
                category = "Unknown"
                gene_position = "Unknown"
                gene_name = "Unknown"
                repeat = ["NA", "NA", "NA"]

                for column in line:
                    if re.match("mutation_category", column):
                        assert category == "Unknown"
                        assert len(column.split("=")) == 2
                        category = column.split("=")[1]
                        continue
                    if re.match("frequency", column):
                        assert freq == "Unknown"
                        assert len(column.split("=")) == 2
                        freq = column.split("=")[1]
                        continue
                    if re.match("gene_position", column):
                        assert gene_position == "Unknown"
                        assert len(column.split("=")) == 2
                        gene_position = column.split("=")[1]
                        continue
                    if re.match("gene_name", column):
                        assert gene_name == "Unknown"
                        assert len(column.split("=")) == 2
                        gene_name = column.split("=")[1]
                        continue
                    if re.match("repeat_seq", column):
                        assert repeat[0] == "NA"
                        assert len(column.split("=")) == 2
                        repeat[0] = column.split("=")[1]
                        continue
                    if re.match("repeat_new_copies", column):
                        assert repeat[2] == "NA"
                        assert len(column.split("=")) == 2
                        repeat[2] = column.split("=")[1]
                        continue
                    if re.match("repeat_ref_copies", column):
                        assert repeat[1] == "NA"
                        assert len(column.split("=")) == 2
                        repeat[1] = column.split("=")[1]
                        continue


                master_dict[gd_file]['Mutations'][line[1]]["Frequency"] = freq
                master_dict[gd_file]['Mutations'][line[1]]["Category"] = category

                if repeat == ["NA", "NA", "NA"]:
                    repeat = "NA"  # repeat not found
                else:
                    repeat = "_".join(repeat)  #join to form: sequence_original-copies_new-copies

                master_dict[gd_file]['Mutations'][line[1]]["Repeat(seq_orig-copies_new-copies)"] = repeat

                if gene_position == "Unknown" or gene_name == "Unknown":
                    if category not in ["large_deletion", "small_indel", "mobile_element_insertion"]:
                        assert re.search("intergenic", "\t".join(line)), line

                master_dict[gd_file]['Mutations'][line[1]]["Gene_Name"] = name_replace(gene_name)
                master_dict[gd_file]['Mutations'][line[1]]["Gene_position"] = gene_position


# write summary file out
with open(args.gd + "/" + args.output, "w") as output:
    print_order = ["Type", "Reference", "Position", "Size", "Frequency", "Category", "Gene_Name", "Gene_position", "Repeat(seq_orig-copies_new-copies)", "Mobile Element Name"]
    print>>output, "\t".join(["Time", "Strain", "Replicate"] + print_order )

    for gd_file in master_dict:
        time = master_dict[gd_file]["Time"]
        strain = master_dict[gd_file]["Strain"]
        replicate = master_dict[gd_file]["Replicate"]
        for mutation in master_dict[gd_file]['Mutations']:
            # if master_dict[gd_file]['Mutations'][mutation]["Size"] == "0":
            #     continue
            print>>output, "\t".join([time, strain, replicate] + [master_dict[gd_file]['Mutations'][mutation][x] for x in print_order])
