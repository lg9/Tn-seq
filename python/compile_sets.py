#!/usr/bin/env python
# Purpose: Combine one or more sets of read counts into a single file.
# A set consists of two files, one with all reads and one with q=0 reads
#
# Copyright (c) 2014 University of Washington

#-----------------------------------------------------
# modules
#-----------------------------------------------------
import sys
import optparse
import os.path
import re

#------------------------------------------------------------------------------
# constants
#------------------------------------------------------------------------------
SHOW_TOTALS = False

#------------------------------------------------------------------------------
# functions
#------------------------------------------------------------------------------
def init_options():
    usage = "usage: %prog [options] <input file(s)>"
    parser = optparse.OptionParser(prog=sys.argv[0],
                                   usage=usage,
                                   add_help_option=True)
    parser.add_option("-l", "--ok_locs_file", action="store", type="string", dest="ok_locs_file", default=None,
                      help="file (with header) listing insertion locations (replicon, position, direction) to count " +
                      "(and ignore all other locations); without this option, all locations are counted")
    parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", help="output file")
    opts, args = parser.parse_args()

    if len(args)==0 or opts.outfile is None:
        parser.print_help()
        exit(1)

    return opts, args

def get_ok_locs(ok_locs_file):
    ok_locs = set()
    with open(ok_locs_file, "r") as fh:
        header = fh.readline()
        for line in fh:
            (replicon, pos, strand) = line.rstrip().split("\t")[:3]
            ok_locs.add((replicon, pos, strand))
    return ok_locs

def read_files(infiles, ok_locs_file=None):
    if ok_locs_file:
        ok_locations = get_ok_locs(ok_locs_file)
        print "Acceptable locations loaded: " + str(len(ok_locations))
    totals = dict()
    filereads = dict()
    for infile in infiles:
        print "reading file " + infile
        with open(infile, "r") as fh:
            try:
                header = fh.readline()
                for line in fh:
                    (replicon, pos, strand, readcount) = line.rstrip().split("\t")
                    if ok_locs_file:
                        # (adjust for possible back-end sequencing indication (lower case strand designation))
                        fe_strand = strand    # front-end strand should be upper case
                        if strand == "f":
                            fe_strand = "R"
                        elif strand == "r":
                            fe_strand = "F"
                        # skip the location if it's not one designated to be counted:
                        if (replicon, pos, fe_strand) not in ok_locations: 
                            continue
                    totals[(replicon, pos, strand)] = totals.get((replicon, pos, strand), 0) + float(readcount)
                    filereads[(infile, replicon, pos, strand)] = readcount
            except StopIteration:
                break
    return (totals, filereads)

def write_compiled(totals, filereads, infiles, outfile):
    with open(outfile, "w") as fh:
        if SHOW_TOTALS:
            total_header = "\tAllReads"
        else:
            total_header = ""
        file_list = "\t".join(map(os.path.basename, infiles))
        fh.write("Replicon\tPosition\tDirection\t" + file_list + total_header + "\n")
        for (replicon, pos, strand) in sorted(totals, key=totals.get, reverse=True):
            outline = replicon + "\t" + pos + "\t" + strand
            for infile in infiles:
                numreads = filereads.get((infile, replicon, pos, strand), 0)
                outline += "\t" + str(numreads)
            if SHOW_TOTALS:
                total_val = "\t" + str(totals[(pos, strand)])
            else:
                total_val = ""
            fh.write(outline + total_val + "\n")
    print "total positions tabulated: " + str(len(totals))

#-----------------------------------------------------
# main
#-----------------------------------------------------
def main():
    # Get command line options
    opts, args = init_options()
    infiles = args

    print "Compiling sets of read counts"
    (totals, filereads) = read_files(infiles, opts.ok_locs_file)

    write_compiled(totals, filereads, infiles, opts.outfile)

if __name__ == "__main__":
    main()
