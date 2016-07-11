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
    parser.add_option("-l", "--locations_screen", action="store", type="string", dest="ok_locs_file", default=None,
                      help="file (with header) listing insertion locations (replicon, position, direction) to " +
                      "include in compiled list (and ignore all other locations); without this option, " +
                      "all locations are counted")
    parser.add_option("-o", "--output_file", action="store", type="string", dest="outfile", required=True,
                      help="output file")
    opts, args = parser.parse_args()

    if len(args)==0:
        parser.print_help()
        exit(1)

    return opts, args

def get_ok_locs(ok_locs_file):
    ok_locs = list()
    with open(ok_locs_file, "r") as fh:
        try:
            header = fh.readline()
            for line in fh:
                (replicon, pos, strand) = line.rstrip().split("\t")[:3]
                ok_locs.append((replicon, pos, strand))
        except StopIteration:
            break
    return ok_locs

def read_files(infiles, ok_locs_file=None):
    if ok_locs_file:
        ok_locations = get_ok_locs(ok_locs_file)
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
                        if (replicon, pos, strand) not in ok_locations:
                            # skip the location if it's not a one designated to be counted
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
    (totals, filereads) = read_files(infiles, args.ok_locs_file)

    write_compiled(totals, filereads, infiles, opts.outfile)

if __name__ == "__main__":
    main()
