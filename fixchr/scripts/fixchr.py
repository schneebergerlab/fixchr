#!/usr/bin/env python3

# -*- coding: utf-8 -*-
"""
Created on Wed May 10 13:05:51 2017

@author: goel
"""
import argparse
from pathlib import Path

from fixchr import __version__


def fixchr(args):
    # Check that correct version of python is being used
    # import logging
    import os

    from fixchr.scripts.dotplot import drawdotplot
    from fixchr.scripts.func import (
        checkdir,
        homchr,
        mylogger,
        readcoords,
        readfasta,
        revcomp,
        setlogconfig,
        writefasta,
    )

    ## Define logger
    setlogconfig(args.log)
    # logger = logging.getLogger("fixchr")
    logger = mylogger("fixchr")

    # Set CWD and check if it exists
    args.dir = Path(args.dir)
    if not args.dir.exists():
        args.dir.mkdir(parents=True, exist_ok=True)

    # Check prefix
    if os.sep in args.prefix:
        logger.warning(
            "For specifying output folder use --dir, use --prefix for modifying the output file names. Current --prefix ({}) may result in crashes.".format(
                args.prefix
            )
        )

    prefix = str(args.dir / args.prefix)

    ###################################################################
    # Read alignments and compare lengths with genome fasta
    ###################################################################
    cfin = args.infile.name
    ref = args.ref.name
    qry = args.qry.name
    ftype = args.ftype
    csize = args.csize
    # Set CIGAR FLAG
    if ftype in ["S", "B", "P"]:
        cigar = True
    elif ftype == "T":
        cigar = False

    # coords = readcoords(alpaf, ftype='P', filter=args.f, cigar=args.cigar)
    # Read coords
    logger.info("Reading alignments")
    coords = readcoords(cfin, ftype=ftype, f=True, cigar=cigar)
    coords.to_csv("input_alignments.txt", index=False, header=True, sep="\t")
    # Get dotplot for initial alignments
    refg = readfasta(ref)
    qryg = readfasta(qry)
    rchrs_len = {k: len(v) for k, v in refg.items()}
    qchrs_len = {k: len(v) for k, v in qryg.items()}
    # TODO: Add parsers for AGP files and plot dimensions
    drawdotplot(coords, rchrs_len, qchrs_len, out=f"{prefix}fixchr.input.pdf")
    # drawdotplot(args.sam.name, out=args.o.name, ragp=args.ragp, qagp=args.qagp,
    #             minsize=args.m, width=args.W, height=args.H)
    logger.info("Selecting homologous chromosomes")
    coords2, assigned = homchr(coords, rchrs_len, qchrs_len, csize=csize)
    coords2.sort_values(["aChr", "aStart", "aEnd"], inplace=True)
    coords2.to_csv(
        f"{prefix}fixchr.homologous_alignments.txt", index=False, header=True, sep="\t"
    )
    drawdotplot(coords2, rchrs_len, qchrs_len, out="homologous.pdf")
    logger.info("Checking chromosome direction")
    rv = checkdir(coords2, assigned)
    if len(rv) > 0:
        logger.warning(f"Inverting query chromosomes: {list(rv)}")
    refout = {k: refg[k] for k in assigned.keys()}
    qryout = {k: qryg[k] for k in assigned.values()}
    for k in rv:
        qryout[k] = revcomp(qryout[k])
        coords2.loc[coords.bChr == k, "bStart"] = (
            qchrs_len[k] - coords2.loc[coords.bChr == k, "bStart"] + 1
        )
        coords2.loc[coords.bChr == k, "bEnd"] = (
            qchrs_len[k] - coords2.loc[coords.bChr == k, "bEnd"] + 1
        )
        coords2.loc[coords.bChr == k, ["bStart", "bEnd"]] = coords2.loc[
            coords.bChr == k, ["bEnd", "bStart"]
        ].to_numpy()
    logger.info("Writing output files")
    coords2.to_csv(
        f"{prefix}fixchr.homologous_strand_corrected_alignments.txt",
        index=False,
        header=True,
        sep="\t",
    )
    drawdotplot(
        coords2,
        rchrs_len,
        qchrs_len,
        out=f"{prefix}fixchr.homologous_strand_corrected.pdf",
    )
    writefasta(refout, f"{prefix}fixchr.ref.filtered.fa")
    writefasta(qryout, f"{prefix}fixchr.qry.filtered.fa")
    logger.info("Finished")


# END


def main(cmd):
    parser = argparse.ArgumentParser(
        "Filter and reorient genomes to get homologous chromosomes",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    optional = parser._action_groups.pop()
    required = parser.add_argument_group("Input Files")
    required.add_argument(
        "-c",
        dest="infile",
        help="File containing alignment coordinates",
        type=argparse.FileType("r"),
        required=True,
    )
    required.add_argument(
        "-r",
        dest="ref",
        help="Genome A (which is considered as reference for the alignments). Required for local variation (large indels, CNVs) identification.",
        type=argparse.FileType("r"),
    )
    required.add_argument(
        "-q",
        dest="qry",
        help="Genome B (which is considered as query for the alignments). Required for local variation (large indels, CNVs) identification.",
        type=argparse.FileType("r"),
    )

    other = parser.add_argument_group("Additional arguments")
    other.add_argument(
        "-F",
        dest="ftype",
        help="Input file type. T: Table, S: SAM, B: BAM, P: PAF",
        default="T",
        choices=["T", "S", "B", "P"],
    )
    other.add_argument(
        "-f",
        dest="f",
        help="As a default, low quality and small alignments are filtered out. Use this parameter to use the full list of alignments without any filtering.",
        default=True,
        action="store_false",
    )
    other.add_argument(
        "--contig_size",
        dest="csize",
        help="Minimum contig size. Contigs smaller than this value would be filtered out.",
        default=100000,
        type=int,
    )
    other.add_argument(
        "--dir",
        dest="dir",
        help="path to working directory (if not current directory). All files must be in this directory.",
        action="store",
        default=".",
    )
    other.add_argument(
        "--prefix",
        dest="prefix",
        help="Prefix to add before the output file Names",
        type=str,
        default="",
    )

    optional.add_argument(
        "--log",
        dest="log",
        help="log level",
        type=str,
        default="INFO",
        choices=["DEBUG", "INFO", "WARN"],
    )
    optional.add_argument(
        "--version", action="version", version="{version}".format(version=__version__)
    )
    parser._action_groups.append(optional)
    args = parser.parse_args(cmd)
    fixchr(args)


# END
