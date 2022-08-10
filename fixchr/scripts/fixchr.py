#!/usr/bin/env python3

# -*- coding: utf-8 -*-
"""
Created on Wed May 10 13:05:51 2017

@author: goel
"""
import argparse
from fixchr import __version__


def fixchr(args):
    # Check that correct version of python is being used
    import logging
    import os
    import sys
    from fixchr.func import readcoords, checkdir, revcomp, setlogconfig
    from fixchr.dotplot import drawdotplot

    ## Define logger
    setlogconfig(args.log)
    logger = logging.getLogger("fixchr")

    # Set CWD and check if it exists
    if args.dir is None:
        args.dir = os.getcwd() + os.sep
    else:
        if os.path.isdir(args.dir):
            args.dir = args.dir + os.sep
        else:
            logger.error(args.dir + ' is not a valid folder. Exiting.')
            sys.exit()

    # Check prefix
    if os.sep in args.prefix:
        logger.warning('For specifying output folder use --dir, use --prefix for modifying the output file names. Current --prefix ({}) may result in crashes.'.format(args.prefix))


    # Set CIGAR FLAG
    if args.ftype in ['S', 'B', 'P']:
        args.cigar = True


    ###################################################################
    # Read alignments and compare lengths with genome fasta
    ###################################################################


    # from syri.scripts.func import readfasta
    # import numpy as np




    alpaf = 'GCA_900660825.1.paf'
    albam = 'GCA_900660825.1.sorted.bam'
    ref = 'TAIR10_Filtered.fasta.gz'
    qry = 'GCA_900660825.1_Ath.Ler-0.MPIPZ.v1.0_genomic.fna.gz'

    # coords = readcoords(alpaf, ftype='P', filter=args.f, cigar=args.cigar)
    # Read coords
    coords = readcoords(albam, ftype='B', f=True, cigar=True)
    coords.to_csv("input_alignments.txt", index=False, header=True, sep='\t')
    # Get dotplot for initial alignments
    refg = readfasta(ref)
    qryg = readfasta(qry)
    rchrs_len = {k: len(v) for k, v in refg.items()}
    qchrs_len = {k: len(v) for k, v in qryg.items()}
    # TODO: Add parsers for AGP files and plot dimensions
    drawdotplot(coords, rchrs_len, qchrs_len, out='input.pdf')
    # drawdotplot(args.sam.name, out=args.o.name, ragp=args.ragp, qagp=args.qagp,
    #             minsize=args.m, width=args.W, height=args.H)

    coords2, assigned = homchr(coords, rchrs_len, qchrs_len, csize=100000)
    coords2.sort_values(['aChr', 'aStart', 'aEnd'], inplace=True)
    coords2.to_csv("homologous_chromosome_alignments.txt", index=False, header=True, sep='\t')
    drawdotplot(coords2, rchrs_len, qchrs_len, out='homologous_selected.pdf')

    rv = checkdir(coords2, assigned)
    if len(rv) > 0:
        logger.warning(f"Inverting query chromosomes: {list(rv)}")

    refout = {k: refg[k] for k in assigned.keys()}
    qryout = {k: qryg[k] for k in assigned.values()}
    for k in rv:
        qryout[k] = revcomp(qryout[k])
        coords2.loc[coords.bChr == k, 'bStart'] = qchrs_len[k] - coords2.loc[coords.bChr == k, 'bStart'] + 1
        coords2.loc[coords.bChr == k, 'bEnd'] = qchrs_len[k] - coords2.loc[coords.bChr == k, 'bEnd'] + 1
        coords2.loc[coords.bChr == k, ['bStart', 'bEnd']] = coords2.loc[coords.bChr == k, ['bEnd', 'bStart']].to_numpy()
    coords2.to_csv("homologous_strand_corrected_alignments.txt", index=False, header=True, sep='\t')
    drawdotplot(coords2, rchrs_len, qchrs_len, out='homologous_strand_corrected.pdf')
    writefasta(refout, "ref.filtered.fa")
    writefasta(qryout, "qry.filtered.fa")
    logger.info("Finished")


