#!/usr/bin/env python3
import matplotlib
matplotlib.use('agg')

def samgetchrsize(sam):
    import pysam
    samfile = pysam.AlignmentFile(sam, "r")
    rchrs = samfile.references
    rchrs_len = {rchrs[i]:samfile.lengths[i] for i in range(len(rchrs))}
    qchrs = set()
    qchrs_len = {}
    for r in samfile:
        if r.qname not in qchrs:
            qchrs = set(list(qchrs) + [r.qname])
            qchrs_len[r.qname] = r.infer_read_length() if r.infer_read_length() is not None else len(r.query_sequence)
    return([rchrs_len, qchrs_len])


def readagp(f):
    from collections import defaultdict
    outdict = defaultdict(list)
    with open(f, 'r') as fin:
        for line in fin:
            if line[0] == '#': continue
            line = line.strip().split()
            if line[4] not in 'W': continue
            outdict[line[0]].append([int(line[1]), int(line[2])])
    return(outdict)


# def drawdotplot(sam, out='dotplot.pdf', ragp='', qagp='', minsize=1000, width=12, height=12):
def drawdotplot(coords, rchrs_len, qchrs_len, out='dotplot.pdf', ragp=None, qagp=None, minsize=1000, width=12, height=12):
    '''
    :param al: whole-genome alignment SAM file
    :param rchr_size:
    :param qchr_size:
    :param ragp: Path to reference AGP file
    :param qagp: Path to query AGP file
    :return:
    '''
    import logging
    from collections import deque
    from matplotlib import pyplot as plt
    from matplotlib.patches import Rectangle
    import numpy as np

    # print("reading sam")
    # al = samtocoords(sam)
    # print('done reading')
    # rchrs_len, qchrs_len = samgetchrsize(sam)
    logger = logging.getLogger('drawdotplot')
    al = coords.copy()
    ragpdata = {}
    qagpdata = {}
    if ragp is not None: ragpdata = readagp(ragp.name)
    if qagp is not None: qagpdata = readagp(qagp.name)
    #TODO: Ensure that there are no errors when there are no alignments for a sequence
    alachr = set(np.unique(al.aChr))
    albchr = set(np.unique(al.bChr))
    rchrs = sorted([k for k in rchrs_len.keys() if k in alachr], key=lambda x: rchrs_len[x], reverse=True)
    qchrs = sorted([k for k in qchrs_len.keys() if k in albchr], key=lambda x: qchrs_len[x], reverse=True)
    # qchrs = sorted(qchrs_len.keys(), key=lambda x: qchrs_len[x], reverse=True)
    rcumsum = deque([0])
    qcumsum = deque([0])

    for i in range(1, len(rchrs)):
        cumsum = sum([rchrs_len[rchrs[j]] for j in range(i)])
        al.loc[al['aChr'] == rchrs[i], 'aStart'] += cumsum
        al.loc[al['aChr'] == rchrs[i], 'aEnd'] += cumsum
        if ragp is not None:
            for k in range(len(ragpdata[rchrs[i]])):
                ragpdata[rchrs[i]][k][0] += cumsum
                ragpdata[rchrs[i]][k][1] += cumsum
        rcumsum.append(cumsum)

    for i in range(1, len(qchrs)):
        cumsum = sum([qchrs_len[qchrs[j]] for j in range(i)])
        al.loc[al['bChr'] == qchrs[i], 'bStart'] += cumsum
        al.loc[al['bChr'] == qchrs[i], 'bEnd'] += cumsum
        if qagp is not None:
            for k in range(len(qagpdata[qchrs[i]])):
                qagpdata[qchrs[i]][k][0] += cumsum
                qagpdata[qchrs[i]][k][1] += cumsum
        qcumsum.append(cumsum)

    al_data = deque()
    for row in al.itertuples(index=False):
        if row[4] < 1000 or row[5] < minsize: next
        al_data.append([row[0], row[1]])
        al_data.append([row[2], row[3]])
        if row[8] == 1: al_data.append('r')
        if row[8] == -1: al_data.append('b')
    al_data = list(al_data)
    xticks = [rcumsum[j]+(rchrs_len[rchrs[j]])/2 for j in range(len(rchrs))]
    yticks = [qcumsum[j]+(qchrs_len[qchrs[j]])/2 for j in range(len(qchrs))]

    logger.info('starting drawing')

    figure = plt.figure(figsize=[width, height])
    ax = plt.subplot(1, 1, 1)
    ax.margins(x=0, y=0)
    ax.set_xlim([0, sum([rchrs_len[k] for k in rchrs])])
    ax.set_ylim([0, sum([qchrs_len[k] for k in qchrs])])
    for r in rcumsum:
        ax.axvline(r, linestyle='--', color='black', alpha=1, linewidth=0.2)
    for r in qcumsum:
        ax.axhline(r, linestyle='--', color='black', alpha=1, linewidth=0.2)
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.set_xticklabels(rchrs, rotation=90)
    ax.set_yticklabels(qchrs)
    ax.plot(*al_data)

    for v in ragpdata.values():
        sorted_v = sorted(v, key=lambda x: x[0])
        draw = -1
        for i in range(len(sorted_v)):
            draw *= -1
            # print(draw)
            if draw == 1:
                ax.add_patch(Rectangle((sorted_v[i][0], 0), sorted_v[i][1] - sorted_v[i][0], sum(list(qchrs_len.values())),
                                       edgecolor='lightgrey',
                                       facecolor='lightgrey',
                                       fill=True,
                                       alpha=0.25,
                                       lw=0))
    for v in qagpdata.values():
        sorted_v = sorted(v, key=lambda x: x[0])
        draw = -1
        for i in range(len(sorted_v)):
            draw *= -1
            # print(draw)
            if draw == 1:
                ax.add_patch(Rectangle((0, sorted_v[i][0]), sum(list(qchrs_len.values())), sorted_v[i][1] - sorted_v[i][0],
                                       edgecolor='lightgrey',
                                       facecolor='lightgrey',
                                       fill=True,
                                       alpha=0.25,
                                       lw=0))

    plt.tight_layout()
    plt.savefig(out)
    plt.close()
    logger.info("Finished successfully. Exiting")
# END


def dotplot(args):
    from fixchr.scripts.func import readcoords, readfasta
    cfin = args.infile.name
    ref = args.ref.name
    qry = args.qry.name
    ftype = args.ftype
    fout = args.o.name

    coords = readcoords(cfin, ftype=ftype, f=True, cigar=True)
    coords.to_csv("input_alignments.txt", index=False, header=True, sep='\t')
    # Get dotplot for initial alignments
    refg = readfasta(ref)
    qryg = readfasta(qry)
    rchrs_len = {k: len(v) for k, v in refg.items()}
    qchrs_len = {k: len(v) for k, v in qryg.items()}
    # TODO: Add parsers for AGP files and plot dimensions
    # drawdotplot(coords, rchrs_len, qchrs_len, out='fout')
    drawdotplot(coords, rchrs_len, qchrs_len, out=fout, ragp=args.ragp, qagp=args.qagp, minsize=args.m, width=args.W, height=args.H)
# END


def main(cmd):
    import argparse
    parser = argparse.ArgumentParser("Draw dotplot for whole-genome alignments", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-c", dest="infile", help="File containing alignment coordinates", type=argparse.FileType('r'), required=True)
    parser.add_argument("-r", dest="ref", help="Genome A (which is considered as reference for the alignments). Required for local variation (large indels, CNVs) identification.", type=argparse.FileType('r'))
    parser.add_argument("-q", dest="qry", help="Genome B (which is considered as query for the alignments). Required for local variation (large indels, CNVs) identification.", type=argparse.FileType('r'))
    parser.add_argument('-F', dest="ftype", help="Input file type. T: Table, S: SAM, B: BAM, P: PAF", default="T", choices=['T', 'S', 'B', 'P'])
    parser.add_argument('-o', dest='o', help='Output file name', type=argparse.FileType('w'), default='dotplot.pdf')
    parser.add_argument('--ragp', dest='ragp', help='Reference .agp file describing contig order and orientation', type=argparse.FileType('r'))
    parser.add_argument('--qagp', dest='qagp', help='Query .agp file describing contig order and orientation', type=argparse.FileType('r'))
    parser.add_argument('-m', dest='m', help='Minimum alignment size to plot', type=int, default=1000)
    parser.add_argument('-W', dest='W', help='Width of the output plot', type=int, default=12)
    parser.add_argument('-H', dest='H', help='Width of the output plot', type=int, default=12)
    # TODO: Add parameter to input list of sequences to consider/remove
    args = parser.parse_args()
    # print(args)
    dotplot(args)



