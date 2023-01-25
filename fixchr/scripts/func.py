import logging

def setlogconfig(lg, fin=''):
    """
    :param lg: Log-level
    :return:
    """
    import logging.config
    logging.config.dictConfig({
        'version': 1,
        'disable_existing_loggers': False,
        'formatters': {
            'log_file': {
                'format': "%(asctime)s - %(name)s - %(levelname)s - %(funcName)s:%(lineno)d - %(message)s",
            },
            'stdout': {
                'format': "%(name)s - %(levelname)s - %(message)s",
            },
        },
        'handlers': {
            'stdout': {
                'class': 'logging.StreamHandler',
                'formatter': 'stdout',
                'level': 'WARNING',
            },
            # 'log_file': {
            #     'class': 'logging.FileHandler',
            #     'filename': fin,
            #     'mode': 'a',
            #     'formatter': 'log_file',
            #     'level': lg,
            # },
        },
        'loggers': {
            '': {
                'level': lg,
                'handlers': ['stdout'],
                # 'handlers': ['stdout', 'log_file'],
            },
        },
    })
#END

def mylogger(logname):
    from hometools.classes import CustomFormatter
    import logging
    logger = logging.getLogger(logname)
    handler = logging.StreamHandler()
    handler.setFormatter(CustomFormatter())
    logger.addHandler(handler)
    logging.basicConfig(level=logging.INFO)
    logger.propagate = False
    return logger
# END

class CustomFormatter(logging.Formatter):
    '''
    https://betterstack.com/community/questions/how-to-color-python-logging-output/
    '''

    grey = "\x1b[0;49;90m"
    green = "\x1b[0;49;32m"
    yellow = "\x1b[0;49;93m"
    red = "\x1b[0;49;31m"
    bold_red = "\x1b[0;49;31;21m"
    reset = "\x1b[0m"
    format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s (%(filename)s:%(lineno)d)"

    FORMATS = {
        logging.DEBUG: grey + format + reset,
        logging.INFO: green + format + reset,
        logging.WARNING: yellow + format + reset,
        logging.ERROR: red + format + reset,
        logging.CRITICAL: bold_red + format + reset
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)
# END

def readfasta(f):
    # TODO: This takes too long when used with getchr for large genomes. Try to optimise FASTA reading when the entire genome is not needed.
    from gzip import open as gzopen
    from gzip import BadGzipFile
    from collections import deque
    import sys
    out = {}
    chrid = ''
    chrseq = deque()
    # Test if the file is Gzipped or not
    with gzopen(f, 'rb') as fin:
        try:
            fin.read(1)
            isgzip = True
        except BadGzipFile:
            isgzip = False
    try:
        if isgzip:
            with gzopen(f, 'rb') as fin:
                for line in fin:
                    if b'>' in line:
                        if chrid != '':
                            out[chrid] = ''.join(chrseq)
                            chrid = line.strip().split(b'>')[1].split(b' ')[0].decode()
                            chrseq = deque()
                        else:
                            chrid = line.strip().split(b'>')[1].split(b' ')[0].decode()
                        if chrid in out.keys():
                            sys.exit(" Duplicate chromosome IDs are not accepted. Chromosome ID {} is duplicated. Provided chromosome with unique IDs".format(chrid))
                    else:
                        chrseq.append(line.strip().decode())
        else:
            with open(f, 'r') as fin:
                for line in fin:
                    if '>' in line:
                        if chrid != '':
                            out[chrid] = ''.join(chrseq)
                            chrid = line.strip().split('>')[1].split(' ')[0]
                            chrseq = deque()
                        else:
                            chrid = line.strip().split('>')[1].split(' ')[0]
                        if chrid in out.keys():
                            sys.exit(" Duplicate chromosome IDs are not accepted. Chromosome ID {} is duplicated. Provided chromosome with unique IDs".format(chrid))
                    else:
                        chrseq.append(line.strip())
    except Exception as e:
        raise Exception(e)

    if chrid != '':
        out[chrid] = ''.join(chrseq)
    # TODO: add check for the validation of input fasta files
    return out
# END


def writefasta(fa, f):
    """
    :param fa: dictionary. Keys are chromosome ids. Values are sequence.
    :param f: Output file
    :return:
    Can output .gz file if output file name ends with .gz
    """
    # TODO: write bgzip file
    # from pysam import tabix_compress as gzopen
    from gzip import open as gzopen
    isgzip = f.rsplit('.', 1)[-1] == 'gz'
    # b, nl >> Begin character, new line
    op, openstr, b, nl = (gzopen, 'wb', b'>', b'\n') if isgzip else (open, 'w', '>', '\n')
    with op(f, openstr) as fo:
        for k, v in fa.items():
            if isgzip:
                k = k.encode()
                v = v.encode()
            fo.write(b+k+nl)
            fo.write(nl.join([v[i:i+60] for i in range(0, len(v), 60)]) + nl)
            # for i in range(0, len(v), 60):
            #     fo.write(v[i:(i+60)]+nl)
# END


def seqsize(args):
    fin = args.fasta.name
    out = [(chrom, len(seq)) for chrom, seq in readfasta(fin).items()]
    #TODO: Add output file
    for i in out:
        print(i[0], i[1], sep="\t")
    if args.t:
        print("Genome_length", sum([i[1] for i in out]), sep="\t")
# END


def samtocoords(f):
    import sys
    import logging
    from pandas import DataFrame
    from collections import deque
    logger = logging.getLogger('SAM reader')
    rc = {}  # Referece chromosomes
    rcs = {}  # Selected chromosomes
    al = deque()  # Individual alignment
    try:
        with open(f, 'r') as fin:
            for l in fin:
                if l[:3] == '@SQ':
                    c, s = 0, 0
                    for h in l.strip().split()[1:]:
                        h = h.split(':')
                        if h[0] == 'SN': c = h[1]
                        if h[0] == 'LN': s = int(h[1])
                    rcs[c] = s
                    continue
                elif l[0] == '@':
                    continue

                l = l.split('\t')[:6]
                # if l[1] == '2064': break
                if l[2] == '*':
                    logger.warning(l[
                                       0] + ' do not align with any reference sequence and cannot be analysed. Remove all unplaced scaffolds and contigs from the assemblies.')  # Skip rows corresponding to non-mapping sequences (contigs/scaffolds)
                    continue

                if 'M' in l[5]:
                    logger.error(
                        'Incorrect CIGAR string found. CIGAR string can only have I/D/H/S/X/=. CIGAR STRING: ' + l[5])
                    sys.exit()
                cgt = [[int(j[0]), j[1]] for j in [i.split(';') for i in
                                                   l[5].replace('S', ';S,').replace('H', ';H,').replace('=',
                                                                                                        ';=,').replace(
                                                       'X', ';X,').replace('I', ';I,').replace('D', ';D,').split(',')[
                                                   :-1]]]
                if len(cgt) > 2:
                    if True in [True if i[1] in ['S', 'H'] else False for i in cgt[1:-1]]:
                        logger.error(
                            f"Incorrect CIGAR string found. Clipped bases inside alignment. H/S can only be in the terminal. CIGAR STRING: {cgt}")
                        sys.exit()

                bf = '{:012b}'.format(int(l[1]))

                rs = int(l[3])
                re = rs - 1 + sum([i[0] for i in cgt if i[1] in ['X', '=', 'D']])

                if bf[7] == '0':  # forward alignment
                    if cgt[0][1] == '=':
                        qs = 1
                    elif cgt[0][1] in ['S', 'H']:
                        qs = cgt[0][0] + 1
                    else:
                        print('ERROR: CIGAR string starting with non-matching base')
                    qe = qs - 1 + sum([i[0] for i in cgt if i[1] in ['X', '=', 'I']])
                elif bf[7] == '1':  # inverted alignment
                    if cgt[-1][1] == '=':
                        qs = 1
                    elif cgt[-1][1] in ['S', 'H']:
                        qs = cgt[-1][0] + 1
                    else:
                        print('ERROR: CIGAR string starting with non-matching base')
                    qe = qs - 1 + sum([i[0] for i in cgt if i[1] in ['X', '=', 'I']])
                    qs, qe = qe, qs

                al.append([
                    rs,
                    re,
                    qs,
                    qe,
                    abs(re - rs) + 1,
                    abs(qs - qe) + 1,
                    format((sum([i[0] for i in cgt if i[1] == '=']) / sum(
                        [i[0] for i in cgt if i[1] in ['=', 'X', 'I', 'D']])) * 100, '.2f'),
                    1,
                    1 if bf[7] == '0' else -1,
                    l[2],
                    l[0],
                    "".join([str(i[0]) + i[1] for i in cgt if i[1] in ['=', 'X', 'I', 'D']])
                ])
                rcs[l[2]] = 1
            rcs = list(rcs.keys())
            for k in list(rc.keys()):
                if k not in rcs: logger.warning(l[
                                                    0] + ' do not align with any query sequence and cannot be analysed. Remove all unplaced scaffolds and contigs from the assemblies.')
    except Exception as e:
        logger.error('Error in reading SAM file: ' + str(e))
        sys.exit()
    al = DataFrame(list(al))
    al[6] = al[6].astype('float')
    al.sort_values([9, 0, 1, 2, 3, 10], inplace=True, ascending=True)
    al.index = range(len(al.index))
    return al
# END


def readsambam(fin, ftype='B'):
    import pysam
    import logging
    import sys
    import numpy as np
    import pandas as pd
    logger = logging.getLogger('Reading BAM/SAM file')
    try:
        if ftype == 'B':
            findata = pysam.AlignmentFile(fin, 'rb')
        elif ftype == 'S':
            return samtocoords(fin)
        else:
            raise ValueError("Wrong parameter")
    except ValueError as e:
        logger.error("Error in opening BAM/SAM file. " + str(e))
        sys.exit()
    except OSError as e:
        logger.error("Error in reading input file." + str(e))
        sys.exit()
    except Exception as e:
        logger.error("Unexpected error in opening BAM/SAM file. " + str(e))
        sys.exit()

    try:
        qry_prim = {}
        ref_prim = {}
        cgdict = {1: 'I', 2: 'D', 7: '=', 8: 'X'}
        coords = {}
        index = 0
        for aln in findata:
            index += 1
            ## Check whether every sequence has at least one primary alignment
            if aln.reference_name is not None:
                if aln.reference_name not in ref_prim.keys():
                    ref_prim[aln.reference_name] = False
            if aln.query_name not in qry_prim.keys():
                qry_prim[aln.query_name] = False
            if aln.reference_name is not None:
                if not ref_prim[aln.reference_name]:
                    if aln.flag < 256:
                        ref_prim[aln.reference_name] = True
            if not qry_prim[aln.query_name]:
                if aln.flag < 256:
                    qry_prim[aln.query_name] = True

            ## Pass non-alinging chromosomes
            if aln.cigarstring is None:
                logger.warning(aln.query_name + ' do not align with any reference chromosome and cannot be analysed')
                continue

            ## Check CIGAR:
            if False in [False if i[0] not in [1, 2, 4, 5, 7, 8] else True for i in aln.cigartuples]:
                logger.error(
                    "Incorrect CIGAR string found. CIGAR string can only have I/D/H/S/X/=. CIGAR STRING: " + str(
                        aln.cigarstring))
                sys.exit()
            if len(aln.cigartuples) > 2:
                if True in [True if i[0] in [4, 5] else False for i in aln.cigartuples[1:-1]]:
                    logger.error(
                        "Incorrect CIGAR string found. Clipped bases inside alignment. H/S can only be in the terminal. CIGAR STRING: " + aln.cigarstring)
                    sys.exit()

            ## Parse information from the aln object
            astart = aln.reference_start + 1
            aend = aln.reference_end
            is_inv = True if np.binary_repr(aln.flag, 12)[7] == '1' else False
            if not is_inv:
                if aln.cigartuples[0][0] in [4, 5]:
                    bstart = aln.cigartuples[0][1] + 1
                else:
                    bstart = 1
                bend = bstart + aln.query_alignment_length - 1
            else:
                if aln.cigartuples[-1][0] in [4, 5]:
                    bend = aln.cigartuples[-1][1] + 1
                else:
                    bend = 1
                bstart = bend + aln.query_alignment_length - 1
            alen = abs(aend - astart) + 1
            blen = abs(bend - bstart) + 1
            iden = format((sum([i[1] for i in aln.cigartuples if i[0] == 7]) / sum(
                [i[1] for i in aln.cigartuples if i[0] in [1, 2, 7, 8]])) * 100, '.2f')
            adir = 1
            bdir = -1 if is_inv else 1
            achr = aln.reference_name
            bchr = aln.query_name
            cg = "".join([str(i[1]) + cgdict[i[0]] for i in aln.cigartuples if i[0] not in [4, 5]])
            coords[index] = [astart, aend, bstart, bend, alen, blen, iden, adir, bdir, achr, bchr, cg]

        ## Give warning for chromosomes which do not have any primary alignment
        for k, v in ref_prim.items():
            if not v:
                logger.warning(
                    'No primary alignment found for reference sequence ' + k + '. This could mean that the entire chromosome ' + k + ' is reapeated.')
        for k, v in qry_prim.items():
            if not v:
                logger.warning(
                    'No primary alignment found for query sequence ' + k + '. This could mean that the entire chromosome ' + k + ' is reapeated.')

        ## Return alignments
        coords = pd.DataFrame.from_dict(coords, orient='index')
        coords.sort_values([9, 0, 1, 2, 3, 10], inplace=True, ascending=True)
        coords.index = range(len(coords.index))
        coords[6] = coords[6].astype('float')
        return coords
    except Exception as e:
        logger.error("Error in reading BAM/SAM file. " + str(e))
        sys.exit()
# END


def readpaf(paf):
    import logging
    import sys
    from collections import deque
    import pandas as pd
    coords = deque()
    logger = logging.getLogger('Reading BAM/SAM file')
    try:
        with open(paf, 'r') as fin:
            for line in fin:
                line = line.strip().split()
                astart = int(line[7]) + 1
                aend = int(line[8])
                adir = 1
                bdir = 1 if line[4] == '+' else -1
                bstart = int(line[2]) + 1 if bdir == 1 else int(line[3])
                bend = int(line[3]) if bdir == 1 else int(line[2]) + 1
                alen = abs(aend - astart) + 1
                blen = abs(bend - bstart) + 1 if bdir == 1 else bstart - bend + 1
                cg = [i.split(":")[-1] for i in line[12:] if i[:2] == 'cg']
                if len(cg) != 1:
                    logger.error("CIGAR string is not present in PAF at line {}. Exiting.".format("\t".join(line)))
                    sys.exit()
                cg = cg[0]
                ## Check CIGAR:
                if not all([True if i[1] in {'I', 'D', 'H', 'S', 'X', '='} else False for i in cgtpl(cg)]):
                    logger.error(
                        "Incorrect CIGAR string found. CIGAR string can only have I/D/H/S/X/=. CIGAR STRING: " + str(
                            cg))
                    sys.exit()
                if len(cgtpl(cg)) > 2:
                    if any([True if i[1] in {'H', 'S'} else False for i in cgtpl(cg)]):
                        logger.error(
                            "Incorrect CIGAR string found. Clipped bases inside alignment. H/S can only be in the terminal. CIGAR STRING: " + str(
                                cg))
                        sys.exit()

                iden = round((sum([int(i[0]) for i in cgtpl(cg) if i[1] == '=']) / sum(
                    [int(i[0]) for i in cgtpl(cg) if i[1] in {'=', 'X', 'D', 'I'}])) * 100, 2)
                achr = line[5]
                bchr = line[0]
                coords.append([astart, aend, bstart, bend, alen, blen, iden, adir, bdir, achr, bchr, cg])
        coords = pd.DataFrame(coords)
        coords.sort_values([9, 0, 1, 2, 3, 10], inplace=True, ascending=True)
        coords.index = range(len(coords.index))
        coords[6] = coords[6].astype('float')
        return coords
    except FileNotFoundError:
        logger.error("Cannot open {} file. Exiting".format(paf))
        sys.exit()
    except ValueError as e:
        logger.error("Error in reading PAF: {}. Exiting".format(e))
        sys.exit()
# END


def readcoords(coordsfin, ftype, f, cigar=False):
    import logging
    import pandas as pd
    import numpy as np
    import sys
    logger = logging.getLogger('Reading Coords')
    logger.debug(ftype)
    if ftype == 'T':
        logger.info("Reading input from .tsv file")
        try:
            coords = pd.read_table(coordsfin, header=None)
        except pd.errors.ParserError:
            coords = pd.read_table(coordsfin, header=None, engine="python")
        except Exception as e:
            logger.error("Error in reading the alignment file. " + str(e))
            sys.exit()
    elif ftype == 'S':
        logger.info("Reading input from SAM file")
        try:
            coords = readsambam(coordsfin, ftype='S')
        except Exception as e:
            logger.error("Error in reading the alignment file. " + str(e))
            sys.exit()
    elif ftype == 'B':
        logger.info("Reading input from BAM file")
        try:
            coords = readsambam(coordsfin, ftype='B')
        except Exception as e:
            logger.error("Error in reading the alignment file" + str(e))
            sys.exit()
    elif ftype == 'P':
        logger.info("Reading input from PAF file")
        try:
            coords = readpaf(coordsfin)
        except Exception as e:
            logger.error("Error in reading the alignment file" + str(e))
            sys.exit()
    else:
        logger.error("Incorrect alignment file type specified.")
        sys.exit()

    if not cigar:
        if coords.shape[1] >= 12:
            coords = coords.iloc[:, 0:11]
        coords.columns = ["aStart", "aEnd", "bStart", "bEnd", "aLen", "bLen", "iden", "aDir", "bDir", "aChr", "bChr"]
    else:
        if coords.shape[1] > 12:
            coords = coords.iloc[:, 0:12]
        coords.columns = ["aStart", "aEnd", "bStart", "bEnd", "aLen", "bLen", "iden", "aDir", "bDir", "aChr", "bChr",
                          'cigar']

    # Sanity check input file
    try:
        coords.aStart = coords.aStart.astype('int')
    except ValueError:
        logger.error('astart is not int')
        sys.exit()

    try:
        coords.aEnd = coords.aEnd.astype('int')
    except ValueError:
        logger.error('aend is not int')
        sys.exit()

    try:
        coords.bStart = coords.bStart.astype('int')
    except ValueError:
        logger.error('bstart is not int')
        sys.exit()

    try:
        coords.bEnd = coords.bEnd.astype('int')
    except ValueError:
        logger.error('abend is not int')
        sys.exit()

    try:
        coords.aLen = coords.aLen.astype('int')
    except ValueError:
        logger.error('alen is not int')
        sys.exit()

    try:
        coords.bLen = coords.bLen.astype('int')
    except ValueError:
        logger.error('blen is not int')
        sys.exit()

    try:
        coords.iden = coords.iden.astype('float')
    except ValueError:
        logger.error('iden is not float')
        sys.exit()

    try:
        coords.aDir = coords.aDir.astype('int')
    except ValueError:
        logger.error('aDir is not int')
        sys.exit()

    if any(coords.aDir != 1):
        logger.error('aDir can only have values 1')
        sys.exit()

    try:
        coords.bDir = coords.bDir.astype('int')
    except ValueError:
        logger.error('bDir is not int')
        sys.exit()

    for i in coords.bDir:
        if i not in [1, -1]:
            logger.error('bDir can only have values 1/-1')
            sys.exit()

    try:
        coords.aChr = coords.aChr.astype(str)
    except:
        logger.error('aChr is not string')
        sys.exit()

    try:
        coords.bChr = coords.bChr.astype(str)
    except:
        logger.error('bChr is not string')
        sys.exit()

    # Filter small alignments
    if f:
        coords = coords.loc[coords.iden > 90]
        coords = coords.loc[(coords.aLen > 100) & (coords.bLen > 100)]

    ## check for bstart > bend when bdir is -1
    check = np.unique(coords.loc[coords.bDir == -1, 'bStart'] > coords.loc[coords.bDir == -1, 'bEnd'])
    if len(check) > 1:
        logger.error(
            'Inconsistent start and end position for inverted alignment in query genome. For inverted alignments, either all bstart < bend or all bend > bstart')
        sys.exit()
    elif len(check) == 0:
        logger.info('No Inverted alignments present.')
    elif check[0]:
        pass
    else:
        logger.info('For inverted alignments, bstart was less than bend. Swapping them.')
        coords.loc[coords.bDir == -1, 'bStart'] = coords.loc[coords.bDir == -1, 'bStart'] + coords.loc[
            coords.bDir == -1, 'bEnd']
        coords.loc[coords.bDir == -1, 'bEnd'] = coords.loc[coords.bDir == -1, 'bStart'] - coords.loc[
            coords.bDir == -1, 'bEnd']
        coords.loc[coords.bDir == -1, 'bStart'] = coords.loc[coords.bDir == -1, 'bStart'] - coords.loc[
            coords.bDir == -1, 'bEnd']
    coords.sort_values(['aChr', 'aStart', 'aEnd', 'bChr', 'bStart', 'bEnd'], inplace=True)
    return coords
# END

inttr = lambda x: [int(x[0]), x[1]]
def cgtpl(cg, to_int=False):
    """
    Takes a cigar string as input and returns a cigar tuple
    """
    for i in "MIDNSHPX=":
        cg = cg.replace(i, ';'+i+',')
    if to_int:
        return [inttr(i.split(';')) for i in cg.split(',')[:-1]]
    else:
        return [i.split(';') for i in cg.split(',')[:-1]]
# END


def mergeRanges(ranges):
    """
    Take a 2D numpy array, with each row as a range and return merged ranges
    i.e. ranges which are overlapping would be combined.
    :param ranges:
    :return:
    """
    from collections import deque
    import numpy as np
    if len(ranges) < 2:
        return ranges
    for i in ranges:
        if i[0] > i[1]:
            i[1], i[0] = i[0], i[1]
    ranges = ranges[ranges[:, 0].argsort()]
    min_value = ranges[0, 0]
    max_value = ranges[0, 1]
    out_range = deque()
    for i in ranges[1:]:
        if i[0] > max_value:
            out_range.append([min_value, max_value])
            min_value = i[0]
            max_value = i[1]
        elif i[1] > max_value:
            max_value = i[1]
    out_range.append([min_value, max_value])
    return np.array(out_range)
# END


def homchr(coords, rchrs_len, qchrs_len, csize=100000):
    '''
    Selects alignments corresponding to homologous chromosomes.
    When the chromosome IDs are different, then the smaller contigs are filtered
    out.
    Homologous chromosomes are selected based on the number of bases aligned
    between them
    :param coords: Input alignments (pandas.core.frame.DataFrame)
    :param rchrs_len: Reference chromosome lengths (dict)
    :param qchrs_len: Query chromosome lengths (dict)
    :param csize: Contig length threshold (int)
    :return:
    '''
    from collections import defaultdict
    import numpy as np
    # import logging
    logger = mylogger("homchr")
    if np.unique(coords.aChr).tolist() != np.unique(coords.bChr).tolist():
        logger.warning('Chromosomes IDs do not match in the two genomes')
        # Filter contigs smaller than csize
        logger.debug('Removing smaller contigs.')
        sr = [k for k, v in rchrs_len.items() if v < csize]
        sq = [k for k, v in qchrs_len.items() if v < csize]
        coords2 = coords.loc[~coords.aChr.isin(sr)].copy()
        coords2 = coords2.loc[~coords2.bChr.isin(sq)]
    else:
        coords2 = coords
    chromrmap = defaultdict(dict)
    chromqmap = defaultdict(dict)
    for i in np.unique(coords2.bChr):
        for j in np.unique(coords2.aChr):
            a = np.array(coords2.loc[(coords2.bChr == i) & (coords2.aChr == j), ["aStart", "aEnd"]])
            a = mergeRanges(a)
            chromrmap[j][i] = len(a) + (a[:, 1] - a[:, 0]).sum()
            a = np.array(coords2.loc[(coords2.bChr == i) & (coords2.aChr == j), ["bStart", "bEnd"]])
            a = mergeRanges(a)
            chromqmap[i][j] = len(a) + (a[:, 1] - a[:, 0]).sum()
    assigned = {}
    for chrom in chromrmap:
        qmax = max(chromrmap[chrom].items(), key=lambda x: x[1])[0]
        rmax = max(chromqmap[qmax].items(), key=lambda x: x[1])[0]
        if rmax == chrom:
            assigned[chrom] = qmax
    logger.debug('Selecting homologous chromosomes.')
    coords2 = coords2.loc[coords2.aChr.isin(list(assigned.keys())) & coords2.bChr.isin(list(assigned.values()))]
    return coords2, assigned
# END


def checkdir(coords, assigned):
    import numpy as np
    import logging
    from collections import deque
    logger = logging.getLogger("checkdir")
    al = coords.copy()
    rv = deque()
    achrs = np.unique(al.aChr).tolist()
    for achr in achrs:
        dir_range = mergeRanges(np.array(al.loc[(al.aChr == achr) & (al.bChr == assigned[achr]) & (al.bDir == 1), ["aStart", "aEnd"]]))
        dir_len = len(dir_range) + (dir_range[:, 1] - dir_range[:, 0]).sum()
        inv_range = mergeRanges(np.array(al.loc[(al.aChr == achr) & (al.bChr == assigned[achr]) & (al.bDir == -1), ["aStart", "aEnd"]]))
        inv_len = len(inv_range) + (inv_range[:, 1] - inv_range[:, 0]).sum()
        if inv_len > dir_len:
            rv.append(assigned[achr])
    return rv
# END


def revcomp(seq):
    assert type(seq) == str
    old = 'ACGTRYKMBDHVacgtrykmbdhv'
    rev = 'TGCAYRMKVHDBtgcayrmkvhdb'
    tab = str.maketrans(old, rev)
    return seq.translate(tab)[::-1]
# END

