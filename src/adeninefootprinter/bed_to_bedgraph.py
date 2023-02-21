import numpy as np
from .trainmodel import load_ref_file


def bed_to_bedgraph(bed: str, bedgraph: str, ref: str, name: str):
    """
    Calculate the feature/block density in BED files and store in BEDGRAPH
    :param bed: [File1,File2,...] Name of BED file(s)
    :param bedgraph: Name of BEDGRAPH
    :param ref: Reference Dictionary File
    :param name: Name of the BEDGRAPH track
    :return:
    """
    assert isinstance(bed, str)
    assert isinstance(bedgraph, str)
    assert isinstance(ref, str)
    assert isinstance(name, str)

    ref_at_array = load_ref_file(ref)
    aa_read_loc = {c: np.zeros((2, len(a))) for c, a in ref_at_array.items()}
    for fn in bed.split(','):
        fn = fn.strip()
        print('Loading BED file ', fn)
        with open(fn) as filep:
            next(filep)
            for line in filep:
                ent = line.split()
                chrom = ent[0]
                start = int(ent[1])
                end = int(ent[2])
                bstarts = list(map(int, ent[-1].split(',')))
                bsizes = list(map(int, ent[-2].split(',')))
                tarr = aa_read_loc[chrom]
                tarr[1, start: end] += 1
                for s, l in zip(bstarts, bsizes):
                    tarr[0, start + s: start + s + l] += 1

    header = 'track type=bedGraph name="{}" description="{}" visibility=full color=0,0,200\n'.format(name, name)
    tmp = [header]

    for chrom, acc in aa_read_loc.items():
        acc = acc[0, :] / (acc[1, :] + 1e-6)
        acc = np.round(acc * 100, decimals=1)
        diff = np.where(np.diff(acc) != 0)[0]
        p = 0
        for d in diff:
            tmp.append('{chrom}\t{start}\t{end}\t{acc}\n'.format(chrom=chrom, start=p, end=d + 1, acc=acc[p]))
            p = d + 1
        tmp.append('{chrom}\t{start}\t{end}\t{acc}\n'.format(chrom=chrom, start=p, end=len(acc), acc=acc[p]))
    print('Writing BEDGRAPH file', bedgraph)
    with open(bedgraph, 'w') as filep:
        filep.write(''.join(tmp))