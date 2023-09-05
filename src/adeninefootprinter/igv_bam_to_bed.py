import pysam

from .trainmodel import load_ref_file


def igv_bam_to_bed(bam: str, acc: str, nuc: str):
    """
    Convert BAM files after prediction to BED format
    :param bam: BAM file name
    :param acc: File name of BED file for accessible regions
    :param nuc: File name of BED file for nucleosome regions
    :return:
    """
    assert isinstance(bam, str)
    assert isinstance(acc, str)
    assert isinstance(nuc, str)
    print('Start Converting ', bam)
    with pysam.AlignmentFile(bam, 'rb') as inbam:
        refs = list(inbam.header.references)
        refs.sort()
        feat_chrom = {r: [] for r in refs}

        for read in inbam:
            chrom = read.reference_name
            cig = read.cigartuples
            name = read.query_name
            start = read.reference_start
            tmp_feat = {0: [], 8: []}
            p = start
            cur_tp = -1
            cur_s = p
            cm6a = 0
            for tp, num in cig:
                if tp != 1:
                    if cur_tp == -1:
                        cur_tp = tp
                        cur_s = p
                    else:
                        if tp != cur_tp:
                            if cur_tp != 2:
                                tmp_feat[cur_tp].append((cur_s, p - 1))
                            cur_tp = tp
                            cur_s = p
                    p += num
                else:
                    cm6a += 1
            if cur_tp != 2:
                tmp_feat[cur_tp].append((cur_s, p - 1))

            feat_chrom[chrom].append(
                (name, start, p, tmp_feat))
    print('Writing Accessible Regions in ', acc)
    content = ['track name=Accessible description="Accessible Region in Reads" useScore=0']
    for chrom, feats in feat_chrom.items():
        for name, start, end, tmp_feat in feats:
            bstart = []
            bsize = []
            if len(tmp_feat[0]) == 0:
                continue
            for i, (s, e) in enumerate(tmp_feat[0]):
                if i == 0:
                    if s != start:
                        bstart.append(0)
                        bsize.append(0)
                bstart.append(s - start)
                bsize.append(e - s + 1)
                if i == len(tmp_feat[0]) - 1:
                    if e != end - 1:
                        bstart.append(end - start - 1)
                        bsize.append(1)
            content.append(
                '{chrom}\t{start}\t{end}\t{name}\t0\t.\t{start}\t{end}\t0\t{cblock}\t{bsizes}\t{bstarts}'.format(
                    chrom=chrom, name=name, start=start, end=end, cblock=len(bsize),
                    bsizes=','.join(map(str, bsize)), bstarts=','.join(map(str, bstart))
                ))
    with open(acc, 'w') as filep:
        filep.write('\n'.join(content))
    print('Writing Nucleosome Regions in ', nuc)
    content = ['track name=Nucleosome description="Nucleosome Region in Reads" useScore=0']
    for chrom, feats in feat_chrom.items():
        for name, start, end, tmp_feat in feats:
            bstart = []
            bsize = []
            if len(tmp_feat[8]) == 0:
                continue
            for i, (s, e) in enumerate(tmp_feat[8]):
                if i == 0:
                    if s != start:
                        bstart.append(0)
                        bsize.append(0)
                bstart.append(s - start)
                bsize.append(e - s + 1)
                if i == len(tmp_feat[8]) - 1:
                    if e != end - 1:
                        bstart.append(end - start - 1)
                        bsize.append(1)
            content.append(
                '{chrom}\t{start}\t{end}\t{name}\t0\t.\t{start}\t{end}\t0\t{cblock}\t{bsizes}\t{bstarts}'.format(
                    chrom=chrom, name=name, start=start, end=end, cblock=len(bsize),
                    bsizes=','.join(map(str, bsize)), bstarts=','.join(map(str, bstart))
                ))
    with open(nuc, 'w') as filep:
        filep.write('\n'.join(content))
    print('Finished Conversion to BED files')
