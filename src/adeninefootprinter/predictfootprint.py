import pickle
import numpy as np
import json
import copy
import subprocess
import pysam

from .trainmodel import load_ref_file
from .trainmodel import calc_local_m6A_adj_pvalue


def load_model_parameter(model: str):
    """
    Load Model paramters from model json file
    :param model: json file for model parameters
    :return: linear_model, window_size, adj_pval_accessible, adj_pval_nucleosome
    """
    assert isinstance(model, str)
    model_params = {}
    with open(model) as filep:
        tmp = json.load(filep)
    for at, (slope, inter) in tmp['Model'].items():
        at = int(at)
        slope = float(slope)
        inter = float(inter)
        model_params[at] = (slope, inter)
    window = int(tmp['Window'])
    qval_acc = float(tmp['qval_acc'])
    qval_nuc = float(tmp['qval_nuc'])
    print('Model Parameters Loaded')
    print('Window Size ', window)
    print('Adj. p-value for Accessible Region: ', qval_acc)
    print('Adj. p-value for Nucleosome Region: ', qval_nuc)
    return model_params, window, qval_acc, qval_nuc


def predict_footprint_from_bed(bed: str, model: str, prefix: str, ref: str,
                               header_file=None, filter_acc=True, max_amb_fuse=150):
    assert isinstance(bed, str)
    assert isinstance(model, str)
    assert isinstance(prefix, str)
    assert isinstance(ref, str)

    ref_dict = load_ref_file(ref=ref)
    model_params, window, qval_acc, qval_nuc = load_model_parameter(model=model)
    assert window > 2
    assert 0 < qval_acc <= 1
    assert 0 < qval_nuc <= 1
    assert max(model_params) <= window
    assert min(model_params) >= 0

    win_unit = np.ones(window)
    empty_read = 0
    no_acc_read = 0
    all_read = 0
    print('Start Prediction of Nucleosome Footprints and Accessible Regions')
    # header info
    if header_file:
        with pysam.AlignmentFile(header_file, 'rb') as filep:
            header = filep.header
    else:
        header = {
            'HD': {'VN': '1.0'},
            'SQ': []
        }
        for chrom, arr in sorted(ref_dict.items()):
            header['SQ'].append({'LN':len(arr), 'SN': chrom})
    output = '{}.predict.bam'.format(prefix)
    with pysam.AlignmentFile(output, 'wb', header=header) as obam:
        for fn in bed.split(','):
            fn = fn.strip()
            print('Loading ', fn)
            with open(fn) as filep:
                for line in filep:
                    if line.startswith('track'):
                        continue
                    ent = line.split()
                    if len(ent) < 12:
                        continue
                    all_read += 1
                    chrom = ent[0]
                    start = int(ent[1])
                    end = int(ent[2])
                    qname = ent[3]
                    readlength = end - start
                    if end - start < window:
                        continue
                    m6a = int(ent[9])
                    mloc = map(int, ent[11].split(','))
                    tmp = ref_dict[chrom][start: end]
                    readat = np.sum(tmp)
                    if readat < 2:
                        continue
                    readave = m6a / readat * 100
                    tmp_m = np.zeros_like(tmp, dtype=int)
                    tmp = np.convolve(tmp, win_unit, 'same')
                    tmp = tmp.astype(dtype=int)
                    for p in mloc:
                        tmp_m[p] = 1
                    mloc = copy.deepcopy(tmp_m)
                    tmp_m = np.convolve(tmp_m, win_unit, 'same')
                    q_acc, q_nuc = calc_local_m6A_adj_pvalue(
                        m6a_array=tmp_m, at_array=tmp, readave=readave, local_model=model_params
                    )
                    is_acc = q_acc < qval_acc
                    is_nuc = q_nuc < qval_nuc
                    # expand the boundaries
                    is_acc = np.convolve(is_acc, win_unit, 'same') > 0
                    is_nuc = np.convolve(is_nuc, win_unit, 'same') > 0

                    tmpr = np.zeros(readlength, dtype=np.int) + 2  # low-confidence: del
                    tmpr[is_acc] = 0  # Accessible region: match
                    tmpr[is_nuc] = 8  # nucleosome region: sub; boundary between accessible and nucleosome: nucleosome
                    if filter_acc:
                        # Skip reads with no High-confidence accessible regions
                        if np.sum(tmpr == 0) == 0:
                            no_acc_read += 1
                            continue
                    new_read = pysam.AlignedSegment()
                    new_read.query_name = qname
                    new_read.reference_start = start
                    new_read.reference_name = chrom
                    new_read.mapping_quality = 50
                    new_read.query_qualities = [90] * readlength
                    # tidy cigar data
                    edgepos = np.where(np.diff(tmpr) != 0)[0]
                    # Filter read with no feature
                    if len(edgepos) == 0:
                        empty_read += 1
                        continue

                    # combine ambiguous region which attached to one nucleosome end
                    pos0 = 0
                    for idx, pos1 in enumerate(edgepos):
                        tp = tmpr[pos1]
                        if tp != 2:
                            pos0 = pos1 + 1
                            continue
                        if idx == 0:
                            if tmpr[pos1 + 1] == 8 and pos1 <= max_amb_fuse:
                                tmpr[:pos1 + 1] = 8
                        else:
                            if tmpr[pos0 - 1] == 8 or tmpr[pos1 + 1] == 8:
                                if pos1 - pos0 + 1 <= max_amb_fuse:
                                    tmpr[pos0: pos1 + 1] = 8
                        pos0 = pos1 + 1
                    else:
                        if tmpr[pos0] == 2:
                            if tmpr[pos0 - 1] == 8:
                                if len(tmpr) - pos0 <= max_amb_fuse:
                                    tmpr[pos0:] = 8

                    tmpr[mloc] += 1  # label m6A

                    # generate cigar
                    new_cig = []
                    pos0 = 0
                    tp = tmpr[0]
                    for pos1 in np.where(np.diff(tmpr) != 0)[0]:
                        pos1 += 1
                        if tp % 2 == 0:
                            new_cig.append((tp, pos1 - pos0))
                        else:
                            for n in range(pos1 - pos0):
                                new_cig.append((tp - 1, 1))
                                new_cig.append((1, 1))
                        tp = tmpr[pos1]
                        pos0 = pos1
                    if tp % 2 == 0:
                        new_cig.append((tp, len(tmpr) - pos0))
                    else:
                        for n in range(len(tmpr) - pos0):
                            new_cig.append((tp - 1, 1))
                            new_cig.append((1, 1))

                    # merge the 1 bp region due to m6A sites

                    if len(new_cig) == 1:
                        new_read.cigartuples = new_cig
                    else:
                        new_cig2 = []
                        pre = new_cig[0]
                        is_merge = False
                        for idx, c in enumerate(new_cig[1:]):
                            if pre[0] == c[0]:
                                new_cig2.append((pre[0], int(pre[1] + c[1])))
                                is_merge = True
                                pre = c
                            else:
                                if is_merge:
                                    pre = c
                                    is_merge = False
                                else:
                                    new_cig2.append(pre)
                                    pre = c
                        if not is_merge:
                            new_cig2.append(pre)
                        new_read.cigartuples = new_cig2
                    obam.write(new_read)
    print('Prediction Completed. Results stored in ', output)
    if filter_acc:
        print('Analyzed {} Reads, {} Reads were filtered due to no High-quality accessible regions'.format(all_read,
                                                                                                           empty_read + no_acc_read))
    else:
        print('Analyzed {} Reads, {} Reads were filtered due to no High-quality regions'.format(all_read, empty_read))
    output2 = output.replace('.bam', '.sort.bam')
    subprocess.run(['samtools', 'sort', '-o', output2, output], check=True)
    subprocess.run(['samtools', 'index', output2], check=True)
    print('Sorted Results stored in ', output2)


def predict_footprint_from_bam(bam: str, model: str, prefix: str, ref: str,
                               filter_acc=True, max_amb_fuse=150):
    assert isinstance(bam, str)
    assert isinstance(model, str)
    assert isinstance(prefix, str)
    assert isinstance(ref, str)

    ref_dict = load_ref_file(ref=ref)
    model_params, window, qval_acc, qval_nuc = load_model_parameter(model=model)
    assert window > 2
    assert 0 < qval_acc <= 1
    assert 0 < qval_nuc <= 1
    assert max(model_params) <= window
    assert min(model_params) >= 0

    win_unit = np.ones(window)
    empty_read = 0
    no_acc_read = 0
    all_read = 0
    print('Start Prediction of Nucleosome Footprints and Accessible Regions in ', bam)
    # header info
    with pysam.AlignmentFile(bam, 'rb') as inbam:
        header = inbam.header
        output = '{}.predict.sort.bam'.format(prefix)
        with pysam.AlignmentFile(output, 'wb', header=header) as obam:
            for read in inbam:
                new_read = copy.deepcopy(read)
                chrom = read.reference_name
                start = read.reference_start
                cigar = read.cigartuples
                readlength = np.sum([num if tp == 0 else 0 for tp, num in cigar])
                end = start + readlength
                all_read += 1
                if readlength + 1 < window:
                    continue
                tmp = ref_dict[chrom][start: end]
                readlength = len(tmp)
                readat = np.sum(tmp)
                if readat < 2:
                    continue
                p = -1
                m6a = 0
                tmp_m = np.zeros_like(tmp, dtype=int)
                for tp, num in cigar:
                    if tp == 0:
                        if num == 0:
                            continue
                        p += num - 1
                    else:
                        p += 1
                        m6a += 1
                        tmp_m[p] = 1
                readave = m6a / readat * 100

                tmp = np.convolve(tmp, win_unit, 'same')
                tmp = tmp.astype(dtype=int)
                mloc = copy.deepcopy(tmp_m)
                tmp_m = np.convolve(tmp_m, win_unit, 'same')
                tmp_m = tmp_m.astype(dtype=int)
                q_acc, q_nuc = calc_local_m6A_adj_pvalue(
                    m6a_array=tmp_m, at_array=tmp, readave=readave, local_model=model_params
                )
                is_acc = q_acc < qval_acc
                is_nuc = q_nuc < qval_nuc
                # expand the boundaries
                is_acc = np.convolve(is_acc, win_unit, 'same') > 0
                is_nuc = np.convolve(is_nuc, win_unit, 'same') > 0

                tmpr = np.zeros_like(is_acc, dtype=np.int) + 2  # low-confidence: del
                tmpr[is_acc] = 0  # Accessible region: match
                tmpr[is_nuc] = 8  # nucleosome region: sub; boundary between accessible and nucleosome: nucleosome
                if filter_acc:
                    # Skip reads with no High-confidence accessible regions
                    if np.sum(tmpr < 1) == 0:
                        no_acc_read += 1
                        continue
                # tidy cigar data
                edgepos = np.where(np.diff(tmpr) != 0)[0]
                # Filter read with no feature
                if len(edgepos) == 0:
                    empty_read += 1
                    continue

                # combine ambiguous region which attached to one nucleosome end
                pos0 = 0
                for idx, pos1 in enumerate(edgepos):
                    tp = tmpr[pos1]
                    if tp != 2:
                        pos0 = pos1 + 1
                        continue
                    if idx == 0:
                        if tmpr[pos1 + 1] == 8 and pos1 <= max_amb_fuse:
                            tmpr[:pos1 + 1] = 8
                    else:
                        if tmpr[pos0 - 1] == 8 or tmpr[pos1 + 1] == 8:
                            if pos1 - pos0 + 1 <= max_amb_fuse:
                                tmpr[pos0: pos1 + 1] = 8
                    pos0 = pos1 + 1
                else:
                    if tmpr[pos0] == 2:
                        if tmpr[pos0 - 1] == 8:
                            if len(tmpr) - pos0 <= max_amb_fuse:
                                tmpr[pos0:] = 8

                tmpr[mloc > 0] += 1  # label m6A

                # generate cigar
                new_cig = []
                pos0 = 0
                tp = tmpr[0]
                for pos1 in np.where(np.diff(tmpr) != 0)[0]:
                    pos1 += 1
                    if tp % 2 == 0:
                        new_cig.append((tp, pos1 - pos0))
                    else:
                        for n in range(pos1 - pos0):
                            new_cig.append((tp - 1, 1))
                            new_cig.append((1, 1))
                    tp = tmpr[pos1]
                    pos0 = pos1
                if tp % 2 == 0:
                    new_cig.append((tp, len(tmpr) - pos0))
                else:
                    for n in range(len(tmpr) - pos0):
                        new_cig.append((tp - 1, 1))
                        new_cig.append((1, 1))

                # merge the 1 bp region due to m6A sites

                if len(new_cig) == 1:
                    new_read.cigartuples = new_cig
                else:
                    new_cig2 = []
                    pre = new_cig[0]
                    is_merge = False
                    for idx, c in enumerate(new_cig[1:]):
                        if pre[0] == c[0]:
                            new_cig2.append((pre[0], int(pre[1] + c[1])))
                            is_merge = True
                            pre = c
                        else:
                            if is_merge:
                                pre = c
                                is_merge = False
                            else:
                                new_cig2.append(pre)
                                pre = c
                    if not is_merge:
                        new_cig2.append(pre)
                    new_read.cigartuples = new_cig2
                obam.write(new_read)
        print('Prediction Completed. Results stored in ', output)
        if filter_acc:
            print('Analyzed {} Reads, {} Reads were filtered due to no High-quality accessible regions'.format(
                all_read, empty_read + no_acc_read))
        else:
            print('Analyzed {} Reads, {} Reads were filtered due to no High-quality regions'.format(all_read, empty_read))
        subprocess.run(['samtools', 'index', output], check=True)


def predict_footprint(bam: str, bed: str, model: str, prefix: str, ref: str,
                      header_file=None, filter_acc=True, max_amb_fuse=150):
    assert isinstance(model, str)
    assert isinstance(prefix, str)
    assert isinstance(ref, str)
    if bam:
        assert isinstance(bam, str)
        predict_footprint_from_bam(bam=bam, model=model, prefix=prefix, ref=ref,
                                   filter_acc=filter_acc, max_amb_fuse=max_amb_fuse)
    elif bed:
        assert isinstance(bed, str)
        predict_footprint_from_bed(bed=bed, model=model, prefix=prefix, ref=ref, header_file=header_file,
                                   filter_acc=filter_acc, max_amb_fuse=max_amb_fuse)
    else:
        print('Please enter the file for prediction!')

