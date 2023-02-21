import pickle
import numpy as np
import json
import pandas as pd
import statsmodels.api as sm
from scipy.stats import binom
from statsmodels.stats.multitest import fdrcorrection


# minimal number of windows for Model calculation
min_window_count = 10
# min quantile for OLS to filter valid read m6A
min_quantile = .05
# max quantile for OLS to filter valid read m6A
max_quantile = .99
# minimal valid read average-local average point for OLS
min_ols_point = 3


def load_ref_file(ref: str):
    """
    Load Reference Dictionary
    :param ref: Name of the reference index generated by footprint index
    :return: dictionary of arrays for AT locations
    """
    assert isinstance(ref, str)
    with open(ref, 'rb') as filep:
        ref_dict = pickle.load(filep)
    return ref_dict


def load_bed_for_local_m6a(bed_list: str, ref_dict: dict, window: int, prefix: str):
    """
    Load BED files and calculate the average local m6A at given read average m6A (rounded to 1%)
    :param bed_list: list of BED files, separated by ,
    :param ref_dict: reference dictionary for AT locations
    :param window: Local window size
    :param prefix: prefix of the output files
    :return: Pandas DataFrame of local m6A information
    """
    assert isinstance(bed_list, str)
    assert isinstance(ref_dict, dict)
    assert isinstance(window, int)
    assert window > 3
    assert isinstance(prefix, str)
    print('Start Calculation of Local Methylation Model')
    bed_list = bed_list.split(',')
    win_unit = np.ones(window, dtype=np.int)
    local_m6a = {}
    for fn in bed_list:
        fn = fn.strip()
        print('Loading BED file {}'.format(fn))
        with open(fn) as filep:
            for line in filep:
                if line.startswith('track'):
                    continue
                ent = line.split()
                if len(ent) < 12:
                    continue
                chrom = ent[0]
                start = int(ent[1])
                end = int(ent[2])
                m6a = int(ent[9])
                mloc = map(int, ent[11].split(','))
                tmp = ref_dict[chrom][start: end]
                if end - start < window:
                    continue
                readat = np.sum(tmp)
                if readat < 3:
                    continue
                readm6a = m6a / readat * 100
                readm6a = int(round(readm6a, 0))
                local_m6a.setdefault(readm6a, {})
                tmp_m = np.zeros_like(tmp, dtype=int)
                tmp = np.convolve(tmp, win_unit, 'valid')
                tmp = tmp.astype(dtype=int)
                for p in mloc:
                    tmp_m[p] = 1
                tmp_m = np.convolve(tmp_m, win_unit, 'valid')
                for p in np.where(tmp > 2)[0]:
                    at = tmp[p]
                    m6a = tmp_m[p]
                    local_m6a[readm6a].setdefault(at, [0, 0])
                    local_m6a[readm6a][at][0] += m6a
                    local_m6a[readm6a][at][1] += 1
        print('Finished Loading {}'.format(fn))
    tmp_pd = []
    for readm6a, rdict in local_m6a.items():
        for at, ct in rdict.items():
            # skip rare windows to exclude strange results
            if ct[1] < min_window_count:
                continue
            ave = ct[0] / ct[1] / at * 100
            tmp_pd.append({'ATContext': at,
                           'AveLocalMethylation': ave,
                           'ReadMethylation': readm6a,
                           'NumWindow': ct[1]})
    tmp_pd = pd.DataFrame(tmp_pd)
    info_name = '{}.Window={}.csv'.format(prefix, window)
    print('Writing local m6A info in {}'.format(info_name))
    tmp_pd.to_csv(info_name)
    print('Finished Loading BED Files')
    return tmp_pd


def calc_ols_model(local_m6a: pd.DataFrame, window: int, prefix: str):
    """
    Calculate the local m6A model using OLS
    :param local_m6a: DataFrame of local m6A from load_bed_for_local_m6a
    :param window: Window size
    :param prefix: Prefix for output file
    :return: Dictionary of Model parameters
    """
    assert isinstance(local_m6a, pd.DataFrame)
    assert 'AveLocalMethylation' in local_m6a
    assert 'ATContext' in local_m6a
    assert 'ReadMethylation' in local_m6a
    assert 'NumWindow' in local_m6a
    assert isinstance(window, int)
    assert isinstance(prefix, str)

    # Obtain read m6A range for OLS
    # The windows for each m6A level (rounded to 1%) are counted
    # The read m6A from min_percentile to max_percentile
    print('Start Model Setup')
    tmp = []
    print('Start Setup Read m6A threshold')
    print('Min Quantile: {}, Max Quantile: {}'.format(min_quantile, max_quantile))
    total = np.sum(local_m6a['NumWindow'])
    for ent in local_m6a.groupby('ReadMethylation')['NumWindow'].agg(np.sum).items():
        tmp.append((ent[0], ent[1]))
    tmp.sort()
    min_count = total * min_quantile
    max_count = total * max_quantile
    count = 0
    min_m6a = -1
    max_m6a = -1
    for read_m6a, ct in tmp:
        count += ct
        if min_m6a < 0 and count >= min_count:
            min_m6a = read_m6a
        if count >= max_count:
            max_m6a = read_m6a
            break
    print('Min Read m6A: {} Max Read m6A: {}'.format(min_m6a, max_m6a))
    print('Start OLS Model')
    # Exclude Extreme Average Read conditions
    local_m6a = local_m6a[(local_m6a['ReadMethylation'] >= min_m6a) &
                          (local_m6a['ReadMethylation'] <= max_m6a)]
    tmp_model = {'Model': {}, 'Window': window, 'MinReadm6A': min_m6a, 'MaxReadm6A': max_m6a}
    model_info = []
    for at in range(3, window + 1):
        tmp = local_m6a[local_m6a['ATContext'] == at]
        if len(tmp) < min_ols_point:
            continue
        X = tmp['ReadMethylation']
        Y = tmp['AveLocalMethylation']
        X = sm.add_constant(X)
        model = sm.OLS(Y, X)
        results = model.fit()
        const, slope = results.params
        r2 = results.rsquared
        model_info.append({'ATContext': at, 'Slope': slope,
                           'Intercept': const, 'R2': r2})
        tmp_model['Model'][at] = (slope, const)
    ols_name = '{}.OLS.csv'.format(prefix)
    model_info = pd.DataFrame(model_info)
    model_info.to_csv(ols_name)
    print('OLS table saved in {}'.format(ols_name))
    print('Model Setup Finished.')
    return tmp_model


def calc_local_m6A_adj_pvalue(local_model: dict, m6a_array: np.ndarray,
                              at_array: np.ndarray, readave: float):
    """
    Calculate the adjusted p-value for a single read.
    :param local_model: Linear model between Local m6A and Read m6A
    :param m6a_array: Local m6A count in windows
    :param at_array: Local AT context in windows
    :param readave: Average Read m6A (%)
    :return: adj_pvalue_accessible and adj_pvalue_nucleosome
    """
    adj_ra = np.array([local_model[x][0] * readave / 100 + local_model[x][1] / 100
                       if x in local_model else np.nan for x in at_array])
    pval_acc = binom.sf(m6a_array, at_array, adj_ra) + binom.pmf(m6a_array, at_array, adj_ra)
    np.nan_to_num(pval_acc, nan=1, copy=False)
    _, qval_acc = fdrcorrection(pval_acc, 0.05, method='p')
    pval_nuc = binom.cdf(m6a_array, at_array, adj_ra)
    np.nan_to_num(pval_nuc, nan=1, copy=False)
    _, qval_nuc = fdrcorrection(pval_nuc, 0.05, method='p')
    return qval_acc, qval_nuc


def calc_adj_pvalue(bed_list: str, model_params: dict, prefix: str, ref_dict: dict, fp_rate: float):
    """
    Calculate the adjusted p-value threshold to predict Nucleosome and accessible region using gDNA control
    :param bed_list: BED files with m6A sites in gDNA control
    :param model_params: Linear Model between Local m6A and read m6A
    :param prefix: Prefix of the output data
    :param ref_dict: Reference Dictionary for AT locations
    :param fp_rate: False-positive Nucleosome and Accessible regions in gDNA for adj. pvalues
    :return: None
    """
    assert isinstance(bed_list, str)
    assert isinstance(model_params, dict)
    assert 'Window' in model_params
    assert 'Model' in model_params
    assert isinstance(prefix, str)
    assert isinstance(ref_dict, dict)

    qval_acc = []
    qval_nuc = []

    window = model_params['Window']
    win_unit = np.ones(window)
    local_model = model_params['Model']
    print('Start False-Positive Estimation in gDNA control for Adjusted P-value thresholds')
    print('Loading BED files for gDNA samples')
    for fn in bed_list.split(','):
        fn = fn.strip()
        print('Loading BED file {}'.format(fn))
        with open(fn) as filep:
            for line in filep:
                if line.startswith('track'):
                    continue
                ent = line.split()
                if len(ent) < 12:
                    continue
                chrom = ent[0]
                start = int(ent[1])
                end = int(ent[2])
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
                tmp_m = np.convolve(tmp_m, win_unit, 'same')
                q_acc, q_nuc = calc_local_m6A_adj_pvalue(
                    m6a_array=tmp_m, at_array=tmp, readave=readave, local_model=local_model
                )
                qval_acc.append(q_acc)
                qval_nuc.append(q_nuc)
    print('Finished Loading BED files')
    qval_acc = np.concatenate(qval_acc)
    qval_nuc = np.concatenate(qval_nuc)
    q_acc = np.quantile(qval_acc, fp_rate)
    q_nuc = np.quantile(qval_nuc, fp_rate)
    q_acc = np.round(q_acc, 3)
    q_nuc = np.round(q_nuc, 3)
    print('Adj. p-value for Accessible Regions: {:.3f}'.format(q_acc))
    print('Adj. p-value for Nucleosome Regions: {:.3f}'.format(q_nuc))
    model_params['qval_acc'] = q_acc
    model_params['qval_nuc'] = q_nuc
    model_params['fp_rate'] = fp_rate
    model_fn = '{}.model.json'.format(prefix)
    print('Model Parameters saved in {}'.format(model_fn))
    with open(model_fn, 'w') as filep:
        json.dump(model_params, filep, indent=2)


def train_model(bed_list: str, window: int, prefix: str, fprate: float, ref: str):
    assert isinstance(bed_list, str)
    assert isinstance(window, int)
    assert window > 2
    assert isinstance(prefix, str)
    assert isinstance(fprate, float)
    assert 0 < fprate <= 1
    assert isinstance(ref, str)
    print('Start Model Training')
    ref_dict = load_ref_file(ref=ref)
    local_m6a = load_bed_for_local_m6a(bed_list=bed_list,
                                       ref_dict=ref_dict,
                                       window=window,
                                       prefix=prefix)
    model_params = calc_ols_model(local_m6a=local_m6a,
                                  window=window,
                                  prefix=prefix)
    calc_adj_pvalue(bed_list=bed_list,
                    model_params=model_params,
                    prefix=prefix,
                    ref_dict=ref_dict,
                    fp_rate=fprate)
    print('Training Finished')










