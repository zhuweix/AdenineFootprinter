__all__ = ['footprinter', 'indexreference', 'trainmodel', 'predictfootprint']
__version__ = '1.0'
__author__ = 'Zhuwei Xu'

import os
import re
import copy
import pickle
import subprocess
import sys
import pysam
import numpy as np
from Bio.Seq import Seq
from Bio import SeqIO
from statsmodels.stats.multitest import fdrcorrection
import statsmodels.api as sm
from scipy.stats import binom
