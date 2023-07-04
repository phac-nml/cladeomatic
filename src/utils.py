import os.path
import shutil
import sys
import time
import psutil
import pandas as pd
import numpy as np
import fastparquet as fp
import tables
from numba import jit
from numba.typed import List
import pyarrow.parquet as pq
import re
from src.constants import MIN_FILE_SIZE


def get_file_length(f):
    return int(os.popen(f'wc -l {f}').read().split()[0])

def get_file_header(f):
    return str(os.popen(f'head -n1 {f}').read())

def get_file_footer(f):
    return str(os.popen(f'tail -n1 {f}').read())

def is_matrix_valid(f):
    num_lines = get_file_length(f)
    footer = get_file_footer(f).split("\t")
    if num_lines == len(footer):
        return True
    return False

def is_file_ok(f):
    status = True
    if not os.path.isfile(f):
        status = False
    elif get_file_length(f) < 2:
        status = False
    elif os.path.getsize(f) < MIN_FILE_SIZE:
        status = False

    return status

