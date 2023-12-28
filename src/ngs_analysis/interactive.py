"""Common imports for interactive work
"""
from collections import Counter
from glob2 import glob
import logging
import os
import re
import sys

from natsort import natsorted, natsort_keygen
import numpy as np
import pandas as pd

import warnings
from math import ceil
from functools import partial
from pandas import IndexSlice as pdx
from tqdm.auto import tqdm
import IPython
from IPython.display import display, Image

IPython.get_ipython().run_line_magic('load_ext', 'autoreload')
IPython.get_ipython().run_line_magic('autoreload', '2')
