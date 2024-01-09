from glob import glob
from natsort import natsorted
import numpy as np


def assert_unique(df, *cols):
    """Each argument can be a single column or a list of columns.
    """
    a = df.shape[0]
    for col in cols:
        b = ~df[col].duplicated()
        if a != b.sum():
            counts = df[col].value_counts()
            raise ValueError(
                f'{b.sum()} / {a} entries are unique for column(s) {col}, '
                f'the most frequent duplicate is {counts.index[0]} ({counts.iloc[0]} entries)'
            )
    return df


def nglob(x):
    return natsorted(glob(x))


def sample_at_coverage(coverage, *xs, seed=0):
    rs = np.random.RandomState(seed)
    assert len(set(len(x) for x in xs)) == 1, 'same length'
    n = int(len(xs[0]) * coverage)
    arr = []
    for x in xs:
        arr += [rs.choice(x, size=n)]
    return arr
