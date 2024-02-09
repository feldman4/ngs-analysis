from glob import glob
from natsort import natsorted
import numpy as np
import pandas as pd
import re
import os
from string import Formatter


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


def csv_frame(files_or_search, progress=lambda x: x, add_file=None, file_pat=None,  
              include_cols=None, exclude_cols=None, ignore_missing=True, 
              ignore_empty=True, ignore_index=True, sort=False, format='auto', **kwargs):
    """Convenience function, pass either a list of files or a 
    glob wildcard search term.

    :param files_or_search: a single search string or a list of files; a search
        string can use format fields enclosed by braces
    :param progress: a progress function, such as tqdm
    :param add_file: filenames will be stored in this column
    :param file_pat: a regex with capture groups applied to filenames
    :param include_cols: columns to include
    :param exclude_cols: columns to exclude
    :param ignore_missing: if given a list of files, ignore the missing ones
    :param ignore_empty: ignore empty files
    :param ignore_index: passed to `pd.concat`
    :param sort: passed to `pd.concat`
    """
    import parse
    reader = pd.read_csv
    if format == 'auto' and files_or_search.endswith('.pq'):
        reader = pd.read_parquet
    elif format in ('parquet', 'pq'):
        reader = pd.read_parquet


    def read(f):
        try:
            df = reader(f, **kwargs)
        except pd.errors.EmptyDataError:
            return None
        if add_file is not None:
            df[add_file] = f
        if include_cols is not None:
            keep = [x for x in df.columns if re.match(include_cols, x)]
            df = df[keep]
        if exclude_cols is not None:
            keep = [x for x in df.columns if not re.match(exclude_cols, x)]
            df = df[keep]
        if file_pat is not None:
            match = re.match(f'.*?{file_pat}.*', f)
            if match is None:
                raise ValueError(f'{file_pat} failed to match {f}')
            if match.groupdict():
                for k,v in match.groupdict().items():
                    df[k] = v
            else:
                if add_file is None:
                    raise ValueError(f'must provide `add_file` or named groups in {file_pat}')
                first = match.groups()[0]
                df[add_file] = first
        return df
    
    # set up `files` and `extra_fields` depending on if `files_or_search` is a list, format string,
    # or regular glob string
    extra_fields = {}
    if isinstance(files_or_search, str):
        if '{' in files_or_search:
            search = files_or_search
            fieldnames = [fname for _, fname, _, _ in Formatter().parse(search) if fname is not None]
            search_glob = search.format('*', **{x: '*' for x in fieldnames})
            files = nglob(search_glob)
            extra_fields.update({f: parse.parse(search, f).named for f in files})
        else:
            files = natsorted(glob(files_or_search))
    else:
        files = files_or_search
        if ignore_missing:
            files = [f for f in files if os.path.exists(f) or 'gs://'    in f]
    
    arr = []
    for f in progress(files):
        df = read(f)
        if df is None:
            if ignore_empty:
                continue
            else:
                raise ValueError(f'Empty table: {f}')
        arr += [df.assign(**extra_fields.get(f, {}))]

    return pd.concat(arr, ignore_index=ignore_index, sort=sort)


def assign_format(df, **kwargs):
    """Apply a format string to each row.
    df.pipe(assign_format, new_col='{a}.{b}')
    """
    for k, v in kwargs.items():
        df = df.assign(
            **{k: lambda x: x.apply(lambda row: v.format(**row), axis=1)})
    return df
