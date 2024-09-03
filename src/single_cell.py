from __future__ import annotations
import os
import re
import signal
from decimal import Decimal
from datetime import date, datetime, time, timedelta
from itertools import chain, islice
from packaging import version
from pathlib import Path
from textwrap import fill
from typing import Any, Callable, Dict, ItemsView, Iterable, KeysView, \
    Literal, Mapping, Sequence, ValuesView, Union


class ignore_sigint:
    """
    Ignore Ctrl + C when importing certain modules, to avoid errors due to
    incomplete imports.
    """
    def __enter__(self):
        signal.signal(signal.SIGINT, signal.SIG_IGN)
    
    def __exit__(self, *_):
        signal.signal(signal.SIGINT, signal.default_int_handler)


with ignore_sigint():
    import h5py
    import numpy as np
    import polars as pl
    from scipy.sparse import csr_array, csc_array, csr_matrix, csc_matrix, \
        hstack, vstack
    from scipy.stats import rankdata

from .utils import bonferroni, bincount, check_bounds, check_dtype, \
    check_R_variable_name, check_type, check_types, cython_inline, \
    cython_type, fdr, filter_columns, generate_palette, getnnz, is_integer, \
    plural, prange, sparse_matrix_vector_op, Timer, to_tuple

Color = Union[str, float, np.floating,
              tuple[Union[int, np.integer], Union[int, np.integer],
                    Union[int, np.integer]],
              tuple[Union[int, np.integer], Union[int, np.integer],
                    Union[int, np.integer], Union[int, np.integer]],
              tuple[Union[float, np.floating], Union[float, np.floating],
                    Union[float, np.floating]],
              tuple[Union[float, np.floating], Union[float, np.floating],
                    Union[float, np.floating], Union[float, np.floating]]]
Indexer = Union[int, np.integer, str, slice,
                np.ndarray[1, Union[np.integer, np.bool_]], pl.Series,
                list[Union[int, np.integer, str, bool, np.bool_]]]
Scalar = Union[str, int, float, Decimal, date, time, datetime, timedelta, bool,
               bytes]
NestedScalarOrArrayDict = \
    Dict[str, Union[str, int, np.integer, float, np.floating, bool, np.bool_,
         np.ndarray[Any, Any], 'NestedScalarOrArrayDict']]
SingleCellColumn = \
    Union[str, pl.Expr, pl.Series, np.ndarray,
          Callable[['SingleCell'], Union[pl.Series, np.ndarray]]]
PseudobulkColumn = \
    Union[str, pl.Expr, pl.Series, np.ndarray,
          Callable[['Pseudobulk', str], Union[pl.Series, np.ndarray]]]


class SingleCell:
    """
    A lightweight alternative to AnnData for representing single-cell data.
    
    Has slots for:
    - X: a scipy sparse array of counts per cell and gene
    - obs: a polars DataFrame of cell metadata
    - var: a polars DataFrame of gene metadata
    - obsm: a dictionary of NumPy arrays and polars DataFrames of cell metadata
    - varm: a dictionary of NumPy arrays and polars DataFrames of gene metadata
    - uns: a dictionary of scalars (strings, numbers or Booleans) or NumPy
           arrays, or nested dictionaries thereof
    as well as `obs_names` and `var_names`, aliases for obs[:, 0] and
    var[:, 0].
    
    Why is X a sparse array rather than matrix? Aside from being more modern
    and consistent with np.array, it's also faster, as explained at
    github.com/scipy/scipy/blob/2aee5efcbe3720f41fe55f336f492ae0acbecdee/scipy/
    sparse/_base.py#L1333-L1337.
    """
    # noinspection PyUnresolvedReferences
    def __init__(self,
                 X: csr_array | csc_array | csr_matrix | csc_matrix |
                    'AnnData' | str | Path,
                 obs: pl.DataFrame | None = None,
                 var: pl.DataFrame | None = None,
                 obsm: dict[str, np.ndarray[2, Any] | pl.DataFrame] |
                       None = None,
                 varm: dict[str, np.ndarray[2, Any] | pl.DataFrame] |
                       None = None,
                 uns: NestedScalarOrArrayDict | None = None,
                 *,
                 X_key: str | None = None,
                 assay: str | None = None,
                 obs_columns: str | Iterable[str] = None,
                 var_columns: str | Iterable[str] = None,
                 num_threads: int | np.integer | None = 1) -> None:
        """
        Load a SingleCell dataset from a file, or create one from an in-memory
        AnnData object or count matrix + metadata.
        
        The supported file types are:
        - AnnData (.h5ad)
        - 10x (.h5 or .mtx)
        - Seurat (.rds) - requires the ryp Python-R bridge
        - SingleCellExperiment (.rds) - requires the ryp Python-R bridge
        
        By default, when an AnnData file, AnnData object, Seurat file, or
        SingleCellExperiment file contains both raw and normalized counts, only
        the raw counts will be loaded. To load normalized counts instead, use
        the `X_key` argument.
        
        To create a SingleCell dataset from an in-memory Seurat or
        SingleCellExperiment object, use `SingleCell.from_seurat()` or
        `SingleCell.from_sce()`.
        
        Args:
            X: the data as a sparse array or matrix (with rows = cells,
               columns = genes), AnnData object, AnnData .h5ad file, 10x .h5
               or .mtx.gz file, or Seurat or SingleCellExperiment .rds file. If
               `X` is a 10x .mtx.gz file, barcodes.tsv.gz and features.tsv.gz
               are assumed to be in the same directory, unless custom paths to
               these files are specified via the `obs` and/or `var` arguments.
            obs: a polars DataFrame of metadata for each cell (row of X), or
                 if X is a 10x .mtx.gz file, an optional filename for
                 cell-level metadata (which is otherwise assumed to be at
                 barcodes.tsv.gz in the same directory as the .mtx.gz file)
            var: a polars DataFrame of metadata for each gene (column of X), or
                 if X is a 10x .mtx.gz file, an optional filename for
                 gene-level metadata (which is otherwise assumed to be at
                 features.tsv.gz in the same directory as the .mtx.gz file)
            obsm: a dictionary of NumPy arrays and polars DataFrames of
                  metadata for each cell. Keys must be strings.
            varm: a dictionary of NumPy arrays and polars DataFrames of
                  metadata for each gene. Keys must be strings.
            uns: a dictionary of unstructured metadata. Keys must be strings;
                 values can be scalars (strings, numbers or Booleans), NumPy
                 arrays, or nested dictionaries thereof.
            X_key: when X is an AnnData object or .h5ad or .rds filename, which
                   location within the object or file to use as X.
                   - If X is an AnnData object, the count matrix to use as X.
                     If None, defaults to `self.layers['UMIs']` or `self.raw.X`
                     if present, otherwise `self.X`.
                   - If X is an .h5ad filename, the name of the key in the
                     .h5ad file to use as X. If None, defaults to
                     `'layers/UMIs'` or `'raw/X'` if present, otherwise `'X'`.
                     Tip: use `SingleCell.ls(h5ad_file)` to see the structure
                     of an .h5ad file without loading it, to figure out which
                     key to use as `X_key`.
                   - If X is a Seurat .rds filename, the slot within the active
                     assay (or the assay specified by the `assay` argument, if
                     not None) to use as X. If None, defaults to `'counts'`.
                     If available, set to `'data'` to load the normalized
                     counts, or `'scale.data'` to load the normalized and
                     scaled counts.
                   - If X is a SingleCellExperiment .rds filename, the element
                     within `sce_object@assays@data` to use as `X`. If None,
                     defaults to `'counts'`. If available, set to `'logcounts'`
                     to load the normalized counts.
            assay: if X is a Seurat .rds filename, the name of the assay within
                   the Seurat object to load data from. Defaults to
                   `seurat_object@active_assay` (usually `'RNA'`).
            obs_columns: if X is an .h5ad filename, the columns of obs to load.
                         If not specified, load all columns. Specifying only a
                         subset of columns can speed up reading. Not supported
                         for .h5 files, since they only have a single obs
                         column (`'barcodes'`), nor for Seurat and
                         SingleCellExperiment .rds files, since .rds files do
                         not support partial loading.
            var_columns: if X is an .h5ad or .h5 filename, the columns of var
                         to load. If not specified, load all columns.
                         Specifying only a subset of columns can speed up
                         reading. Not supported for Seurat and
                         SingleCellExperiment .rds files, since .rds files do
                         not support partial loading.
            num_threads: the number of threads to use when reading .h5ad and
                         .h5 files; set `num_threads=None` to use all available
                         cores (as determined by `os.cpu_count()`)
        
        Note:
            Both ordered and unordered categorical columns of obs and var will
            be loaded as polars Enums rather than polars Categoricals. This is
            because polars Categoricals use a shared numerical encoding across
            columns, so their codes are not [0, 1, 2, ...] like pandas
            categoricals and polars Enums are. Using Categoricals leads to a
            large overhead (~25%) when loading obs from an .h5ad file, for
            example.
        
        Note:
            SingleCell does not support dense matrices, which are highly
            memory-inefficient for single-cell data. Passing a NumPy array as
            the `X` argument will give an error; if for some reason your data
            has been improperly stored as a dense matrix, convert it to a
            sparse matrix first with `csr_matrix(numpy_array)`). However, when
            loading from disk or converting from other formats, dense matrices
            will be automatically converted to sparse matrices, to avoid giving
            an error when loading or converting.
        """
        type_string = str(type(X))
        is_anndata = type_string.startswith("<class 'anndata")
        if is_anndata:
            with ignore_sigint():
                from anndata import AnnData
            if not isinstance(X, AnnData):
                is_anndata = False
        is_filename = isinstance(X, (str, Path))
        if is_filename:
            X = str(X)
            is_h5ad = X.endswith('.h5ad')
            is_h5 = X.endswith('.h5')
            is_hdf5 = is_h5ad or is_h5
        else:
            is_h5ad = False
            is_hdf5 = False
        if is_h5ad:
            if obs_columns is not None:
                obs_columns = to_tuple(obs_columns)
                if len(obs_columns) == 0:
                    error_message = 'obs_columns is empty'
                    raise ValueError(error_message)
                check_types(obs_columns, 'obs_columns', str, 'strings')
        else:
            if obs_columns is not None:
                error_message = (
                    'obs_columns can only be specified when loading an .h5ad '
                    'file')
                raise ValueError(error_message)
        if is_hdf5:
            if var_columns is not None:
                var_columns = to_tuple(var_columns)
                if len(var_columns) == 0:
                    error_message = 'var_columns is empty'
                    raise ValueError(error_message)
                check_types(var_columns, 'var_columns', str, 'strings')
            if num_threads is None:
                num_threads = os.cpu_count()
            else:
                check_type(num_threads, 'num_threads', int,
                           'a positive integer')
                check_bounds(num_threads, 'num_threads', 1)
        else:
            if var_columns is not None:
                error_message = (
                    'var_columns can only be specified when loading an .h5ad '
                    'or .h5 file')
                raise ValueError(error_message)
            if num_threads != 1:
                error_message = (
                    'num_threads can only be specified when loading an .h5ad '
                    'or .h5 file')
                raise ValueError(error_message)
        if isinstance(X, (csr_array, csc_array, csr_matrix, csc_matrix)):
            for prop, prop_name in (X_key, 'X_key'), (assay, 'assay'):
                if prop is not None:
                    error_message = (
                        f'when X is a sparse array or matrix, {prop_name} '
                        f'must be None')
                    raise ValueError(error_message)
            check_type(obs, 'obs', pl.DataFrame, 'a polars DataFrame')
            check_type(var, 'var', pl.DataFrame, 'a polars DataFrame')
            obsm = {} if obsm is None else obsm.copy()
            varm = {} if varm is None else varm.copy()
            uns = {} if uns is None else uns.copy()
            for field, field_name in (obsm, 'obsm'), (varm, 'varm'):
                for key, value in field.items():
                    if not isinstance(key, str):
                        error_message = (
                            f'all keys of {field_name} must be strings, but '
                            f'it contains a key of type '
                            f'{type(key).__name__!r}')
                        raise TypeError(error_message)
                    if isinstance(value, np.ndarray):
                        if value.ndim != 2:
                            error_message = (
                                f'all values of {field_name} must be 2D NumPy '
                                f'arrays or polars DataFrames, but '
                                f'{field_name}[{key!r}] is a {value.ndim:,}D '
                                f'NumPy array')
                            raise ValueError(error_message)
                    elif not isinstance(value, pl.DataFrame):
                        error_message = (
                            f'all values of {field_name} must be NumPy '
                            f'arrays or polars DataFrames, but {field_name}'
                            f'[{key!r}] has type {type(value).__name__!r}')
                        raise TypeError(error_message)
                valid_uns_types = str, int, np.integer, float, np.floating, \
                    bool, np.bool_, np.ndarray
                for description, value in SingleCell._iter_uns(uns):
                    if not isinstance(value, valid_uns_types):
                        error_message = (
                            f'all values of uns must be scalars (strings, '
                            f'numbers or Booleans) or NumPy arrays, or nested '
                            f'dictionaries thereof, but {description} has '
                            f'type {type(value).__name__!r}')
                        raise TypeError(error_message)
            if isinstance(X, csr_matrix):
                X = csr_array(X)
            if isinstance(X, csc_matrix):
                X = csc_array(X)
            self._X = X
            self._obs = obs
            self._var = var
            self._obsm = obsm
            self._varm = varm
            self._uns = uns
        elif is_filename:
            filename = os.path.expanduser(X)
            # noinspection PyUnboundLocalVariable
            if is_h5ad:
                if not os.path.exists(filename):
                    error_message = f'.h5ad file {X!r} does not exist'
                    raise FileNotFoundError(error_message)
                for prop, prop_name in (obs, 'obs'), (var, 'var'), \
                        (obsm, 'obsm'), (varm, 'varm'), (uns, 'uns'), \
                        (assay, 'assay'):
                    if prop is not None:
                        error_message = (
                            f'when loading an .h5ad file, {prop_name} must be '
                            f'None')
                        raise ValueError(error_message)
                # See anndata.readthedocs.io/en/latest/fileformat-prose.html
                # for the AnnData on-disk format specification
                with h5py.File(filename) as h5ad_file:
                    # Load obs and var
                    self._obs = SingleCell._read_h5ad_dataframe(
                        h5ad_file, 'obs', columns=obs_columns,
                        num_threads=num_threads)
                    self._var = SingleCell._read_h5ad_dataframe(
                        h5ad_file, 'var', columns=var_columns,
                        num_threads=num_threads)
                    # Load obsm
                    if 'obsm' in h5ad_file:
                        obsm = h5ad_file['obsm']
                        self._obsm = {
                            key: value[:]
                                 if isinstance(value, h5py.Dataset) else
                                 SingleCell._read_h5ad_dataframe(
                                     h5ad_file, f'obsm/{key}',
                                     num_threads=num_threads)
                            for key, value in obsm.items()}
                    else:
                        self._obsm = {}
                    # Load varm
                    if 'varm' in h5ad_file:
                        varm = h5ad_file['varm']
                        self._varm = {
                            key: value[:]
                                 if isinstance(value, h5py.Dataset) else
                                 SingleCell._read_h5ad_dataframe(
                                     h5ad_file, f'varm/{key}',
                                     num_threads=num_threads)
                            for key, value in varm.items()}
                    else:
                        self._varm = {}
                    # Load uns
                    if 'uns' in h5ad_file:
                        self._uns = SingleCell._read_uns(h5ad_file['uns'])
                    else:
                        self._uns = {}
                    # Load X
                    if X_key is None:
                        has_layers_UMIs = 'layers/UMIs' in h5ad_file
                        has_raw_X = 'raw/X' in h5ad_file
                        if has_layers_UMIs and has_raw_X:
                            error_message = (
                                "both layers['UMIs'] and raw.X are present; "
                                "this should never happen in well-formed "
                                "AnnData files")
                            raise ValueError(error_message)
                        X_key = 'layers/UMIs' if has_layers_UMIs else \
                            'raw/X' if has_raw_X else 'X'
                    else:
                        check_type(X_key, 'X_key', str, 'a string')
                        if X_key not in h5ad_file:
                            error_message = (
                                f'X_key {X_key!r} is not present in the .h5ad '
                                f'file')
                            raise ValueError(error_message)
                    X = h5ad_file[X_key]
                    matrix_class = X.attrs['encoding-type'] \
                        if 'encoding-type' in X.attrs else \
                        X.attrs['h5sparse_format'] + '_matrix'
                    if matrix_class == 'csr_matrix':
                        array_class = csr_array
                    elif matrix_class == 'csc_matrix':
                        array_class = csc_array
                    else:
                        error_message = (
                            f"X has unsupported encoding-type "
                            f"{matrix_class!r}, but should be 'csr_matrix' or "
                            f"'csc_matrix' (csr is preferable for speed)")
                        raise ValueError(error_message)
                    self._X = array_class((
                        self._read_dataset(X['data'], num_threads),
                        self._read_dataset(X['indices'], num_threads),
                        self._read_dataset(X['indptr'], num_threads)),
                        shape=X.attrs['shape'] if 'shape' in X.attrs else
                              X.attrs['h5sparse_shape'])
            else:
                # noinspection PyUnboundLocalVariable
                if is_h5:
                    if not os.path.exists(filename):
                        error_message = f'.h5 file {X!r} does not exist'
                        raise FileNotFoundError(error_message)
                    for prop, prop_name in (obs, 'obs'), (var, 'var'), \
                            (obsm, 'obsm'), (varm, 'varm'), (uns, 'uns'), \
                            (X_key, 'X_key'), (assay, 'assay'):
                        if prop is not None:
                            error_message = (
                                f'when loading an .h5 file, {prop_name} must '
                                f'be None')
                            raise ValueError(error_message)
                    with h5py.File(filename) as h5_file:
                        matrix = h5_file['matrix']
                        features = matrix['features']
                        self._obs = pl.Series('barcodes',
                                              matrix['barcodes'][:])\
                            .cast(pl.String)\
                            .to_frame()
                        var_columns = \
                            ['name', 'id', 'feature_type', 'genome'] + \
                            [column for column in
                             ('pattern', 'read', 'sequence')
                             if column in features]
                        self._var = pl.DataFrame([
                            pl.Series(column, features[column][:])
                            .cast(pl.String) for column in var_columns])
                        self._X = csr_array((
                            self._read_dataset(matrix['data'], num_threads),
                            self._read_dataset(matrix['indices'], num_threads),
                            self._read_dataset(matrix['indptr'], num_threads)),
                            shape=matrix['shape'][:][::-1])
                        self._obsm = {}
                        self._varm = {}
                        self._uns = {}
                elif X.endswith('.mtx.gz'):
                    if not os.path.exists(filename):
                        error_message = f'10x file {X} does not exist'
                        raise FileNotFoundError(error_message)
                    for prop, prop_name, prop_description in (
                            (obs, 'obs',
                             'a barcodes.tsv.gz file of cell-level metadata'),
                            (var, 'var',
                             'a features.tsv.gz file of gene-level metadata')):
                        if prop is not None and not \
                                isinstance(prop, (str, Path)):
                            error_message = (
                                f'when loading a 10x .mtx.gz file, '
                                f'{prop_name} must be None or the path to '
                                f'{prop_description}')
                            raise TypeError(error_message)
                    for prop, prop_name in \
                            (obsm, 'obsm'), (varm, 'varm'), (uns, 'uns'), \
                            (X_key, 'X_key'), (assay, 'assay'):
                        if prop is not None:
                            error_message = (
                                f'when loading an .h5ad file, {prop_name} '
                                f'must be None')
                            raise ValueError(error_message)
                    from scipy.io import mmread
                    self._X = csr_array(mmread(X).T.tocsr())
                    self._obs = pl.read_csv(
                        f'{os.path.dirname(X)}/barcodes.tsv.gz'
                        if obs is None else obs,
                        has_header=False, new_columns=['cell'])
                    self._var = pl.read_csv(
                        f'{os.path.dirname(X)}/features.tsv.gz'
                        if var is None else var,
                        has_header=False, new_columns=['gene'])
                    self._obsm = {}
                    self._varm = {}
                    self._uns = {}
                elif X.endswith('.rds'):
                    if not os.path.exists(filename):
                        error_message = f'.rds file {X} does not exist'
                        raise FileNotFoundError(error_message)
                    for prop, prop_name in (obs, 'obs'), (var, 'var'), \
                            (obsm, 'obsm'), (varm, 'varm'), (uns, 'uns'):
                        if prop is not None:
                            error_message = (
                                f'when loading an .rds file, {prop_name} '
                                f'must be None')
                            raise ValueError(error_message)
                    from ryp import r, to_py, to_r
                    r(f'.SingleCell.object = readRDS({X!r})')
                    try:
                        if X_key is None:
                            X_key = 'counts'
                        classes = to_py('class(.SingleCell.object)',
                                        squeeze=False)
                        if len(classes) == 1:
                            if classes[0] == 'Seurat':
                                r('suppressPackageStartupMessages('
                                  'library(SeuratObject))')
                                self._X, self._obs, self._var, self._obsm, \
                                    self._uns = SingleCell._from_seurat(
                                        '.SingleCell.object', assay=assay,
                                        slot=X_key, constructor=True)
                                self._varm = {}
                            elif classes[0] == 'SingleCellExperiment':
                                if assay is not None:
                                    error_message = (
                                        f'when loading a SingleCellExperiment '
                                        f'.rds file, assay must be None')
                                    raise ValueError(error_message)
                                r('suppressPackageStartupMessages('
                                  'library(SingleCellExperiment))')
                                self._X, self._obs, self._var, self._obsm, \
                                    self._uns = SingleCell._from_sce(
                                        '.SingleCell.object', slot=X_key,
                                        constructor=True)
                                self._varm = {}
                            else:
                                error_message = (
                                    f'the R object loaded from {X} must be a '
                                    f'Seurat or SingleCellExperiment object, '
                                    f'but has class {classes[0]!r}')
                                raise TypeError(error_message)
                        elif len(classes) == 0:
                            error_message = (
                                f'the R object loaded from {X} must be a '
                                f'Seurat or SingleCellExperiment object, but '
                                f'has no class')
                            raise TypeError(error_message)
                        else:
                            classes_string = \
                                ', '.join(f'{c!r}' for c in classes[:-1])
                            error_message = (
                                f'the R object loaded from {X} must be a '
                                f'Seurat object, but has classes '
                                f'{classes_string} and {classes[-1]!r}')
                            raise TypeError(error_message)
                    finally:
                        r('rm(.SingleCell.object)')
                else:
                    error_message = (
                        f'X is a filename with unknown extension '
                        f'{".".join(X.split(".")[1:])}; it must be .h5ad '
                        f'(AnnData), .rds (Seurat) or .mtx.gz (10x)')
                    raise ValueError(error_message)
                if obs_columns is not None:
                    self._obs = self._obs.select(obs_columns)
                if var_columns is not None:
                    self._var = self._var.select(var_columns)
        elif is_anndata:
            for prop, prop_name in (obs, 'obs'), (var, 'var'), \
                    (obsm, 'obsm'), (varm, 'varm'), (uns, 'uns'), \
                    (assay, 'assay'):
                if prop is not None:
                    error_message = (
                        f'when initializing a SingleCell dataset from an '
                        f'AnnData object, {prop_name} must be None')
                    raise ValueError(error_message)
            if X_key is None:
                has_layers_UMIs = 'UMIs' in X._layers
                has_raw_X = hasattr(X._raw, '_X')
                if has_layers_UMIs and has_raw_X:
                    error_message = (
                        "both layers['UMIs'] and raw.X are present; this "
                        "should never happen in well-formed AnnData objects")
                    raise ValueError(error_message)
                counts = X._layers['UMIs'] if has_layers_UMIs else \
                    X._raw._X if has_raw_X else X._X
                if not isinstance(counts, (
                        csr_array, csc_array, csr_matrix, csc_matrix)):
                    error_message = (
                        f'to initialize a SingleCell dataset from an AnnData '
                        f'object, its X must be a csr_array, csc_array, '
                        f'csr_matrix, or csc_matrix, but it has type '
                        f'{type(counts).__name__!r}. Either convert X to a '
                        f'csr_array or csc_array, or specify a custom X via '
                        f'the X_key argument')
                    raise TypeError(error_message)
            else:
                check_type(X_key, 'X_key',
                           (csr_array, csc_array, csr_matrix, csc_matrix),
                           'a csr_array, csc_array, csr_matrix, or csc_matrix')
                counts = X_key
            self._X = counts if isinstance(counts, (csr_array, csc_array)) \
                else csr_array(counts) if isinstance(counts, csr_matrix) else \
                    csc_matrix(counts)
            for attr in '_obs', '_var':
                df = getattr(X, attr)
                if df.index.name is None:
                    df = df.rename_axis('_index')  # for consistency with .h5ad
                # Convert Categoricals with string categories to polars Enums,
                # and Categoricals with other types of categories (e.g.
                # integers) to non-categorical columns of the corresponding
                # polars dtype (e.g. pl.Int64) since polars only supports
                # string categories
                schema_overrides = {}
                cast_dict = {}
                for column, dtype in \
                        df.dtypes[df.dtypes == 'category'].items():
                    categories_dtype = dtype.categories.dtype
                    if categories_dtype == object:
                        schema_overrides[column] = pl.Enum(dtype.categories)
                    else:
                        cast_dict[column] = categories_dtype
                setattr(self, attr, pl.from_pandas(
                    df.astype(cast_dict), schema_overrides=schema_overrides,
                    include_index=True))
            self._obsm: dict[str, np.ndarray] = dict(X._obsm)
            self._varm: dict[str, np.ndarray] = dict(X._varm)
            self._uns = dict(X._uns)
        else:
            error_message = (
                f'X must be a csc_array, csr_array, csc_matrix, or '
                f'csr_matrix, an AnnData object, or an .h5ad (AnnData), .rds '
                f'(Seurat) or .mtx.gz (10x) filename, but has type '
                f'{type(X).__name__!r}')
            raise TypeError(error_message)
        # Check that shapes match
        num_cells, num_genes = self._X.shape
        if len(self._obs) == 0:
            error_message = 'len(obs) is 0: no cells remain'
            raise ValueError(error_message)
        if len(self._var) == 0:
            error_message = 'len(var) is 0: no genes remain'
            raise ValueError(error_message)
        if len(self._obs) != num_cells:
            error_message = (
                f'len(obs) is {len(self._obs):,}, but X.shape[0] is '
                f'{num_cells:,}')
            raise ValueError(error_message)
        if len(self._var) != num_genes:
            error_message = (
                f'len(var) is {len(self._var):,}, but X.shape[1] is '
                f'{num_genes:,}')
            raise ValueError(error_message)
        for key, value in self._obsm.items():
            if len(value) != num_cells:
                error_message = (
                    f'len(obsm[{key!r}]) is {len(value):,}, but X.shape[0] is '
                    f'{num_cells:,}')
                raise ValueError(error_message)
        for key, value in self._varm.items():
            if len(value) != num_genes:
                error_message = (
                    f'len(varm[{key!r}]) is {len(value):,}, but X.shape[0] is '
                    f'{num_genes:,}')
                raise ValueError(error_message)
        # Set `uns['normalized']` and `uns['QCed']` to False if not set yet;
        # if set and not a Boolean, move to `uns['_normalized']`/`uns['_QCed']`
        for key in 'normalized', 'QCed':
            if key in self._uns:
                if not isinstance(self._uns[key], bool):
                    import warnings
                    new_key = f'_{key}'
                    warning_message = (
                        f'uns[{key!r}] already exists and is not Boolean; '
                        f'moving it to uns[{new_key!r}]')
                    warnings.warn(warning_message)
                    self._uns[new_key] = self._uns[key]
                    self._uns[key] = False
            else:
                self._uns[key] = False
    
    @property
    def X(self) -> csr_array | csc_array:
        return self._X
    
    @X.setter
    def X(self, X: csr_array | csc_array | csr_matrix | csc_matrix) -> None:
        if isinstance(X, (csr_array, csc_array)):
            pass
        elif isinstance(X, csr_matrix):
            X = csr_array(X)
        elif isinstance(X, csc_matrix):
            X = csc_array(X)
        else:
            error_message = (
                f'new X must be a csr_array, csc_array, csr_matrix, or '
                f'csc_matrix, but has type {type(X).__name__!r}')
            raise TypeError(error_message)
        if X.shape != self._X.shape:
            error_message = (
                f'new X is {X.shape[0]:,} × {X.shape[1]:,}, but old X is '
                f'{self._X.shape[0]:,} × {self._X.shape[1]:,}')
            raise ValueError(error_message)
        self._X = X
    
    @property
    def obs(self) -> pl.DataFrame:
        return self._obs
    
    @obs.setter
    def obs(self, obs: pl.DataFrame) -> None:
        check_type(obs, 'obs', pl.DataFrame, 'a polars DataFrame')
        if len(obs) != len(self._obs):
            error_message = (
                f'new obs has length {len(obs):,}, but old obs has length '
                f'{len(self._obs):,}')
            raise ValueError(error_message)
        self._obs = obs

    @property
    def var(self) -> pl.DataFrame:
        return self._var
    
    @var.setter
    def var(self, var: pl.DataFrame) -> None:
        check_type(var, 'var', pl.DataFrame, 'a polars DataFrame')
        if len(var) != len(self._var):
            error_message = (
                f'new var has length {len(var):,}, but old var has length '
                f'{len(self._var):,}')
            raise ValueError(error_message)
        self._var = var
    
    @property
    def obsm(self) -> dict[str, np.ndarray[2, Any] | pl.DataFrame]:
        return self._obsm
    
    @obsm.setter
    def obsm(self, obsm: dict[str, np.ndarray[2, Any] | pl.DataFrame]) -> None:
        num_cells = len(self)
        for key, value in obsm.items():
            if not isinstance(key, str):
                error_message = (
                    f'all keys of obsm must be strings, but new obsm contains '
                    f'a key of type {type(key).__name__!r}')
                raise TypeError(error_message)
            if isinstance(value, np.ndarray):
                if value.ndim != 2:
                    error_message = (
                        f'all values of obsm must be 2D NumPy arrays or '
                        f'polars DataFrames, but new obsm[{key!r}] is a '
                        f'{value.ndim:,}D NumPy array')
                    raise ValueError(error_message)
            elif not isinstance(value, pl.DataFrame):
                error_message = (
                    f'all values of obsm must be NumPy arrays or polars '
                    f'DataFrames, but new obsm[{key!r}] has type '
                    f'{type(value).__name__!r}')
                raise TypeError(error_message)
            if len(obsm) != num_cells:
                error_message = (
                    f'the length of new obsm[{key!r}] is {len(value):,}, but '
                    f'X.shape[0] is {self._X.shape[0]:,}')
                raise ValueError(error_message)
        self._obsm = obsm.copy()
    
    @property
    def varm(self) -> dict[str, np.ndarray[2, Any] | pl.DataFrame]:
        return self._varm
    
    @varm.setter
    def varm(self, varm: dict[str, np.ndarray[2, Any] | pl.DataFrame]) -> None:
        num_cells = len(self)
        for key, value in varm.items():
            if not isinstance(key, str):
                error_message = (
                    f'all keys of varm must be strings, but new varm '
                    f'contains a key of type {type(key).__name__!r}')
                raise TypeError(error_message)
            if isinstance(value, np.ndarray):
                if value.ndim != 2:
                    error_message = (
                        f'all values of varm must be 2D NumPy arrays or '
                        f'polars DataFrames, but new varm[{key!r}] is a '
                        f'{value.ndim:,}D NumPy array')
                    raise ValueError(error_message)
            elif not isinstance(value, pl.DataFrame):
                error_message = (
                    f'all values of varm must be NumPy arrays or polars '
                    f'DataFrames, but new varm[{key!r}] has type '
                    f'{type(value).__name__!r}')
                raise TypeError(error_message)
            if len(varm) != num_cells:
                error_message = (
                    f'the length of new varm[{key!r}] is {len(value):,}, but '
                    f'X.shape[0] is {self._X.shape[0]:,}')
                raise ValueError(error_message)
        self._varm = varm.copy()
    
    @property
    def uns(self) -> NestedScalarOrArrayDict:
        return self._uns
    
    @uns.setter
    def uns(self, uns: NestedScalarOrArrayDict) -> None:
        valid_uns_types = str, int, np.integer, float, np.floating, \
            bool, np.bool_, np.ndarray
        for description, value in SingleCell._iter_uns(uns):
            if not isinstance(value, valid_uns_types):
                error_message = (
                    f'all values of uns must be scalars (strings, numbers or '
                    f'Booleans) or NumPy arrays, or nested dictionaries '
                    f'thereof, but {description} has type '
                    f'{type(value).__name__!r}')
                raise TypeError(error_message)
        self._uns = uns.copy()
    
    @staticmethod
    def _iter_uns(uns: NestedScalarOrArrayDict, *, prefix: str = 'uns') -> \
            Iterable[tuple[str, str | int | float | np.integer | np.floating |
                                bool | np.bool_ | np.ndarray[Any, Any]]]:
        """
        Recurse through uns, yielding tuples of a string describing each key
        (e.g. "uns['a']['b']") and the corresponding value.
        
        Args:
            uns: an uns dictionary

        Yields:
            Length-2 tuples where the first element is a string describing each
            key, and the second element is the corresponding value.
        """
        for key, value in uns.items():
            key = f'{prefix}[{key!r}]'
            if isinstance(value, dict):
                SingleCell._iter_uns(value, prefix=key)
            else:
                yield key, value
    
    @staticmethod
    def _read_uns(uns_group: h5py.Group) -> NestedScalarOrArrayDict:
        """
        Recursively load uns from an .h5ad file.
        
        Args:
            uns_group: uns as an h5py.Group

        Returns:
            The loaded uns.
        """
        return {key: SingleCell._read_uns(value)
                     if isinstance(value, h5py.Group) else
                     (pl.Series(value[:]).cast(pl.String).to_numpy()
                      if value.shape else value[()].decode('utf-8'))
                     if value.dtype == object else
                     (value[:] if value.shape else value[()].item())
                for key, value in uns_group.items()}
    
    @staticmethod
    def _save_uns(uns: NestedScalarOrArrayDict,
                  uns_group: h5py.Group,
                  h5ad_file: h5py.File) -> None:
        """
        Recursively save uns to an .h5ad file.
        
        Args:
            uns: an uns dictionary
            uns_group: uns as an h5py.Group
            h5ad_file: an `h5py.File` open in write mode
        """
        uns_group.attrs['encoding-type'] = 'dict'
        uns_group.attrs['encoding-version'] = '0.1.0'
        for key, value in uns.items():
            if isinstance(value, dict):
                SingleCell._save_uns(value, uns_group.create_group(key),
                                     h5ad_file)
            else:
                dataset = uns_group.create_dataset(key, data=value)
                dataset.attrs['encoding-type'] = \
                    ('string-array' if value.dtype == object else 'array') \
                    if isinstance(value, np.ndarray) else \
                    'string' if isinstance(value, str) else 'numeric-scalar'
                dataset.attrs['encoding-version'] = '0.2.0'
    
    @property
    def obs_names(self) -> pl.Series:
        return self._obs[:, 0]
    
    @property
    def var_names(self) -> pl.Series:
        return self._var[:, 0]
    
    def set_obs_names(self, column: str) -> SingleCell:
        """
        Sets a column as the new first column of obs, i.e. the obs_names.
        
        Args:
            column: the column name in obs; must have String, Categorical, or
                    Enum data type

        Returns:
            A new SingleCell dataset with `column` as the first column of obs.
            If `column` is already the first column, return this dataset
            unchanged.
        """
        obs = self._obs
        check_type(column, 'column', str, 'a string')
        if column == obs.columns[0]:
            return self
        if column not in obs:
            error_message = f'{column!r} is not a column of obs'
            raise ValueError(error_message)
        check_dtype(obs[column], f'obs[{column!r}]',
                    (pl.String, pl.Categorical, pl.Enum))
        # noinspection PyTypeChecker
        return SingleCell(X=self._X,
                          obs=obs.select(column, pl.exclude(column)),
                          var=self._var, obsm=self._obsm, varm=self._varm,
                          uns=self._uns)
    
    def set_var_names(self, column: str) -> SingleCell:
        """
        Sets a column as the new first column of var, i.e. the var_names.
        
        Args:
            column: the column name in var; must have String, Categorical, or
                    Enum data type

        Returns:
            A new SingleCell dataset with `column` as the first column of var.
            If `column` is already the first column, return this dataset
            unchanged.
        """
        var = self._var
        check_type(column, 'column', str, 'a string')
        if column == var.columns[0]:
            return self
        if column not in var:
            error_message = f'{column!r} is not a column of var'
            raise ValueError(error_message)
        check_dtype(self._var[column], f'var[{column!r}]',
                    (pl.String, pl.Categorical, pl.Enum))
        # noinspection PyTypeChecker
        return SingleCell(X=self._X, obs=self._obs,
                          var=var.select(column, pl.exclude(column)),
                          obsm=self._obsm, varm=self._varm, uns=self._uns)
    
    @staticmethod
    def _read_datasets(datasets: Sequence[h5py.Dataset],
                       num_threads: int | np.integer) -> \
            dict[str, np.ndarray[1, Any]]:
        """
        Read a sequence of HDF5 datasets into a dictionary of 1D NumPy arrays.
        Assume all are from the same file (this is not checked).
        
        Args:
            datasets: a sequence of `h5py.Dataset` objects to read
            num_threads: the number of threads to use when reading; if >1,
                         spawn multiple processes and read into shared memory

        Returns:
            A dictionary of NumPy arrays with the contents of the datasets; the
            keys are taken from each dataset's `name` attribute.
        """
        if len(datasets) == 0:
            return {}
        import multiprocessing
        # noinspection PyUnresolvedReferences
        from multiprocessing.sharedctypes import _new_value
        import resource
        
        # Increase the maximum number of possible file descriptors this process
        # can use, since dataframes with hundreds of columns can easily exhaust
        # the common default limit of 1024 file descriptors (`ulimit -n`)
        
        soft_limit, hard_limit = resource.getrlimit(resource.RLIMIT_NOFILE)
        if soft_limit < hard_limit:
            resource.setrlimit(resource.RLIMIT_NOFILE,
                               (hard_limit, hard_limit))
        
        # Allocate shared memory for each dataset. Use _new_value() instead of
        # multiprocessing.Array() to avoid the memset at github.com/python/
        # cpython/blob/main/Lib/multiprocessing/sharedctypes.py#L62.
        
        # noinspection PyTypeChecker
        buffers = {dataset.name: _new_value(
            len(dataset) * np.ctypeslib.as_ctypes_type(dataset.dtype))
            for dataset in datasets}
        
        # Spawn num_threads processes: the first loads the first len(dataset) /
        # num_threads elements of each dataset, the second loads the next
        # len(dataset) / num_threads, etc. Because the chunks loaded by each
        # process are non-overlapping, there's no need to lock.
        
        filename = datasets[0].file.filename
        
        def read_dataset_chunks(thread_index: int) -> None:
            try:
                with h5py.File(filename) as h5ad_file:
                    for dataset_name, buffer in buffers.items():
                        chunk_size = len(buffer) // num_threads
                        start = thread_index * chunk_size
                        end = len(buffer) \
                            if thread_index == num_threads - 1 else \
                            start + chunk_size
                        chunk = np.s_[start:end]
                        dataset = h5ad_file[dataset_name]
                        dataset.read_direct(np.frombuffer(
                            buffer, dtype=dataset.dtype),
                            source_sel=chunk, dest_sel=chunk)
            except KeyboardInterrupt:
                pass  # do not print KeyboardInterrupt tracebacks, just return
        
        processes = []
        for thread_index in range(num_threads):
            process = multiprocessing.Process(
                target=read_dataset_chunks, args=(thread_index,))
            processes.append(process)
            process.start()
        for process in processes:
            process.join()
        # Wrap the shared memory in NumPy arrays
        arrays = {
            dataset_name: np.frombuffer(buffer, dtype=dataset.dtype)
            for (dataset_name, buffer), dataset in
            zip(buffers.items(), datasets)}
        return arrays
    
    @staticmethod
    def _read_dataset(dataset: h5py.Dataset,
                      num_threads: int | np.integer,
                      preloaded_datasets: dict[str, h5py.Dataset] |
                                          None = None) -> np.ndarray[1, Any]:
        """
        Read an HDF5 dataset into a 1D NumPy array.
        
        Args:
            dataset: the `h5py.Dataset` to read
            num_threads: the number of threads to use for reading; if >1, spawn
                         multiple processes and read into shared memory, unless
                         - the array is small enough (heuristically, under 10k
                           elements) that it's probably faster to read serially
                         - the array is so small that there's less than one
                           element to read per thread
                         - the array has `dtype=object` (not compatible with
                           shared memory)
            preloaded_datasets: a dictionary of preloaded datasets, or None to
                                always load from scratch. If specified and
                                `dataset.name` is in `preloaded_datasets`,
                                just return `preloaded_datasets[dataset.name]`.
        
        Returns:
            A 1D NumPy array with the contents of the dataset.
        """
        if preloaded_datasets is not None and \
                dataset.name in preloaded_datasets:
            return preloaded_datasets[dataset.name]
        dtype = dataset.dtype
        min_size = max(10_000, num_threads)
        if num_threads == 1 or dtype == object or dataset.size < min_size:
            return dataset[:]
        else:
            return SingleCell._read_datasets([dataset], num_threads)\
                [dataset.name]
    
    @staticmethod
    def _preload_datasets(group: h5py.Group,
                          num_threads: int | np.integer = 1) -> \
            dict[str, np.ndarray[1, Any]]:
        """
        Given a group from an .h5ad file, preload all datasets inside it,
        except for those where:
        - the array is small enough (heuristically, under 10k elements) that
          it's probably faster to read serially
        - the array is so small that there's less than one element to read per
          thread the array has `dtype=object` (not compatible with shared
          memory)
        
        Args:
            group: an `h5py.Group` to preload
            num_threads: the number of threads to use when preloading; if
                         `num_threads == 1`, do not preload

        Returns:
            A (possibly empty) dictionary of preloaded datasets.
        """
        if num_threads == 1:
            return {}
        datasets = []
        min_size = max(10_000, num_threads)
        group.visititems(
            lambda name, node: datasets.append(node)
            if isinstance(node, h5py.Dataset) and node.dtype != object
            and node.size >= min_size else None)
        return SingleCell._read_datasets(datasets, num_threads)
    
    @staticmethod
    def _read_h5ad_dataframe(h5ad_file: h5py.File,
                             key: str,
                             columns: str | Sequence[str] | None = None,
                             num_threads: int | np.integer = 1) -> \
            pl.DataFrame:
        """
        Load obs or var from an .h5ad file as a polars DataFrame.
        
        Args:
            h5ad_file: an `h5py.File` open in read mode
            key: the key to load as a DataFrame, e.g. `'obs'` or `'var'`
            columns: the column(s) of the DataFrame to load; the index column
                     is always loaded as the first column, regardless of
                     whether it is specified here, and then the remaining
                     columns are loaded in the order specified
            num_threads: the number of threads to use when reading
        
        Returns:
            A polars DataFrame of the data in h5ad_file[key].
        """
        group = h5ad_file[key]
        # Special case: the entire obs or var may rarely be a single NumPy
        # structured array (dtype=void)
        if isinstance(group, h5py.Dataset) and \
                np.issubdtype(group.dtype, np.void):
            data = pl.from_numpy(group[:])
            data = data.with_columns(pl.col(pl.Binary).cast(pl.String))
            return data
        preloaded_datasets = SingleCell._preload_datasets(group, num_threads)
        data = {}
        if columns is None:
            columns = group.attrs['column-order']
        else:
            columns = [column for column in to_tuple(columns)
                       if column != group.attrs['_index']]
            for column in columns:
                if column not in group.attrs['column-order']:
                    error_message = f'{column!r} is not a column of {key}'
                    raise ValueError(error_message)
        for column in chain((group.attrs['_index'],), columns):
            value = group[column]
            encoding_type = value.attrs.get('encoding-type')
            if encoding_type == 'categorical' or (
                    isinstance(value, h5py.Group) and all(
                    key == 'categories' or key == 'codes'
                    for key in value.keys())) or 'categories' in value.attrs:
                # Sometimes, the categories are stored in a different place
                # which is pointed to by value.attrs['categories']
                if 'categories' in value.attrs:
                    category_object = h5ad_file[value.attrs['categories']]
                    category_encoding_type = None
                    # noinspection PyTypeChecker
                    codes = SingleCell._read_dataset(
                        value, num_threads, preloaded_datasets)
                else:
                    category_object = value['categories']
                    category_encoding_type = \
                        category_object.attrs.get('encoding-type')
                    codes = SingleCell._read_dataset(
                        value['codes'], num_threads, preloaded_datasets)
                # Sometimes, the categories are themselves nullable
                # integer or Boolean arrays
                if category_encoding_type == 'nullable-integer' or \
                        category_encoding_type == 'nullable-boolean' or (
                        isinstance(category_object, h5py.Group) and all(
                        key == 'values' or key == 'mask'
                        for key in category_object.keys())):
                    data[column] = pl.Series(SingleCell._read_dataset(
                        category_object['values'], num_threads,
                        preloaded_datasets)[codes])
                    mask = pl.Series(SingleCell._read_dataset(
                        category_object['mask'], num_threads,
                        preloaded_datasets)[codes] | (codes == -1))
                    has_missing = mask.any()
                    if has_missing:
                        data[column] = data[column].set(mask, None)
                    continue
                # noinspection PyTypeChecker
                categories = SingleCell._read_dataset(
                    category_object, num_threads, preloaded_datasets)
                mask = pl.Series(codes == -1)
                has_missing = mask.any()
                # polars does not (as of version 0.20.2) support Categoricals
                # or Enums with non-string categories, so if the categories are
                # not strings, just map the codes to the categories.
                if category_encoding_type == 'array' or (
                        isinstance(category_object, h5py.Dataset) and
                        category_object.dtype != object):
                    data[column] = pl.Series(categories[codes],
                                             nan_to_null=True)
                    if has_missing:
                        data[column] = data[column].set(mask, None)
                elif category_encoding_type == 'string-array' or (
                        isinstance(category_object, h5py.Dataset) and
                        category_object.dtype == object):
                    if has_missing:
                        codes[mask] = 0
                    data[column] = pl.Series(codes, dtype=pl.UInt32)
                    if has_missing:
                        data[column] = data[column].set(mask, None)
                    # noinspection PyUnresolvedReferences
                    data[column] = data[column].cast(
                        pl.Enum(pl.Series(categories).cast(pl.String)))
                else:
                    encoding = \
                        f'encoding-type {category_encoding_type!r}' \
                        if category_encoding_type is not None else \
                            'encoding'
                    error_message = (
                        f'{column!r} column of {key!r} is a categorical '
                        f'with unsupported {encoding}')
                    raise ValueError(error_message)
            elif encoding_type == 'nullable-integer' or \
                    encoding_type == 'nullable-boolean' or (
                    isinstance(value, h5py.Group) and all(
                    key == 'values' or key == 'mask' for key in value.keys())):
                values = SingleCell._read_dataset(
                    value['values'], num_threads, preloaded_datasets)
                mask = SingleCell._read_dataset(
                    value['mask'], num_threads, preloaded_datasets)
                data[column] = pl.Series(values).set(pl.Series(mask), None)
            elif encoding_type == 'array' or (
                    isinstance(value, h5py.Dataset) and value.dtype != object):
                data[column] = pl.Series(SingleCell._read_dataset(
                    value, num_threads, preloaded_datasets), nan_to_null=True)
            elif encoding_type == 'string-array' or (
                    isinstance(value, h5py.Dataset) and value.dtype == object):
                data[column] = SingleCell._read_dataset(
                    value, num_threads, preloaded_datasets)
            else:
                encoding = f'encoding-type {encoding_type!r}' \
                    if encoding_type is not None else 'encoding'
                error_message = \
                    f'{column!r} column of {key!r} has unsupported {encoding}'
                raise ValueError(error_message)
        # NumPy doesn't support encoding object-dtyped string arrays as UTF-8,
        # so do the conversion in polars instead
        data = pl.DataFrame(data)\
            .with_columns(pl.col(pl.Binary).cast(pl.String))
        return data
    
    @staticmethod
    def read_obs(h5ad_file: h5py.File | str | Path,
                 columns: str | Iterable[str] | None = None,
                 num_threads: int | np.integer | None = 1) -> pl.DataFrame:
        """
        Load just obs from an .h5ad file as a polars DataFrame.
        
        Args:
            h5ad_file: an .h5ad filename
            columns: the column(s) of obs to load; if None, load all columns
            num_threads: the number of threads to use when reading; set
                         `num_threads=None` to use all available cores (as
                         determined by `os.cpu_count()`)
    
        Returns:
            A polars DataFrame of the data in obs.
        """
        check_type(h5ad_file, 'h5ad_file', (str, Path),
                   'a string or pathlib.Path')
        h5ad_file = str(h5ad_file)
        if not h5ad_file.endswith('.h5ad'):
            error_message = f".h5ad file {h5ad_file!r} must end with '.h5ad'"
            raise ValueError(error_message)
        filename = os.path.expanduser(h5ad_file)
        if not os.path.exists(filename):
            error_message = f'.h5ad file {h5ad_file!r} does not exist'
            raise FileNotFoundError(error_message)
        if columns is not None:
            columns = to_tuple(columns)
            if len(columns) == 0:
                error_message = 'no columns were specified'
                raise ValueError(error_message)
            check_types(columns, 'columns', str, 'strings')
        if num_threads is None:
            num_threads = os.cpu_count()
        else:
            check_type(num_threads, 'num_threads', int, 'a positive integer')
            check_bounds(num_threads, 'num_threads', 1)
        with h5py.File(filename) as f:
            return SingleCell._read_h5ad_dataframe(
                f, 'obs', columns=columns, num_threads=num_threads)
    
    @staticmethod
    def read_var(h5ad_file: str | Path,
                 columns: str | Iterable[str] | None = None,
                 num_threads: int | np.integer | None = 1) -> pl.DataFrame:
        """
        Load just var from an .h5ad file as a polars DataFrame.
        
        Args:
            h5ad_file: an .h5ad filename
            columns: the column(s) of var to load; if None, load all columns
            num_threads: the number of threads to use when reading; set
                         `num_threads=None` to use all available cores (as
                         determined by `os.cpu_count()`)
    
        Returns:
            A polars DataFrame of the data in var.
        """
        check_type(h5ad_file, 'h5ad_file', (str, Path),
                   'a string or pathlib.Path')
        h5ad_file = str(h5ad_file)
        if not h5ad_file.endswith('.h5ad'):
            error_message = f".h5ad file {h5ad_file!r} must end with '.h5ad'"
            raise ValueError(error_message)
        filename = os.path.expanduser(h5ad_file)
        if not os.path.exists(filename):
            error_message = f'.h5ad file {h5ad_file!r} does not exist'
            raise FileNotFoundError(error_message)
        if columns is not None:
            columns = to_tuple(columns)
            if len(columns) == 0:
                error_message = 'no columns were specified'
                raise ValueError(error_message)
            check_types(columns, 'columns', str, 'strings')
        if num_threads is None:
            num_threads = os.cpu_count()
        else:
            check_type(num_threads, 'num_threads', int, 'a positive integer')
            check_bounds(num_threads, 'num_threads', 1)
        with h5py.File(filename) as f:
            return SingleCell._read_h5ad_dataframe(
                f, 'var', columns=columns, num_threads=num_threads)
    
    @staticmethod
    def read_obsm(h5ad_file: str | Path,
                  keys: str | Iterable[str] | None = None,
                  num_threads: int | np.integer | None = 1) -> \
            dict[str, np.ndarray[2, Any] | pl.DataFrame]:
        """
        Load just obsm from an .h5ad file as a polars DataFrame.
        
        Args:
            h5ad_file: an .h5ad filename
            keys: the keys(s) of obsm to load; if None, load all keys
            num_threads: the number of threads to use when reading DataFrame
                         keys of obsm; set `num_threads=None` to use all
                         available cores (as determined by `os.cpu_count()`)
        
        Returns:
            A dictionary of NumPy arrays and polars DataFrames of the data in
            obsm.
        """
        check_type(h5ad_file, 'h5ad_file', (str, Path),
                   'a string or pathlib.Path')
        h5ad_file = str(h5ad_file)
        if not h5ad_file.endswith('.h5ad'):
            error_message = f".h5ad file {h5ad_file!r} must end with '.h5ad'"
            raise ValueError(error_message)
        filename = os.path.expanduser(h5ad_file)
        if not os.path.exists(filename):
            error_message = f'.h5ad file {h5ad_file!r} does not exist'
            raise FileNotFoundError(error_message)
        if keys is not None:
            keys = to_tuple(keys)
            if len(keys) == 0:
                error_message = 'no keys were specified'
                raise ValueError(error_message)
            check_types(keys, 'keys', str, 'strings')
        with h5py.File(filename) as f:
            if 'obsm' in f:
                obsm = f['obsm']
                if keys is None:
                    return {key: value[:]
                                 if isinstance(value, h5py.Dataset) else
                                 SingleCell._read_h5ad_dataframe(
                                     f, f'obsm/{key}', num_threads=num_threads)
                            for key, value in obsm.items()}
                else:
                    for key_index, key in enumerate(keys):
                        if key not in obsm:
                            error_message = (
                                f'keys[{key_index}] is {key!r}, which is not '
                                f'a key of obsm')
                            raise ValueError(error_message)
                    return {key: obsm[key][:]
                                 if isinstance(obsm[key], h5py.Dataset) else
                                 SingleCell._read_h5ad_dataframe(
                                    f, f'obsm/{key}', num_threads=num_threads)
                            for key in keys}
            else:
                if keys is not None:
                    error_message = 'keys was specified, but obsm is empty'
                    raise ValueError(error_message)
                return {}
         
    @staticmethod
    def read_varm(h5ad_file: str | Path,
                  keys: str | Iterable[str] | None = None,
                  num_threads: int | np.integer | None = 1) -> \
            dict[str, np.ndarray[2, Any] | pl.DataFrame]:
        """
        Load just varm from an .h5ad file as a polars DataFrame.
        
        Args:
            h5ad_file: an .h5ad filename
            keys: the keys(s) of varm to load; if None, load all keys
            num_threads: the number of threads to use when reading DataFrame
                         keys of varm; set `num_threads=None` to use all
                         available cores (as determined by `os.cpu_count()`)
        
        Returns:
            A dictionary of NumPy arrays and polars DataFrames of the data in
            varm.
        """
        check_type(h5ad_file, 'h5ad_file', (str, Path),
                   'a string or pathlib.Path')
        h5ad_file = str(h5ad_file)
        if not h5ad_file.endswith('.h5ad'):
            error_message = f".h5ad file {h5ad_file!r} must end with '.h5ad'"
            raise ValueError(error_message)
        filename = os.path.expanduser(h5ad_file)
        if not os.path.exists(filename):
            error_message = f'.h5ad file {h5ad_file!r} does not exist'
            raise FileNotFoundError(error_message)
        if keys is not None:
            keys = to_tuple(keys)
            if len(keys) == 0:
                error_message = 'no keys were specified'
                raise ValueError(error_message)
            check_types(keys, 'keys', str, 'strings')
        with h5py.File(filename) as f:
            if 'varm' in f:
                varm = f['varm']
                if keys is None:
                    return {key: value[:]
                                 if isinstance(value, h5py.Dataset) else
                                 SingleCell._read_h5ad_dataframe(
                                     f, f'varm/{key}', num_threads=num_threads)
                            for key, value in varm.items()}
                else:
                    for key_index, key in enumerate(keys):
                        if key not in varm:
                            error_message = (
                                f'keys[{key_index}] is {key!r}, which is not '
                                f'a key of varm')
                            raise ValueError(error_message)
                    return {key: varm[key][:]
                                 if isinstance(varm[key], h5py.Dataset) else
                                 SingleCell._read_h5ad_dataframe(
                                    f, f'varm/{key}', num_threads=num_threads)
                            for key in keys}
            else:
                if keys is not None:
                    error_message = 'keys was specified, but varm is empty'
                    raise ValueError(error_message)
                return {}
    
    @staticmethod
    def read_uns(h5ad_file: str | Path) -> NestedScalarOrArrayDict:
        """
        Load just uns from an .h5ad file as a dictionary.
        
        Args:
            h5ad_file: an .h5ad filename
        
        Returns:
            A dictionary of the data in uns.
        """
        check_type(h5ad_file, 'h5ad_file', (str, Path),
                   'a string or pathlib.Path')
        h5ad_file = str(h5ad_file)
        if not h5ad_file.endswith('.h5ad'):
            error_message = f".h5ad file {h5ad_file!r} must end with '.h5ad'"
            raise ValueError(error_message)
        filename = os.path.expanduser(h5ad_file)
        if not os.path.exists(filename):
            error_message = f'.h5ad file {h5ad_file!r} does not exist'
            raise FileNotFoundError(error_message)
        with h5py.File(filename) as f:
            if 'uns' in f:
                return SingleCell._read_uns(f['uns'])
            else:
                return {}
    
    @staticmethod
    def _print_matrix_info(X: h5py.Group | h5py.Dataset, X_name: str) -> None:
        """
        Given a key of an .h5ad file representing a sparse or dense matrix,
        print its shape, data type and (if sparse) numbr of non-zero elements.
        
        Args:
            X: the key in the .h5ad file representing the matrix, as a Group or
               Dataset object
            X_name: the name of the key
        """
        is_sparse = isinstance(X, h5py.Group)
        if is_sparse:
            data = X['data']
            shape = X.attrs['shape'] if 'shape' in X.attrs else \
                X.attrs['h5sparse_shape']
            dtype = str(data.dtype)
            nnz = data.shape[0]
            print(f'{X_name}: {shape[0]:,} × {shape[1]:,} sparse matrix with '
                  f'{nnz:,} non-zero elements, data type {dtype!r}, and '
                  f'first non-zero element = {data[0]:.6g}')
        else:
            shape = X.shape
            dtype = str(X.dtype)
            print(f'{X_name}: {shape[0]:,} × {shape[1]:,} dense matrix with '
                  f'data type {dtype!r} and first non-zero element = '
                  f'{X[0, 0]:.6g}')
    
    @staticmethod
    def ls(h5ad_file: str | Path) -> None:
        """
        Print the fields in an .h5ad file. This can be useful e.g. when
        deciding which count matrix to load via the `X_key` argument to
        `SingleCell()`.
        
        Args:
            h5ad_file: an .h5ad filename
        """
        check_type(h5ad_file, 'h5ad_file', (str, Path),
                   'a string or pathlib.Path')
        h5ad_file = str(h5ad_file)
        if not h5ad_file.endswith('.h5ad'):
            error_message = f".h5ad file {h5ad_file!r} must end with '.h5ad'"
            raise ValueError(error_message)
        filename = os.path.expanduser(h5ad_file)
        if not os.path.exists(filename):
            error_message = f'.h5ad file {h5ad_file!r} does not exist'
            raise FileNotFoundError(error_message)
        terminal_width = os.get_terminal_size().columns
        attrs = 'obs', 'var', 'obsm', 'varm', 'obsp', 'varp', 'uns'
        with h5py.File(filename) as f:
            # X
            SingleCell._print_matrix_info(f['X'], 'X')
            # layers
            if 'layers' in f:
                layers = f['layers']
                if len(layers) > 0:
                    for layer_name, layer in layers.items():
                        SingleCell._print_matrix_info(
                            layer, f'layers[{layer_name!r}]')
            # obs, var, obsm, varm, obsp, varp, uns
            for attr in attrs:
                if attr in f:
                    entries = f[attr]
                    if (attr == 'obs' or attr == 'var') and \
                            isinstance(entries, h5py.Dataset) and \
                            np.issubdtype(entries.dtype, np.void):
                        entries = entries.dtype.fields
                    if len(entries) > 0:
                        print(fill(f'{attr}: {", ".join(entries)}',
                                   width=terminal_width,
                                   subsequent_indent=' ' * (len(attr) + 2)))
            # raw
            if 'raw' in f:
                raw = f['raw']
                if len(raw) > 0:
                    print('raw:')
                    if 'X' in raw:
                        SingleCell._print_matrix_info(raw['X'], '    X')
                    if 'layers' in raw:
                        layers = raw['layers']
                        if len(layers) > 0:
                            for layer_name, layer in layers.items():
                                SingleCell._print_matrix_info(
                                    layer, f'    layers[{layer_name!r}]')
                    for attr in attrs:
                        if attr in raw:
                            entries = raw[attr]
                            if (attr == 'obs' or attr == 'var') and \
                                    isinstance(entries, h5py.Dataset) and \
                                    np.issubdtype(entries.dtype, np.void):
                                entries = entries.dtype.fields
                            if len(entries) > 0:
                                print(fill(f'    {attr}: {", ".join(entries)}',
                                           width=terminal_width,
                                           subsequent_indent=' ' * (
                                                   len(attr) + 6)))
    
    def __eq__(self, other: SingleCell) -> bool:
        """
        Test for equality with another SingleCell dataset.
        
        Args:
            other: the other SingleCell dataset to test for equality with

        Returns:
            Whether the two SingleCell datasets are identical.
        """
        if not isinstance(other, SingleCell):
            error_message = (
                f'the left-hand operand of `==` is a SingleCell dataset, but '
                f'the right-hand operand has type {type(other).__name__!r}')
            raise TypeError(error_message)
        # noinspection PyUnresolvedReferences
        return self._obs.equals(other._obs) and \
               self._var.equals(other._var) and \
               self._obsm.keys() == other._obsm.keys() and \
               self._varm.keys() == other._varm.keys() and \
            all(type(other._obsm[key]) is type(value) and
                (np.array_equal(other._obsm[key], value, equal_nan=True)
                 if isinstance(value, np.ndarray) else
                 other._obsm[key].equals(value))
                for key, value in self._obsm.items()) and \
            all(type(other._varm[key]) is type(value) and
                (np.array_equal(other._varm[key], value, equal_nan=True)
                 if isinstance(value, np.ndarray) else
                 other._varm[key].equals(value))
                for key, value in self._varm.items()) and \
            SingleCell._eq_uns(self._uns, other._uns) and \
            self._X.nnz == other._X.nnz and not (self._X != other._X).nnz
    
    @staticmethod
    def _eq_uns(uns: NestedScalarOrArrayDict,
                other_uns: NestedScalarOrArrayDict) -> bool:
        """
        Test whether two uns are equal.
        
        Args:
            uns: an uns
            other_uns: another uns

        Returns:
            Whether the two uns are equal.
        """
        return uns.keys() == other_uns.keys() and all(
            isinstance(value, dict) and isinstance(other_value, dict) and
            SingleCell._eq_uns(value, other_value) or
            isinstance(value, np.ndarray) and
            isinstance(other_value, np.ndarray) and
            np.array_equal(value, other_value, equal_nan=True) or
            not isinstance(other_value, (dict, np.ndarray)) and
            value == other_value
            for (key, value), (other_key, other_value) in
            zip(uns.items(), other_uns.items()))
    
    @staticmethod
    def _getitem_error(item: Indexer) -> None:
        """
        Raise an error if the indexer is invalid.
        
        Args:
            item: the indexer
        """
        types = tuple(type(elem).__name__ for elem in to_tuple(item))
        if len(types) == 1:
            types = types[0]
        error_message = (
            f'SingleCell indices must be cells, a length-1 tuple of (cells,), '
            f'or a length-2 tuple of (cells, genes). Cells and genes must '
            f'each be a string or integer; a slice of strings or integers; or '
            f'a list, NumPy array, or polars Series of strings, integers, or '
            f'Booleans. You indexed with: {types}.')
        raise ValueError(error_message)
    
    @staticmethod
    def _getitem_by_string(df: pl.DataFrame, string: str) -> int:
        """
        Get the index where df[:, 0] == string, raising an error if no rows or
        multiple rows match.
        
        Args:
            df: a DataFrame (obs or var)
            string: the string to find the index of in the first column of df

        Returns:
            The integer index of the string within the first column of df.
        """
        first_column = df.columns[0]
        try:
            return df\
                .select(pl.int_range(pl.len(), dtype=pl.Int32)
                        .alias('_SingleCell_getitem'), first_column)\
                .row(by_predicate=pl.col(first_column) == string)\
                [0]
        except pl.exceptions.NoRowsReturnedError:
            raise KeyError(string)
    
    @staticmethod
    def _getitem_process(item: Indexer, index: int, df: pl.DataFrame) -> \
            list[int] | slice | pl.Series:
        """
        Process an element of an item passed to __getitem__().
        
        Args:
            item: the item
            index: the index of the element to process
            df: the DataFrame (obs or var) to process the element with respect
                to

        Returns:
            A new indexer indicating the rows/columns to index.
        """
        subitem = item[index]
        if is_integer(subitem):
            return [subitem]
        elif isinstance(subitem, str):
            return [SingleCell._getitem_by_string(df, subitem)]
        elif isinstance(subitem, slice):
            start = subitem.start
            stop = subitem.stop
            step = subitem.step
            if isinstance(start, str):
                start = SingleCell._getitem_by_string(df, start)
            elif start is not None and not is_integer(start):
                SingleCell._getitem_error(item)
            if isinstance(stop, str):
                stop = SingleCell._getitem_by_string(df, stop)
            elif stop is not None and not is_integer(stop):
                SingleCell._getitem_error(item)
            if step is not None and not is_integer(step):
                SingleCell._getitem_error(item)
            return slice(start, stop, step)
        elif isinstance(subitem, (list, np.ndarray, pl.Series)):
            if not isinstance(subitem, pl.Series):
                subitem = pl.Series(subitem)
            if subitem.is_null().any():
                error_message = 'your indexer contains missing values'
                raise ValueError(error_message)
            if subitem.dtype == pl.String or subitem.dtype == \
                    pl.Categorical or subitem.dtype == pl.Enum:
                indices = subitem\
                    .to_frame(df.columns[0])\
                    .join(df.with_columns(_SingleCell_index=pl.int_range(
                              pl.len(), dtype=pl.Int32)),
                          on=df.columns[0], how='left')\
                    ['_SingleCell_index']
                if indices.null_count():
                    error_message = subitem.filter(indices.is_null())[0]
                    raise KeyError(error_message)
                return indices
            elif subitem.dtype.is_integer() or subitem.dtype == pl.Boolean:
                return subitem
            else:
                SingleCell._getitem_error(item)
        else:
            SingleCell._getitem_error(item)
            
    def __getitem__(self, item: Indexer | tuple[Indexer, Indexer]) -> \
            SingleCell:
        """
        Subset to specific cell(s) and/or gene(s).
        
        Index with a tuple of `(cells, genes)`. If `cells` and `genes` are
        integers, arrays/lists/slices of integers, or arrays/lists of Booleans,
        the result will be a SingleCell dataset subset to `X[cells, genes]`,
        `obs[cells]`, `var[genes]`, `obsm[cells]`, and `varm[genes]`. However,
        `cells` and/or `genes` can instead be strings (or arrays or slices of
        strings), in which case they refer to the first column of obs
        (`obs_names`) and/or var (`var_names`), respectively.
        
        Examples:
        - Subset to one cell, for all genes:
          sc['CGAATTGGTGACAGGT-L8TX_210916_01_B05-1131590416']
          sc[2]
        - Subset to one gene, for all cells:
          sc[:, 'APOE']
          sc[:, 13196]
        - Subset to one cell and one gene:
          sc['CGAATTGGTGACAGGT-L8TX_210916_01_B05-1131590416', 'APOE']
          sc[2, 13196]
        - Subset to a range of cells and genes:
          sc['CGAATTGGTGACAGGT-L8TX_210916_01_B05-1131590416':
             'CCCTCTCAGCAGCCTC-L8TX_211007_01_A09-1135034522',
             'APOE':'TREM2']
          sc[2:6, 13196:34268]
        - Subset to specific cells and genes:
          sc[['CGAATTGGTGACAGGT-L8TX_210916_01_B05-1131590416',
              'CCCTCTCAGCAGCCTC-L8TX_211007_01_A09-1135034522']]
          sc[:, pl.Series(['APOE', 'TREM2'])]
          sc[['CGAATTGGTGACAGGT-L8TX_210916_01_B05-1131590416',
              'CCCTCTCAGCAGCCTC-L8TX_211007_01_A09-1135034522'],
              np.array(['APOE', 'TREM2'])]
        
        Args:
            item: the item to index with
        
        Returns:
            A new SingleCell dataset subset to the specified cells and/or
            genes.
        """
        if not isinstance(item, (int, str, slice, tuple, list,
                                 np.ndarray, pl.Series)):
            error_message = (
                f'SingleCell datasets must be indexed with an integer, '
                f'string, slice, tuple, list, NumPy array, or polars Series, '
                f'but you tried to index with an object of type '
                f'{type(item).__name__!r}')
            raise TypeError(error_message)
        if isinstance(item, tuple):
            if not 1 <= len(item) <= 2:
                self._getitem_error(item)
        else:
            item = item,
        rows = self._getitem_process(item, 0, self._obs)
        rows_are_Series = isinstance(rows, pl.Series)
        if rows_are_Series:
            boolean_Series = rows.dtype == pl.Boolean
            obs = self._obs.filter(rows) if boolean_Series else self._obs[rows]
        else:
            boolean_Series = False
            obs = self._obs[rows]
        if self._obsm:
            rows_NumPy = rows.to_numpy() if rows_are_Series else rows
            obsm = {key: (value.filter(rows) if boolean_Series else
                          value[rows]) if isinstance(value, pl.DataFrame) else
                         value[rows_NumPy]
                    for key, value in self._obsm.items()}
        else:
            obsm = {}
        if len(item) == 1:
            return SingleCell(X=self._X[rows], obs=obs, var=self._var,
                              obsm=obsm, varm=self._varm, uns=self._uns)
        cols = self._getitem_process(item, 1, self._var)
        cols_are_Series = isinstance(cols, pl.Series)
        if cols_are_Series:
            boolean_Series = cols.dtype == pl.Boolean
            var = self._var.filter(cols) if boolean_Series else self._var[cols]
        else:
            boolean_Series = False
            var = self._var[cols]
        if self._varm:
            cols_NumPy = cols.to_numpy() if cols_are_Series else cols
            varm = {key: (value.filter(cols) if boolean_Series else
                          value[cols]) if isinstance(value, pl.DataFrame) else
                         value[cols_NumPy]
                    for key, value in self._varm.items()}
        else:
            varm = {}
        X = self._X[rows, cols] \
            if isinstance(rows, slice) or isinstance(cols, slice) else \
            self._X[np.ix_(rows, cols)]
        return SingleCell(X=X, obs=obs, var=var, obsm=obsm, varm=varm,
                          uns=self._uns)
    
    def cell(self, cell: str) -> np.ndarray[1, Any]:
        """
        Get the row of X corresponding to a single cell, based on the cell's
        name in obs_names.
        
        Args:
            cell: the name of the cell in obs_names
        
        Returns:
            The corresponding row of X, as a dense 1D NumPy array with zeros
            included.
        """
        row_index = SingleCell._getitem_by_string(self._obs, cell)
        return self._X[[row_index]].toarray().squeeze()
    
    def gene(self, gene: str) -> np.ndarray[1, Any]:
        """
        Get the column of X corresponding to a single gene, based on the gene's
        name in var_names.
        
        Args:
            gene: the name of the gene in var_names
        
        Returns:
            The corresponding column of X, as a dense 1D NumPy array with zeros
            included.
        """
        column_index = SingleCell._getitem_by_string(self._var, gene)
        return self._X[:, [column_index]].toarray().squeeze()
    
    def __len__(self) -> int:
        """
        Get the number of cells in this SingleCell dataset.
        
        Returns:
            The number of cells.
        """
        return self._X.shape[0]
       
    def __repr__(self) -> str:
        """
        Get a string representation of this SingleCell dataset.
        
        Returns:
            A string summarizing the dataset.
        """
        descr = (
            f'SingleCell dataset with {len(self._obs):,} '
            f'{plural("cell", len(self._obs))} (obs), {len(self._var):,} '
            f'{plural("gene", len(self._var))} (var), and {self._X.nnz:,} '
            f'non-zero {"entries" if self._X.nnz != 1 else "entry"} (X)')
        terminal_width = os.get_terminal_size().columns
        for attr in 'obs', 'var', 'obsm', 'varm', 'uns':
            entries = getattr(self, attr).columns \
                if attr == 'obs' or attr == 'var' else getattr(self, attr)
            if len(entries) > 0:
                descr += '\n' + fill(
                    f'    {attr}: {", ".join(entries)}',
                    width=terminal_width,
                    subsequent_indent=' ' * (len(attr) + 6))
        return descr
    
    @property
    def shape(self) -> tuple[int, int]:
        """
        Get the shape of this SingleCell dataset.
        
        Returns:
            A length-2 tuple where the first element is the number of cells,
            and the second is the number of genes.
        """
        return self._X.shape
        
    @staticmethod
    def _write_h5ad_dataframe(h5ad_file: h5py.File,
                              df: pl.DataFrame,
                              key: str,
                              preserve_strings: bool) -> None:
        """
        Write obs or var to an .h5ad file.
        
        Args:
            h5ad_file: an `h5py.File` open in write mode
            df: the DataFrame to write, e.g. obs or var
            key: the name of the DataFrame, e.g. `'obs'` or `'var'`
            preserve_strings: if False, encode string columns with duplicate
                              values as Enums to save space; if True, preserve
                              these columns as string columns
        """
        # Create a group for the data frame and add top-level metadata
        group = h5ad_file.create_group(key)
        group.attrs['_index'] = df.columns[0]
        group.attrs['column-order'] = df.columns[1:]
        group.attrs['encoding-type'] = 'dataframe'
        group.attrs['encoding-version'] = '0.2.0'
        for column in df:
            dtype = column.dtype
            if dtype == pl.String:
                if column.null_count() or \
                        not preserve_strings and column.is_duplicated().any():
                    column = column\
                        .cast(pl.Enum(column.unique(maintain_order=True)
                                      .drop_nulls()))
                    dtype = column.dtype
                else:
                    dataset = group.create_dataset(column.name,
                                                   data=column.to_numpy())
                    dataset.attrs['encoding-type'] = 'string-array'
                    dataset.attrs['encoding-version'] = '0.2.0'
                    continue
            if dtype == pl.Enum or dtype == pl.Categorical:
                is_Enum = dtype == pl.Enum
                subgroup = group.create_group(column.name)
                subgroup.attrs['encoding-type'] = 'categorical'
                subgroup.attrs['encoding-version'] = '0.2.0'
                subgroup.attrs['ordered'] = is_Enum
                categories = column.cat.get_categories()
                if not is_Enum:
                    column = column.cast(pl.Enum(categories))
                codes = column.to_physical().fill_null(-1)
                subgroup.create_dataset('codes', data=codes.to_numpy())
                if len(categories) == 0:
                    subgroup.create_dataset('categories', shape=(0,),
                                            dtype=h5py.special_dtype(vlen=str))
                else:
                    subgroup.create_dataset('categories',
                                            data=categories.to_numpy())
            elif dtype.is_float():
                # Nullable floats are not supported, so convert null to NaN
                dataset = group.create_dataset(
                    column.name, data=column.fill_null(np.nan).to_numpy())
                dataset.attrs['encoding-type'] = 'array'
                dataset.attrs['encoding-version'] = '0.2.0'
            elif dtype == pl.Boolean or dtype.is_integer():
                is_boolean = dtype == pl.Boolean
                if column.null_count():
                    # Store as nullable integer/Boolean
                    subgroup = group.create_group(column.name)
                    subgroup.attrs['encoding-type'] = \
                        f'nullable-{"boolean" if is_boolean else "integer"}'
                    subgroup.attrs['encoding-version'] = '0.1.0'
                    subgroup.create_dataset(
                        'values',
                        data=column.fill_null(False if is_boolean else 1)
                        .to_numpy())
                    subgroup.create_dataset(
                        'mask', data=column.is_null().to_numpy())
                else:
                    # Store as regular integer/Boolean
                    dataset = group.create_dataset(column.name,
                                                   data=column.to_numpy())
                    dataset.attrs['encoding-type'] = 'array'
                    dataset.attrs['encoding-version'] = '0.2.0'
            else:
                error_message = \
                    f'internal error: unsupported data type {dtype!r}'
                raise TypeError(error_message)
    
    def save(self,
             filename: str | Path,
             *,
             overwrite: bool = False,
             preserve_strings: bool = False,
             sce: bool = False) -> None:
        """
        Save this SingleCell dataset to an AnnData .h5ad file, 10x .h5 or
        .mtx.gz file, or Seurat or SingleCellExperiment .rds file.
        
        Args:
            filename: an AnnData .h5ad file, 10x .h5 or .mtx.gz file, or Seurat
                      or SingleCellExperiment .rds file to save to. File format
                      will be inferred from the extension. If the extension is
                      .rds, the `sce` argument will determine whether to save
                      to a Seurat or SingleCellExperiment object.
                      - When saving to a 10x .h5 file, `obs['barcodes']`,
                        `var['feature_type']`, `var['genome']`, `var['id']`,
                        and `var['name']` must all exist. Only X and these
                        columns will be saved, along with whichever of
                        `var['pattern']`, `var['read']`, and `var['sequence']`
                        exist. All of these columns (if they exist) must be
                        String, Categorical or Enum.
                      - When saving to a 10x .mtx.gz file, barcodes.tsv.gz and
                        features.tsv.gz will be created in the same directory.
                        Only X, obs and var will be saved.
                      - When saving to a Seurat .rds file, to match the
                        requirements of Seurat objects, the `'X_'` prefix
                        (often used by Scanpy) will be removed from each key of
                        obsm where it is present (e.g. `'X_umap'` will become
                        `'umap'`).
                      - When swaving to a Seurat .rds file, Seurat will add
                       `'orig.ident'`, `'nCount_RNA'` and `'nFeature_RNA'` as
                        gene-level metadata by default; you can disable the
                        calculation of the latter two columns with:
                        
                        ```python
                        from ryp import r
                        r('options(Seurat.object.assay.calcn = FALSE)')
                        ```
                        
                      - When saving to a Seurat .rds file, non-string keys of
                        uns will not be saved.
                      - When saving to a Seurat or SingleCellExperiment .rds
                        file, varm will not be saved.
            overwrite: if False, raises an error if (any of) the file(s) exist;
                       if True, overwrites them
            preserve_strings: if False, encode string columns with duplicate
                              values as Enums to save space, when saving to
                              AnnData .h5ad or Seurat or SingleCellExperiment
                              .rds; if True, preserve these columns as string
                              columns. (Regardless of the value of
                              `preserve_strings`, String columns with `null`
                              values will be encoded as Enums when saving to
                              .h5ad, since the .h5ad format cannot represent
                              them otherwise.)
            sce: if True and the extension of filename is .rds, save to a
                 SingleCellExperiment object instead of a Seurat object
        """
        # Check inputs
        check_type(filename, 'filename', (str, Path),
                   'a string or pathlib.Path')
        filename = str(filename)
        is_anndata = filename.endswith('.h5ad')
        is_h5 = filename.endswith('.h5')
        is_mtx = filename.endswith('.mtx.gz')
        is_rds = filename.endswith('.rds')
        if not (is_anndata or is_h5 or is_mtx or is_rds):
            error_message = (
                f"filename {filename!r} does not end with '.h5ad', '.h5', "
                f"'.mtx.gz', or '.rds'")
            raise ValueError(error_message)
        filename_expanduser = os.path.expanduser(filename)
        if not overwrite and os.path.exists(filename_expanduser):
            error_message = (
                f'filename {filename!r} already exists; set overwrite=True '
                f'to overwrite')
            raise FileExistsError(error_message)
        check_type(sce, 'sce', bool, 'Boolean')
        if sce and not is_rds:
            error_message = 'sce can only be True when saving to an .rds file'
            raise ValueError(error_message)
        # Raise an error if obs, var, or DataFrame keys of obsm or varm contain
        # columns with unsupported data types (anything but float, int, String,
        # Categorical, Enum, Boolean)
        valid_polars_dtypes = pl.FLOAT_DTYPES | pl.INTEGER_DTYPES | \
                              {pl.String, pl.Categorical, pl.Enum, pl.Boolean}
        for df, df_name in (self._obs, 'obs'), (self._var, 'var'):
            for column, dtype in df.schema.items():
                if dtype.base_type() not in valid_polars_dtypes:
                    error_message = (
                        f'{df_name}[{column!r}] has the data type '
                        f'{dtype.base_type()!r}, which is not supported when '
                        f'saving')
                    raise TypeError(error_message)
        for field, field_name in (self._obsm, 'obsm'), (self._varm, 'varm'):
            for key, value in field.items():
                if not isinstance(value, pl.DataFrame):
                    continue
                for column, dtype in value.schema.items():
                    if dtype.base_type() not in valid_polars_dtypes:
                        error_message = (
                            f'{field}[{key!r}][{column!r}] has the data type '
                            f'{dtype.base_type()!r}, which is not supported '
                            f'when saving')
                        raise TypeError(error_message)
        # Raise an error if obsm, varm or uns contain NumPy arrays with
        # unsupported data types (datetime64, timedelta64, unstructured void).
        # Do not specifically check `dtype=object` to avoid extra overhead.
        for field, field_name in \
                (self._obsm, 'obsm'), (self._varm, 'varm'), (self._uns, 'uns'):
            for key, value in field.items():
                if not isinstance(value, np.ndarray):
                    continue
                if value.dtype.type == np.void and value.dtype.names is None:
                    error_message = (
                        f'{field_name}[{key!r}] is an unstructured void '
                        f'array, which is not supported when saving')
                    raise TypeError(error_message)
                elif value.dtype == np.datetime64:
                    error_message = (
                        f'{field_name}[{key!r}] is a datetime64 array, which '
                        f'is not supported when saving')
                    raise TypeError(error_message)
                elif value.dtype == np.timedelta64:
                    error_message = (
                        f'{field_name}[{key!r}] is a timedelta64 array, which '
                        f'is not supported when saving')
                    raise TypeError(error_message)
        # Save, depending on the file extension
        if is_anndata:
            try:
                with h5py.File(filename_expanduser, 'w') as h5ad_file:
                    # Add top-level metadata
                    h5ad_file.attrs['encoding-type'] = 'anndata'
                    h5ad_file.attrs['encoding-version'] = '0.1.0'
                    # Save obs and var
                    SingleCell._write_h5ad_dataframe(
                        h5ad_file, self._obs, 'obs', preserve_strings)
                    SingleCell._write_h5ad_dataframe(
                        h5ad_file, self._var, 'var', preserve_strings)
                    # Save obsm
                    if self._obsm:
                        obsm = h5ad_file.create_group('obsm')
                        obsm.attrs['encoding-type'] = 'dict'
                        obsm.attrs['encoding-version'] = '0.1.0'
                        for key, value in self._obsm.items():
                            if isinstance(value, pl.DataFrame):
                                SingleCell._write_h5ad_dataframe(
                                    h5ad_file, value, f'obsm/{key}',
                                    preserve_strings)
                            else:
                                obsm.create_dataset(key, data=value)
                    # Save varm
                    if self._varm:
                        varm = h5ad_file.create_group('varm')
                        varm.attrs['encoding-type'] = 'dict'
                        varm.attrs['encoding-version'] = '0.1.0'
                        for key, value in self._varm.items():
                            if isinstance(value, pl.DataFrame):
                                SingleCell._write_h5ad_dataframe(
                                    h5ad_file, value, f'varm/{key}',
                                    preserve_strings)
                            else:
                                varm.create_dataset(key, data=value)
                    # Save uns
                    if self._uns:
                        SingleCell._save_uns(self._uns,
                                             h5ad_file.create_group('uns'),
                                             h5ad_file)
                    # Save X
                    X = h5ad_file.create_group('X')
                    X.attrs['encoding-type'] = 'csr_matrix' \
                        if isinstance(self._X, csr_array) else 'csc_matrix'
                    X.attrs['encoding-version'] = '0.1.0'
                    X.attrs['shape'] = self._X.shape
                    X.create_dataset('data', data=self._X.data)
                    X.create_dataset('indices', data=self._X.indices)
                    X.create_dataset('indptr', data=self._X.indptr)
            except:
                if os.path.exists(filename_expanduser):
                    os.unlink(filename_expanduser)
                raise
        elif is_h5:
            obs_columns = ['barcodes']
            var_columns = ['feature_type', 'genome', 'id', 'name']
            for columns, df, df_name in (obs_columns, self._obs, 'obs'), \
                    (var_columns, self._var, 'var'):
                for column in columns:
                    if column not in df:
                        error_message = (
                            f'{column!r} was not found in {df_name}, but is a '
                            f'required column when saving to a 10x .h5 file')
                        raise ValueError(error_message)
                    check_dtype(df[column], f'{df_name}[{column!r}]',
                                (pl.String, pl.Categorical, pl.Enum))
            all_tag_keys = ['genome']
            for column in 'pattern', 'read', 'sequence':
                if column in self._var:
                    check_dtype(self._var[column], f'var[{column!r}]',
                                (pl.String, pl.Categorical, pl.Enum))
                    var_columns.append(column)
                    all_tag_keys.append(column)
            try:
                with h5py.File(filename_expanduser, 'w') as h5_file:
                    matrix = h5_file.create_group('matrix')
                    matrix.create_dataset('barcodes',
                                          data=self._obs[:, 0].to_numpy())
                    matrix.create_dataset('data', data=self._X.data)
                    features = matrix.create_group('features')
                    matrix.create_dataset('indices', data=self._X.indices)
                    matrix.create_dataset('indptr', data=self._X.indptr)
                    matrix.create_dataset('shape', data=self._X.shape[::-1])
                    features.create_dataset('_all_tag_keys', data=all_tag_keys)
                    for column in var_columns:
                        features.create_dataset(
                            column, data=self._var[column].to_numpy())
            except:
                if os.path.exists(filename_expanduser):
                    os.unlink(filename_expanduser)
                raise
        elif is_mtx:
            from scipy.io import mmwrite
            barcode_filename = os.path.join(
                os.path.dirname(filename_expanduser), 'barcodes.tsv.gz')
            feature_filename = os.path.join(
                os.path.dirname(filename_expanduser), 'features.tsv.gz')
            if not overwrite:
                for ancillary_filename in barcode_filename, feature_filename:
                    if os.path.exists(ancillary_filename):
                        error_message = (
                            f'{ancillary_filename!r} already exists; set '
                            f'overwrite=True to overwrite')
                        raise FileExistsError(error_message)
            try:
                mmwrite(filename_expanduser, self._X.T)
                self._obs.write_csv(barcode_filename, include_header=False)
                self._var.write_csv(feature_filename, include_header=False)
            except:
                if os.path.exists(filename_expanduser):
                    os.unlink(filename_expanduser)
                if os.path.exists(barcode_filename):
                    os.unlink(barcode_filename)
                if os.path.exists(feature_filename):
                    os.unlink(feature_filename)
                raise
        else:
            from ryp import r
            if preserve_strings:
                sc = self
            else:
                # use .pipe(lambda df: df.filter(...) if len(df) > 0 else df)
                # instead of .filter(...) to work around
                # github.com/pola-rs/polars/issues/18170
                enumify = lambda df: df.cast({
                    row[0]: pl.Enum(row[1]) for row in df
                    .select(pl.selectors.string()
                            .unique(maintain_order=True)
                            .implode()
                            .list.drop_nulls())
                    .unpivot()
                    .pipe(lambda df: df.filter(pl.col.value.list.len() ==
                                               len(df))
                          if len(df) > 0 else df)
                    .rows()})
                sc = SingleCell(X=self._X, obs=enumify(self._obs),
                                var=enumify(self._var), obsm=self._obsm,
                                uns=self._uns)
            if sce:
                sc.to_sce('.SingleCell.object')
            else:
                sc.to_seurat('.SingleCell.object')
            try:
                r(f'saveRDS(.SingleCell.object, {filename_expanduser!r})')
            except:
                if os.path.exists(filename_expanduser):
                    os.unlink(filename_expanduser)
                raise
            finally:
                r('rm(.SingleCell.object)')
    
    def _get_column(self,
                    obs_or_var_name: Literal['obs', 'var'],
                    column: SingleCellColumn,
                    variable_name: str,
                    dtypes: pl.datatypes.classes.DataTypeClass | str |
                            tuple[pl.datatypes.classes.DataTypeClass | str,
                            ...],
                    *,
                    QC_column: pl.Series | None = None,
                    allow_missing: bool = False,
                    allow_null: bool = False,
                    custom_error: str | None = None) -> pl.Series | None:
        """
        Get a column of the same length as obs/var, or None if the column is
        missing from obs/var and `allow_missing=True`.
        
        Args:
            obs_or_var_name: the name of the DataFrame the column is with
                             respect to, i.e. `'obs'` or `'var'`
            column: a string naming a column of obs/var, a polars expression
                    that evaluates to a single column when applied to obs/var,
                    a polars Series or 1D NumPy array of the same length as
                    obs/var, or a function that takes in `self` and returns a
                    polars Series or 1D NumPy array of the same length as
                    obs/var
            variable_name: the name of the variable corresponding to `column`
            dtypes: the required dtype(s) of the column
            QC_column: an optional column of cells passing QC. If specified,
                       the presence of null values will only raise an error for
                       cells passing QC. Has no effect when `allow_null=True`.
            allow_missing: whether to allow `column` to be a string missing
                           from obs/var, returning None in this case
            allow_null: whether to allow `column` to contain null values
            custom_error: a custom error message for when `column` is a string
                          and is not found in obs/var, and
                          `allow_missing=False`; use `{}` as a placeholder for
                          the name of the column
        
        Returns:
            A polars Series of the same length as obs/var, or None if the
            column is missing from obs/var and `allow_missing=True`.
        """
        obs_or_var = self._obs if obs_or_var_name == 'obs' else self._var
        if isinstance(column, str):
            if column in obs_or_var:
                column = obs_or_var[column]
            elif allow_missing:
                return None
            else:
                error_message = \
                    f'{column!r} is not a column of {obs_or_var_name}' \
                    if custom_error is None else \
                    custom_error.format(f'{column!r}')
                raise ValueError(error_message)
        elif isinstance(column, pl.Expr):
            column = obs_or_var.select(column)
            if column.width > 1:
                error_message = (
                    f'{variable_name} is a polars expression that expands to '
                    f'{column.width:,} columns rather than 1')
                raise ValueError(error_message)
            column = column.to_series()
        elif isinstance(column, pl.Series):
            if len(column) != len(obs_or_var):
                error_message = (
                    f'{variable_name} is a polars Series of length '
                    f'{len(column):,}, which differs from the length of '
                    f'{obs_or_var_name} ({len(obs_or_var):,})')
                raise ValueError(error_message)
        elif isinstance(column, np.ndarray):
            if len(column) != len(obs_or_var):
                error_message = (
                    f'{variable_name} is a NumPy array of length '
                    f'{len(column):,}, which differs from the length of '
                    f'{obs_or_var_name} ({len(obs_or_var):,})')
                raise ValueError(error_message)
            column = pl.Series(variable_name, column)
        elif callable(column):
            column = column(self)
            if isinstance(column, np.ndarray):
                if column.ndim != 1:
                    error_message = (
                        f'{variable_name} is a function that returns a '
                        f'{column.ndim:,}D NumPy array, but must return a '
                        f'polars Series or 1D NumPy array')
                    raise ValueError(error_message)
                column = pl.Series(variable_name, column)
            elif not isinstance(column, pl.Series):
                error_message = (
                    f'{variable_name} is a function that returns a variable '
                    f'of type {type(column).__name__}, but must return a '
                    f'polars Series or 1D NumPy array')
                raise TypeError(error_message)
            if len(column) != len(obs_or_var):
                error_message = (
                    f'{variable_name} is a function that returns a column of '
                    f'length {len(column):,}, which differs from the length '
                    f'of {obs_or_var_name} ({len(obs_or_var):,})')
                raise ValueError(error_message)
        else:
            error_message = (
                f'{variable_name} must be a string column name, a polars '
                f'expression, a polars Series, a 1D NumPy array, or a '
                f'function that returns a polars Series or 1D NumPy array '
                f'when applied to this SingleCell dataset, but has type '
                f'{type(column).__name__!r}')
            raise TypeError(error_message)
        check_dtype(column, variable_name, dtypes)
        if not allow_null:
            if QC_column is None:
                null_count = column.null_count()
                if null_count > 0:
                    error_message = (
                        f'{variable_name} contains {null_count:,} '
                        f'{plural("null value", null_count)}, but must not '
                        f'contain any')
                    raise ValueError(error_message)
            else:
                null_count = (column.is_null() & QC_column).sum()
                if null_count > 0:
                    error_message = (
                        f'{variable_name} contains {null_count:,} '
                        f'{plural("null value", null_count)} for cells '
                        f'passing QC, but must not contain any')
                    raise ValueError(error_message)
        return column
    
    @staticmethod
    def _get_columns(obs_or_var_name: Literal['obs', 'var'],
                     datasets: Sequence[SingleCell],
                     columns: SingleCellColumn | None |
                              Sequence[SingleCellColumn | None],
                     variable_name: str,
                     dtypes: pl.datatypes.classes.DataTypeClass | str |
                             tuple[pl.datatypes.classes.DataTypeClass | str,
                                   ...],
                     *,
                     QC_columns: list[pl.Series | None] = None,
                     allow_None: bool = True,
                     allow_missing: bool = False,
                     allow_null: bool = False,
                     custom_error: str | None = None) -> \
            list[pl.Series | None]:
        """
        Get a column of the same length as obs/var from each dataset.
        
        Args:
            obs_or_var_name: the name of the DataFrame the column is with
                             respect to, i.e. `'obs'` or `'var'`
            datasets: a sequence of SingleCell datasets
            columns: a string naming a column of obs/var, a polars expression
                     that evaluates to a single column when applied to obs/var,
                     a polars Series or 1D NumPy array of the same length as
                     obs/var, or a function that takes in `self` and returns a
                     polars Series or 1D NumPy array of the same length as
                     obs/var. Or, a Sequence of these, one per dataset in
                     `datasets`. May also be None (or a Sequence containing
                     None) if `allow_None=True`.
            variable_name: the name of the variable corresponding to `columns`
            dtypes: the required dtype(s) of the columns
            QC_columns: an optional column of cells passing QC for each
                        dataset. If not None for a given dataset, the presence
                        of null values for that dataset will only raise an
                        error for cells passing QC. Has no effect when
                        `allow_null=True`.
            allow_None: whether to allow `columns` or its elements to be None
            allow_missing: whether to allow `columns` to be a string (or
                           contain strings) missing from certain datasets'
                           obs/var, returning None for these datasets
            allow_null: whether to allow `columns` to contain null values
            custom_error: a custom error message for when `column` is a string
                          and is not found in obs/var, and
                          `allow_missing=False`; use `{}` as a placeholder for
                          the name of the column
        
        Returns:
            A list of polars Series of the same length as `datasets`, where
            each Series has the same length as the corresponding dataset's
            obs/var. Or, if `columns` is None (or if some elements are None) or
            missing from obs/var (when `allow_missing=True`), a list of None
            (or where the corresponding elements are None).
        """
        if columns is None:
            if not allow_None:
                error_message = f'{variable_name} is None'
                raise TypeError(error_message)
            return [None] * len(datasets)
        if isinstance(columns, Sequence) and not isinstance(columns, str):
            if len(columns) != len(datasets):
                error_message = (
                    f'{variable_name} has length {len(columns):,}, but you '
                    f'specified {len(datasets):,} datasets')
                raise ValueError(error_message)
            if not allow_None and any(column is None for column in columns):
                error_message = \
                    f'{variable_name} contains an element that is None'
                raise TypeError(error_message)
            if QC_columns is None:
                return [dataset._get_column(
                    obs_or_var_name=obs_or_var_name, column=column,
                    variable_name=variable_name, dtypes=dtypes,
                    allow_null=allow_null, allow_missing=allow_missing,
                    custom_error=custom_error)
                    if column is not None else None
                    for dataset, column in zip(datasets, columns)]
            else:
                return [dataset._get_column(
                    obs_or_var_name=obs_or_var_name, column=column,
                    variable_name=variable_name, dtypes=dtypes,
                    QC_column=QC_column, allow_null=allow_null,
                    allow_missing=allow_missing, custom_error=custom_error)
                    if column is not None else None
                    for dataset, column, QC_column in
                    zip(datasets, columns, QC_columns)]
        else:
            if QC_columns is None:
                return [dataset._get_column(
                    obs_or_var_name=obs_or_var_name, column=columns,
                    variable_name=variable_name, dtypes=dtypes,
                    allow_null=allow_null, allow_missing=allow_missing,
                    custom_error=custom_error) for dataset in datasets]
            else:
                return [dataset._get_column(
                    obs_or_var_name=obs_or_var_name, column=columns,
                    variable_name=variable_name, dtypes=dtypes,
                    QC_column=QC_column, allow_null=allow_null,
                    allow_missing=allow_missing, custom_error=custom_error)
                    for dataset, QC_column in zip(datasets, QC_columns)]
    
    # noinspection PyUnresolvedReferences
    def to_anndata(self, *, QC_column: str | None = 'passed_QC') -> 'AnnData':
        """
        Converts this SingleCell dataset to an AnnData object.
        
        Make sure to remove cells failing QC with `filter_obs(QC_column)`
        first, or specify `subset=True` in `qc()`. Alternatively, to include
        cells failing QC in the AnnData object, set `QC_column` to None.
        
        Note that there is no `from_anndata()`; simply do
        `SingleCell(anndata_object)` to initialize a SingleCell dataset from an
        in-memory AnnData object.
        
        Args:
            QC_column: if not None, give an error if this column is present in
                       obs and not all cells pass QC
        
        Returns:
            An AnnData object. For AnnData versions older than 0.11.0, which
            do not support csr_array/csc_array, counts will be converted to
            csr_matrix/csc_matrix.
        """
        with ignore_sigint():
            import anndata
            import pandas as pd
            import pyarrow as pa
        valid_dtypes = pl.FLOAT_DTYPES | pl.INTEGER_DTYPES | \
                       {pl.String, pl.Categorical, pl.Enum, pl.Boolean}
        for df, df_name in (self._obs, 'obs'), (self._var, 'var'):
            for column, dtype in df.schema.items():
                if dtype.base_type() not in valid_dtypes:
                    error_message = (
                        f'{df_name}[{column!r}] has the data type '
                        f'{dtype.base_type()!r}, which is not supported by '
                        f'AnnData')
                    raise TypeError(error_message)
        if QC_column is not None:
            check_type(QC_column, 'QC_column', str, 'a string')
            if QC_column in self._obs:
                QCed_cells = self._obs[QC_column]
                check_dtype(QCed_cells, f'obs[{QC_column!r}]',
                            pl.Boolean)
                if QCed_cells.null_count() or not QCed_cells.all():
                    error_message = (
                        f'not all cells pass QC; remove cells failing QC with '
                        f'filter_obs({QC_column!r}) or by specifying '
                        f'subset=True in qc(), or set QC_column=None to '
                        f'include them in the AnnData object')
                    raise ValueError(error_message)
        type_mapping = {
            pa.int8(): pd.Int8Dtype(), pa.int16(): pd.Int16Dtype(),
            pa.int32(): pd.Int32Dtype(), pa.int64(): pd.Int64Dtype(),
            pa.uint8(): pd.UInt8Dtype(), pa.uint16(): pd.UInt16Dtype(),
            pa.uint32(): pd.UInt32Dtype(), pa.uint64(): pd.UInt64Dtype(),
            pa.string(): pd.StringDtype(storage='pyarrow'),
            pa.bool_(): pd.BooleanDtype()}
        to_pandas = lambda df: df\
            .to_pandas(split_blocks=True, types_mapper=type_mapping.get)\
            .set_index(df.columns[0])
        return anndata.AnnData(
            X=self._X if version.parse(anndata.__version__) >=
                         version.parse('0.11.0') else
              csr_matrix(self._X) if isinstance(self._X, csr_array)
              else csc_matrix(self._X),
            obs=to_pandas(self._obs), var=to_pandas(self._var),
            obsm=self._obsm, varm=self._varm, uns=self._uns)
    
    @staticmethod
    def _from_seurat(seurat_object_name: str,
                     *,
                     assay: str | None,
                     slot: str,
                     constructor: bool) -> \
            tuple[csr_array | csc_array, pl.DataFrame, pl.DataFrame,
                  dict[str, np.ndarray[2, Any]], NestedScalarOrArrayDict]:
        """
        Create a SingleCell dataset from an in-memory Seurat object loaded with
        the ryp Python-R bridge. Used by `__init__()` and `from_seurat()`.
        
        Args:
            seurat_object_name: the name of the Seurat object in the ryp R
                                workspace
            assay: the name of the assay within the Seurat object to load data
                   from; if None, defaults to `seurat_object@active.assay`
            slot: the slot within the active assay (or the assay specified by
                  the `assay` argument, if not None) to use as X. Set to
                  `'data'` to load the normalized counts, or `'scale.data'` to
                  load the normalized and scaled counts. If dense, will be
                  automatically converted to a sparse array.
            constructor: whether this method is being called from the
                         constructor or from `from_seurat()`
        
        Returns:
            A length-5 tuple of (X, obs, var, obsm, uns).
        """
        from ryp import r, to_py
        if assay is None:
            assay = to_py(f'{seurat_object_name}@active.assay')
        elif to_py(f'{seurat_object_name}@{assay}') is None:
            error_message = (
                f'assay {assay!r} does not exist in the Seurat object '
                f'{seurat_object_name!r}')
            raise ValueError(error_message)
        # If Seurat v5, merge layers if necessary, and use $slot
        # instead of @slot for X and meta.data instead of
        # meta.features for var
        if to_py(f'inherits({seurat_object_name}@assays${assay}, "Assay5")'):
            if not to_py(f'"{slot}" %in% names({seurat_object_name}@assays$'
                         f'{assay}@layers)'):
                error_message = (
                    f'slot {slot!r} does not exist in '
                    f'{seurat_object_name}@assays${assay}@layers')
                raise ValueError(error_message)
            if to_py(f'length({seurat_object_name}@assays${assay}@'
                     f'layers)') > 1:
                r(f'{seurat_object_name}@assays${assay} = '
                  f'JoinLayers({seurat_object_name}@assays${assay}, "{slot}")')
            X_slot = f'{seurat_object_name}@assays${assay}${slot}'
            var = to_py(f'{seurat_object_name}@assays${assay}@meta.data')
        else:
            # unlike v5 objects, v3 objects indicate the absence of a slot with
            # a 0 x 0 matrix
            if not (to_py(f'"{slot}" %in% slotNames({seurat_object_name}@'
                         f'assays${assay})') and
                    to_py(f'prod(dim({seurat_object_name}@assays${assay}@'
                          f'{slot}))') > 0):
                error_message = (
                    f'slot {slot!r} does not exist in '
                    f'{seurat_object_name}@assays${assay}')
                raise ValueError(error_message)
            X_slot = f'{seurat_object_name}@assays${assay}@{slot}'
            var = to_py(f'{seurat_object_name}@assays${assay}@meta.features')
        X_classes = tuple(to_py(f'class({X_slot})', squeeze=False))
        if X_classes == ('dgCMatrix',):
            X = to_py(X_slot).T
        elif X_classes == ('matrix', 'array'):
            X = csr_array(to_py(X_slot, format='numpy').T)
        else:
            error_message = (
                f'{X_slot} must be a dgCMatrix (column-oriented sparse '
                f'matrix) or matrix, but has ')
            if len(X_classes) == 0:
                error_message += 'no classes'
            elif len(X_classes) == 1:
                error_message += f'the class {X_classes[0]!r}'
            else:
                error_message += (
                    f'the classes '
                    f'{", ".join(f"{c!r}" for c in X_classes[:-1])} and '
                    f'{X_classes[-1]}')
            error_message += (
                f'; consider setting {"X_key" if constructor else "slot"} to '
                f'a different value than {slot!r}')
            raise TypeError(error_message)
        obs_key = f'{seurat_object_name}@meta.data'
        obs = to_py(obs_key, index='_index' if to_py(
            f'"cell" %in% {obs_key}') else 'cell')
        if var is None:
            var = to_py(f'rownames({seurat_object_name}@assays${assay}@'
                        f'{slot}').to_frame('gene')
        reductions = to_py(f'names({seurat_object_name}@reductions)')
        obsm = {key: to_py(
            f'{seurat_object_name}@reductions${key}@cell.embeddings',
            format='numpy') for key in reductions
            if not to_py(f'is.null({seurat_object_name}@reductions${key})')} \
            if reductions is not None else {}
        obs = obs.cast({column.name: pl.Enum(column.cat.get_categories())
                        for column in obs.select(pl.col(pl.Categorical))})
        var = var.cast({column.name: pl.Enum(column.cat.get_categories())
                        for column in var.select(pl.col(pl.Categorical))})
        uns = to_py(f'{seurat_object_name}@misc')
        return X, obs, var, obsm, uns
    
    @staticmethod
    def from_seurat(seurat_object_name: str, *,
                    assay: str | None = None,
                    slot: str = 'counts') -> SingleCell:
        """
        Create a SingleCell dataset from a Seurat object that has already been
        loaded into memory via the ryp Python-R bridge. To load a Seurat object
        from disk, use e.g. `SingleCell('filename.rds')`.
        
        Args:
            seurat_object_name: the name of the Seurat object in the ryp R
                                workspace
            assay: the name of the assay within the Seurat object to load data
                   from; if None, defaults to seurat_object@active.assay
            slot: the slot within the active assay (or the assay specified by
                  the `assay` argument, if not None) to use as X. Defaults to
                  `'counts'`. Set to `'data'` to load the normalized counts,
                  or `'scale.data'` to load the normalized and scaled counts.
                  If dense, will be automatically converted to a sparse array.
        
        Returns:
            The corresponding SingleCell dataset.
        """
        from ryp import to_py
        check_type(seurat_object_name, 'seurat_object_name', str, 'a string')
        check_R_variable_name(seurat_object_name, 'seurat_object_name')
        if assay is not None:
            check_type(assay, 'assay', str, 'a string')
        check_type(slot, 'slot', str, 'a string')
        if not to_py(f'inherits({seurat_object_name}, "Seurat")'):
            classes = to_py(f'class({seurat_object_name})', squeeze=False)
            error_message = (
                f'the R object named by seurat_object_name, '
                f'{seurat_object_name}, must be a Seurat object, but has ')
            if len(classes) == 0:
                error_message += 'no classes'
            elif len(classes) == 1:
                error_message += f'the class {classes[0]!r}'
            else:
                error_message += (
                    f'the classes '
                    f'{", ".join(f"{c!r}" for c in classes[:-1])} and '
                    f'{classes[-1]!r}')
            raise TypeError(error_message)
        X, obs, var, obsm, uns = \
            SingleCell._from_seurat(seurat_object_name, assay=assay, slot=slot,
                                    constructor=False)
        return SingleCell(X=X, obs=obs, var=var, obsm=obsm, uns=uns)
    
    def to_seurat(self,
                  seurat_object_name: str,
                  *,
                  QC_column: str | None = 'passed_QC') -> None:
        """
        Convert this SingleCell dataset to a Seurat object (version 3, not
        version 5) in the ryp R workspace. varm, DataFrame keys of obsm, and
        non-string keys of uns will not be converted.
        
        Make sure to remove cells failing QC with `filter_obs(QC_column)`
        first, or specify `subset=True` in `qc()`. Alternatively, to include
        cells failing QC in the Seurat object, set `QC_column` to None.
        
        When converting to Seurat, to match the requirements of Seurat objects,
        the `'X_'` prefix (often used by Scanpy) will be removed from each key
        of obsm where it is present (e.g. `'X_umap'` will become `'umap'`).
        Seurat will also add `'orig.ident'`, `'nCount_RNA'` and
        `'nFeature_RNA'` as gene-level metadata by default; you can disable
        the calculation of the latter two columns with:
        
        ```python
        from ryp import r
        r('options(Seurat.object.assay.calcn = FALSE)')
        ```
        
        Args:
            seurat_object_name: the name of the R variable to assign the Seurat
                                object to
            QC_column: if not None, give an error if this column is present in
                       obs and not all cells pass QC
        """
        from ryp import r, to_py, to_r
        r('suppressPackageStartupMessages(library(SeuratObject))')
        if self._X.nnz > 2_147_483_647:
            error_message = (
                f'X has {self._X.nnz:,} non-zero elements, more than '
                f'INT32_MAX (2,147,483,647), the maximum supported in R')
            raise ValueError(error_message)
        check_type(seurat_object_name, 'seurat_object_name', str, 'a string')
        check_R_variable_name(seurat_object_name, 'seurat_object_name')
        if QC_column is not None:
            check_type(QC_column, 'QC_column', str, 'a string')
            if QC_column in self._obs:
                QCed_cells = self._obs[QC_column]
                check_dtype(QCed_cells, f'obs[{QC_column!r}]',
                            pl.Boolean)
                if QCed_cells.null_count() or not QCed_cells.all():
                    error_message = (
                        f'not all cells pass QC; remove cells failing QC with '
                        f'filter_obs({QC_column!r}) or by specifying '
                        f'subset=True in qc(), or set QC_column=None to '
                        f'include them in the Seurat object')
                    raise ValueError(error_message)
        valid_dtypes = pl.FLOAT_DTYPES | pl.INTEGER_DTYPES | \
                       {pl.String, pl.Categorical, pl.Enum, pl.Boolean}
        for df, df_name in (self._obs, 'obs'), (self._var, 'var'):
            for column, dtype in df.schema.items():
                if dtype.base_type() not in valid_dtypes:
                    error_message = (
                        f'{df_name}[{column!r}] has the data type '
                        f'{dtype.base_type()!r}, which is not supported when '
                        f'converting to a Seurat object')
                    raise TypeError(error_message)
        is_string = self.var_names.dtype == pl.String
        num_with_underscores = self.var_names.str.contains('_').sum() \
            if is_string else \
            self.var_names.cat.get_categories().str.contains('_').sum()
        if num_with_underscores:
            var_names_expression = f'pl.col.{self.var_names.name}' \
                if self.var_names.name.isidentifier() else \
                f'pl.col({self.var_names.name!r})'
            error_message = (
                f"var_names contains {num_with_underscores:,}"
                f"{'' if is_string else ' unique'} gene "
                f"{plural('name', num_with_underscores)} with "
                f"underscores, which are not supported by Seurat; Seurat "
                f"recommends changing the underscores to dashes, which you "
                f"can do with .with_columns_var({var_names_expression}"
                f"{'' if is_string else '.cast(pl.String)'}"
                f".str.replace_all('_', '-'))")
            raise ValueError(error_message)
        try:
            to_r(self._X.T, '.SingleCell.X.T', rownames=self.var_names,
                 colnames=self.obs_names)
            try:
                to_r(self._obs.drop(self.obs_names.name), '.SingleCell.obs',
                     rownames=self.obs_names)
                try:
                    r(f'{seurat_object_name} = CreateSeuratObject('
                      'CreateAssayObject(.SingleCell.X.T), '
                      'meta.data=.SingleCell.obs)')
                    # Reverse the column name-mangling introduced by Seurat
                    r('.SingleCell.cols_to_ignore = '
                      'c("orig.ident", "nCount_RNA", "nFeature_RNA")')
                    try:
                        r(f'names({seurat_object_name}@meta.data)['
                          f'!names({seurat_object_name}@meta.data) %in% '
                          f'.SingleCell.cols_to_ignore] = '
                          f'names(.SingleCell.obs)[!names(.SingleCell.obs) '
                          f'%in% .SingleCell.cols_to_ignore]')
                    finally:
                        r('rm(.SingleCell.cols_to_ignore)')
                finally:
                    r('rm(.SingleCell.obs)')
            finally:
                r('rm(.SingleCell.X.T)')
            to_r(self._var.drop(self.var_names.name), '.SingleCell.var',
                 rownames=self.var_names)
            try:
                r(f'{seurat_object_name}@assays$RNA@meta.features = '
                  f'.SingleCell.var')
            finally:
                r('rm(.SingleCell.var)')
            if self._obsm:
                for key, value in self._obsm.items():
                    if isinstance(value, pl.DataFrame):
                        continue
                    # Remove the initial X_ from the reduction name and suffix
                    # with `'_'` when creating the key, like
                    # mojaveazure.github.io/seurat-disk/reference/Convert.html
                    key = key.removeprefix('X_')
                    to_r(value, '.SingleCell.value', rownames=self.obs_names,
                         colnames=[f'{key}_{i}'
                                   for i in range(1, value.shape[1] + 1)])
                    try:
                        r(f'{seurat_object_name}@reductions${key} = '
                          f'CreateDimReducObject(.SingleCell.value, '
                          f'key="{key}_", assay="RNA")')
                    finally:
                        r('rm(.SingleCell.value)')
            if self._uns:
                to_r({key: value for key, value in self._uns.items()
                      if isinstance(value, str)}, '.SingleCell.uns')
                try:
                    r(f'{seurat_object_name}@misc = .SingleCell.uns')
                finally:
                    r('rm(.SingleCell.uns)')
        except:
            if to_py(f'exists({seurat_object_name!r})'):
                r(f'rm({seurat_object_name})')
            raise
    
    @staticmethod
    def _get_DFrame(dframe: str, *, index: str) -> pl.DataFrame:
        """
        Convert a DFrame in the ryp R workspace to a Python DataFrame, raising
        an error if the DFrame contains any nested data stuctures.
        
        Args:
            dframe: the name of the DFrame object in the ryp R workspace
            index: the name of the column containing the rownames in the output
                   DataFrame; if this column name is already present in the
                   DFrame, fall back to calling the rownames column `'_index'`

        Returns:
            A polars DataFrame containing the data in `dframe`.
        """
        from ryp import to_py
        df = to_py(f'{dframe}@listData')
        if not all(isinstance(value, pl.Series) for value in df.values()):
            error_message = (
                f'{dframe} contains nested data; unnest before converting to '
                f'a SingleCell dataset')
            raise ValueError(error_message)
        if index in df.keys():
            index = '_index'
        df = pl.DataFrame({index: to_py(f'{dframe}@rownames')} | df)
        return df
    
    @staticmethod
    def _from_sce(sce_object_name: str,
                  *,
                  slot: str,
                  constructor: bool) -> \
            tuple[csr_array | csc_array, pl.DataFrame, pl.DataFrame,
                  dict[str, np.ndarray[2, Any]], NestedScalarOrArrayDict]:
        """
        Create a SingleCell dataset from an in-memory SingleCellExperiment
        object loaded with the ryp Python-R bridge. Used by `__init__()` and
        `from_sce()`.
        
        Args:
            sce_object_name: the name of the SingleCellExperiment object in the
                             ryp R workspace
            slot: the element within `sce_object@assays@data` to use as `X`.
                  Set to `'counts'` to load raw counts, or `'logcounts'` to
                  load the normalized counts if available. If dense, will be
                  automatically converted to a sparse array.
            constructor: whether this method is being called from the
                         constructor or from `from_sce()`
        
        Returns:
            A length-5 tuple of (X, obs, var, obsm, uns).
        """
        from ryp import to_py
        X_slot = f'{sce_object_name}@assays@data${slot}'
        X_classes = tuple(to_py(f'class({X_slot})', squeeze=False))
        if X_classes == ('dgCMatrix',):
            X = to_py(X_slot).T
        elif X_classes == ('matrix', 'array'):
            X = csr_array(to_py(X_slot, format='numpy').T)
        else:
            error_message = (
                f'{X_slot} must be a dgCMatrix (column-oriented sparse '
                f'matrix) or matrix, but has ')
            if len(X_classes) == 0:
                error_message += 'no classes'
            elif len(X_classes) == 1:
                error_message += f'the class {X_classes[0]!r}'
            else:
                error_message += (
                    f'the classes '
                    f'{", ".join(f"{c!r}" for c in X_classes[:-1])} and '
                    f'{X_classes[-1]}')
            error_message += (
                f'; consider setting {"X_key" if constructor else "slot"} to '
                f'a different value than {slot!r}')
            raise TypeError(error_message)
        obs = SingleCell._get_DFrame(f'colData({sce_object_name})',
                                     index='cell')
        var = SingleCell._get_DFrame(f'rowData({sce_object_name})',
                                     index='gene')
        obs = obs.cast({column.name: pl.Enum(column.cat.get_categories())
                        for column in obs.select(pl.col(pl.Categorical))})
        var = var.cast({column.name: pl.Enum(column.cat.get_categories())
                        for column in var.select(pl.col(pl.Categorical))})
        obsm = to_py(f'reducedDims({sce_object_name})@listData',
                     format='numpy')
        uns = to_py(f'{sce_object_name}@metadata', format='numpy')
        return X, obs, var, obsm, uns
    
    @staticmethod
    def from_sce(sce_object_name: str,
                 *,
                 slot: str = 'counts') -> SingleCell:
        """
        Create a SingleCell dataset from a SingleCellExperiment object that has
        already been loaded into memory via the ryp Python-R bridge. To load a
        SingleCellExperiment object from disk, use e.g.
        `SingleCell('filename.rds')`.
        
        Args:
            sce_object_name: the name of the SingleCellExperiment object in the
                             ryp R workspace
            slot: the element within `{sce_object_name}@assays@data` to use as
                  `X`. Defaults to `'counts'`. If available, set to
                  `'logcounts'` to load the normalized counts. If dense, will
                  be automatically converted to a sparse array.
        
        Returns:
            The corresponding SingleCell dataset.
        """
        from ryp import to_py
        check_type(sce_object_name, 'sce_object_name', str, 'a string')
        check_R_variable_name(sce_object_name, 'sce_object_name')
        check_type(slot, 'slot', str, 'a string')
        if not to_py(f'inherits({sce_object_name}, "SingleCellExperiment")'):
            classes = to_py(f'class({sce_object_name})', squeeze=False)
            error_message = (
                f'the R object named by sce_object_name, {sce_object_name}, '
                f'must be a SingleCellExperiment object, but has ')
            if len(classes) == 0:
                error_message += 'no classes'
            elif len(classes) == 1:
                error_message += f'the class {classes[0]!r}'
            else:
                error_message += (
                    f'the classes '
                    f'{", ".join(f"{c!r}" for c in classes[:-1])} and '
                    f'{classes[-1]!r}')
            raise TypeError(error_message)
        X, obs, var, obsm, uns = \
            SingleCell._from_sce(sce_object_name, slot=slot, constructor=False)
        return SingleCell(X=X, obs=obs, var=var, obsm=obsm, uns=uns)
        
    def to_sce(self,
               sce_object_name: str,
               *,
               QC_column: str | None = 'passed_QC') -> None:
        """
        Convert this SingleCell dataset to a SingleCellExperiment object in the
        ryp R workspace. varm and DataFrame keys of obsm will not be converted.
        
        Make sure to remove cells failing QC with `filter_obs(QC_column)`
        first, or specify `subset=True` in `qc()`. Alternatively, to include
        cells failing QC in the SingleCellExperiment object, set `QC_column` to
        None.
        
        Args:
            sce_object_name: the name of the R variable to assign the
                             SingleCellExperiment object to
            QC_column: if not None, give an error if this column is present in
                       obs and not all cells pass QC
        """
        from ryp import r, to_py, to_r
        r('suppressPackageStartupMessages(library(SingleCellExperiment))')
        if self._X.nnz > 2_147_483_647:
            error_message = (
                f'X has {self._X.nnz:,} non-zero elements, more than '
                f'INT32_MAX (2,147,483,647), the maximum supported in R')
            raise ValueError(error_message)
        check_type(sce_object_name, 'sce_object_name', str, 'a string')
        check_R_variable_name(sce_object_name, 'sce_object_name')
        if QC_column is not None:
            check_type(QC_column, 'QC_column', str, 'a string')
            if QC_column in self._obs:
                QCed_cells = self._obs[QC_column]
                check_dtype(QCed_cells, f'obs[{QC_column!r}]',
                            pl.Boolean)
                if QCed_cells.null_count() or not QCed_cells.all():
                    error_message = (
                        f'not all cells pass QC; remove cells failing QC with '
                        f'filter_obs({QC_column!r}) or by specifying '
                        f'subset=True in qc(), or set QC_column=None to '
                        f'include them in the SingleCellExperiment object')
                    raise ValueError(error_message)
        valid_dtypes = pl.FLOAT_DTYPES | pl.INTEGER_DTYPES | \
                       {pl.String, pl.Categorical, pl.Enum, pl.Boolean}
        for df, df_name in (self._obs, 'obs'), (self._var, 'var'):
            for column, dtype in df.schema.items():
                if dtype.base_type() not in valid_dtypes:
                    error_message = (
                        f'{df_name}[{column!r}] has the data type '
                        f'{dtype.base_type()!r}, which is not supported when '
                        f'converting to a SingleCellExperiment object')
                    raise TypeError(error_message)
        try:
            to_r(self._X.T, '.SingleCell.X.T', rownames=self.var_names,
                 colnames=self.obs_names)
            try:
                to_r(self._obs.drop(self.obs_names.name), '.SingleCell.obs',
                     rownames=self.obs_names)
                try:
                    to_r(self._var.drop(self.var_names.name),
                         '.SingleCell.var', rownames=self.var_names)
                    try:
                        r(f'{sce_object_name} = SingleCellExperiment('
                          f'assays = list(counts = .SingleCell.X.T), '
                          f'colData = S4Vectors::DataFrame('
                          f'    .SingleCell.obs, check.names=FALSE), '
                          f'rowData = S4Vectors::DataFrame('
                          f'    .SingleCell.var, check.names=FALSE))')
                    finally:
                        r('rm(.SingleCell.var)')
                finally:
                    r('rm(.SingleCell.obs)')
            finally:
                r('rm(.SingleCell.X.T)')
            if self._obsm:
                for key, value in self._obsm.items():
                    if isinstance(value, pl.DataFrame):
                        continue
                    to_r(value, '.SingleCell.value', rownames=self.obs_names,
                             colnames=[f'{key}_{i}'
                                       for i in range(1, value.shape[1] + 1)])
                    try:
                        r(f'reducedDim({sce_object_name}, {key!r}) = '
                          f'.SingleCell.value')
                    finally:
                        r('rm(.SingleCell.value)')
            if self._uns:
                to_r(self._uns, '.SingleCell.uns')
                try:
                    r(f'{sce_object_name}@metadata = .SingleCell.uns')
                finally:
                    r('rm(.SingleCell.uns)')
        except:
            if to_py(f'exists({sce_object_name!r})'):
                r(f'rm({sce_object_name})')
            raise
        
    def copy(self, deep: bool = False) -> SingleCell:
        """
        Make a deep (if `deep=True`) or shallow copy of this SingleCell
        dataset.
        
        Returns:
            A copy of the SingleCell dataset. Since polars DataFrames are
            immutable, obs and var will always point to the same underlying
            data as the original. The only difference when deep=True is that X
            will point to a fresh copy of the data, rather than the same data.
            Watch out: when deep=False, any modifications to X will modify both
            copies!
        """
        check_type(deep, 'deep', bool, 'Boolean')
        # noinspection PyTypeChecker
        return SingleCell(X=self._X.copy() if deep else self._X, obs=self._obs,
                          var=self._var, obsm=self._obsm, varm=self._varm,
                          uns=self._uns)
    
    def concat_obs(self,
                   datasets: SingleCell,
                   *more_datasets: SingleCell,
                   flexible: bool = False) -> SingleCell:
        """
        Concatenate the cells of multiple SingleCell datasets.
        
        By default, all datasets must have the same var, varm and uns. They
        must also have the same columns in obs and the same keys in obsm, with
        the same data types.
        
        Conversely, if `flexible=True`, subset to genes present in all datasets
        (according to the first column of var, i.e. `var_names`) before
        concatenating. Subset to columns of var and keys of varm and uns that
        are identical in all datasets after this subsetting. Also, subset to
        columns of obs and keys of obsm that are present in all datasets, and
        have the same data types. All datasets' `obs_names` must have the same
        name and dtype, and similarly for `var_names`.
        
        The one exception to the obs "same data type" rule: if a column is Enum
        in some datasets and Categorical in others, or Enum in all datasets but
        with different categories in each dataset, that column will be retained
        as an Enum column (with the union of the categories) in the
        concatenated obs.
        
        Args:
            datasets: one or more SingleCell datasets to concatenate with this
                      one
            *more_datasets: additional SingleCell datasets to concatenate with
                            this one, specified as positional arguments
            flexible: whether to subset to genes, columns of obs and var, and
                      keys of obsm, varm and uns common to all datasets before
                      concatenating, rather than raising an error on any
                      mismatches
        
        Returns:
            The concatenated SingleCell dataset.
        """
        # Check inputs
        if isinstance(datasets, SingleCell):
            datasets = datasets,
        datasets = (self,) + datasets + more_datasets
        if len(datasets) == 1:
            error_message = \
                'need at least one other SingleCell dataset to concatenate'
            raise ValueError(error_message)
        check_types(datasets[1:], 'datasets', SingleCell,
                    'SingleCell datasets')
        check_type(flexible, 'flexible', bool, 'Boolean')
        # Perform either flexible or non-flexible concatenation
        if flexible:
            # Check that `obs_names` and `var_names` have the same name and
            # data type across all datasets
            obs_names_name = self.obs_names.name
            if not all(dataset.obs_names.name == obs_names_name
                       for dataset in datasets[1:]):
                error_message = (
                    'not all SingleCell datasets have the same name for the '
                    'first column of obs (the obs_names column)')
                raise ValueError(error_message)
            var_names_name = self.var_names.name
            if not all(dataset.var_names.name == var_names_name
                       for dataset in datasets[1:]):
                error_message = (
                    'not all SingleCell datasets have the same name for the '
                    'first column of var (the var_names column)')
                raise ValueError(error_message)
            obs_names_dtype = self.obs_names.dtype
            if not all(dataset.obs_names.dtype == obs_names_dtype
                       for dataset in datasets[1:]):
                error_message = (
                    'not all SingleCell datasets have the same data type for '
                    'the first column of obs (the obs_names column)')
                raise TypeError(error_message)
            var_names_dtype = self.var_names.dtype
            if not all(dataset.var_names.dtype == var_names_dtype
                       for dataset in datasets[1:]):
                error_message = (
                    'not all SingleCell datasets have the same data type for '
                    'the first column of var (the var_names column)')
                raise TypeError(error_message)
            # Subset to genes in common across all datasets
            genes_in_common = self.var_names\
                .filter(self.var_names
                        .is_in(pl.concat([dataset.var_names
                                          for dataset in datasets[1:]])))
            if len(genes_in_common) == 0:
                error_message = \
                    'no genes are shared across all SingleCell datasets'
                raise ValueError(error_message)
            datasets = [dataset[:, genes_in_common] for dataset in datasets]
            # Subset to columns of var and keys of varm and uns that are
            # identical in all datasets after this subsetting
            var_columns_in_common = [
                column.name for column in datasets[0]._var[:, 1:]
                if all(column.name in dataset._var and
                       dataset._var[column.name].equals(column)
                       for dataset in datasets[1:])]
            varm = self._varm
            # noinspection PyUnresolvedReferences
            varm_keys_in_common = [
                key for key, value in varm.items()
                if all(key in dataset._varm and
                       type(varm[key]) is type(dataset._varm[key]) and
                       dataset._varm[key].dtype == value.dtype and
                       (np.array_equal(dataset._varm[key], value,
                                       equal_nan=True)
                        if isinstance(varm[key], np.ndarray) else
                        dataset._varm[key].equals(varm[key]))
                       for dataset in datasets[1:])]
            # noinspection PyTypeChecker,PyUnresolvedReferences
            uns_keys_in_common = [
                key for key, value in self._uns.items()
                if isinstance(value, dict) and
                   all(isinstance(dataset._uns[key], dict) and
                       SingleCell._eq_uns(value, dataset._uns[key])
                       for dataset in datasets[1:]) or
                   isinstance(value, np.ndarray) and
                   all(isinstance(dataset._uns[key], np.ndarray) and
                       np.array_equal(dataset._uns[key], value, equal_nan=True)
                       for dataset in datasets[1:]) or
                   all(not isinstance(dataset._uns[key], (dict, np.ndarray))
                       and dataset._uns[key] == value
                       for dataset in datasets[1:])]
            for dataset in datasets:
                dataset._var = dataset._var.select(dataset.var_names,
                                                   *var_columns_in_common)
                dataset._varm = {key: dataset._varm[key]
                                  for key in varm_keys_in_common}
                dataset._uns = {key: dataset._uns[key]
                                for key in uns_keys_in_common}
            # Subset to columns of obs and keys of obsm that are present in all
            # datasets, and have the same data types. Also include columns of
            # obs that are Enum in some datasets and Categorical in others, or
            # Enum in all datasets but with different categories in each
            # dataset; cast these to Categorical.
            obs_mismatched_categoricals = {
                column for column, dtype in self._obs[:, 1:]
                .select(pl.col(pl.Categorical, pl.Enum)).schema.items()
                if all(column in dataset._obs and
                       dataset._obs[column].dtype in (pl.Categorical, pl.Enum)
                       for dataset in datasets[1:]) and
                   not all(dataset._obs[column].dtype == dtype
                           for dataset in datasets[1:])}
            obs_columns_in_common = [
                column
                for column, dtype in islice(self._obs.schema.items(), 1, None)
                if column in obs_mismatched_categoricals or
                   all(column in dataset._obs and
                       dataset._obs[column].dtype == dtype
                       for dataset in datasets[1:])]
            cast_dict = {column: pl.Enum(
                pl.concat([dataset._obs[column].cat.get_categories()
                           for dataset in datasets])
                .unique(maintain_order=True))
                for column in obs_mismatched_categoricals}
            for dataset in datasets:
                dataset._obs = dataset._obs\
                    .select(dataset.obs_names, *obs_columns_in_common)\
                    .cast(cast_dict)
            obsm_keys_in_common = [
                key for key, value in self._obsm.items()
                if all(key in dataset._obsm and
                       type(dataset._obsm[key]) is type(value) and
                       dataset._obsm[key].dtype == value.dtype
                       for dataset in datasets[1:])]
            for dataset in datasets:
                dataset._obsm = {key: dataset._obsm[key]
                                 for key in obsm_keys_in_common}
        else:  # non-flexible
            # Check that all var, varm and uns are identical
            var = self._var
            for dataset in datasets[1:]:
                if not dataset._var.equals(var):
                    error_message = (
                        'all SingleCell datasets must have the same var, '
                        'unless flexible=True')
                    raise ValueError(error_message)
            varm = self._varm
            for dataset in datasets[1:]:
                if dataset._varm.keys() != varm.keys() or \
                        any(type(dataset._varm[key]) is not type(value) or
                            dataset._varm[key].dtype != value.dtype
                            for key, value in varm.items()) or \
                        any((dataset._varm[key] != value).any()
                            for key, value in varm.items()):
                    error_message = (
                        'all SingleCell datasets must have the same varm, '
                        'unless flexible=True')
                    raise ValueError(error_message)
            for dataset in datasets[1:]:
                if not SingleCell._eq_uns(self._uns, dataset._uns):
                    error_message = (
                        'all SingleCell datasets must have the same uns, '
                        'unless flexible=True')
                    raise ValueError(error_message)
            # Check that all obs have the same columns and data types
            schema = self._obs.schema
            for dataset in datasets[1:]:
                if dataset._obs.schema != schema:
                    error_message = (
                        'all SingleCell datasets must have the same columns '
                        'in obs, with the same data types, unless '
                        'flexible=True')
                    raise ValueError(error_message)
            # Check that all obsm have the same keys and data types
            obsm = self._obsm
            for dataset in datasets[1:]:
                if dataset._obsm.keys() != obsm.keys() or any(
                        type(dataset._obsm[key]) is not type(value) or
                        dataset._obsm[key].dtype != value.dtype
                        for key, value in obsm.items()):
                    error_message = (
                        'all SingleCell datasets must have the same keys in '
                        'obsm, with the same data types, unless flexible=True')
                    raise ValueError(error_message)
        # Concatenate; output should be CSR when there's a mix of inputs
        X = vstack([dataset._X for dataset in datasets],
                   format='csr' if any(isinstance(dataset._X, csr_array)
                                       for dataset in datasets) else 'csc')
        obs = pl.concat([dataset._obs for dataset in datasets])
        # noinspection PyTypeChecker
        obsm = {key: np.concatenate([dataset._obsm[key]
                                     for dataset in datasets])
                if isinstance(value, np.ndarray) else
                pl.concat([dataset._obsm[key] for dataset in datasets])
                for key, value in self._obsm.items()}
        return SingleCell(X=X, obs=obs, var=datasets[0]._var, obsm=obsm,
                          varm=datasets[0]._varm, uns=datasets[0]._uns)

    def concat_var(self,
                   datasets: SingleCell,
                   *more_datasets: SingleCell,
                   flexible: bool = False) -> SingleCell:
        """
        Concatenate the genes of multiple SingleCell datasets.
        
        By default, all datasets must have the same obs, obsm and uns. They
        must also have the same columns in var and the same keys in varm, with
        the same data types.
        
        Conversely, if `flexible=True`, subset to genes present in all datasets
        (according to the first column of obs, i.e. `obs_names`) before
        concatenating. Subset to columns of obs and keys of obsm and uns that
        are identical in all datasets after this subsetting. Also, subset to
        columns of var and keys of varm that are present in all datasets, and
        have the same data types. All datasets' `var_names` must have the same
        name and dtype, and similarly for `obs_names`.
        
        The one exception to the var "same data type" rule: if a column is Enum
        in some datasets and Categorical in others, or Enum in all datasets but
        with different categories in each dataset, that column will be retained
        as an Enum column (with the union of the categories) in the
        concatenated var.
        
        Args:
            datasets: one or more SingleCell datasets to concatenate with this
                      one
            *more_datasets: additional SingleCell datasets to concatenate with
                            this one, specified as positional arguments
            flexible: whether to subset to genes, columns of var and obs, and
                      keys of varm, obsm and uns common to all datasets before
                      concatenating, rather than raising an error on any
                      mismatches
        
        Returns:
            The concatenated SingleCell dataset.
        """
        # Check inputs
        if isinstance(datasets, SingleCell):
            datasets = datasets,
        datasets = (self,) + datasets + more_datasets
        if len(datasets) == 1:
            error_message = \
                'need at least one other SingleCell dataset to concatenate'
            raise ValueError(error_message)
        check_types(datasets[1:], 'datasets', SingleCell,
                    'SingleCell datasets')
        check_type(flexible, 'flexible', bool, 'Boolean')
        # Perform either flexible or non-flexible concatenation
        if flexible:
            # Check that `var_names` and `obs_names` have the same name and
            # data type across all datasets
            var_names_name = self.var_names.name
            if not all(dataset.var_names.name == var_names_name
                       for dataset in datasets[1:]):
                error_message = (
                    'not all SingleCell datasets have the same name for the '
                    'first column of var (the var_names column)')
                raise ValueError(error_message)
            obs_names_name = self.obs_names.name
            if not all(dataset.obs_names.name == obs_names_name
                       for dataset in datasets[1:]):
                error_message = (
                    'not all SingleCell datasets have the same name for the '
                    'first column of obs (the obs_names column)')
                raise ValueError(error_message)
            var_names_dtype = self.var_names.dtype
            if not all(dataset.var_names.dtype == var_names_dtype
                       for dataset in datasets[1:]):
                error_message = (
                    'not all SingleCell datasets have the same data type for '
                    'the first column of var (the var_names column)')
                raise TypeError(error_message)
            obs_names_dtype = self.obs_names.dtype
            if not all(dataset.obs_names.dtype == obs_names_dtype
                       for dataset in datasets[1:]):
                error_message = (
                    'not all SingleCell datasets have the same data type for '
                    'the first column of obs (the obs_names column)')
                raise TypeError(error_message)
            # Subset to genes in common across all datasets
            genes_in_common = self.obs_names\
                .filter(self.obs_names
                        .is_in(pl.concat([dataset.obs_names
                                          for dataset in datasets[1:]])))
            if len(genes_in_common) == 0:
                error_message = \
                    'no genes are shared across all SingleCell datasets'
                raise ValueError(error_message)
            datasets = [dataset[:, genes_in_common] for dataset in datasets]
            # Subset to columns of obs and keys of obsm and uns that are
            # identical in all datasets after this subsetting
            obs_columns_in_common = [
                column.name for column in datasets[0]._obs[:, 1:]
                if all(column.name in dataset._obs and
                       dataset._obs[column.name].equals(column)
                       for dataset in datasets[1:])]
            obsm = self._obsm
            # noinspection PyUnresolvedReferences
            obsm_keys_in_common = [
                key for key, value in obsm.items()
                if all(key in dataset._obsm and
                       type(obsm[key]) is type(dataset._obsm[key]) and
                       dataset._obsm[key].dtype == value.dtype and
                       (np.array_equal(dataset._obsm[key], value,
                                       equal_nan=True)
                        if isinstance(obsm[key], np.ndarray) else
                        dataset._obsm[key].equals(obsm[key]))
                       for dataset in datasets[1:])]
            # noinspection PyTypeChecker,PyUnresolvedReferences
            uns_keys_in_common = [
                key for key, value in self._uns.items()
                if isinstance(value, dict) and
                   all(isinstance(dataset._uns[key], dict) and
                       SingleCell._eq_uns(value, dataset._uns[key])
                       for dataset in datasets[1:]) or
                   isinstance(value, np.ndarray) and
                   all(isinstance(dataset._uns[key], np.ndarray) and
                       np.array_equal(dataset._uns[key], value, equal_nan=True)
                       for dataset in datasets[1:]) or
                   all(not isinstance(dataset._uns[key], (dict, np.ndarray))
                       and value == dataset._uns[key]
                       for dataset in datasets[1:])]
            for dataset in datasets:
                dataset._obs = dataset._obs.select(dataset.obs_names,
                                                   *obs_columns_in_common)
                dataset._obsm = {key: dataset._obsm[key]
                                  for key in obsm_keys_in_common}
                dataset._uns = {key: dataset._uns[key]
                                for key in uns_keys_in_common}
            # Subset to columns of var and keys of varm that are present in all
            # datasets, and have the same data types. Also include columns of
            # var that are Enum in some datasets and Categorical in others, or
            # Enum in all datasets but with different categories in each
            # dataset; cast these to Categorical.
            var_mismatched_categoricals = {
                column for column, dtype in self._var[:, 1:]
                .select(pl.col(pl.Categorical, pl.Enum)).schema.items()
                if all(column in dataset._var and
                       dataset._var[column].dtype in (pl.Categorical, pl.Enum)
                       for dataset in datasets[1:]) and
                   not all(dataset._var[column].dtype == dtype
                           for dataset in datasets[1:])}
            var_columns_in_common = [
                column
                for column, dtype in islice(self._var.schema.items(), 1, None)
                if column in var_mismatched_categoricals or
                   all(column in dataset._var and
                       dataset._var[column].dtype == dtype
                       for dataset in datasets[1:])]
            cast_dict = {column: pl.Enum(
                pl.concat([dataset._var[column].cat.get_categories()
                           for dataset in datasets])
                .unique(maintain_order=True))
                for column in var_mismatched_categoricals}
            for dataset in datasets:
                dataset._var = dataset._var\
                    .select(dataset.var_names, *var_columns_in_common)\
                    .cast(cast_dict)
            varm_keys_in_common = [
                key for key, value in self._varm.items()
                if all(key in dataset._varm and
                       type(dataset._varm[key]) is type(value) and
                       dataset._varm[key].dtype == value.dtype
                       for dataset in datasets[1:])]
            for dataset in datasets:
                dataset._varm = {key: dataset._varm[key]
                                  for key in varm_keys_in_common}
        else:  # non-flexible
            # Check that all obs, obsm and uns are identical
            obs = self._obs
            for dataset in datasets[1:]:
                if not dataset._obs.equals(obs):
                    error_message = (
                        'all SingleCell datasets must have the same obs, '
                        'unless flexible=True')
                    raise ValueError(error_message)
            obsm = self._obsm
            for dataset in datasets[1:]:
                if dataset._obsm.keys() != obsm.keys() or \
                        any(type(dataset._obsm[key]) is not type(value) or
                            dataset._obsm[key].dtype != value.dtype
                            for key, value in obsm.items()) or \
                        any((dataset._obsm[key] != value).any()
                            for key, value in obsm.items()):
                    error_message = (
                        'all SingleCell datasets must have the same obsm, '
                        'unless flexible=True')
                    raise ValueError(error_message)
            for dataset in datasets[1:]:
                if not SingleCell._eq_uns(self._uns, dataset._uns):
                    error_message = (
                        'all SingleCell datasets must have the same uns, '
                        'unless flexible=True')
                    raise ValueError(error_message)
            # Check that all var have the same columns and data types
            schema = self._var.schema
            for dataset in datasets[1:]:
                if dataset._var.schema != schema:
                    error_message = (
                        'all SingleCell datasets must have the same columns '
                        'in var, with the same data types, unless '
                        'flexible=True')
                    raise ValueError(error_message)
            # Check that all varm have the same keys and data types
            varm = self._varm
            for dataset in datasets[1:]:
                if dataset._varm.keys() != varm.keys() or any(
                        type(dataset._varm[key]) is not type(value) or
                        dataset._varm[key].dtype != value.dtype
                        for key, value in varm.items()):
                    error_message = (
                        'all SingleCell datasets must have the same keys in '
                        'varm, with the same data types, unless flexible=True')
                    raise ValueError(error_message)
        # Concatenate; output should be CSR when there's a mix of inputs
        X = hstack([dataset._X for dataset in datasets],
                   format='csr' if any(isinstance(dataset._X, csr_array)
                                       for dataset in datasets) else 'csc')
        var = pl.concat([dataset._var for dataset in datasets])
        # noinspection PyTypeChecker
        varm = {key: np.concatenate([dataset._varm[key]
                                     for dataset in datasets])
                if isinstance(value, np.ndarray) else
                pl.concat([dataset._varm[key] for dataset in datasets])
                for key, value in self._varm.items()}
        return SingleCell(X=X, obs=datasets[0]._obs, var=var,
                          obsm=datasets[0]._obsm, varm=varm,
                          uns=datasets[0]._uns)

    def split_by_obs(self,
                     column: SingleCellColumn,
                     *,
                     QC_column: SingleCellColumn | None = 'passed_QC',
                     sort: bool = False) -> tuple[SingleCell, ...]:
        """
        The opposite of concat_obs(): splits a SingleCell dataset into a tuple
        of SingleCell datasets, one per unique value of a column of obs.

        Args:
            column: a String, Categorical or Enum column of obs to split by.
                    Can be a column name, a polars expression, a polars Series,
                    a 1D NumPy array, or a function that takes in this
                    SingleCell dataset and returns a polars Series or 1D NumPy
                    array. Can contain null entries: the corresponding cells
                    will not be included in the result.
            QC_column: an optional Boolean column of obs indicating which cells
                       passed QC. Can be a column name, a polars expression, a
                       polars Series, a 1D NumPy array, or a function that
                       takes in this SingleCell dataset and returns a polars
                       Series or 1D NumPy array. Set to None to include all
                       cells. Cells failing QC will not be selected when
                       subsampling.
            sort: if True, sort the SingleCell datasets in the returned tuple
                  in decreasing size. If False, sort in order of each value's
                  first appearance in `column`.
        
        Returns:
            A tuple of SingleCell datasets, one per unique value of `column`.
        """
        if QC_column is not None:
            QC_column = self._get_column(
                'obs', QC_column, 'QC_column', pl.Boolean,
                allow_missing=QC_column == 'passed_QC')
        column = self._get_column('obs', column, 'column',
                                  (pl.String, pl.Categorical, pl.Enum),
                                  QC_column=QC_column, allow_null=True)
        check_type(sort, 'sort', pl.Boolean, 'Boolean')
        values = (column.value_counts(sort=True).to_series().drop_nulls()
                  if sort else column.unique(maintain_order=True))
        if QC_column is None:
            return tuple(self.filter_obs(column == value) for value in values)
        else:
            return tuple(self.filter_obs(column.eq(value) & QC_column)
                         for value in values)
    
    def split_by_var(self,
                     column: SingleCellColumn,
                     *,
                     sort: bool = False) -> tuple[SingleCell, ...]:
        """
        The opposite of concat_var(): splits a SingleCell dataset into a tuple
        of SingleCell datasets, one per unique value of a column of var.

        Args:
            column: a String, Categorical or Enum column of var to split by.
                    Can be a column name, a polars expression, a polars Series,
                    a 1D NumPy array, or a function that takes in this
                    SingleCell dataset and returns a polars Series or 1D NumPy
                    array. Can contain null entries: the corresponding genes
                    will not be included in the result.
            sort: if True, sort the SingleCell datasets in the returned tuple
                  in decreasing size. If False, sort in order of each value's
                  first appearance in `column`.
        
        Returns:
            A tuple of SingleCell datasets, one per unique value of `column`.
        """
        column = self._get_column('var', column, 'column',
                                  (pl.String, pl.Categorical, pl.Enum),
                                  allow_null=True)
        check_type(sort, 'sort', pl.Boolean, 'Boolean')
        return tuple(self.filter_var(column == value) for value in
                     (column.value_counts(sort=True).to_series().drop_nulls()
                      if sort else column.unique(maintain_order=True)))
    
    def tocsr(self) -> SingleCell:
        """
        Make a copy of this SingleCell dataset, converting X to a csr_array.
        Raise an error if X is already a csr_array.
        
        Returns:
            A copy of this SingleCell dataset, with X as a csr_array.
        """
        if isinstance(self._X, csr_array):
            error_message = 'X is already a csr_array'
            raise TypeError(error_message)
        return SingleCell(X=self._X.tocsr(), obs=self._obs, var=self._var,
                          obsm=self._obsm, varm=self._varm, uns=self._uns)
    
    def tocsc(self) -> SingleCell:
        """
        Make a copy of this SingleCell dataset, converting X to a csc_array.
        Raise an error if X is already a csc_array.
        
        This function is provided for completeness, but csc_array is a far
        better format for cell-wise operations like pseudobulking.
        
        Returns:
            A copy of this SingleCell dataset, with X as a csc_array.
        """
        if isinstance(self._X, csc_array):
            error_message = 'X is already a csc_array'
            raise TypeError(error_message)
        return SingleCell(X=self._X.tocsc(), obs=self._obs, var=self._var,
                          obsm=self._obsm, varm=self._varm, uns=self._uns)
    
    def filter_obs(self,
                   *predicates: pl.Expr | pl.Series | str |
                                Iterable[pl.Expr | pl.Series | str] | bool |
                                list[bool] | np.ndarray[1, np.bool_],
                   **constraints: Any) -> SingleCell:
        """
        Equivalent to `df.filter()` from polars, but applied to both obs/obsm
        and X.
        
        Args:
            *predicates: one or more column names, expressions that evaluate to
                         Boolean Series, Boolean Series, lists of Booleans,
                         and/or 1D Boolean NumPy arrays
            **constraints: column filters: `name=value` filters to cells
                           where the column named `name` has the value `value`
        
        Returns:
            A new SingleCell dataset filtered to cells passing all the
            Boolean filters in `predicates` and `constraints`.
        """
        obs = self._obs\
            .with_columns(_SingleCell_index=pl.int_range(pl.len(),
                                                         dtype=pl.Int32))\
            .filter(*predicates, **constraints)
        mask = obs['_SingleCell_index'].to_numpy()
        return SingleCell(X=self._X[mask],
                          obs=obs.drop('_SingleCell_index'), var=self._var,
                          obsm={key: value[mask]
                                for key, value in self._obsm.items()},
                          varm=self._varm, uns=self._uns)
    
    def filter_var(self,
                   *predicates: pl.Expr | pl.Series | str |
                                Iterable[pl.Expr | pl.Series | str] | bool |
                                list[bool] | np.ndarray[1, np.bool_],
                   **constraints: Any) -> SingleCell:
        """
        Equivalent to `df.filter()` from polars, but applied to both var/varm
        and X.
        
        Args:
            *predicates: one or more column names, expressions that evaluate to
                         Boolean Series, Boolean Series, lists of Booleans,
                         and/or 1D Boolean NumPy arrays
            **constraints: column filters: `name=value` filters to genes
                           where the column named `name` has the value `value`
        
        Returns:
            A new SingleCell dataset filtered to genes passing all the
            Boolean filters in `predicates` and `constraints`.
        """
        var = self._var\
            .with_columns(_SingleCell_index=pl.int_range(pl.len(),
                                                         dtype=pl.Int32))\
            .filter(*predicates, **constraints)
        return SingleCell(X=self._X[:, var['_SingleCell_index'].to_numpy()],
                          obs=self._obs, var=var.drop('_SingleCell_index'),
                          obsm=self._obsm, varm=self._varm, uns=self._uns)
    
    def select_obs(self,
                   *exprs: Scalar | pl.Expr | pl.Series |
                           Iterable[Scalar | pl.Expr | pl.Series],
                   **named_exprs: Scalar | pl.Expr | pl.Series) -> SingleCell:
        """
        Equivalent to `df.select()` from polars, but applied to obs. obs_names
        will be automatically included as the first column, if not included
        explicitly.
        
        Args:
            *exprs: column(s) to select, specified as positional arguments.
                    Accepts expression input. Strings are parsed as column
                    names, other non-expression inputs are parsed as literals.
            **named_exprs: additional columns to select, specified as keyword
                           arguments. The columns will be renamed to the
                           keyword used.
        
        Returns:
            A new SingleCell dataset with
            `obs=obs.select(*exprs, **named_exprs)`, and obs_names as the first
            column unless already included explicitly.
        """
        obs = self._obs.select(*exprs, **named_exprs)
        if self.obs_names.name not in obs:
            obs = obs.select(self.obs_names, pl.all())
        return SingleCell(X=self._X, obs=obs, var=self._var, obsm=self._obsm,
                          varm=self._varm, uns=self._uns)
    
    def select_var(self,
                   *exprs: Scalar | pl.Expr | pl.Series |
                           Iterable[Scalar | pl.Expr | pl.Series],
                   **named_exprs: Scalar | pl.Expr | pl.Series) -> SingleCell:
        """
        Equivalent to `df.select()` from polars, but applied to var. var_names
        will be automatically included as the first column, if not included
        explicitly.
        
        Args:
            *exprs: column(s) to select, specified as positional arguments.
                    Accepts expression input. Strings are parsed as column
                    names, other non-expression inputs are parsed as literals.
            **named_exprs: additional columns to select, specified as keyword
                           arguments. The columns will be renamed to the
                           keyword used.
        
        Returns:
            A new SingleCell dataset with
            `var=var.select(*exprs, **named_exprs)`, and var_names as the first
            column unless already included explicitly.
        """
        var = self._var.select(*exprs, **named_exprs)
        if self.var_names.name not in var:
            var = var.select(self.var_names, pl.all())
        return SingleCell(X=self._X, obs=self._obs, var=var, obsm=self._obsm,
                          varm=self._varm, uns=self._uns)

    def select_obsm(self, keys: str | Iterable[str], *more_keys: str) -> \
            SingleCell:
        """
        Subsets obsm to the specified key(s).
        
        Args:
            keys: key(s) to select
            *more_keys: additional keys to select, specified as positional
                        arguments
        
        Returns:
            A new SingleCell dataset with obsm subset to the specified key(s).
        """
        keys = to_tuple(keys) + more_keys
        check_types(keys, 'keys', str, 'strings')
        for key in keys:
            if key not in self._obsm:
                error_message = \
                    f'tried to select {key!r}, which is not a key of obsm'
                raise ValueError(error_message)
        return SingleCell(X=self._X, obs=self._obs, var=self._var,
                          obsm={key: value for key, value in self._obsm.items()
                                if key in keys},
                          varm=self._varm, uns=self._uns)
    
    def select_varm(self, keys: str | Iterable[str], *more_keys: str) -> \
            SingleCell:
        """
        Subsets varm to the specified key(s).
        
        Args:
            keys: key(s) to select
            *more_keys: additional keys to select, specified as positional
                        arguments
        
        Returns:
            A new SingleCell dataset with varm subset to the specified key(s).
        """
        keys = to_tuple(keys) + more_keys
        check_types(keys, 'keys', str, 'strings')
        for key in keys:
            if key not in self._varm:
                error_message = \
                    f'tried to select {key!r}, which is not a key of varm'
                raise ValueError(error_message)
        return SingleCell(X=self._X, obs=self._obs, var=self._var,
                          obsm=self._obsm,
                          varm={key: value for key, value in self._varm.items()
                                if key in keys},
                          uns=self._uns)
    
    def select_uns(self, keys: str | Iterable[str], *more_keys: str) -> \
            SingleCell:
        """
        Subsets uns to the specified key(s).
        
        Args:
            keys: key(s) to select
            *more_keys: additional keys to select, specified as positional
                        arguments
        
        Returns:
            A new SingleCell dataset with uns subset to the specified key(s).
        """
        keys = to_tuple(keys) + more_keys
        check_types(keys, 'keys', str, 'strings')
        for key in keys:
            if key not in self._uns:
                error_message = \
                    f'tried to select {key!r}, which is not a key of uns'
                raise ValueError(error_message)
        return SingleCell(X=self._X, obs=self._obs, var=self._var,
                          obsm=self._obsm, varm=self._varm,
                          uns={key: value for key, value in self._uns.items()
                                if key in keys})
    
    def with_columns_obs(self,
                         *exprs: Scalar | pl.Expr | pl.Series |
                                 Iterable[Scalar | pl.Expr | pl.Series],
                         **named_exprs: Scalar | pl.Expr | pl.Series) -> \
            SingleCell:
        """
        Equivalent to `df.with_columns()` from polars, but applied to obs.
        
        Args:
            *exprs: column(s) to add, specified as positional arguments.
                    Accepts expression input. Strings are parsed as column
                    names, other non-expression inputs are parsed as literals.
            **named_exprs: additional columns to add, specified as keyword
                           arguments. The columns will be renamed to the
                           keyword used.
        
        Returns:
            A new SingleCell dataset with
            `obs=obs.with_columns(*exprs, **named_exprs)`.
        """
        # noinspection PyTypeChecker
        return SingleCell(X=self._X,
                          obs=self._obs.with_columns(*exprs, **named_exprs),
                          var=self._var, obsm=self._obsm, varm=self._varm,
                          uns=self._uns)
    
    def with_columns_var(self,
                         *exprs: Scalar | pl.Expr | pl.Series |
                                 Iterable[Scalar | pl.Expr | pl.Series],
                         **named_exprs: Scalar | pl.Expr | pl.Series) -> \
            SingleCell:
        """
        Equivalent to `df.with_columns()` from polars, but applied to var.
        
        Args:
            *exprs: column(s) to add, specified as positional arguments.
                    Accepts expression input. Strings are parsed as column
                    names, other non-expression inputs are parsed as literals.
            **named_exprs: additional columns to add, specified as keyword
                           arguments. The columns will be renamed to the
                           keyword used.
        
        Returns:
            A new SingleCell dataset with
            `var=var.with_columns(*exprs, **named_exprs)`.
        """
        # noinspection PyTypeChecker
        return SingleCell(X=self._X, obs=self._obs,
                          var=self._var.with_columns(*exprs, **named_exprs),
                          obsm=self._obsm, varm=self._varm, uns=self._uns)
    
    def with_obsm(self,
                  obsm: dict[str, np.ndarray[2, Any] | pl.DataFrame] = {},
                  **more_obsm: np.ndarray[2, Any]) -> SingleCell:
        """
        Adds one or more keys to obsm, overwriting existing keys with the same
        names if present.
        
        Args:
            obsm: a dictionary of keys to add to (or overwrite in) obsm
            **more_obsm: additional keys to add to (or overwrite in) obsm,
                         specified as keyword arguments

        Returns:
            A new SingleCell dataset with the new key(s) added to or
            overwritten in obsm.
        """
        check_type(obsm, 'obsm', dict, 'a dictionary')
        for key, value in obsm.items():
            if not isinstance(key, str):
                error_message = (
                    f'all keys of obsm must be strings, but it contains a key '
                    f'of type {type(key).__name__!r}')
                raise TypeError(error_message)
        obsm |= more_obsm
        if len(obsm) == 0:
            error_message = \
                'obsm is empty and no keyword arguments were specified'
            raise ValueError(error_message)
        for key, value in obsm.items():
            if isinstance(value, np.ndarray):
                if value.ndim != 2:
                    error_message = (
                        f'all values of obsm must be 2D NumPy arrays or '
                        f'polars DataFrames, but new obsm[{key!r}] is a '
                        f'{value.ndim:,}D NumPy array')
                    raise ValueError(error_message)
            elif not isinstance(value, pl.DataFrame):
                error_message = (
                    f'all values of obsm must be NumPy arrays or polars '
                    f'DataFrames, but new obsm[{key!r}] has type '
                    f'{type(value).__name__!r}')
                raise TypeError(error_message)
            if len(value) != self._X.shape[0]:
                error_message = (
                    f'len(obsm[{key!r}]) is {len(value):,}, but X.shape[0] is '
                    f'{self._X.shape[0]:,}')
                raise ValueError(error_message)
        return SingleCell(X=self._X, obs=self._obs, var=self._var,
                          obsm=self._obsm | obsm, varm=self._varm,
                          uns=self._uns)
    
    def with_varm(self,
                  varm: dict[str, np.ndarray[2, Any] | pl.DataFrame] = {},
                  **more_varm: np.ndarray[2, Any]) -> SingleCell:
        """
        Adds one or more keys to varm, overwriting existing keys with the same
        names if present.
        
        Args:
            varm: a dictionary of keys to add to (or overwrite in) varm
            **more_varm: additional keys to add to (or overwrite in) varm,
                         specified as keyword arguments

        Returns:
            A new SingleCell dataset with the new key(s) added to or
            overwritten in varm.
        """
        check_type(varm, 'varm', dict, 'a dictionary')
        for key, value in varm.items():
            if not isinstance(key, str):
                error_message = (
                    f'all keys of varm must be strings, but it contains a key '
                    f'of type {type(key).__name__!r}')
                raise TypeError(error_message)
        varm |= more_varm
        if len(varm) == 0:
            error_message = \
                'varm is empty and no keyword arguments were specified'
            raise ValueError(error_message)
        for key, value in varm.items():
            if isinstance(value, np.ndarray):
                if value.ndim != 2:
                    error_message = (
                        f'all values of varm must be 2D NumPy arrays or '
                        f'polars DataFrames, but new varm[{key!r}] is a '
                        f'{value.ndim:,}D NumPy array')
                    raise ValueError(error_message)
            elif not isinstance(value, pl.DataFrame):
                error_message = (
                    f'all values of varm must be NumPy arrays or polars '
                    f'DataFrames, but new varm[{key!r}] has type '
                    f'{type(value).__name__!r}')
                raise TypeError(error_message)
            if len(value) != self._X.shape[0]:
                error_message = (
                    f'len(varm[{key!r}]) is {len(value):,}, but X.shape[0] is '
                    f'{self._X.shape[0]:,}')
                raise ValueError(error_message)
        return SingleCell(X=self._X, obs=self._obs, var=self._var,
                          obsm=self._obsm, varm=self._varm | varm,
                          uns=self._uns)
    
    def with_uns(self,
                 uns: dict[str, NestedScalarOrArrayDict] = {},
                  **more_uns: NestedScalarOrArrayDict) -> SingleCell:
        """
        Adds one or more keys to uns, overwriting existing keys with the same
        names if present.
        
        Args:
            uns: a dictionary of keys to add to (or overwrite in) uns
            **more_uns: additional keys to add to (or overwrite in) uns,
                        specified as keyword arguments

        Returns:
            A new SingleCell dataset with the new key(s) added to or
            overwritten in uns.
        """
        check_type(uns, 'uns', dict, 'a dictionary')
        for key, value in uns.items():
            if not isinstance(key, str):
                error_message = (
                    f'all keys of uns must be strings, but it contains a key '
                    f'of type {type(key).__name__!r}')
                raise TypeError(error_message)
        uns |= more_uns
        if len(uns) == 0:
            error_message = \
                'uns is empty and no keyword arguments were specified'
            raise ValueError(error_message)
        valid_uns_types = str, int, np.integer, float, np.floating, \
            bool, np.bool_, np.ndarray
        for description, value in SingleCell._iter_uns(uns):
            if not isinstance(value, valid_uns_types):
                error_message = (
                    f'all values of uns must be scalars (strings, numbers or '
                    f'Booleans) or NumPy arrays, or nested dictionaries '
                    f'thereof, but {description} has type '
                    f'{type(value).__name__!r}')
                raise TypeError(error_message)
        return SingleCell(X=self._X, obs=self._obs, var=self._var,
                          obsm=self._obsm, varm=self._varm,
                          uns=self._uns | uns)
    
    def drop_obs(self,
                 columns: pl.type_aliases.ColumnNameOrSelector |
                          Iterable[pl.type_aliases.ColumnNameOrSelector],
                 *more_columns: pl.type_aliases.ColumnNameOrSelector) -> \
            SingleCell:
        """
        Create a new SingleCell dataset with `columns` and `more_columns`
        removed from obs.
        
        Args:
            columns: columns(s) to drop
            *more_columns: additional columns to drop, specified as
                              positional arguments
        
        Returns:
            A new SingleCell dataset with the column(s) removed.
        """
        columns = to_tuple(columns) + more_columns
        return SingleCell(X=self._X, obs=self._obs.drop(columns),
                          var=self._var, obsm=self._obsm, varm=self._varm,
                          uns=self._uns)

    def drop_var(self,
                 columns: pl.type_aliases.ColumnNameOrSelector |
                          Iterable[pl.type_aliases.ColumnNameOrSelector],
                 *more_columns: pl.type_aliases.ColumnNameOrSelector) -> \
            SingleCell:
        """
        Create a new SingleCell dataset with `columns` and `more_columns`
        removed from var.
        
        Args:
            columns: columns(s) to drop
            *more_columns: additional columns to drop, specified as
                           positional arguments
        
        Returns:
            A new SingleCell dataset with the column(s) removed.
        """
        columns = to_tuple(columns) + more_columns
        return SingleCell(X=self._X, obs=self._obs,
                          var=self._var.drop(columns), obsm=self._obsm,
                          varm=self._varm, uns=self._uns)
    
    def drop_obsm(self, keys: str | Iterable[str], *more_keys: str) -> \
            SingleCell:
        """
        Create a new SingleCell dataset with `keys` and `more_keys` removed
        from obsm.
        
        Args:
            keys: key(s) to drop
            *more_keys: additional keys to drop, specified as positional
                        arguments
        
        Returns:
            A new SingleCell dataset with the specified key(s) removed from
            obsm.
        """
        keys = to_tuple(keys) + more_keys
        check_types(keys, 'keys', str, 'strings')
        for key in keys:
            if key not in self._obsm:
                error_message = \
                    f'tried to drop {key!r}, which is not a key of obsm'
                raise ValueError(error_message)
        return SingleCell(X=self._X, obs=self._obs, var=self._var,
                          obsm={key: value for key, value in self._obsm.items()
                                if key not in keys},
                          varm=self._varm, uns=self._uns)
    
    def drop_varm(self, keys: str | Iterable[str], *more_keys: str) -> \
            SingleCell:
        """
        Create a new SingleCell dataset with `keys` and `more_keys` removed
        from varm.
        
        Args:
            keys: key(s) to drop
            *more_keys: additional keys to drop, specified as positional
                        arguments
        
        Returns:
            A new SingleCell dataset with the specified key(s) removed from
            varm.
        """
        keys = to_tuple(keys) + more_keys
        check_types(keys, 'keys', str, 'strings')
        for key in keys:
            if key not in self._varm:
                error_message = \
                    f'tried to drop {key!r}, which is not a key of varm'
                raise ValueError(error_message)
        return SingleCell(X=self._X, obs=self._obs, var=self._var,
                          obsm=self._obsm,
                          varm={key: value for key, value in self._varm.items()
                                if key not in keys},
                          uns=self._uns)
    
    def drop_uns(self, keys: str | Iterable[str], *more_keys: str) -> \
            SingleCell:
        """
        Create a new SingleCell dataset with `keys` and `more_keys` removed
        from uns.
        
        Args:
            keys: key(s) to drop
            *more_keys: additional keys to drop, specified as positional
                        arguments
        
        Returns:
            A new SingleCell dataset with the specified key(s) removed from
            uns.
        """
        keys = to_tuple(keys) + more_keys
        check_types(keys, 'keys', str, 'strings')
        for key in keys:
            if key not in self._uns:
                error_message = \
                    f'tried to drop {key!r}, which is not a key of uns'
                raise ValueError(error_message)
        return SingleCell(X=self._X, obs=self._obs, var=self._var,
                          obsm=self._obsm, varm=self._varm,
                          uns={key: value for key, value in self._uns.items()
                                if key not in keys})
    
    def rename_obs(self, mapping: dict[str, str] | Callable[[str], str]) -> \
            SingleCell:
        """
        Create a new SingleCell dataset with column(s) of obs renamed.
        
        Rename column(s) of obs.
        
        Args:
            mapping: the renaming to apply, either as a dictionary with the old
                     names as keys and the new names as values, or a function
                     that takes an old name and returns a new name

        Returns:
            A new SingleCell dataset with the column(s) of obs renamed.
        """
        return SingleCell(X=self._X, obs=self._obs.rename(mapping),
                          var=self._var, obsm=self._obsm, varm=self._varm,
                          uns=self._uns)
    
    def rename_var(self, mapping: dict[str, str] | Callable[[str], str]) -> \
            SingleCell:
        """
        Create a new SingleCell dataset with column(s) of var renamed.
        
        Args:
            mapping: the renaming to apply, either as a dictionary with the old
                     names as keys and the new names as values, or a function
                     that takes an old name and returns a new name

        Returns:
            A new SingleCell dataset with the column(s) of var renamed.
        """
        return SingleCell(X=self._X, obs=self._obs,
                          var=self._var.rename(mapping), obsm=self._obsm,
                          varm=self._varm, uns=self._uns)
    
    def rename_obsm(self, mapping: dict[str, str] | Callable[[str], str]) -> \
            SingleCell:
        """
        Create a new SingleCell dataset with key(s) of obsm renamed.
        
        Args:
            mapping: the renaming to apply, either as a dictionary with the old
                     names as keys and the new names as values, or a function
                     that takes an old name and returns a new name

        Returns:
            A new SingleCell dataset with the key(s) of obsm renamed.
        """
        check_types(mapping.keys(), 'mapping.keys()', str, 'strings')
        check_types(mapping.values(), 'mapping.values()', str, 'strings')
        if isinstance(mapping, dict):
            for key, new_key in mapping.items():
                if key not in self._obsm:
                    error_message = \
                        f'tried to rename {key!r}, which is not a key of obsm'
                    raise ValueError(error_message)
                if new_key in self._obsm:
                    error_message = (
                        f'tried to rename obsm[{key!r}] to obsm[{new_key!r}], '
                        f'but obsm[{new_key!r}] already exists')
                    raise ValueError(error_message)
            obsm = {mapping.get(key, key): value
                    for key, value in self._obsm.items()}
        elif isinstance(mapping, Callable):
            obsm = {}
            for key, value in self._obsm.items():
                new_key = mapping(key)
                if not isinstance(new_key, str):
                    error_message = (
                        f'tried to rename obsm[{key!r}] to a non-string value '
                        f'of type {type(new_key).__name__!r}')
                    raise TypeError(error_message)
                if new_key in self._obsm:
                    error_message = (
                        f'tried to rename obsm[{key!r}] to obsm[{new_key!r}], '
                        f'but obsm[{new_key!r}] already exists')
                    raise ValueError(error_message)
                obsm[new_key] = value
        else:
            error_message = (
                f'mapping must be a dictionary or function, but has type '
                f'{type(mapping).__name__!r}')
            raise TypeError(error_message)
        return SingleCell(X=self._X, obs=self._obs, var=self._var, obsm=obsm,
                          varm=self._varm, uns=self._uns)
    
    def rename_varm(self, mapping: dict[str, str] | Callable[[str], str]) -> \
            SingleCell:
        """
        Create a new SingleCell dataset with key(s) of varm renamed.
        
        Args:
            mapping: the renaming to apply, either as a dictionary with the old
                     names as keys and the new names as values, or a function
                     that takes an old name and returns a new name

        Returns:
            A new SingleCell dataset with the key(s) of varm renamed.
        """
        check_types(mapping.keys(), 'mapping.keys()', str, 'strings')
        check_types(mapping.values(), 'mapping.values()', str, 'strings')
        if isinstance(mapping, dict):
            for key, new_key in mapping.items():
                if key not in self._varm:
                    error_message = \
                        f'tried to rename {key!r}, which is not a key of varm'
                    raise ValueError(error_message)
                if new_key in self._varm:
                    error_message = (
                        f'tried to rename varm[{key!r}] to varm[{new_key!r}], '
                        f'but varm[{new_key!r}] already exists')
                    raise ValueError(error_message)
            varm = {mapping.get(key, key): value
                    for key, value in self._varm.items()}
        elif isinstance(mapping, Callable):
            varm = {}
            for key, value in self._varm.items():
                new_key = mapping(key)
                if not isinstance(new_key, str):
                    error_message = (
                        f'tried to rename varm[{key!r}] to a non-string value '
                        f'of type {type(new_key).__name__!r}')
                    raise TypeError(error_message)
                if new_key in self._varm:
                    error_message = (
                        f'tried to rename varm[{key!r}] to varm[{new_key!r}], '
                        f'but varm[{new_key!r}] already exists')
                    raise ValueError(error_message)
                varm[new_key] = value
        else:
            error_message = (
                f'mapping must be a dictionary or function, but has type '
                f'{type(mapping).__name__!r}')
            raise TypeError(error_message)
        return SingleCell(X=self._X, obs=self._obs, var=self._var,
                          obsm=self._obsm, varm=varm, uns=self._uns)
    
    def rename_uns(self, mapping: dict[str, str] | Callable[[str], str]) -> \
            SingleCell:
        """
        Create a new SingleCell dataset with key(s) of uns renamed.
        
        Args:
            mapping: the renaming to apply, either as a dictionary with the old
                     names as keys and the new names as values, or a function
                     that takes an old name and returns a new name

        Returns:
            A new SingleCell dataset with the key(s) of uns renamed.
        """
        check_types(mapping.keys(), 'mapping.keys()', str, 'strings')
        check_types(mapping.values(), 'mapping.values()', str, 'strings')
        if isinstance(mapping, dict):
            for key, new_key in mapping.items():
                if key not in self._uns:
                    error_message = \
                        f'tried to rename {key!r}, which is not a key of uns'
                    raise ValueError(error_message)
                if new_key in self._uns:
                    error_message = (
                        f'tried to rename uns[{key!r}] to uns[{new_key!r}], '
                        f'but uns[{new_key!r}] already exists')
                    raise ValueError(error_message)
            uns = {mapping.get(key, key): value
                   for key, value in self._uns.items()}
        elif isinstance(mapping, Callable):
            uns = {}
            for key, value in self._uns.items():
                new_key = mapping(key)
                if not isinstance(new_key, str):
                    error_message = (
                        f'tried to rename uns[{key!r}] to a non-string value '
                        f'of type {type(new_key).__name__!r}')
                    raise TypeError(error_message)
                if new_key in self._uns:
                    error_message = (
                        f'tried to rename uns[{key!r}] to uns[{new_key!r}], '
                        f'but uns[{new_key!r}] already exists')
                    raise ValueError(error_message)
                uns[new_key] = value
        else:
            error_message = (
                f'mapping must be a dictionary or function, but has type '
                f'{type(mapping).__name__!r}')
            raise TypeError(error_message)
        return SingleCell(X=self._X, obs=self._obs, var=self._var,
                          obsm=self._obsm, varm=self._varm, uns=uns)
    
    def cast_X(self, dtype: np._typing.DTypeLike) -> SingleCell:
        """
        Cast X to the specified data type.
        
        Args:
            dtype: a NumPy data type

        Returns:
            A new SingleCell dataset with X cast to the specified data type.
        """
        return SingleCell(X=self._X.astype(dtype),
                          obs=self._obs, var=self._var, obsm=self._obsm,
                          varm=self._varm, uns=self._uns)
    
    def cast_obs(self,
                 dtypes: Mapping[pl.type_aliases.ColumnNameOrSelector |
                                 pl.type_aliases.PolarsDataType,
                                 pl.type_aliases.PolarsDataType] |
                         pl.type_aliases.PolarsDataType,
                 *,
                 strict: bool = True) -> SingleCell:
        """
        Cast column(s) of obs to the specified data type(s).
        
        Args:
            dtypes: a mapping of column names (or selectors) to data types, or
                    a single data type to which all columns will be cast
            strict: whether to raise an error if a cast could not be done (for
                    instance, due to numerical overflow)

        Returns:
            A new SingleCell dataset with column(s) of obs cast to the
            specified data type(s).
        """
        return SingleCell(X=self._X, obs=self._obs.cast(dtypes, strict=strict),
                          var=self._var, obsm=self._obsm, varm=self._varm,
                          uns=self._uns)
    
    def cast_var(self,
                 dtypes: Mapping[pl.type_aliases.ColumnNameOrSelector |
                                 pl.type_aliases.PolarsDataType,
                                 pl.type_aliases.PolarsDataType] |
                         pl.type_aliases.PolarsDataType,
                 *,
                 strict: bool = True) -> SingleCell:
        """
        Cast column(s) of var to the specified data type(s).
        
        Args:
            dtypes: a mapping of column names (or selectors) to data types, or
                    a single data type to which all columns will be cast
            strict: whether to raise an error if a cast could not be done (for
                    instance, due to numerical overflow)

        Returns:
            A new SingleCell dataset with column(s) of var cast to the
            specified data type(s).
        """
        return SingleCell(X=self._X, obs=self._obs,
                          var=self._var.cast(dtypes, strict=strict),
                          obsm=self._obsm, varm=self._varm, uns=self._uns)
    
    def join_obs(self,
                 other: pl.DataFrame,
                 on: str | pl.Expr | Sequence[str | pl.Expr] | None = None,
                 *,
                 left_on: str | pl.Expr | Sequence[str | pl.Expr] |
                          None = None,
                 right_on: str | pl.Expr | Sequence[str | pl.Expr] |
                           None = None,
                 suffix: str = '_right',
                 validate: Literal['m:m', 'm:1', '1:m', '1:1'] = 'm:m',
                 join_nulls: bool = False,
                 coalesce: bool = True) -> SingleCell:
        """
        Left join obs with another DataFrame.
        
        Args:
            other: a polars DataFrame to join obs with
            on: the name(s) of the join column(s) in both DataFrames
            left_on: the name(s) of the join column(s) in obs
            right_on: the name(s) of the join column(s) in `other`
            suffix: a suffix to append to columns with a duplicate name
            validate: checks whether the join is of the specified type. Can be:
                      - 'm:m' (many-to-many): the default, no checks performed.
                      - '1:1' (one-to-one): check that none of the values in
                        the join column(s) appear more than once in obs or more
                        than once in `other`.
                      - '1:m' (one-to-many): check that none of the values in
                        the join column(s) appear more than once in obs.
                      - 'm:1' (many-to-one): check that none of the values in
                        the join column(s) appear more than once in `other`.
            join_nulls: whether to include null as a valid value to join on.
                        By default, null values will never produce matches.
            coalesce: if True, coalesce each of the pairs of join columns
                      (the columns in `on` or `left_on`/`right_on`) from obs
                      and `other` into a single column, filling missing values
                      from one with the corresponding values from the other.
                      If False, include both as separate columns, adding
                      `suffix` to the join columns from `other`.
        
        Returns:
            A new SingleCell dataset with the columns from `other` joined to
            obs.
        
        Note:
            If a column of `on`, `left_on` or `right_on` is Enum in obs and
            Categorical in `other` (or vice versa), or Enum in both but with
            different categories in each, that pair of columns will be
            automatically cast to a common Enum data type (with the union of
            the categories) before joining.
        """
        # noinspection PyTypeChecker
        check_type(other, 'other', pl.DataFrame, 'a polars DataFrame')
        left = self._obs
        right = other
        if on is None:
            if left_on is None and right_on is None:
                error_message = (
                    "either 'on' or both of 'left_on' and 'right_on' must be "
                    "specified")
                raise ValueError(error_message)
            elif left_on is None:
                error_message = \
                    'right_on is specified, so left_on must be specified'
                raise ValueError(error_message)
            elif right_on is None:
                error_message = \
                    'left_on is specified, so right_on must be specified'
                raise ValueError(error_message)
            left_columns = left.select(left_on)
            right_columns = right.select(right_on)
        else:
            if left_on is not None:
                error_message = "'on' is specified, so 'left_on' must be None"
                raise ValueError(error_message)
            if right_on is not None:
                error_message = "'on' is specified, so 'right_on' must be None"
                raise ValueError(error_message)
            left_columns = left.select(on)
            right_columns = right.select(on)
        left_cast_dict = {}
        right_cast_dict = {}
        for left_column, right_column in zip(left_columns, right_columns):
            left_dtype = left_column.dtype
            right_dtype = right_column.dtype
            if left_dtype == right_dtype:
                continue
            if (left_dtype == pl.Enum or left_dtype == pl.Categorical) and (
                    right_dtype == pl.Enum or right_dtype == pl.Categorical):
                common_dtype = \
                    pl.Enum(pl.concat([left_column.cat.get_categories(),
                                       right_column.cat.get_categories()])
                            .unique(maintain_order=True))
                left_cast_dict[left_column.name] = common_dtype
                right_cast_dict[right_column.name] = common_dtype
            else:
                error_message = (
                    f'obs[{left_column.name!r}] has data type '
                    f'{left_dtype.base_type()!r}, but '
                    f'other[{right_column.name!r}] has data type '
                    f'{right_dtype.base_type()!r}')
                raise TypeError(error_message)
        if left_cast_dict is not None:
            left = left.cast(left_cast_dict)
            right = right.cast(right_cast_dict)
        obs = left.join(right, on=on, how='left', left_on=left_on,
                        right_on=right_on, suffix=suffix, validate=validate,
                        join_nulls=join_nulls, coalesce=coalesce)
        if len(obs) > len(self):
            other_on = to_tuple(right_on if right_on is not None else on)
            assert other.select(other_on).is_duplicated().any()
            duplicate_column = other_on[0] if len(other_on) == 1 else \
                next(column for column in other_on
                     if other[column].is_duplicated().any())
            error_message = (
                f'other[{duplicate_column!r}] contains duplicate values, so '
                f'it must be deduplicated before being joined on')
            raise ValueError(error_message)
        return SingleCell(X=self._X, obs=obs, var=self._var, obsm=self._obsm,
                          varm=self._varm, uns=self._uns)
    
    def join_var(self,
                 other: pl.DataFrame,
                 on: str | pl.Expr | Sequence[str | pl.Expr] | None = None,
                 *,
                 left_on: str | pl.Expr | Sequence[str | pl.Expr] |
                          None = None,
                 right_on: str | pl.Expr | Sequence[str | pl.Expr] |
                           None = None,
                 suffix: str = '_right',
                 validate: Literal['m:m', 'm:1', '1:m', '1:1'] = 'm:m',
                 join_nulls: bool = False,
                 coalesce: bool = True) -> SingleCell:
        """
        Join var with another DataFrame.
        
        Args:
            other: a polars DataFrame to join var with
            on: the name(s) of the join column(s) in both DataFrames
            left_on: the name(s) of the join column(s) in var
            right_on: the name(s) of the join column(s) in `other`
            suffix: a suffix to append to columns with a duplicate name
            validate: checks whether the join is of the specified type. Can be:
                      - 'm:m' (many-to-many): the default, no checks performed.
                      - '1:1' (one-to-one): check that none of the values in
                        the join column(s) appear more than once in var or more
                        than once in `other`.
                      - '1:m' (one-to-many): check that none of the values in
                        the join column(s) appear more than once in var.
                      - 'm:1' (many-to-one): check that none of the values in
                        the join column(s) appear more than once in `other`.
            join_nulls: whether to include null as a valid value to join on.
                        By default, null values will never produce matches.
            coalesce: if True, coalesce each of the pairs of join columns
                      (the columns in `on` or `left_on`/`right_on`) from obs
                      and `other` into a single column, filling missing values
                      from one with the corresponding values from the other.
                      If False, include both as separate columns, adding
                      `suffix` to the join columns from `other`.
        
        Returns:
            A new SingleCell dataset with the columns from `other` joined to
            var.
        
        Note:
            If a column of `on`, `left_on` or `right_on` is Enum in obs and
            Categorical in `other` (or vice versa), or Enum in both but with
            different categories in each, that pair of columns will be
            automatically cast to a common Enum data type (with the union of
            the categories) before joining.
        """
        check_type(other, 'other', pl.DataFrame, 'a polars DataFrame')
        left = self._var
        right = other
        if on is None:
            if left_on is None and right_on is None:
                error_message = (
                    "either 'on' or both of 'left_on' and 'right_on' must be "
                    "specified")
                raise ValueError(error_message)
            elif left_on is None:
                error_message = \
                    'right_on is specified, so left_on must be specified'
                raise ValueError(error_message)
            elif right_on is None:
                error_message = \
                    'left_on is specified, so right_on must be specified'
                raise ValueError(error_message)
            left_columns = left.select(left_on)
            right_columns = right.select(right_on)
        else:
            if left_on is not None:
                error_message = "'on' is specified, so 'left_on' must be None"
                raise ValueError(error_message)
            if right_on is not None:
                error_message = "'on' is specified, so 'right_on' must be None"
                raise ValueError(error_message)
            left_columns = left.select(on)
            right_columns = right.select(on)
        left_cast_dict = {}
        right_cast_dict = {}
        for left_column, right_column in zip(left_columns, right_columns):
            left_dtype = left_column.dtype
            right_dtype = right_column.dtype
            if left_dtype == right_dtype:
                continue
            if (left_dtype == pl.Enum or left_dtype == pl.Categorical) and (
                    right_dtype == pl.Enum or right_dtype == pl.Categorical):
                common_dtype = \
                    pl.Enum(pl.concat([left_column.cat.get_categories(),
                                       right_column.cat.get_categories()])
                            .unique(maintain_order=True))
                left_cast_dict[left_column.name] = common_dtype
                right_cast_dict[right_column.name] = common_dtype
            else:
                error_message = (
                    f'var[{left_column.name!r}] has data type '
                    f'{left_dtype.base_type()!r}, but '
                    f'other[{right_column.name!r}] has data type '
                    f'{right_dtype.base_type()!r}')
                raise TypeError(error_message)
        if left_cast_dict is not None:
            left = left.cast(left_cast_dict)
            right = right.cast(right_cast_dict)
        # noinspection PyTypeChecker
        var = left.join(right, on=on, how='left', left_on=left_on,
                        right_on=right_on, suffix=suffix, validate=validate,
                        join_nulls=join_nulls, coalesce=coalesce)
        if len(var) > len(self):
            other_on = to_tuple(right_on if right_on is not None else on)
            assert other.select(other_on).is_duplicated().any()
            duplicate_column = other_on[0] if len(other_on) == 1 else \
                next(column for column in other_on
                     if other[column].is_duplicated().any())
            error_message = (
                f'other[{duplicate_column!r}] contains duplicate values, so '
                f'it must be deduplicated before being joined on')
            raise ValueError(error_message)
        return SingleCell(X=self._X, obs=self._obs, var=var, obsm=self._obsm,
                          varm=self._varm, uns=self._uns)
    
    def peek_obs(self, row: int = 0) -> None:
        """
        Print a row of obs (the first row, by default) with each column on its
        own line.
        
        Args:
            row: the index of the row to print
        """
        check_type(row, 'row', int, 'an integer')
        with pl.Config(tbl_rows=-1):
            print(self._obs[row].unpivot(variable_name='column'))
    
    def peek_var(self, row: int = 0) -> None:
        """
        Print a row of var (the first row, by default) with each column on its
        own line.
        
        Args:
            row: the index of the row to print
        """
        check_type(row, 'row', int, 'an integer')
        with pl.Config(tbl_rows=-1):
            print(self._var[row].unpivot(variable_name='column'))
    
    def subsample_obs(self,
                      n: int | np.integer | None = None,
                      *,
                      fraction: int | float | np.integer | np.floating |
                                None = None,
                      QC_column: SingleCellColumn | None = 'passed_QC',
                      by_column: SingleCellColumn | None = None,
                      subsample_column: str | None = None,
                      seed: int | np.integer = 0,
                      overwrite: bool = False) -> SingleCell:
        """
        Subsample a specific number or fraction of cells.
        
        Args:
            n: the number of cells to return; mutually exclusive with
               `fraction`
            fraction: the fraction of cells to return; mutually exclusive with
                      `n`
            QC_column: an optional Boolean column of obs indicating which cells
                       passed QC. Can be a column name, a polars expression, a
                       polars Series, a 1D NumPy array, or a function that
                       takes in this SingleCell dataset and returns a polars
                       Series or 1D NumPy array. Set to None to include all
                       cells. Cells failing QC will not be selected when
                       subsampling, and will not count towards the denominator
                       of `fraction`; QC_column will not appear in the returned
                       SingleCell object, since it would be redundant.
            by_column: an optional String, Categorical, Enum, or integer column
                       of obs to subsample by. Can be a column name, a polars
                       expression, a polars Series, a 1D NumPy array, or a
                       function that takes in this SingleCell dataset and
                       returns a polars Series or 1D NumPy array. Specifying
                       `by_column` ensures that the same fraction of cells with
                       each value of `by_column` are subsampled. When combined
                       with `n`, to make sure the total number of samples is
                       exactly `n`, some of the smallest groups may be
                       oversampled by one element, or some of the largest
                       groups may be undersampled by one element. Can contain
                       null entries: the corresponding cells will not be
                       included in the result.
            subsample_column: an optional name of a Boolean column to add to
                              obs indicating the subsampled genes; if None,
                              subset to these genes instead
            seed: the random seed to use when subsampling
            overwrite: if True, overwrite `subsample_column` if already present
                       in obs, instead of raising an error. Must be False when
                       `subsample_column` is None.
        
        Returns:
            A new SingleCell dataset subset to the subsampled cells, or if
            `subsample_column` is not None, the full dataset with
            `subsample_column` added to obs.
        """
        if n is not None:
            check_type(n, 'n', int, 'a positive integer')
            check_bounds(n, 'n', 1)
        elif fraction is not None:
            check_type(fraction, 'fraction', float,
                       'a floating-point number between 0 and 1')
            check_bounds(fraction, 'fraction', 0, 1, left_open=True,
                         right_open=True)
        else:
            error_message = 'one of n and fraction must be specified'
            raise ValueError(error_message)
        if n is not None and fraction is not None:
            error_message = 'only one of n and fraction must be specified'
            raise ValueError(error_message)
        if QC_column is not None:
            QC_column = self._get_column(
                'obs', QC_column, 'QC_column', pl.Boolean,
                allow_missing=QC_column == 'passed_QC')
        check_type(overwrite, 'overwrite', bool, 'Boolean')
        if subsample_column is not None:
            check_type(subsample_column, 'subsample_column', str, 'a string')
            if not overwrite and subsample_column in self._obs:
                error_message = (
                    f'subsample_column {subsample_column!r} is already a '
                    f'column of obs; did you already run subsample_obs()? Set '
                    f'overwrite=True to overwrite.')
                raise ValueError(error_message)
        elif overwrite:
            error_message = \
                'overwrite must be False when subsample_column is None'
            raise ValueError(error_message)
        check_type(seed, 'seed', int, 'an integer')
        if by_column is not None:
            by_column = self._get_column(
                'obs', by_column, 'by_column',
                (pl.String, pl.Categorical, pl.Enum, 'integer'),
                QC_column=QC_column, allow_null=True)
            if QC_column is not None:
                by_column = by_column.filter(QC_column)
            by_frame = by_column.to_frame()
            by_name = by_column.name
            if n is not None:
                # Get a vector of the number of elements to sample per group.
                # The total sample size should exactly match the original n; if
                # necessary, oversample the smallest groups or undersample the
                # largest groups to make this happen.
                group_counts = by_frame\
                    .group_by(by_name)\
                    .agg(pl.len(), n=(n * pl.len() / len(by_column))
                                     .round().cast(pl.Int32))\
                    .drop_nulls(by_name)
                diff = n - group_counts['n'].sum()
                if diff != 0:
                    group_counts = group_counts\
                        .sort('len', descending=diff < 0)\
                        .with_columns(n=pl.col.n +
                                        pl.int_range(pl.len(), dtype=pl.Int32)
                                        .lt(abs(diff)).cast(pl.Int32) *
                                        pl.lit(diff).sign())
                selected = by_frame\
                    .join(group_counts, on=by_name)\
                    .select(pl.int_range(pl.len(), dtype=pl.Int32)
                            .shuffle(seed=seed)
                            .over(by_name)
                            .lt(pl.col.n))\
                    .to_series()
            else:
                selected = by_frame\
                    .select(pl.int_range(pl.len(), dtype=pl.Int32)
                            .shuffle(seed=seed)
                            .over(by_name)
                            .lt((fraction * pl.len().over(by_name)).round()))\
                    .to_series()
        elif QC_column is not None:
            selected = pl.int_range(QC_column.sum(), dtype=pl.Int32, 
                                    eager=True)\
                .shuffle(seed=seed)\
                .lt(n if fraction is None else (fraction * pl.len()).round())
        else:
            selected = self._obs\
                .select(pl.int_range(pl.len(), dtype=pl.Int32)
                        .shuffle(seed=seed)
                        .lt(n if fraction is None else
                            (fraction * pl.len()).round()))\
                .to_series()
        if QC_column is not None:
            # Back-project from QCed cells to all cells, filling with nulls
            selected = pl.when(QC_column)\
                .then(selected.gather(QC_column.cum_sum().cast(pl.Int32) - 1))
        sc = self.filter_obs(selected) if subsample_column is None else \
            self.with_columns_obs(selected.alias(subsample_column))
        if QC_column is not None:
            # noinspection PyTypeChecker
            sc._obs = sc._obs.drop(QC_column.name)
        return sc
    
    def subsample_var(self,
                      n: int | np.integer | None = None,
                      *,
                      fraction: int | float | np.integer | np.floating |
                                None = None,
                      by_column: SingleCellColumn | None = None,
                      subsample_column: str | None = None,
                      seed: int | np.integer = 0, 
                      overwrite: bool = False) -> SingleCell:
        """
        Subsample a specific number or fraction of genes.
        
        Args:
            n: the number of genes to return; mutually exclusive with
               `fraction`
            fraction: the fraction of genes to return; mutually exclusive with
                      `n`
            by_column: an optional String, Categorical, Enum, or integer column
                       of var to subsample by. Can be a column name, a polars
                       expression, a polars Series, a 1D NumPy array, or a
                       function that takes in this SingleCell dataset and
                       returns a polars Series or 1D NumPy array. Specifying
                       `by_column` ensures that the same fraction of genes with
                       each value of `by_column` are subsampled. When combined
                       with `n`, to make sure the total number of samples is
                       exactly `n`, some of the smallest groups may be
                       oversampled by one element, or some of the largest
                       groups may be undersampled by one element. Can contain
                       null entries: the corresponding genes will not be
                       included in the result.
            subsample_column: an optional name of a Boolean column to add to
                              var indicating the subsampled genes; if None,
                              subset to these genes instead
            seed: the random seed to use when subsampling
            overwrite: if True, overwrite `subsample_column` if already present
                       in var, instead of raising an error. Must be False when
                       `subsample_column` is None.
        
        Returns:
            A new SingleCell dataset subset to the subsampled genes, or if
            `subsample_column` is not None, the full dataset with
            `subsample_column` added to var.
        """
        if n is not None:
            check_type(n, 'n', int, 'a positive integer')
            check_bounds(n, 'n', 1)
        elif fraction is not None:
            check_type(fraction, 'fraction', float,
                       'a floating-point number between 0 and 1')
            check_bounds(fraction, 'fraction', 0, 1, left_open=True,
                         right_open=True)
        else:
            error_message = 'one of n and fraction must be specified'
            raise ValueError(error_message)
        if n is not None and fraction is not None:
            error_message = 'only one of n and fraction must be specified'
            raise ValueError(error_message)
        check_type(overwrite, 'overwrite', bool, 'Boolean')
        if subsample_column is not None:
            check_type(subsample_column, 'subsample_column', str, 'a string')
            if not overwrite and subsample_column in self._var:
                error_message = (
                    f'subsample_column {subsample_column!r} is already a '
                    f'column of var; did you already run subsample_var()? Set '
                    f'overwrite=True to overwrite.')
                raise ValueError(error_message)
        elif overwrite:
            error_message = \
                'overwrite must be False when subsample_column is None'
            raise ValueError(error_message)
        check_type(seed, 'seed', int, 'an integer')
        if by_column is not None:
            by_column = self._get_column(
                'var', by_column, 'by_column',
                (pl.String, pl.Categorical, pl.Enum, 'integer'),
                allow_null=True)
            by_frame = by_column.to_frame()
            by_name = by_column.name
            if n is not None:
                # Get a vector of the number of elements to sample per group.
                # The total sample size should exactly match the original n; if
                # necessary, oversample the smallest groups or undersample the
                # largest groups to make this happen.
                group_counts = by_frame\
                    .group_by(by_name)\
                    .agg(pl.len(), n=(n * pl.len() / len(by_column))
                                     .round().cast(pl.Int32))\
                    .drop_nulls(by_name)
                diff = n - group_counts['n'].sum()
                if diff != 0:
                    group_counts = group_counts\
                        .sort('len', descending=diff < 0)\
                        .with_columns(n=pl.col.n +
                                        pl.int_range(pl.len(), dtype=pl.Int32)
                                        .lt(abs(diff)).cast(pl.Int32) *
                                        pl.lit(diff).sign())
                selected = by_frame\
                    .join(group_counts, on=by_name)\
                    .select(pl.int_range(pl.len(), dtype=pl.Int32)
                            .shuffle(seed=seed)
                            .over(by_name)
                            .lt(pl.col.n))
            else:
                selected = by_frame\
                    .select(pl.int_range(pl.len(), dtype=pl.Int32)
                            .shuffle(seed=seed)
                            .over(by_name)
                            .lt((fraction * pl.len().over(by_name)).round()))
        else:
            selected = self._var\
                .select(pl.int_range(pl.len(), dtype=pl.Int32)
                        .shuffle(seed=seed)
                        .lt(n if fraction is None else
                            (fraction * pl.len()).round()))
        selected = selected.to_series()
        return self.filter_var(selected) if subsample_column is None else \
            self.with_columns_var(selected.alias(subsample_column))
    
    def pipe(self,
             function: Callable[[SingleCell, ...], Any],
             *args: Any,
             **kwargs: Any) -> Any:
        """
        Apply a function to a SingleCell dataset.
        
        Args:
            function: the function to apply
            *args: the positional arguments to the function
            **kwargs: the keyword arguments to the function

        Returns:
            function(self, *args, **kwargs)
        """
        return function(self, *args, **kwargs)
    
    def pipe_X(self,
               function: Callable[[csr_array | csc_array, ...],
                                  csr_array | csc_array],
               *args: Any,
               **kwargs: Any) -> SingleCell:
        """
        Apply a function to a SingleCell dataset's X.
        
        Args:
            function: the function to apply to X
            *args: the positional arguments to the function
            **kwargs: the keyword arguments to the function

        Returns:
            A new SingleCell dataset where the function has been applied to X.
        """
        return SingleCell(X=function(self._X, *args, **kwargs), obs=self._obs,
                          var=self._var, obsm=self._obsm, varm=self._varm,
                          uns=self._uns)
    
    def pipe_obs(self,
                 function: Callable[[pl.DataFrame, ...], pl.DataFrame],
                 *args: Any,
                 **kwargs: Any) -> SingleCell:
        """
        Apply a function to a SingleCell dataset's obs.
        
        Args:
            function: the function to apply to obs
            *args: the positional arguments to the function
            **kwargs: the keyword arguments to the function

        Returns:
            A new SingleCell dataset where the function has been applied to
            obs.
        """
        return SingleCell(X=self._X, obs=function(self._obs, *args, **kwargs),
                          var=self._var, obsm=self._obsm, varm=self._varm,
                          uns=self._uns)
    
    def pipe_var(self,
                 function: Callable[[pl.DataFrame, ...], pl.DataFrame],
                 *args: Any,
                 **kwargs: Any) -> SingleCell:
        """
        Apply a function to a SingleCell dataset's var.
        
        Args:
            function: the function to apply to var
            *args: the positional arguments to the function
            **kwargs: the keyword arguments to the function

        Returns:
            A new SingleCell dataset where the function has been applied to
            var.
        """
        return SingleCell(X=self._X, obs=self._obs,
                          var=function(self._var, *args, **kwargs),
                          obsm=self._obsm, varm=self._varm, uns=self._uns)
    
    def pipe_obsm(self,
                  function: Callable[[dict[str, np.ndarray[2, Any] |
                                                pl.DataFrame], ...],
                                     dict[str, np.ndarray[2, Any] |
                                               pl.DataFrame]],
                  *args: Any,
                  **kwargs: Any) -> SingleCell:
        """
        Apply a function to a SingleCell dataset's obsm.
        
        Args:
            function: the function to apply to obsm
            *args: the positional arguments to the function
            **kwargs: the keyword arguments to the function

        Returns:
            A new SingleCell dataset where the function has been applied to
            obsm.
        """
        return SingleCell(X=self._X, obs=self._obs, var=self._var,
                          obsm=function(self._obsm, *args, **kwargs),
                          varm=self._varm, uns=self._uns)
    
    def pipe_varm(self,
                  function: Callable[[dict[str, np.ndarray[2, Any] |
                                                pl.DataFrame], ...],
                                     dict[str, np.ndarray[2, Any] |
                                               pl.DataFrame]],
                  *args: Any,
                  **kwargs: Any) -> SingleCell:
        """
        Apply a function to a SingleCell dataset's varm.
        
        Args:
            function: the function to apply to varm
            *args: the positional arguments to the function
            **kwargs: the keyword arguments to the function

        Returns:
            A new SingleCell dataset where the function has been applied to
            varm.
        """
        return SingleCell(X=self._X, obs=self._obs, var=self._var,
                          obsm=self._obsm,
                          varm=function(self._varm, *args, **kwargs),
                          uns=self._uns)
    
    def pipe_uns(self,
                 function: Callable[[NestedScalarOrArrayDict, ...],
                                    NestedScalarOrArrayDict],
                 *args: Any,
                 **kwargs: Any) -> SingleCell:
        """
        Apply a function to a SingleCell dataset's uns.
        
        Args:
            function: the function to apply to uns
            *args: the positional arguments to the function
            **kwargs: the keyword arguments to the function

        Returns:
            A new SingleCell dataset where the function has been applied to
            uns.
        """
        return SingleCell(X=self._X, obs=self._obs, var=self._var,
                          obsm=self._obsm, varm=self._varm,
                          uns=function(self._uns, *args, **kwargs))
    
    def qc(self,
           custom_filter: SingleCellColumn | None = None,
           *,
           subset: bool = True,
           QC_column: str = 'passed_QC',
           max_mito_fraction: int | float | np.integer | np.floating |
                              None = 0.1,
           min_genes: int | np.integer | None = 100,
           MALAT1_filter: bool = True,
           num_threads: int | np.integer | None = 1,
           allow_float: bool = False,
           overwrite: bool = False,
           verbose: bool = True) -> SingleCell:
        """
        Adds a Boolean column to obs indicating which cells passed quality
        control (QC), or subsets to these cells if `subset=True`. By default,
        filters to non-doublet cells with a cell-type confidence of ≥90%, ≤10%
        mitochondrial reads, and ≥100 genes detected. Raises an error if any
        gene names appear more than once in var_names (they can be deduplicated
        with `make_var_names_unique()`).
        
        Args:
            custom_filter: an optional Boolean column of obs containing a
                           filter to apply on top of the other QC filters; True
                           elements will be kept. Can be a column name, a
                           polars expression, a polars Series, a 1D NumPy
                           array, or a function that takes in this SingleCell
                           dataset and returns a polars Series or 1D NumPy
                           array.
            subset: whether to subset to cells passing QC, instead of merely
                    adding a `QC_column` to obs. This will roughly double
                    memory usage, but speed up subsequent operations.
            QC_column: the name of a Boolean column to add to obs indicating
                       which cells passed QC, if `subset=False`. Gives an error
                       if obs already has a column with this name, unless
                       `overwrite=True`.
            max_mito_fraction: if not None, filter to cells with <= this
                               fraction of mitochondrial counts. The default
                               value of 10% is in between Seurat's recommended
                               value of 5%, and Scanpy's recommendation to not
                               filter on mitochondrial counts at all, at least
                               initially.
            min_genes: if not None, filter to cells with >= this many genes
                       detected (with non-zero count). The default of 100
                       matches Scanpy's recommended value, while Seurat
                       recommends a minimum of 200.
            MALAT1_filter: if True, filter out cells with 0 expression of the
                           nuclear-expressed lncRNA MALAT1, which likely
                           represent empty droplets or poor-quality cells
                           (biorxiv.org/content/10.1101/2024.07.14.603469v1).
                           There must be exactly one gene in obs_names with the
                           name MALAT1 or Malat1 to use this filter.
            num_threads: the number of threads to use when filtering based on
                         mitochondrial counts (the bottleneck of this
                         function). Set `num_threads=None` to use all available
                         cores (as determined by `os.cpu_count()`).
                         `num_threads` must be 1 (the default) when
                         `max_mito_fraction` is None and `MALAT1_filter` is
                         False, since `num_threads` is only used for the
                         mitochondrial count and MALAT1 filters.
            allow_float: if False, raise an error if `X.dtype` is
                         floating-point (suggesting the user may not be using
                         the raw counts); if True, disable this sanity check
            overwrite: if True, overwrite `QC_column` if already present in
                       obs, instead of raising an error. Must be False when
                       `subset=True`.
            verbose: whether to print how many cells were filtered out at each
                     step of the QC process
        
        Returns:
            A new SingleCell dataset with `QC_column` added to obs, or subset
            to QCed cells if `subset=True`, and `uns['QCed']` set to True.
        """
        X = self._X
        # Check that `self` is not already QCed
        if self._uns['QCed']:
            error_message = (
                "uns['QCed'] is True; did you already run qc()? Set "
                "uns['QCed'] = False or run with_uns(QCed=False) to bypass "
                "this check.")
            raise ValueError(error_message)
        # Check inputs
        if self.var_names.n_unique() < len(self.var):
            error_message = (
                'var_names contains duplicates; deduplicate with '
                'make_var_names_unique()')
            raise ValueError(error_message)
        if custom_filter is not None:
            custom_filter = self._get_column(
                'obs', custom_filter, 'custom_filter', pl.Boolean)
        check_type(subset, 'subset', bool, 'Boolean')
        check_type(overwrite, 'overwrite', bool, 'Boolean')
        if not subset:
            check_type(QC_column, 'QC_column', str, 'a string')
            if not overwrite and QC_column in self._obs:
                error_message = (
                    f'QC_column {QC_column!r} is already a column of obs; did '
                    f'you already run qc()? Set overwrite=True to overwrite.')
                raise ValueError(error_message)
        elif overwrite:
            error_message = 'overwrite must be False when subset is True'
            raise ValueError(error_message)
        if max_mito_fraction is not None:
            check_type(max_mito_fraction, 'max_mito_fraction', (int, float),
                       'a number between 0 and 1, inclusive')
            check_bounds(max_mito_fraction, 'max_mito_fraction', 0, 1)
        if min_genes is not None:
            check_type(min_genes, 'min_genes', int, 'a non-negative integer')
            check_bounds(min_genes, 'min_genes', 0)
        check_type(MALAT1_filter, 'MALAT1_filter', bool, 'Boolean')
        if num_threads is None:
            num_threads = os.cpu_count()
        else:
            check_type(num_threads, 'num_threads', int, 'a positive integer')
            check_bounds(num_threads, 'num_threads', 1)
        if num_threads != 1 and max_mito_fraction is None:
            error_message = (
                'num_threads must be 1 (the default) when max_mito_fraction '
                'is None, since num_threads is only used for the '
                'mitochondrial count filter')
            raise ValueError(error_message)
        check_type(allow_float, 'allow_float', bool, 'Boolean')
        check_type(verbose, 'verbose', bool, 'Boolean')
        # If `allow_float` is False, raise an error if `X` is floating-point
        if not allow_float and np.issubdtype(X.dtype, np.floating):
            error_message = (
                f'qc() requires raw counts but X.dtype is {str(X.dtype)!r}, '
                f'a floating-point data type. If you are sure that all values '
                f'are raw integer counts, i.e. that (X.data == '
                f'X.data.astype(int)).all(), then set allow_float=True.')
            raise TypeError(error_message)
        # Apply the custom filter, if specified
        if verbose:
            print(f'Starting with {len(self):,} cells.')
        mask = None
        if custom_filter is not None:
            if verbose:
                print('Applying the custom filter...')
            mask = custom_filter
            if verbose:
                print(f'{mask.sum():,} cells remain after applying the custom '
                      f'filter.')
        # Filter to cells with ≤ `100 * max_mito_fraction`% mitochondrial
        # counts, if `max_mito_fraction` was specified
        if max_mito_fraction is not None:
            if verbose:
                print(f'Filtering to cells with ≤{100 * max_mito_fraction}% '
                      f'mitochondrial counts...')
            var_names = self.var_names
            if var_names.dtype != pl.String:
                var_names = var_names.cast(pl.String)
            mt_genes = var_names.str.to_uppercase().str.starts_with('MT-')
            if not mt_genes.any():
                error_message = (
                    'no genes are mitochondrial (start with "MT-" or "mt-"); '
                    'this may happen if your var_names are Ensembl IDs (ENSG) '
                    'rather than gene symbols (in which case you should set '
                    'the gene symbols as the var_names with set_var_names()), '
                    'or if mitochondrial genes have already been filtered out '
                    '(in which case you can set max_mito_fraction to None)')
                raise ValueError(error_message)
            mito_mask = np.empty(X.shape[0], dtype=bool)
            prange_import = 'from cython.parallel cimport prange' \
                if num_threads > 1 else ''
            sum_variable_type = X.dtype
            if sum_variable_type == np.float32:
                sum_variable_type = float
            if isinstance(X, csr_array):
                cython_inline(rf'''
                    {prange_import}
                    def mito_mask(
                            const {cython_type(X.dtype)}[::1] data,
                            const {cython_type(X.indices.dtype)}[::1] indices,
                            const {cython_type(X.indptr.dtype)}[::1] indptr,
                            char[::1] mt_genes,
                            const double max_mito_fraction,
                            char[::1] mito_mask,
                            const unsigned num_threads):
                        cdef int row, col
                        cdef {cython_type(sum_variable_type)} row_sum, mt_sum
                        for row in \
                                {prange('indptr.shape[0] - 1', num_threads)}:
                            row_sum = mt_sum = 0
                            for col in range(indptr[row], indptr[row + 1]):
                                row_sum = row_sum + data[col]
                                if mt_genes[indices[col]]:
                                    mt_sum = mt_sum + data[col]
                            mito_mask[row] = (<double> mt_sum / row_sum) <= \
                                             max_mito_fraction
                        ''')['mito_mask'](
                            data=X.data, indices=X.indices, indptr=X.indptr,
                            mt_genes=mt_genes.to_numpy(),
                            max_mito_fraction=max_mito_fraction,
                            mito_mask=mito_mask, num_threads=num_threads)
            else:
                row_sums = np.zeros(X.shape[0], dtype=sum_variable_type)
                mt_sums = np.zeros(X.shape[0], dtype=sum_variable_type)
                cython_inline(rf'''
                    {prange_import}
                    def mito_mask(
                            const {cython_type(X.dtype)}[::1] data,
                            const {cython_type(X.indices.dtype)}[::1] indices,
                            const {cython_type(X.indptr.dtype)}[::1] indptr,
                            char[::1] mt_genes,
                            const double max_mito_fraction,
                            {cython_type(sum_variable_type)}[::1] row_sums,
                            {cython_type(sum_variable_type)}[::1] mt_sums,
                            char[::1] mito_mask,
                            const unsigned num_threads):
                        cdef int row, col, i
                        for col in \
                                {prange('indptr.shape[0] - 1', num_threads)}:
                            if mt_genes[col]:
                                for i in range(indptr[col], indptr[col + 1]):
                                    row_sums[indices[i]] += data[i]
                                    mt_sums[indices[i]] += data[i]
                            else:
                                for i in range(indptr[col], indptr[col + 1]):
                                    row_sums[indices[i]] += data[i]
                        for row in {prange('mito_mask.shape[0]', num_threads)}:
                            mito_mask[row] = \
                                (<double> mt_sums[row] / row_sums[row]) <= \
                                max_mito_fraction
                        ''')['mito_mask'](
                            data=X.data, indices=X.indices, indptr=X.indptr,
                            mt_genes=mt_genes.to_numpy(),
                            max_mito_fraction=max_mito_fraction,
                            row_sums=row_sums, mt_sums=mt_sums,
                            mito_mask=mito_mask, num_threads=num_threads)
            mito_mask = pl.Series(mito_mask)
            if not mito_mask.any():
                error_message = (
                    f'no cells remain after filtering to cells with '
                    f'≤{100 * max_mito_fraction}% mitochondrial counts')
                raise ValueError(error_message)
            if mask is None:
                mask = mito_mask
            else:
                mask &= mito_mask
            if verbose:
                print(f'{mask.sum():,} cells remain after filtering to cells '
                      f'with ≤{100 * max_mito_fraction}% mitochondrial '
                      f'counts.')
        # Filter to cells with ≥ `min_genes` genes detected, if specified
        if min_genes is not None:
            if verbose:
                print(f'Filtering to cells with ≥{min_genes:,} genes '
                      f'detected (with non-zero count)...')
            gene_mask = pl.Series(getnnz(X, axis=1, num_threads=num_threads) >=
                                  min_genes)
            if not gene_mask.any():
                error_message = (
                    f'no cells remain after filtering to cells with '
                    f'≥{min_genes:,} genes detected')
                raise ValueError(error_message)
            if mask is None:
                mask = gene_mask
            else:
                mask &= gene_mask
            if verbose:
                print(f'{mask.sum():,} cells remain after filtering to cells '
                      f'with ≥{min_genes:,} genes detected.')
        # Filter to cells with non-zero MALAT1 expression, if `MALAT1_filter`
        # is True
        if MALAT1_filter:
            if verbose:
                print(f'Filtering to cells with non-zero MALAT1 expression...')
            MALAT1_index = self._var\
                .select(pl.arg_where(pl.col(self.var_names.name).eq('MALAT1') |
                                     pl.col(self.var_names.name).eq('Malat1')))
            if len(MALAT1_index) == 0:
                error_message = (
                    f"neither 'MALAT1' nor 'Malat1' was found in var_names; "
                    f"this may happen if your var_names are Ensembl IDs "
                    f"(ENSG) rather than gene symbols (in which case you "
                    f"should set the gene symbols as the var_names with "
                    f"set_var_names()). Alternatively, set "
                    f"MALAT1_filter=False to disable filtering on MALAT1 "
                    f"expression.")
                raise ValueError(error_message)
            if len(MALAT1_index) == 2:
                error_message = (
                    "both 'MALAT1' and 'Malat1' were found in var_names; if "
                    "this is intentional, rename one of them before running "
                    "qc(), or set MALAT1_filter=False to disable filtering "
                    "on MALAT1 expression")
                raise ValueError(error_message)
            MALAT1_index = MALAT1_index.item()
            # The code below is a faster version of:
            # MALAT1_mask = (X[:, [MALAT1_index]] != 0).toarray().squeeze()
            MALAT1_mask = np.zeros(X.shape[0], dtype=bool)
            if isinstance(X, csr_array):
                prange_import = 'from cython.parallel cimport prange' \
                    if num_threads > 1 else ''
                # Don't use binary search because CSR indices may not be sorted
                cython_inline(rf'''
                    {prange_import}
                    def get_MALAT1_mask_csr(
                            const {cython_type(X.dtype)}[::1] data,
                            const {cython_type(X.indices.dtype)}[::1] indices,
                            const {cython_type(X.indptr.dtype)}[::1] indptr,
                            const int MALAT1_index,
                            char[::1] MALAT1_mask,
                            const unsigned num_threads):
                        cdef int row, col
                        for row in \
                                {prange('indptr.shape[0] - 1', num_threads)}:
                            for col in range(indptr[row], indptr[row + 1]):
                                if indices[col] == MALAT1_index:
                                    MALAT1_mask[row] = True
                                    break
                ''')['get_MALAT1_mask_csr'](
                    data=X.data, indices=X.indices, indptr=X.indptr,
                    MALAT1_index=MALAT1_index, MALAT1_mask=MALAT1_mask,
                    num_threads=num_threads)
            else:
                start = X.indptr[MALAT1_index]
                end = X.indptr[MALAT1_index + 1]
                MALAT1_mask[X.indices[start:end]] = True
            MALAT1_mask = pl.Series(MALAT1_mask)
            if mask is None:
                mask = MALAT1_mask
            else:
                mask &= MALAT1_mask
            if verbose:
                print(f'{mask.sum():,} cells remain after filtering to cells '
                      f'with non-zero MALAT1 expression.')
        # Add the mask of QCed cells as a column, or subset if `subset=True`
        if mask is None:
            error_message = 'no QC filters were specified'
            raise ValueError(error_message)
        if subset:
            if verbose:
                print(f'Subsetting to cells passing QC...')
            sc = self.filter_obs(mask)
        else:
            if verbose:
                print(f'Adding a Boolean column, obs[{QC_column!r}], '
                      f'indicating which cells passed QC...')
            sc = SingleCell(X=X, obs=self._obs.with_columns(
                pl.lit(mask).alias(QC_column)), var=self._var, obsm=self._obsm,
                              varm=self._varm, uns=self._uns)
        sc._uns['QCed'] = True
        return sc
    
    def make_obs_names_unique(self, separator: str = '-') -> SingleCell:
        """
        Make obs names unique by appending `'-1'` to the second occurence of
        a given obs name, `'-2'` to the third occurrence, and so on, where
        `'-'` can be switched to a different string via the `separator`
        argument. Raises an error if any obs_names already contain `separator`.
        
        Args:
            separator: the string connecting the original obs name and the
                       integer suffix

        Returns:
            A new SingleCell dataset with the obs names made unique.
        """
        check_type(separator, 'separator', str, 'a string')
        if self.obs_names.str.contains(separator).any():
            error_message = (
                f'some obs_names already contain the separator {separator!r}; '
                f'did you already run make_obs_names_unique()? If not, set '
                f'the separator argument to a different string.')
            raise ValueError(error_message)
        obs_names = pl.col(self.obs_names.name)
        num_times_seen = pl.int_range(pl.len(), dtype=pl.Int32).over(obs_names)
        return SingleCell(X=self._X,
                          obs=self._obs.with_columns(
                              pl.when(num_times_seen > 0)
                              .then(obs_names + separator +
                                    num_times_seen.cast(str))
                              .otherwise(obs_names)),
                          var=self._var, obsm=self._obsm, varm=self._varm,
                          uns=self._uns)
    
    def make_var_names_unique(self, separator: str = '-') -> SingleCell:
        """
        Make var names unique by appending `'-1'` to the second occurence of
        a given var name, `'-2'` to the third occurrence, and so on, where
        `'-'` can be switched to a different string via the `separator`
        argument. Raises an error if any var_names already contain `separator`.
        
        Args:
            separator: the string connecting the original var name and the
                       integer suffix

        Returns:
            A new SingleCell dataset with the var names made unique.
        """
        var_names = pl.col(self.var_names.name)
        num_times_seen = pl.int_range(pl.len(), dtype=pl.Int32).over(var_names)
        return SingleCell(X=self._X,
                          obs=self._obs,
                          var=self._var.with_columns(
                              pl.when(num_times_seen > 0)
                              .then(var_names + separator +
                                    num_times_seen.cast(str))
                              .otherwise(var_names)),
                          obsm=self._obsm, varm=self._varm, uns=self._uns)
    
    def get_sample_covariates(self,
                              ID_column: SingleCellColumn,
                              *,
                              QC_column: SingleCellColumn |
                                         None = 'passed_QC') -> pl.DataFrame:
        """
        Get a DataFrame of sample-level covariates, i.e. the columns of obs
        that are the same for all cells within each sample.
        
        Args:
            ID_column: a column of obs containing sample IDs. Can be a column
                       name, a polars expression, a polars Series, a 1D NumPy
                       array, or a function that takes in this SingleCell
                       dataset and returns a polars Series or 1D NumPy array.
            QC_column: an optional Boolean column of obs indicating which cells
                       passed QC. Can be a column name, a polars expression, a
                       polars Series, a 1D NumPy array, or a function that
                       takes in this SingleCell dataset and returns a polars
                       Series or 1D NumPy array. Set to None to include all
                       cells. Cells failing QC will be ignored.
        
        Returns:
            A DataFrame of the sample-level covariates, with ID_column (sorted)
            as the first column.
        """
        if QC_column is not None:
            QC_column = self._get_column(
                'obs', QC_column, 'QC_column', pl.Boolean,
                allow_missing=QC_column == 'passed_QC')
        ID_column = self._get_column('obs', ID_column, 'ID_column',
                                     (pl.String, pl.Categorical, pl.Enum,
                                      'integer'), QC_column=QC_column)
        ID_column_name = ID_column.name
        obs = self._obs
        if QC_column is not None:
            obs = obs.filter(QC_column)
            ID_column = ID_column.filter(QC_column)
        return obs\
            .select(ID_column,
                    *obs
                    .group_by(ID_column)
                    .n_unique()
                    .pipe(filter_columns,
                          (pl.exclude(ID_column_name)
                           if ID_column_name in obs else pl.all()).max().eq(1))
                    .columns)\
            .unique(ID_column_name)\
            .sort(ID_column_name)
    
    def pseudobulk(self,
                   ID_column: SingleCellColumn,
                   cell_type_column: SingleCellColumn,
                   *,
                   QC_column: SingleCellColumn | None = 'passed_QC',
                   additional_obs: pl.DataFrame | None = None,
                   sort_genes: bool = False,
                   num_threads: int | np.integer | None = 1,
                   verbose: bool = True) -> Pseudobulk:
        """
        Pseudobulks a single-cell dataset with sample ID and cell type columns,
        after filtering to cells passing QC according to `QC_column`. Returns a
        Pseudobulk dataset.
        
        You can run this function multiple times at different cell type
        resolutions by setting a different cell_type_column each time.
        
        Args:
            ID_column: a column of obs containing sample IDs. Can be a column
                       name, a polars expression, a polars Series, a 1D NumPy
                       array, or a function that takes in this SingleCell
                       dataset and returns a polars Series or 1D NumPy array.
            cell_type_column: a column of obs containing cell-type labels. Can
                              be a column name, a polars expression, a polars
                              Series, a 1D NumPy array, or a function that
                              takes in this SingleCell dataset and returns a
                              polars Series or 1D NumPy array.
            QC_column: an optional Boolean column of obs indicating which cells
                       passed QC. Can be a column name, a polars expression, a
                       polars Series, a 1D NumPy array, or a function that
                       takes in this SingleCell dataset and returns a polars
                       Series or 1D NumPy array. Set to None to include all
                       cells. Cells failing QC will be excluded from the
                       pseudobulk.
            additional_obs: an optional DataFrame of additional sample-level
                            covariates, which will be joined to the
                            pseudobulk's obs for each cell type
            sort_genes: whether to sort genes in alphabetical order in the
                        pseudobulk
            num_threads: the number of threads to use when pseudobulking;
                         parallelism happens across {sample, cell type} pairs
                         (or just samples, if `cell_type_column` is None).
                         Set `num_threads=None` to use all available cores
                         (as determined by `os.cpu_count()`). For count
                         matrices stored in the usual CSR format,
                         parallelization takes place across cell types and
                         samples, so specifying more cores than the number of
                         cell type-sample pairs may not improve performance.
            verbose: whether to print a warning when X is a csc_array;
                     pseudobulking may be much slower in this case
    
        Returns:
            A Pseudobulk dataset with X (counts), obs (metadata per sample),
            and var (metadata per gene) fields, each dicts across cell types.
            The columns of each cell type's obs will be:
            - 'ID' (a renamed version of `ID_column`)
            - 'num_cells' (the number of cells for that sample and cell type)
            followed by whichever columns of the SingleCell dataset's obs are
            constant across samples.
        """
        X = self._X
        # Check that `self` is QCed and not normalized
        if not self._uns['QCed']:
            error_message = (
                "uns['QCed'] is False; did you forget to run qc()? Set "
                "uns['QCed'] = True or run .with_uns(QCed=True) to bypass "
                "this check.")
            raise ValueError(error_message)
        if self._uns['normalized']:
            error_message = (
                "uns['normalized'] is True; did you already run normalize()?")
            raise ValueError(error_message)
        # Check inputs
        if QC_column is not None:
            QC_column = self._get_column(
                'obs', QC_column, 'QC_column', pl.Boolean,
                allow_missing=QC_column == 'passed_QC')
        ID_column = self._get_column('obs', ID_column, 'ID_column',
                                     (pl.String, pl.Categorical, pl.Enum,
                                      'integer'), QC_column=QC_column)
        ID_column_name = ID_column.name
        cell_type_column = \
            self._get_column('obs', cell_type_column, 'cell_type_column',
                             (pl.String, pl.Categorical, pl.Enum),
                             QC_column=QC_column)
        cell_type_column_name = cell_type_column.name
        num_cells_column_name = 'num_cells'
        for column_description, column_name in ('ID_column', ID_column_name), \
                ('cell_type_column', cell_type_column_name):
            if column_name == num_cells_column_name:
                error_message = (
                    f'{column_description} has the name '
                    f'{num_cells_column_name!r}, which conflicts with the '
                    f'name of the column to be added to the Pseudobulk '
                    f'dataset containing the number of cells of each cell '
                    f'type')
                raise ValueError(error_message)
        check_type(sort_genes, 'sort_genes', bool, 'Boolean')
        if num_threads is None:
            num_threads = os.cpu_count()
        else:
            check_type(num_threads, 'num_threads', int, 'a positive integer')
            check_bounds(num_threads, 'num_threads', 1)
        check_type(verbose, 'verbose', bool, 'Boolean')
        if additional_obs is not None:
            check_type(additional_obs, 'additional_obs', pl.DataFrame,
                       'a polars DataFrame')
            if ID_column_name not in additional_obs:
                error_message = (
                    f'ID_column {ID_column_name!r} is not a column of '
                    f'additional_obs')
                raise ValueError(error_message)
            if ID_column.dtype != additional_obs[ID_column_name].dtype:
                error_message = (
                    f'ID_column {ID_column_name!r} has a different data type '
                    f'in additional_obs than in this SingleCell dataset')
                raise TypeError(error_message)
        # Ensure the first column of var is Categorical, Enum, or String: this
        # is a requirement of the Pseudobulk class. (The first column of obs
        # must be as well, but this will always be true by construction, since
        # it will always be the sample ID.)
        if self.var_names.dtype not in (pl.Categorical, pl.Enum, pl.String):
            error_message = (
                f'the first column of var (var_names) has data type '
                f'{self.obs_names.dtype!r}, but must be Categorical, Enum, '
                f'or String')
            raise ValueError(error_message)
        # Get the row indices in sc.X that will be pseudobulked across for each
        # group (cell type-sample pair)
        # noinspection PyUnboundLocalVariable
        groups = (pl.DataFrame((cell_type_column, ID_column))
                  .with_columns(
                      _SingleCell_group_indices=pl.int_range(pl.len(),
                                                             dtype=pl.Int32))
                  if QC_column is None else
                  pl.DataFrame((cell_type_column, ID_column, QC_column))
                  .with_columns(
                      _SingleCell_group_indices=pl.int_range(pl.len(),
                                                             dtype=pl.Int32))
                  .filter(QC_column.name))\
            .group_by(cell_type_column_name, ID_column_name,
                      maintain_order=True)\
            .agg('_SingleCell_group_indices',
                 pl.len().alias(num_cells_column_name))\
            .sort(cell_type_column_name, ID_column_name)
        # Pseudobulk, storing the result in a preallocated NumPy array
        result = np.zeros((len(groups), X.shape[1]), dtype='int32', order='F')
        prange_import = \
            'from cython.parallel cimport prange' if num_threads > 1 else ''
        if isinstance(X, csr_array):
            group_indices = \
                groups['_SingleCell_group_indices'].explode().to_numpy()
            group_ends = groups[num_cells_column_name].cum_sum().to_numpy()
            cython_inline(rf'''
                {prange_import}
                def groupby_sum_csr(
                        const {cython_type(X.dtype)}[::1] data,
                        const {cython_type(X.indices.dtype)}[::1] indices,
                        const {cython_type(X.indptr.dtype)}[::1] indptr,
                        const int[::1] group_indices,
                        const unsigned[::1] group_ends,
                        int[::1, :] result,
                        const unsigned num_threads):
                    cdef int num_groups = group_ends.shape[0]
                    cdef int group, row, gene
                    cdef unsigned cell
                    # For each group (cell type-sample pair)...
                    for group in {prange('num_groups', num_threads)}:
                        # For each cell within this group...
                        for cell in range(
                                0 if group == 0 else group_ends[group - 1],
                                group_ends[group]):
                            # Get this cell's row index in the sparse matrix
                            row = group_indices[cell]
                            # For each gene (column) that's non-zero for this
                            # cell...
                            for gene in range(indptr[row], indptr[row + 1]):
                                # Add the value at this cell and gene to the
                                # total for this group and gene
                                result[group, indices[gene]] += int(data[gene])
                ''')['groupby_sum_csr'](
                    data=X.data, indices=X.indices, indptr=X.indptr,
                    group_indices=group_indices, group_ends=group_ends,
                    result=result, num_threads=num_threads)
        else:
            if verbose:
                print('Warning: X is a csc_array rather than a csr_array, '
                      'so pseudobulking may be thousands of times slower. '
                      'If you have enough memory, call .tocsr() on your '
                      'SingleCell dataset before pseudobulking. Set '
                      'verbose=False to silence this warning.')
            group_indices = pl.int_range(X.shape[0], dtype=pl.Int32,
                                         eager=True)\
                .to_frame('_SingleCell_group_indices')\
                .join(groups
                      .select('_SingleCell_group_indices',
                              _SingleCell_index=pl.int_range(pl.len(),
                                                             dtype=pl.Int32))
                      .explode('_SingleCell_group_indices'),
                      on='_SingleCell_group_indices', how='left')\
                ['_SingleCell_index']
            if QC_column is not None:
                group_indices = group_indices.fill_null(-1)
                continue_if_missing = 'if group == -1: continue'
            else:
                continue_if_missing = ''
            group_indices = group_indices.to_numpy()
            cython_inline(f'''
                {prange_import}
                def groupby_sum_csc(
                        const {cython_type(X.dtype)}[::1] data,
                        const {cython_type(X.indices.dtype)}[::1] indices,
                        const {cython_type(X.indptr.dtype)}[::1] indptr,
                        const int[::1] group_indices,
                        int[::1, :] result,
                        const unsigned num_threads):
                    cdef int gene, cell, group
                    # For each gene (column of the sparse matrix)...
                    for gene in {prange('result.shape[1]', num_threads)}:
                        # For each cell (row) that's non-zero for this gene...
                        for cell in range(indptr[gene], indptr[gene + 1]):
                            # Get the group index for this cell (-1 if it
                            # failed QC)
                            group = group_indices[indices[cell]]
                            {continue_if_missing}
                            # Add the value at this cell and gene to the total
                            # for this group and gene
                            result[group, gene] += <int> data[cell]
                    ''')['groupby_sum_csc'](
                        data=X.data, indices=X.indices, indptr=X.indptr,
                        group_indices=group_indices, result=result,
                        num_threads=num_threads)
        # Sort genes, if `sort_genes=True`
        cell_type_var = self._var
        if sort_genes:
            result = result[:, cell_type_var[:, 0].arg_sort()]
            cell_type_var = cell_type_var.sort(cell_type_var.columns[0])
        # Break up the results by cell type
        sample_covariates = self.get_sample_covariates(ID_column,
                                                       QC_column=QC_column)
        X, obs, var = {}, {}, {}
        start_index = 0
        for cell_type, count in groups[cell_type_column_name]\
                .value_counts(sort=True).iter_rows():
            end_index = start_index + count
            X[cell_type] = result[start_index:end_index]
            obs[cell_type] = groups.lazy()\
                .select(ID_column_name, num_cells_column_name)\
                .slice(start_index, count)\
                .join(sample_covariates.lazy(), on=ID_column_name, how='left')\
                .pipe(lambda df: df.join(additional_obs.lazy(),
                                         on=ID_column_name, how='left')
                      if additional_obs is not None else df)\
                .rename({ID_column_name: 'ID'})\
                .pipe(lambda df: df if QC_column is None else
                                 df.drop(QC_column.name))\
                .collect()
            var[cell_type] = cell_type_var
            start_index = end_index
        return Pseudobulk(X=X, obs=obs, var=var)
    
    def hvg(self,
            *others: SingleCell,
            QC_column: SingleCellColumn | None |
                       Sequence[SingleCellColumn | None] = 'passed_QC',
            batch_column: SingleCellColumn | None |
                          Sequence[SingleCellColumn | None] = None,
            num_genes: int | np.integer = 2000,
            min_cells: int | np.integer = 3,
            flavor: Literal['seurat_v3', 'seurat_v3_paper'] = 'seurat_v3',
            span: int | float | np.integer | np.floating = 0.3,
            hvg_column: str = 'highly_variable',
            rank_column: str = 'highly_variable_rank',
            num_threads: int | np.integer | None = 1,
            allow_float: bool = False,
            overwrite: bool = False) -> SingleCell | tuple[SingleCell, ...]:
        """
        Select highly variable genes using Seurat's algorithm. Operates on
        raw counts.
        
        By default, uses the same approach as Scanpy's
        `scanpy.pp.highly_variable_genes` function with the `flavor` argument
        set to the non-default value `'seurat_v3'`, and Seurat's
        `FindVariableFeatures` function with the `selection.method` argument
        set to the default value `'vst'`.

        Requires the scikit-misc package; install with:
        pip install --no-deps --no-build-isolation scikit-misc
        
        The general idea is that since genes with higher mean expression tend
        to have higher variance in expression (because they have more non-zero
        values), we want to select genes that have a high variance *relative to
        their mean expression*. Otherwise, we'd only be picking highly
        expressed genes! To correct for the mean-variance relationship, fit a
        LOESS curve fit to the mean-variance trend.
        
        Args:
            others: optional SingleCell datasets to jointly compute highly
                    variable genes across, alongside this one. Each dataset
                    will be treated as a separate batch. If `batch_column` is
                    not None, each dataset AND each distinct value of
                    `batch_column` within each dataset will be treated as a
                    separate batch. Variances will be computed per batch and
                    then aggregated (see `flavor`) across batches.
            QC_column: an optional Boolean column of obs indicating which cells
                       passed QC. Can be a column name, a polars expression, a
                       polars Series, a 1D NumPy array, or a function that
                       takes in this SingleCell dataset and returns a polars
                       Series or 1D NumPy array. Set to None to include all
                       cells. Cells failing QC will be ignored. When `others`
                       is specified, `QC_column` can be a
                       length-`1 + len(others)` sequence of columns,
                       expressions, Series, functions, or None for each
                       dataset (for `self`, followed by each dataset in
                       `others`).
            batch_column: an optional String, Categorical, Enum, or integer
                          column of obs indicating which batch each cell is
                          from. Can be a column name, a polars expression, a
                          polars Series, a 1D NumPy array, or a function that
                          takes in this SingleCell dataset and returns a polars
                          Series or 1D NumPy array. Each batch will be treated
                          as if it were a distinct dataset; this is exactly
                          equivalent to splitting the dataset with
                          `split_by(batch_column)` and then passing each of the
                          resulting datasets to `hvg()`, except that the
                          `min_cells` filter will always be calculated
                          per-dataset rather than per-batch. Variances will be
                          computed per batch and then aggregated (see `flavor`)
                          across batches. Set to None to treat each dataset as
                          having a single batch. When `others` is specified,
                          `batch_column` can be a length-`1 + len(others)`
                          sequence of columns, expressions, Series, functions,
                          or None for each dataset (for `self`, followed by
                          each dataset in `others`).
            num_genes: the number of highly variable genes to select. The
                       default of 2000 matches Seurat and Scanpy's recommended
                       value. Fewer than `num_genes` genes will be selected if
                       not enough genes have non-zero count in >= `min_cells`
                       cells (or when `min_cells` is None, if not enough genes
                       are present).
            min_cells: if not None, filter to genes detected (with non-zero
                       count) in >= this many cells in every dataset, before
                       calculating highly variable genes. The default value of
                       3 matches Seurat and Scanpy's recommended value. Note
                       that genes with zero variance in any dataset will always
                       be filtered out, even if `min_cells` is 0.
            flavor: the highly variable gene algorithm to use. Must be one of
                    `seurat_v3` and `seurat_v3_paper`, both of which match the
                    algorithms with the same name in scanpy. Both algorithms
                    select genes based on two criteria: 1) which genes are
                    ranked as most variable (taking the median of the ranks
                    across batches where the gene is among the top `num_genes`
                    highly variable genes) and 2) the number of batches in
                    which a gene is ranked in among the top `num_genes` in
                    variability. `seurat_v3` ranks genes by 1) and uses 2) to
                    tiebreak, whereas `seurat_v3_paper` ranks genes by 2) and
                    uses 1) to tiebreak. When there is only one batch, both
                    algorithms are the same and only rank based on 1).
            span: the span of the LOESS fit; higher values will lead to more
                  smoothing
            hvg_column: the name of a Boolean column to be added to (each
                        dataset's) var indicating the highly variable genes
            rank_column: the name of an integer column to be added to (each
                         dataset's) var with the rank of each highly variable
                         gene's variance (1 = highest variance, 2 =
                         next-highest, etc.); will be null for non-highly
                         variable genes. In the very unlikely event of ties,
                         the gene that appears first in var will get the lowest
                         rank.
            num_threads: the number of threads to use when finding highly
                         variable genes. Set `num_threads=None` to use all
                         available cores (as determined by `os.cpu_count()`).
            allow_float: if False, raise an error if `X.dtype` is
                         floating-point (suggesting the user may not be using
                         the raw counts); if True, disable this sanity check
            overwrite: if True, overwrite `hvg_column` and/or `rank_column` if
                       already present in var, instead of raising an error
        
        Returns:
            A new SingleCell dataset where var contains an additional Boolean
            column, `hvg_column` (default: 'highly_variable'), indicating the
            `num_genes` most highly variable genes, and `rank_column` (default:
            'highly_variable_rank') indicating the (one-based) rank of each
            highly variable gene's variance. Or, if additional SingleCell
            dataset(s) are specified via the `others` argument, a
            length-`1 + len(others)` tuple of SingleCell datasets with these
            two columns added: `self`, followed by each dataset in `others`.
        
        Note:
            This function may give an incorrect output if the count matrix
            contains explicit zeros (i.e. if `(X.data == 0).any()`): this is
            not checked for, due to speed considerations. In the unlikely event
            that your dataset contains explicit zeros, remove them by running
            `X.eliminate_zeros()` (an in-place operation).
        
        Note:
            This function may not give identical results to Seurat and Scanpy.
            It uses Welford's algorithm to compute the variance, which is more
            numerically stable than Scanpy and Seurat's calculations. If
            multiple genes are tied as the `num_genes`-th most highly variable
            gene in a batch/dataset, this function includes all of them,
            whereas Seurat and Scanpy arbitrarily pick one (or a subset) of
            them. This function also uses the ordering from a stable sort to
            break ties when selecting the final list of highly variable genes,
            instead of the unstable sort used by Seurat and Scanpy.
        """
        # noinspection PyUnresolvedReferences
        from skmisc.loess import loess
        # Check that all elements of `others` are SingleCell datasets
        if others:
            check_types(others, 'others', SingleCell, 'SingleCell datasets')
        datasets = [self] + list(others)
        # Check that all datasets are QCed and not normalized
        if not all(dataset._uns['QCed'] for dataset in datasets):
            suffix = ' for at least one dataset' if others else ''
            error_message = (
                f"uns['QCed'] is False{suffix}; did you forget to run qc()? "
                f"Set uns['QCed'] = True or run with_uns(QCed=True) to bypass "
                f"this check.")
            raise ValueError(error_message)
        if any(dataset._uns['normalized'] for dataset in datasets):
            suffix = ' for at least one dataset' if others else ''
            error_message = (
                f"hvg() requires raw counts but uns['normalized'] is "
                f"True{suffix}; did you already run normalize()?")
            raise ValueError(error_message)
        # Check that there are at least three cells in each dataset (since
        # LOESS seems to need at least three observations to converge)
        if any(len(dataset._obs) < 3 for dataset in datasets):
            suffix = ' for at least one dataset' if others else ''
            error_message = (
                f'there are fewer than three cells{suffix}, so '
                f'highly variable genes cannot be calculated')
            raise ValueError(error_message)
        # Get `QC_column` and `batch_column` from every dataset, if not None
        QC_columns = SingleCell._get_columns(
            'obs', datasets, QC_column, 'QC_column', pl.Boolean,
            allow_missing=True)
        batch_columns = SingleCell._get_columns(
            'obs', datasets, batch_column, 'batch_column',
            (pl.String, pl.Categorical, pl.Enum, 'integer'),
            QC_columns=QC_columns)
        # Check that `num_genes` is a positive integer
        check_type(num_genes, 'num_genes', int, 'a positive integer')
        check_bounds(num_genes, 'num_genes', 1)
        # Check that `min_cells` is a positive integer and at least as large as
        # the number of cells in each dataset
        check_type(min_cells, 'min_cells', int, 'a non-negative integer')
        check_bounds(min_cells, 'min_cells', 0)
        if any(len(dataset._obs) < min_cells for dataset in datasets):
            suffix = ' for at least one dataset' if others else ''
            error_message = (
                f'the number of cells in this dataset ({len(self._obs)}) '
                f'is less than min_cells ({min_cells}){suffix}; increase '
                f'min_cells')
            raise ValueError(error_message)
        # Check that `overwrite` is Boolean
        check_type(overwrite, 'overwrite', bool, 'Boolean')
        # Check that `hvg_column` and `rank_column` are strings and not already
        # in var for any dataset
        for column, column_name in (hvg_column, 'hvg_column'), \
                (rank_column, 'rank_column'):
            check_type(column, column_name, str, 'a string')
            if not overwrite and \
                    any(column in dataset._var for dataset in datasets):
                suffix = ' for at least one dataset' if others else ''
                error_message = (
                    f'{column_name} {column!r} is already a column of '
                    f'var{suffix}; did you already run hvg()? Set '
                    f'overwrite=True to overwrite.')
                raise ValueError(error_message)
        # Check that `num_threads` is a positive integer or None; if None, set
        # to `os.cpu_count()`
        if num_threads is None:
            num_threads = os.cpu_count()
        else:
            check_type(num_threads, 'num_threads', int, 'a positive integer')
            check_bounds(num_threads, 'num_threads', 1)
        # If `allow_float` is False, raise an error if `X` is floating-point
        # for any dataset
        check_type(allow_float, 'allow_float', bool, 'Boolean')
        if not allow_float:
            for dataset in datasets:
                X = dataset._X
                if np.issubdtype(X.dtype, np.floating):
                    error_message = (
                        f'hvg() requires raw counts but X.dtype is '
                        f'{str(X.dtype)!r}, a floating-point data type. If '
                        f'you are sure that all values are raw integer '
                        f'counts, i.e. that (X.data == X.data.astype(int))'
                        f'.all(), then set allow_float=True.')
                    raise TypeError(error_message)
        # Get the universe of genes we'll be considering: those present in any
        # dataset. If there are multiple datasets, also get the indices of
        # these genes in each dataset (with nulls for genes not present in that
        # particular dataset).
        if others:
            # The use of `align_frames` here is a bit wasteful memory-wise,
            # because it creates an identical `'gene'` column for every
            # DataFrame in `genes_and_indices`. Fortunately, it's only one
            # small string column per dataset.
            genes_and_indices = pl.align_frames(*(
                dataset.var[:, 0].to_frame('gene')
                .with_columns(index=pl.int_range(pl.len(), dtype=pl.Int32))
                for dataset in datasets), on='gene')
            # noinspection PyTypeChecker
            genes_in_any_dataset = genes_and_indices[0]['gene']
            # noinspection PyTypeChecker
            dataset_gene_indices = [df['index'] for df in genes_and_indices]
            del genes_and_indices
        else:
            genes_in_any_dataset = self.var_names.rename('gene')
        # Get the batches to calculate variance across (datasets + batches
        # within each dataset)
        if others:
            if batch_column is None:
                # noinspection PyUnboundLocalVariable
                batches = ((dataset._X, QC_column.to_numpy()
                            if QC_column is not None else None, gene_indices)
                           for dataset, QC_column, gene_indices in
                           zip(datasets, QC_columns, dataset_gene_indices))
            else:
                # noinspection PyUnboundLocalVariable
                batches = ((dataset._X,
                            (batch_column.eq(batch) if QC_column is None else
                             batch_column.eq(batch) & QC_column).to_numpy()
                            if batch is not None else
                            (QC_column.to_numpy() if QC_column is not None else
                             None), gene_indices)
                           for dataset, QC_column, batch_column, gene_indices
                           in zip(datasets, QC_columns, batch_columns,
                                  dataset_gene_indices)
                           for batch in ((None,) if batch_column is None else
                                         batch_column.unique()))
        else:
            X = self._X
            batch_column = batch_columns[0]
            if batch_column is None:
                if QC_column is not None and QC_columns[0] is not None:
                    batches = (X, QC_columns[0].to_numpy(), None),
                else:
                    batches = (X, None, None),
            else:
                if QC_column is not None and QC_columns[0] is not None:
                    batches = ((X, (batch_column.eq(batch) & QC_columns[0])
                                   .to_numpy(), None)
                               for batch in batch_column.unique())
                else:
                    batches = ((X, batch_column.eq(batch).to_numpy(), None)
                               for batch in batch_column.unique())
        # Get the variance of each gene in each batch across cells passing QC
        prange_import = \
            'from cython.parallel cimport prange' if num_threads > 1 else ''
        norm_gene_vars = []
        for X, cell_mask, gene_indices in batches:
            num_dataset_genes = X.shape[1]
            mean = np.zeros(num_dataset_genes)
            var = np.zeros(num_dataset_genes)
            nonzero_count = np.zeros(num_dataset_genes, dtype=np.int32)
            is_CSR = isinstance(X, csr_array)
            if is_CSR:
                if cell_mask is None:
                    cell_indices = np.array([], dtype=int)
                    num_cells = X.shape[0]
                else:
                    cell_indices = np.flatnonzero(cell_mask)
                    num_cells = len(cell_indices)
                pranges = {key: prange(key, num_threads, nogil=False)
                           for key in ('num_cells', 'num_dataset_genes',
                                       'indices.shape[0]')}
                cython_inline(rf'''
                    {prange_import}
                    def sparse_mean_var_minor_axis(
                            const {cython_type(X.dtype)}[::1] data,
                            const {cython_type(X.indices.dtype)}[::1] indices,
                            const {cython_type(X.indptr.dtype)}[::1] indptr,
                            const long[::1] cell_indices,
                            const int num_cells,
                            const int num_dataset_genes,
                            double[::1] mean,
                            double[::1] var,
                            int[::1] nonzero_count,
                            const unsigned num_threads):
                        
                        cdef int cell, gene, i, j, zero_count
                        cdef double delta, delta2, zero_fraction, \
                            zero_to_nonzero_ratio
                        cdef double one_over_num_cells = 1.0 / num_cells
                        cdef double one_over_num_cells_minus_one = \
                            1.0 / (num_cells - 1)
                        cdef double num_cells_over_num_cells_minus_1 = \
                            <double> num_cells / (num_cells - 1)
                        
                        {"with nogil" if num_threads > 1 else "if True"}:
                            # Calculate the mean and "unnormalized variance"
                            # (called `M2` in Welford's algorithm) per gene,
                            # across cells with non-zero counts for that gene
                            if cell_indices.shape[0] == 0:
                                for i in {pranges['indices.shape[0]']}:
                                    gene = indices[i]
                                    nonzero_count[gene] += 1
                                    delta = data[i] - mean[gene]
                                    mean[gene] += delta / nonzero_count[gene]
                                    delta2 = data[i] - mean[gene]
                                    var[gene] += delta * delta2
                            else:
                                for i in {pranges['num_cells']}:
                                    cell = cell_indices[i]
                                    for j in range(indptr[cell],
                                                   indptr[cell + 1]):
                                        gene = indices[j]
                                        nonzero_count[gene] += 1
                                        delta = data[j] - mean[gene]
                                        mean[gene] += \
                                            delta / nonzero_count[gene]
                                        delta2 = data[j] - mean[gene]
                                        var[gene] += delta * delta2
                            
                            # Then, calculate the mean and variance across all
                            # cells, by including the contribution from cells
                            # with zero counts
                            for gene in {pranges['num_dataset_genes']}:
                                zero_fraction = \
                                    nonzero_count[gene] * one_over_num_cells
                                mean[gene] *= zero_fraction
                                zero_count = num_cells - nonzero_count[gene]
                                if not zero_count:
                                    var[gene] = \
                                        one_over_num_cells_minus_one * \
                                        var[gene]
                                elif nonzero_count[gene]:
                                    zero_to_nonzero_ratio = \
                                        <double> zero_count / \
                                        nonzero_count[gene]
                                    var[gene] = \
                                        one_over_num_cells_minus_one * \
                                        var[gene] + \
                                        num_cells_over_num_cells_minus_1 * \
                                        zero_to_nonzero_ratio * \
                                        mean[gene] * mean[gene]
                    ''')['sparse_mean_var_minor_axis'](
                        data=X.data, indices=X.indices, indptr=X.indptr,
                        cell_indices=cell_indices, num_cells=num_cells,
                        num_dataset_genes=num_dataset_genes, mean=mean,
                        var=var, nonzero_count=nonzero_count,
                        num_threads=num_threads)
            else:
                if cell_mask is None:
                    cell_mask = np.array([], dtype=bool)
                    num_cells = X.shape[0]
                else:
                    num_cells = cell_mask.sum()
                num_dataset_genes_prange = \
                    prange('num_dataset_genes', num_threads)
                cython_inline(rf'''
                    {prange_import}
                    
                    cdef inline void adjust_for_zero_counts(
                            double[::1] mean,
                            double[::1] var,
                            const int gene,
                            const int num_cells,
                            const int nonzero_count,
                            const double one_over_num_cells,
                            const double one_over_num_cells_minus_one,
                            const double num_cells_over_num_cells_minus_1) \
                            noexcept nogil:
                        
                        cdef int zero_count
                        cdef double zero_fraction, zero_to_nonzero_ratio
                        
                        zero_fraction = nonzero_count * one_over_num_cells
                        mean[gene] *= zero_fraction
                        zero_count = num_cells - nonzero_count
                        if not zero_count:
                            var[gene] = \
                                one_over_num_cells_minus_one * \
                                var[gene]
                        elif nonzero_count:
                            zero_to_nonzero_ratio = \
                                <double> zero_count / nonzero_count
                            var[gene] = \
                                one_over_num_cells_minus_one * \
                                var[gene] + \
                                num_cells_over_num_cells_minus_1 * \
                                zero_to_nonzero_ratio * \
                                mean[gene] * mean[gene]
                    
                    def sparse_mean_var_major_axis(
                            const {cython_type(X.dtype)}[::1] data,
                            const {cython_type(X.indices.dtype)}[::1] indices,
                            const {cython_type(X.indptr.dtype)}[::1] indptr,
                            char[::1] cell_mask,
                            const int num_cells,
                            const int num_dataset_genes,
                            double[::1] mean,
                            double[::1] var,
                            int[::1] nonzero_count,
                            const unsigned num_threads):
                    
                        cdef int gene, gene_index, i
                        cdef double delta, delta2
                        cdef double one_over_num_cells = 1.0 / num_cells
                        cdef double one_over_num_cells_minus_one = \
                            1.0 / (num_cells - 1)
                        cdef double num_cells_over_num_cells_minus_1 = \
                            <double> num_cells / (num_cells - 1)
                        
                        # Calculate the mean and "unnormalized variance"
                        # (called `M2` in Welford's algorithm) per gene,
                        # across cells with non-zero counts for that gene
                        if cell_mask.shape[0] == 0:
                            for gene in {num_dataset_genes_prange}:
                                nonzero_count[gene] = 0
                                for i in range(indptr[gene],
                                               indptr[gene + 1]):
                                    nonzero_count[gene] += 1
                                    delta = data[i] - mean[gene]
                                    mean[gene] += delta / nonzero_count[gene]
                                    delta2 = data[i] - mean[gene]
                                    var[gene] += delta * delta2
                                # Then, calculate the mean and variance
                                # across all cells, by including the
                                # contribution from cells with zero counts
                                adjust_for_zero_counts(
                                    mean, var, gene, num_cells,
                                    nonzero_count[gene], one_over_num_cells,
                                    one_over_num_cells_minus_one,
                                    num_cells_over_num_cells_minus_1)
                        else:
                            # As above, but only include cells where
                            # `cell_mask[indices[i]]` is True
                            for gene in {num_dataset_genes_prange}:
                                nonzero_count[gene] = 0
                                for i in range(indptr[gene],
                                               indptr[gene + 1]):
                                    if cell_mask[indices[i]]:
                                        nonzero_count[gene] += 1
                                        delta = data[i] - mean[gene]
                                        mean[gene] += \
                                            delta / nonzero_count[gene]
                                        delta2 = data[i] - mean[gene]
                                        var[gene] += delta * delta2
                                adjust_for_zero_counts(
                                    mean, var, gene, num_cells,
                                    nonzero_count[gene], one_over_num_cells,
                                    one_over_num_cells_minus_one,
                                    num_cells_over_num_cells_minus_1)
                    ''')['sparse_mean_var_major_axis'](
                        data=X.data, indices=X.indices, indptr=X.indptr,
                        cell_mask=cell_mask, num_cells=num_cells,
                        num_dataset_genes=num_dataset_genes, mean=mean,
                        var=var, nonzero_count=nonzero_count,
                        num_threads=num_threads)
            
            not_constant = var > 0
            y = np.log10(var[not_constant])
            x = np.log10(mean[not_constant])
            model = loess(x, y, span=span)
            model.fit()
            
            estimated_variance = np.empty(num_dataset_genes)
            estimated_variance[not_constant] = model.outputs.fitted_values
            estimated_variance[~not_constant] = 0
            estimated_stddev = np.sqrt(10 ** estimated_variance)
            clip_val = mean + estimated_stddev * np.sqrt(num_cells)
            
            batch_counts_sum = np.zeros(num_dataset_genes)
            squared_batch_counts_sum = np.zeros(num_dataset_genes)
            if is_CSR:
                pranges = {key: prange(key, num_threads) for key in
                           ('indptr.shape[0] - 1', 'cell_indices.shape[0]')}
                # noinspection PyUnboundLocalVariable
                cython_inline(rf'''
                    {prange_import}
                    def clipped_sum(const {cython_type(X.dtype)}[::1] data,
                                    const {cython_type(X.indices.dtype)}[::1]
                                        indices,
                                    const {cython_type(X.indptr.dtype)}[::1]
                                        indptr,
                                    const long[::1] cell_indices,
                                    const double[::1] clip_val,
                                    double[::1] batch_counts_sum,
                                    double[::1] squared_batch_counts_sum,
                                    const unsigned num_threads):
                        cdef int cell, gene, i, j
                        cdef double value
                        
                        if cell_indices.shape[0] == 0:
                            for cell in {pranges['indptr.shape[0] - 1']}:
                                for i in range(indptr[cell], indptr[cell + 1]):
                                    gene = indices[i]
                                    value = data[i]
                                    if value > clip_val[gene]:
                                        value = clip_val[gene]
                                    batch_counts_sum[gene] += value
                                    squared_batch_counts_sum[gene] += \
                                        value ** 2
                        else:
                            for j in {pranges['cell_indices.shape[0]']}:
                                cell = cell_indices[j]
                                for i in range(indptr[cell], indptr[cell + 1]):
                                    gene = indices[i]
                                    value = data[i]
                                    if value > clip_val[gene]:
                                        value = clip_val[gene]
                                    batch_counts_sum[gene] += value
                                    squared_batch_counts_sum[gene] += \
                                        value ** 2
                    ''')['clipped_sum'](
                        data=X.data, indices=X.indices, indptr=X.indptr,
                        cell_indices=cell_indices, clip_val=clip_val,
                        batch_counts_sum=batch_counts_sum,
                        squared_batch_counts_sum=squared_batch_counts_sum,
                        num_threads=num_threads)
            else:
                indptr_prange = prange('indptr.shape[0] - 1', num_threads)
                cython_inline(rf'''
                    {prange_import}
                    def clipped_sum(const {cython_type(X.dtype)}[::1] data,
                                    const {cython_type(X.indices.dtype)}[::1]
                                        indices,
                                    const {cython_type(X.indptr.dtype)}[::1]
                                        indptr,
                                    char[::1] cell_mask,
                                    const double[::1] clip_val,
                                    double[::1] batch_counts_sum,
                                    double[::1] squared_batch_counts_sum,
                                    const unsigned num_threads):
                        cdef int cell, gene, i
                        cdef double value, clip_val_gene
                        
                        if cell_mask.shape[0] == 0:
                            for gene in {indptr_prange}:
                                clip_val_gene = clip_val[gene]
                                for i in range(indptr[gene], indptr[gene + 1]):
                                    value = data[i]
                                    if value > clip_val_gene:
                                        value = clip_val_gene
                                    batch_counts_sum[gene] += value
                                    squared_batch_counts_sum[gene] += \
                                        value ** 2
                        else:
                            for gene in {indptr_prange}:
                                clip_val_gene = clip_val[gene]
                                for i in range(indptr[gene], indptr[gene + 1]):
                                    cell = indices[i]
                                    if cell_mask[cell]:
                                        value = data[i]
                                        if value > clip_val_gene:
                                            value = clip_val_gene
                                        batch_counts_sum[gene] += value
                                        squared_batch_counts_sum[gene] += \
                                            value ** 2
                    ''')['clipped_sum'](
                        data=X.data, indices=X.indices, indptr=X.indptr,
                        cell_mask=cell_mask, clip_val=clip_val,
                        batch_counts_sum=batch_counts_sum,
                        squared_batch_counts_sum=squared_batch_counts_sum,
                        num_threads=num_threads)
            norm_gene_var = pl.Series(
                (1 / ((num_cells - 1) * np.square(estimated_stddev))) *
                ((num_cells * np.square(mean)) + squared_batch_counts_sum -
                 2 * batch_counts_sum * mean))
            # If `min_cells` is non-zero, set variances to null for genes
            # with a non-zero count less than `min_cells`
            if min_cells:
                norm_gene_var = norm_gene_var\
                    .set(pl.Series(nonzero_count < min_cells), None)
            # If there are multiple datasets, `norm_gene_var` is currently with
            # respect to the genes in `dataset.var_names`; map back to the
            # genes in `genes_in_any_dataset`, filling with nulls
            if others:
                norm_gene_var = norm_gene_var[gene_indices]
            norm_gene_vars.append(norm_gene_var)
        
        rank = pl.exclude('gene').rank('min', descending=True)
        final_rank = pl.struct(
            ('median_rank', 'nbatches') if flavor == 'seurat_v3' else
            ('nbatches', 'median_rank')).rank('ordinal')
        # note: the expression for `median_rank` can be replaced by
        # `pl.median_horizontal(pl.exclude('gene'))` once polars implements it
        hvgs = pl.DataFrame([genes_in_any_dataset] + norm_gene_vars)\
            .lazy()\
            .pipe(lambda df: df.drop_nulls(pl.selectors.exclude('gene'))
                             if min_cells or others else df)\
            .with_columns(pl.when(rank <= num_genes).then(rank))\
            .with_columns(nbatches=pl.sum_horizontal(pl.exclude('gene')
                                                     .is_null()),
                          median_rank=pl.concat_list(pl.exclude('gene'))
                                      .explode()
                                      .median().over(pl.int_range(
                                          pl.len(), dtype=pl.Int32)))\
            .select('gene', (final_rank <= num_genes).alias(hvg_column),
                    pl.when(final_rank <= num_genes).then(final_rank)
                    .alias(rank_column))\
            .collect()
        # Return a new SingleCell dataset (or a tuple of datasets, if others
        # is non-empty) containing the highly variable genes
        for dataset_index, dataset in enumerate(datasets):
            new_var = dataset._var\
                .join(hvgs.rename({'gene': dataset.var_names.name}),
                      on=dataset.var_names.name, how='left')\
                .with_columns(pl.col(hvg_column).fill_null(False))
            datasets[dataset_index] = \
                SingleCell(X=dataset._X, obs=dataset._obs, var=new_var,
                           obsm=dataset._obsm, varm=dataset._varm,
                           uns=dataset._uns)
        return tuple(datasets) if others else datasets[0]
    
    def normalize(self,
                  QC_column: SingleCellColumn | None = 'passed_QC',
                  method: Literal['PFlog1pPF', 'log1pPF',
                                  'logCP10k'] = 'PFlog1pPF',
                  num_threads: int | np.integer | None = 1,
                  allow_float: bool = False) -> SingleCell:
        """
        Normalize this SingleCell dataset's counts.
        
        By default, uses the PFlog1pPF method introduced in Booeshaghi et al.
        2022 (biorxiv.org/content/10.1101/2022.05.06.490859v1.full). With
        `method='logCP10k'`, it matches the default settings of Seurat's
        `NormalizeData` function, aside from differences in floating-point
        error.
        
        PFlog1pPF is a three-step process:
        1. Divide each cell's counts by a "size factor", namely the total
        number of counts for that cell, divided by the mean number of counts
        across all cells. Booeshaghi et al. call this process, which performs
        rowwise division of a matrix `X` by the vector
        `X.sum(axis=1) / X.sum(axis=1).mean()`, "proportional fitting" (PF).
        2. Take the logarithm of each entry plus 1, i.e. `log1p()`.
        3. Run an additional round of proportional fitting.
        
        If method='log1pPF', only performs steps 1 and 2 and leaves out step 3.
        Booeshaghi et al. call this method "log1pPF". Ahlmann-Eltze and Huber
        2023 (nature.com/articles/s41592-023-01814-1) recommend this method and
        argue that it outperforms log(CPM) normalization. However, Booeshaghi
        et al. note that log1pPF does not fully normalize for read depth,
        because the log transform of step 2 partially undoes the normalization
        introduced by step 1. This is the reasoning behind their use of step 3:
        to restore full depth normalization. By default, scanpy's
        normalize_total() uses a variation of proportional fitting that divides
        by the median instead of the mean, so it's closest to method='log1pPF'.
        
        If method='logCP10k', uses 10,000 for the denominator of the size
        factors instead of `X.sum(axis=1).mean()`, and leaves out step 3. This
        method is not recommended because it implicitly assumes an
        unrealistically large amount of overdispersion, and performs worse than
        log1pPF and PFlog1pPF in Ahlmann-Eltze and Huber and Booeshaghi et
        al.'s benchmarks. Seurat's NormalizeData() uses logCP10k normalization.
        
        Args:
            QC_column: an optional Boolean column of obs indicating which cells
                       passed QC. Can be a column name, a polars expression, a
                       polars Series, a 1D NumPy array, or a function that
                       takes in this SingleCell dataset and returns a polars
                       Series or 1D NumPy array. Set to None to include all
                       cells. Cells failing QC will still be normalized, but
                       will not count towards the calculation of the mean total
                       count across cells when `method` is `'PFlog1pPF'` or
                       `'log1pPF'`. Has no effect when `method` is
                       `'logCP10k'`.
            method: the normalization method to use (see above)
            num_threads: the number of threads to use when normalizing. Set
                         `num_threads=None` to use all available cores (as
                         determined by `os.cpu_count()`).
            allow_float: if False, raise an error if `X.dtype` is
                         floating-point (suggesting the user may not be using
                         the raw counts); if True, disable this sanity check
        
        Returns:
            A new SingleCell dataset with the normalized counts, and
            `uns['normalized']` set to True.
        """
        # Check that `self` is QCed and not already normalized
        if not self._uns['QCed']:
            error_message = (
                "uns['QCed'] is False; did you forget to run qc()? Set "
                "uns['QCed'] = True or run with_uns(QCed=True) to bypass this "
                "check.")
            raise ValueError(error_message)
        if self._uns['normalized']:
            error_message = \
                "uns['normalized'] is True; did you already run normalize()?"
            raise ValueError(error_message)
        # Get the QC column, if not None
        if QC_column is not None:
            QC_column = self._get_column(
                'obs', QC_column, 'QC_column', pl.Boolean,
                allow_missing=QC_column == 'passed_QC')
        # Check that `method` is one of the three valid methods
        if method not in ('PFlog1pPF', 'log1pPF', 'logCP10k'):
            error_message = \
                "method must be one of 'PFlog1pPF', 'log1pPF', or 'logCP10k'"
            raise ValueError(error_message)
        # Check that `num_threads` is a positive integer or None; if None, set
        # to `os.cpu_count()`
        if num_threads is None:
            num_threads = os.cpu_count()
        else:
            check_type(num_threads, 'num_threads', int, 'a positive integer')
            check_bounds(num_threads, 'num_threads', 1)
        # If `allow_float` is False, raise an error if `X` is floating-point
        X = self._X
        check_type(allow_float, 'allow_float', bool, 'Boolean')
        if not allow_float and np.issubdtype(X.dtype, np.floating):
            error_message = (
                f'normalize() requires raw counts but X.dtype is '
                f'{str(X.dtype)!r}, a floating-point data type. If you are '
                f'sure that all values are raw integer counts, i.e. that '
                f'(X.data == X.data.astype(int)).all(), then set '
                f'allow_float=True.')
            raise TypeError(error_message)
        # Step 1
        rowsums = X.sum(axis=1)
        inverse_size_factors = np.empty_like(rowsums, dtype=float) \
            if np.issubdtype(rowsums.dtype, np.integer) else rowsums
        # Note: QCed cells will have null as the batch, and over() treats
        # null as its own category, so effectively all cells failing QC
        # will be treated as their own batch. This doesn't matter since we
        # never use the counts for these cells anyway.
        np.divide(10_000 if method == 'logCP10k' else
                  rowsums.mean() if QC_column is None else
                  rowsums[QC_column].mean(), rowsums, inverse_size_factors)
        X = sparse_matrix_vector_op(X, '*', inverse_size_factors, axis=0,
                                    return_dtype=float,
                                    num_threads=num_threads)
        # Step 2
        if num_threads == 1:
            np.log1p(X.data, X.data)
        else:
            cython_inline(f'''
                from cython.parallel cimport prange
                from libc.math cimport log1p
                
                def log1p_parallel(double[::1] array,
                                   const unsigned num_threads):
                    cdef int i
                    for i in prange(array.shape[0], nogil=True,
                                    num_threads=num_threads):
                        array[i] = log1p(array[i])
                ''')['log1p_parallel'](X.data, num_threads)
        # Step 3
        if method == 'PFlog1pPF':
            rowsums = X.sum(axis=1)
            inverse_size_factors = rowsums
            np.divide(rowsums.mean() if QC_column is None else
                      rowsums[QC_column].mean(), rowsums, inverse_size_factors)
            sparse_matrix_vector_op(X, '*', inverse_size_factors, axis=0,
                                    inplace=True, return_dtype=float,
                                    num_threads=num_threads)
        sc = SingleCell(X=X, obs=self._obs, var=self._var, obsm=self._obsm,
                        varm=self._varm, uns=self._uns)
        sc._uns['normalized'] = True
        return sc
    
    def PCA(self,
            *others: SingleCell,
            QC_column: SingleCellColumn | None |
                       Sequence[SingleCellColumn | None] = 'passed_QC',
            hvg_column: SingleCellColumn |
                        Sequence[SingleCellColumn] = 'highly_variable',
            PC_key: str = 'PCs',
            num_PCs: int | np.integer = 50,
            seed: int | np.integer = 0,
            overwrite: bool = False,
            verbose: bool = True) -> SingleCell | tuple[SingleCell, ...]:
        """
        Compute principal components using irlba, the package used by Seurat.
        Operates on normalized counts (see `normalize()`).
        
        Install irlba with:
        
        from ryp import r
        r('install.packages("irlba", type="source")')
        
        IMPORTANT: if you already have a copy of irlba from CRAN (e.g.
        installed with Seurat), you will get the error:
        
        RuntimeError: in irlba(X, 50, verbose = FALSE) :
          function 'as_cholmod_sparse' not provided by package 'Matrix'
        
        This error will go away if you install irlba from source as described
        above.
        
        Args:
            others: optional SingleCell datasets to jointly compute principal
                    components across, alongside this one.
            QC_column: an optional Boolean column of obs indicating which cells
                       passed QC. Can be a column name, a polars expression, a
                       polars Series, a 1D NumPy array, or a function that
                       takes in this SingleCell dataset and returns a polars
                       Series or 1D NumPy array. Set to None to include all
                       cells. Cells failing QC will be ignored and have their
                       PCs set to NaN. When `others` is specified, `QC_column`
                       can be a length-`1 + len(others)` sequence of columns,
                       expressions, Series, functions, or None for each dataset
                       (for `self`, followed by each dataset in `others`).
            hvg_column: a Boolean column of var indicating the highly variable
                        genes. Can be a column name, a polars expression, a
                        polars Series, a 1D NumPy array, or a function that
                        takes in this SingleCell dataset and returns a polars
                        Series or 1D NumPy array. Set to None to use all genes.
                        When `others` is specified, `hvg_column`
                        can be a length-`1 + len(others)` sequence of columns,
                        expressions, Series, functions, or None for each
                        dataset (for `self`, followed by each dataset in
                        `others`).
            PC_key: the key of obsm where the principal components will be
                    stored
            num_PCs: the number of top principal components to calculate
            seed: the random seed to use for irlba when computing PCs, via R's
                  set.seed() function
            overwrite: if True, overwrite `PC_key` if already present in obsm,
                       instead of raising an error
            verbose: whether to set the verbose flag in irlba
        
        Returns:
            A new SingleCell dataset where obsm contains an additional key,
            `PC_key` (default: 'PCs'), containing a NumPy array of the top
            `num_PCs` principal components. Or, if additional SingleCell
            dataset(s) are specified via the `others` argument, a
            length-`1 + len(others)` tuple of SingleCell datasets with the PCs
            added: `self`, followed by each dataset in `others`.
        
        Note:
            Unlike Seurat's `RunPCA` function, which requires `ScaleData` to be
            run first, this function does not require the data to be scaled
            beforehand. Instead, it scales the data implicitly. It does this by
            providing the standard deviation and mean of the data to `irlba()`
            via its `scale` and `center` arguments, respectively. This approach
            is much more computationally efficient than explicit scaling, and
            is also taken by Seurat's internal (and currently unused)
            `RunPCA_Sparse` function, which this function is based on.
        """
        from ryp import r, to_py, to_r
        from sklearn.utils.sparsefuncs import mean_variance_axis
        r('suppressPackageStartupMessages(library(irlba))')
        # Check that all elements of `others` are SingleCell datasets
        if others:
            check_types(others, 'others', SingleCell, 'SingleCell datasets')
        datasets = [self] + list(others)
        # Check that all datasets are normalized
        suffix = ' for at least one dataset' if others else ''
        if not all(dataset._uns['normalized'] for dataset in datasets):
            error_message = (
                f"PCA() requires normalized counts but uns['normalized'] is "
                f"False{suffix}; did you forget to run normalize()?")
            raise ValueError(error_message)
        # Raise an error if `X` has an integer data type for any dataset
        for dataset in datasets:
            if np.issubdtype(dataset._X.dtype, np.integer):
                error_message = (
                    f'PCA() requires raw counts, but X.dtype is '
                    f'{str(dataset._X.dtype)!r}, an integer data type'
                    f'{suffix}; did you forget to run normalize() before '
                    f'PCA()?')
                raise TypeError(error_message)
        # Get `QC_column` (if not None) and `hvg_column` from every dataset
        QC_columns = SingleCell._get_columns(
            'obs', datasets, QC_column, 'QC_column', pl.Boolean,
            allow_missing=True)
        hvg_columns = SingleCell._get_columns(
            'var', datasets, hvg_column, 'hvg_column', pl.Boolean,
            allow_None=False,
            custom_error=f'hvg_column {{}} is not a column of var{suffix}; '
                         f'did you forget to run hvg() before PCA()?')
        # Check that `overwrite` is Boolean
        check_type(overwrite, 'overwrite', bool, 'Boolean')
        # Check that `PC_key` is not already in obsm
        check_type(PC_key, 'PC_key', str, 'a string')
        for dataset in datasets:
            if not overwrite and PC_key in dataset._obsm:
                error_message = (
                    f'PC_key {PC_key!r} is already a key of obsm{suffix}; did '
                    f'you already run PCA()? Set overwrite=True to overwrite.')
                raise ValueError(error_message)
        # Check that `num_PCs` is a positive integer
        check_type(num_PCs, 'num_PCs', int, 'a positive integer')
        check_bounds(num_PCs, 'num_PCs', 1)
        # Check that `seed` is an integer
        check_type(seed, 'seed', int, 'an integer')
        # Check that `verbose` is Boolean
        check_type(verbose, 'verbose', bool, 'Boolean')
        # Get the matrix to compute PCA across: a CSC array of counts for
        # highly variable genes (or all genes, if `hvg_column` is None) across
        # cells passing QC. Use X[np.ix_(rows, columns)] as a faster, more
        # memory-efficient alternative to X[rows][:, columns]. Use CSC rather
        # than CSR because irlba has a fast C-based implementation for CSC.
        if others:
            if hvg_column is None:
                genes_in_all_datasets = self.var_names\
                    .filter(self.var_names
                            .is_in(pl.concat([dataset.var_names
                                              for dataset in others])))
            else:
                hvg_in_self = self._var.filter(hvg_columns[0]).to_series()
                genes_in_all_datasets = hvg_in_self\
                    .filter(hvg_in_self
                            .is_in(pl.concat([dataset._var.filter(hvg_col)
                                              .to_series()
                                              for dataset, hvg_col in
                                              zip(others, hvg_columns[1:])])))
            gene_indices = (
                genes_in_all_datasets
                .to_frame()
                .join(dataset._var.with_columns(_SingleCell_index=pl.int_range(
                          pl.len(), dtype=pl.Int32)),
                      left_on=genes_in_all_datasets.name,
                      right_on=dataset.var_names.name, how='left')
                ['_SingleCell_index']
                .to_numpy()
                for dataset in datasets)
            if QC_column is None:
                Xs = [dataset._X[:, genes]
                      for dataset, genes in zip(datasets, gene_indices)]
            else:
                Xs = [dataset._X[np.ix_(QC_col.to_numpy(), genes)]
                      if QC_col is not None else dataset._X[:, genes]
                      for dataset, genes, QC_col in
                      zip(datasets, gene_indices, QC_columns)]
        else:
            if QC_column is None:
                if hvg_column is None:
                    Xs = [dataset._X for dataset in datasets]
                else:
                    Xs = [dataset._X[:, hvg_col.to_numpy()]
                          for dataset, hvg_col in zip(datasets, hvg_columns)]
            else:
                if hvg_column is None:
                    Xs = [dataset._X[QC_col.to_numpy()]
                          if QC_col is not None else dataset._X
                          for dataset, QC_col in zip(datasets, QC_columns)]
                else:
                    Xs = [dataset._X[np.ix_(QC_col.to_numpy(),
                                            hvg_col.to_numpy())]
                          if QC_col is not None else
                          dataset._X[:, hvg_col.to_numpy()]
                          for dataset, QC_col, hvg_col in
                          zip(datasets, QC_columns, hvg_columns)]
        X = vstack(Xs, format='csc')
        num_cells_per_dataset = np.array([X.shape[0] for X in Xs])
        del Xs
        # Check that `num_PCs` is at most the width of this matrix
        check_bounds(num_PCs, 'num_PCs', upper_bound=X.shape[1])
        # Run PCA with irlba (github.com/bwlewis/irlba/blob/master/R/irlba.R)
        # This section is adapted from
        # github.com/satijalab/seurat/blob/master/R/integration.R#L7276-L7317
        # Note: totalvar doesn't seem to be used by irlba, maybe a Seurat bug?
        # Note: mean_variance_axis() can be replaced by Welford's algorithm
        center, feature_var = mean_variance_axis(X, axis=0)
        scale = np.sqrt(feature_var)
        scale.clip(min=1e-8, out=scale)
        to_r(X, '.SingleCell.X')
        try:
            to_r(center, '.SingleCell.center')
            try:
                to_r(scale, '.SingleCell.scale')
                try:
                    r(f'set.seed({seed})')
                    r(f'.SingleCell.PCs = irlba(.SingleCell.X, {num_PCs}, '
                      f'verbose={str(verbose).upper()}, '
                      f'scale=.SingleCell.scale, '
                      f'center=.SingleCell.center)')
                    try:
                        PCs = to_py('.SingleCell.PCs$u', format='numpy') * \
                              to_py('.SingleCell.PCs$d', format='numpy')
                    finally:
                        r('rm(.SingleCell.PCs)')
                finally:
                    r('rm(.SingleCell.scale)')
            finally:
                r('rm(.SingleCell.center)')
        finally:
            r('rm(.SingleCell.X)')
        # Store each dataset's PCs in its obsm
        for dataset_index, (dataset, QC_col, num_cells, end_index) in \
                enumerate(zip(datasets, QC_columns, num_cells_per_dataset,
                              num_cells_per_dataset.cumsum())):
            start_index = end_index - num_cells
            dataset_PCs = PCs[start_index:end_index]
            # If `QC_col` is not None for this dataset, back-project from QCed
            # cells to all cells, filling with NaN
            if QC_col is not None:
                dataset_PCs_QCed = dataset_PCs
                dataset_PCs = np.full((len(dataset),
                                       dataset_PCs_QCed.shape[1]), np.nan)
                dataset_PCs[QC_col.to_numpy()] = dataset_PCs_QCed
            else:
                dataset_PCs = np.ascontiguousarray(dataset_PCs)
            datasets[dataset_index] = SingleCell(
                X=dataset._X, obs=dataset._obs, var=dataset._var,
                obsm=dataset._obsm | {PC_key: dataset_PCs}, varm=self._varm,
                uns=self._uns)
        return tuple(datasets) if others else datasets[0]
        
    def harmonize(self,
                  *others: SingleCell,
                  QC_column: SingleCellColumn | None |
                             Sequence[SingleCellColumn | None] = 'passed_QC',
                  batch_column: SingleCellColumn | None |
                                Sequence[SingleCellColumn | None] = None,
                  PC_key: str = 'PCs',
                  Harmony_key: str = 'Harmony_PCs',
                  num_clusters: int | np.integer | None = None,
                  max_iter_harmony: int | np.integer = 10,
                  max_iter_clustering: int | np.integer | None = 20,
                  block_proportion: int | float | np.integer |
                                    np.floating = 0.05,
                  tol_harmony: int | float | np.integer | np.floating = 1e-4,
                  tol_clustering: int | float | np.integer |
                                  np.floating = 1e-5,
                  ridge_lambda: int | float | np.integer | np.floating = 1,
                  sigma: int | float | np.integer | np.floating = 0.1,
                  theta: int | float | np.integer | np.floating = 2,
                  tau: int | float | np.integer | np.floating = 0,
                  seed: int | np.integer = 0,
                  overwrite: bool = False,
                  verbose: bool = True) -> \
            SingleCell | tuple[SingleCell, ...]:
        """
        Harmonize this SingleCell dataset with other datasets, using Harmony
        (nature.com/articles/s41592-019-0619-0). Harmony was originally written
        in R (github.com/immunogenomics/harmony) but has two Python ports,
        harmony-pytorch (github.com/lilab-bcb/harmony-pytorch), which our
        implementation is based on, and harmonypy
        (github.com/slowkow/harmonypy).
        
        Args:
            others: the other SingleCell datasets to harmonize this one with
            QC_column: an optional Boolean column of obs indicating which cells
                       passed QC. Can be a column name, a polars expression, a
                       polars Series, a 1D NumPy array, or a function that
                       takes in this SingleCell dataset and returns a polars
                       Series or 1D NumPy array. Set to None to include all
                       cells. Cells failing QC will be ignored and have their
                       Harmony embeddings set to NaN. When `others` is
                       specified, `QC_column` can be a length-`1 + len(others)`
                       sequence of columns, expressions, Series, functions, or
                       None for each dataset (for `self`, followed by each
                       dataset in `others`).
            batch_column: an optional String, Categorical, Enum, or integer
                          column of obs indicating which batch each cell is
                          from. Can be a column name, a polars expression, a
                          polars Series, a 1D NumPy array, or a function that
                          takes in this SingleCell dataset and returns a polars
                          Series or 1D NumPy array. Each batch will be treated
                          as if it were a distinct dataset; this is exactly
                          equivalent to splitting the dataset with
                          `split_by(batch_column)` and then passing each of the
                          resulting datasets to `harmonize()`. Set to None to
                          treat each dataset as having a single batch. When
                          `others` is specified, `batch_column` may be a
                          length-`1 + len(others)` sequence of columns,
                          expressions, Series, functions, or None for each
                          dataset (for `self`, followed by each dataset in
                          `others`).
            PC_key: the key of obsm containing the principal components
                    calculated with PCA(), to use as the input to Harmony
            Harmony_key: the key of obsm where the Harmony embeddings will be
                         stored; will be added in-place to both `self` and each
                         of the datasets in `others`!
            num_clusters: the number of clusters used in the Harmony algorithm.
                          If not specified, take the minimum of 100 and
                          floor(number of cells / 30).
            max_iter_harmony: the maximum number of iterations to run Harmony
                              for, if convergence is not achieved. Defaults to
                              10, like the original harmony package,
                              harmony-pytorch, and harmonypy. Set to None to
                              use as many iterations as necessary to achieve
                              convergence.
            max_iter_clustering: the maximum number of iterations to run the
                                 clustering step within each Harmony iteration
                                 for, if convergence is not achieved. Defaults
                                 to 20 iterations, like the original harmony
                                 package and harmonypy; this differs from
                                 the default of 200 iterations used by
                                 harmony-pytorch. Set to None to use as many
                                 iterations as necessary to achieve
                                 convergence.
            block_proportion: the proportion of cells to use in each batch
                              update in the clustering step; must be greater
                              than zero and less than or equal to 1
            tol_harmony: the relative tolerance used to determine whether to
                         stop Harmony before `max_iter_harmony` iterations;
                         must be positive
            tol_clustering: the relative tolerance used to determine whether to
                            stop clustering before `max_iter_clustering`
                            iterations; must be positive
            ridge_lambda: the ridge regression penalty used in the Harmony
                          correction step; must be non-negative
            sigma: the weight of the entropy term in the Harmony objective
                   function; must be non-negative
            theta: the weight of the diversity penalty term in the Harmony
                   objective function; must be non-negative
            tau: the discounting factor on theta; must be non-negative. By
                 default, `tau = 0`, so there is no discounting.
            seed: the random seed for Harmony
            overwrite: if True, overwrite `Harmony_key` if already present in
                       obsm, instead of raising an error
            verbose: whether to print details of the harmonization process
        
        Returns:
            A length-`1 + len(others)` tuple of SingleCell datasets with the
            Harmony embeddings added as obsm[Harmony_key]: `self`, followed by
            each dataset in `others`.
        """
        with ignore_sigint():
            import faiss
        # Check `others`
        if not others:
            error_message = 'others cannot be empty'
            raise ValueError(error_message)
        check_types(others, 'others', SingleCell, 'SingleCell datasets')
        datasets = [self] + list(others)
        # Get `QC_column` and `batch_column` from every dataset, if not None
        QC_columns = SingleCell._get_columns(
            'obs', datasets, QC_column, 'QC_column', pl.Boolean,
            allow_missing=True)
        QC_columns_NumPy = [QC_col.to_numpy() if QC_col is not None else None
                            for QC_col in QC_columns]
        batch_columns = SingleCell._get_columns(
            'obs', datasets, batch_column, 'batch_column',
            (pl.String, pl.Categorical, pl.Enum, 'integer'),
            QC_columns=QC_columns)
        # Check that `PC_key` is a key of obsm for every dataset
        check_type(PC_key, 'PC_key', str, 'a string')
        if not all(PC_key in dataset._obsm for dataset in datasets):
            error_message = (
                f'PC_key {PC_key!r} is not a column of obs for at least one '
                f'dataset; did you forget to run PCA() before harmonize()?')
            raise ValueError(error_message)
        # Check that `overwrite` is Boolean
        check_type(overwrite, 'overwrite', bool, 'Boolean')
        # Check that `Harmony_key` is a string and not already in obsm for any
        # dataset
        check_type(Harmony_key, 'Harmony_key', str, 'a string')
        if not overwrite and \
                any(Harmony_key in dataset._obsm for dataset in datasets):
            error_message = (
                f'Harmony_key {Harmony_key!r} is already a key of obsm for at '
                f'least one dataset; did you already run harmonize()? Set '
                f'overwrite=True to overwrite.')
            raise ValueError(error_message)
        # Check that `num_clusters`, `max_iter_harmony`, and
        # `max_iter_clustering` are None or a positive integer; if either max
        # iter argument is None, set it to `INT32_MAX`
        for parameter, parameter_name in (
                (num_clusters, 'num_clusters'),
                (max_iter_harmony, 'max_iter_harmony'),
                (max_iter_clustering, 'max_iter_clustering')):
            if parameter is not None:
                check_type(parameter, parameter_name, int,
                           'a positive integer')
                check_bounds(parameter, parameter_name, 1)
        if max_iter_harmony is None:
            max_iter_harmony = 2147483647
        if max_iter_clustering is None:
            max_iter_clustering = 2147483647
        # Check that `block_proportion` is a number and that
        # `0 < block_proportion <= 1`
        check_type(block_proportion, 'block_proportion', (int, float),
                   'a number greater than zero and less than or equal to 1')
        check_bounds(block_proportion, 'block_proportion', 0, 1,
                     left_open=True)
        # Check that `tol_harmony` and `tol_clustering` are positive numbers,
        # and that `ridge_lambda`, `sigma`, `theta`, and `tau` are non-negative
        # numbers. If any is an integer, cast it to a float.
        for parameter, parameter_name in (
                (tol_harmony, 'tol_harmony'),
                (tol_clustering, 'tol_clustering')):
            check_type(parameter, parameter_name, (int, float),
                       'a positive number')
            check_bounds(parameter, parameter_name, 0, left_open=True)
        for parameter, parameter_name in (
                (ridge_lambda, 'ridge_lambda'), (sigma, 'sigma'),
                (theta, 'theta'), (tau, 'tau')):
            check_type(parameter, parameter_name, (int, float),
                       'a non-negative number')
            check_bounds(parameter, parameter_name, 0)
        tol_harmony = float(tol_harmony)
        tol_clustering = float(tol_clustering)
        ridge_lambda = float(ridge_lambda)
        sigma = float(sigma)
        theta = float(theta)
        tau = float(tau)
        # Check that `seed` is an integer
        check_type(seed, 'seed', int, 'an integer')
        # Check that `verbose` is Boolean
        check_type(verbose, 'verbose', bool, 'Boolean')
        # Concatenate PCs (`Z`) across datasets; get labels indicating which
        # rows of these concatenated PCs come from each dataset or batch.
        # Check that the PCs are float64 and C-contiguous.
        Z = [dataset._obsm[PC_key] for dataset in datasets]
        for PCs in Z:
            dtype = PCs.dtype
            if dtype != float:
                error_message = (
                    f'obsm[{PC_key!r}].dtype is {dtype!r} for at least one '
                    f'dataset, but must be float64')
                raise TypeError(error_message)
            if not PCs.flags['C_CONTIGUOUS']:
                error_message = (
                    f'obsm[{PC_key!r}].dtype is not C-contiguous for at least '
                    f'one dataset; make it C-contiguous with '
                    f'np.ascontiguousarray(dataset.obsm[{PC_key!r}])')
                raise ValueError(error_message)
        if QC_column is not None:
            Z = [PCs[QCed] if QCed is not None else PCs
                 for PCs, QCed in zip(Z, QC_columns_NumPy)]
        num_cells_per_dataset = np.array(list(map(len, Z)))
        if batch_column is None:
            batch_labels = np.repeat(np.arange(len(num_cells_per_dataset),
                                               dtype=np.uint32),
                                     num_cells_per_dataset)
        else:
            batch_labels = []
            batch_index = 0
            for dataset, QC_col, batch_col in \
                    zip(datasets, QC_columns, batch_columns):
                if batch_col is not None:
                    if QC_col is not None:
                        batch_col = batch_col.filter(QC_col)
                    if batch_col.dtype in (pl.String, pl.Categorical, pl.Enum):
                        if batch_col.dtype != pl.Enum:
                            batch_col = batch_col\
                                .cast(pl.Enum(batch_col.unique().drop_nulls()))
                        batch_col = batch_col.to_physical()
                    batch_labels.append(batch_col.to_numpy() + batch_index)
                    batch_index += batch_col.n_unique()
                else:
                    batch_labels.append(np.full(batch_index,
                                                len(dataset) if QC_col is None
                                                else QC_col.sum()))
                    batch_index += 1
            batch_labels = np.concatenate(batch_labels)
        Z = np.concatenate(Z)
        
        # Run Harmony
        
        cython_functions = cython_inline(rf'''
            cdef extern from "<cmath>" namespace "std" nogil:
                double abs(double x)
                double exp(double x)
                double pow(double base, double exponent)
                double log(double x)
                double sqrt(double x)
            
            from cpython.exc cimport PyErr_CheckSignals
            from scipy.linalg.cython_blas cimport dgemm, dgemv
            
            cdef inline void matrix_multiply(const double[:, ::1] A,
                                             const double[:, ::1] B,
                                             double[:, ::1] C,
                                             const bint transpose_A,
                                             const bint transpose_B,
                                             const double alpha,
                                             const double beta) noexcept nogil:
                # Flip A <-> B and shape[0] <-> shape[1] since our matrices are
                # C-major
                cdef int m, n, k, lda, ldb
                cdef char transA, transB
                if transpose_B:
                    m = B.shape[0]
                    k = B.shape[1]
                    lda = k
                    transA = b'T'
                else:
                    m = B.shape[1]
                    k = B.shape[0]
                    lda = m
                    transA = b'N'
                if transpose_A:
                    n = A.shape[1]
                    ldb = n
                    transB = b'T'
                else:
                    n = A.shape[0]
                    ldb = k
                    transB = b'N'
                dgemm(&transA, &transB, &m, &n, &k, <double *> &alpha,
                      <double *> &B[0, 0], &lda, <double *> &A[0, 0], &ldb,
                      <double *> &beta, &C[0, 0], &m)
            
            cdef inline void matrix_vector_multiply(
                    const double[:, ::1] A,
                    const double[::1] X,
                    double[::1] Y,
                    const bint transpose,
                    const double alpha,
                    const double beta) noexcept nogil:
                # Flip trans since our matrix is C-major
                cdef int m = A.shape[1], n = A.shape[0], incx = 1, incy = 1
                cdef char trans = b'N' if transpose else b'T'
                dgemv(&trans, &m, &n, <double *> &alpha, <double *> &A[0,0],
                      &m, <double *> &X[0], &incx, <double *> &beta, &Y[0],
                      &incy)
            
            cdef inline int rand(long* state) noexcept nogil:
                cdef long x = state[0]
                state[0] = x * 6364136223846793005L + 1442695040888963407L
                cdef int s = (x ^ (x >> 18)) >> 27
                cdef int rot = x >> 59
                return (s >> rot) | (s << ((-rot) & 31))
            
            cdef inline long srand(const long seed) noexcept nogil:
                cdef long state = seed + 1442695040888963407L
                rand(&state)
                return state
            
            cdef inline int randint(const int bound, long* state) \
                    noexcept nogil:
                cdef int r, threshold = -bound % bound
                while True:
                    r = rand(state)
                    if r >= threshold:
                        return r % bound
            
            cdef inline void shuffle_array(int[::1] arr, long* state) \
                    noexcept nogil:
                cdef int i, j, temp
                for i in range(arr.shape[0] - 1, 0, -1):
                    j = randint(i + 1, state)
                    temp = arr[i]
                    arr[i] = arr[j]
                    arr[j] = temp
            
            cdef double compute_objective(double[:, ::1] Z_norm_times_Y_norm,
                                          const double[:, ::1] R,
                                          const double[:, ::1] E,
                                          const double[:, ::1] O,
                                          double[:, ::1] ratio,
                                          const double[::1] theta,
                                          double[::1] theta_times_ratio,
                                          const double sigma,
                                          const int num_cells,
                                          const int num_clusters,
                                          const int num_batches) \
                    noexcept nogil:
                cdef int i, j
                cdef double kmeans_error, entropy_term, diversity_penalty
                kmeans_error = entropy_term = diversity_penalty = 0
                for i in range(num_cells):
                    for j in range(num_clusters):
                        kmeans_error += \
                            R[i, j] * (1 - Z_norm_times_Y_norm[i, j])
                        entropy_term += R[i, j] * log(R[i, j])
                kmeans_error *= 2
                entropy_term *= sigma
                for i in range(num_batches):
                    for j in range(num_clusters):
                        ratio[i, j] = O[i, j] * log(
                            (O[i, j] + 1) / (E[i, j] + 1))
                matrix_vector_multiply(ratio, theta, theta_times_ratio,
                                       transpose=True, alpha=1, beta=0)
                for i in range(num_clusters):
                    diversity_penalty += theta_times_ratio[i]
                diversity_penalty *= sigma
                return kmeans_error + entropy_term + diversity_penalty
            
            def initialize(const double[:, ::1] Z_norm,
                           double[:, ::1] Y_norm,
                           double[:, ::1] Z_norm_times_Y_norm,
                           const int[::1] N_b,
                           double[::1] Pr_b,
                           const unsigned[::1] batch_labels,
                           double[:, ::1] R,
                           double[::1] R_sum,
                           double[:, ::1] E,
                           double[:, ::1] O,
                           const double sigma,
                           double[:, ::1] ratio,
                           double[::1] theta,
                           double[::1] theta_times_ratio,
                           const double tau):
                cdef int num_cells = Z_norm.shape[0]
                cdef int num_batches = E.shape[0]
                cdef int num_clusters = E.shape[1]
                cdef int i, j
                cdef unsigned batch_label
                cdef double norm, objective, base, two_over_sigma = 2 / sigma
                
                # Initialize Pr_b
                for i in range(num_batches):
                    Pr_b[i] = <double> N_b[i] / num_cells
                
                # Initialize R (and R_sum) and O
                matrix_multiply(Z_norm, Y_norm, Z_norm_times_Y_norm,
                                transpose_A=False, transpose_B=True, alpha=1,
                                beta=0)
                for i in range(num_cells):
                    batch_label = batch_labels[i]
                    norm = 0
                    for j in range(num_clusters):
                        R[i, j] = exp(two_over_sigma *
                                      (Z_norm_times_Y_norm[i, j] - 1))
                        norm += R[i, j]
                    norm = 1 / norm
                    for j in range(num_clusters):
                        R[i, j] *= norm
                        O[batch_label, j] += R[i, j]
                        R_sum[j] += R[i, j]
                
                # Initialize E
                for i in range(num_batches):
                    for j in range(num_clusters):
                        E[i, j] = Pr_b[i] * R_sum[j]
                
                # Apply discounting to theta, if specified
                if tau > 0:
                    for i in range(num_batches):
                        base = exp(-N_b[i] / (num_clusters * tau))
                        theta[i] = theta[i] * (1 - base * base)
                
                # Compute and return the initial value of the objective
                # function
                objective = compute_objective(
                    Z_norm_times_Y_norm, R, E, O, ratio, theta,
                    theta_times_ratio, sigma, num_cells, num_clusters,
                    num_batches)
                
                return objective
            
            def clustering(const double[:, ::1] Z_norm,
                           double[:, ::1] Z_norm_in,
                           double[:, ::1] Y_norm,
                           double[:, ::1] Z_norm_times_Y_norm,
                           const double[::1] Pr_b,
                           const unsigned[::1] batch_labels,
                           double[:, ::1] R,
                           double[:, ::1] R_in,
                           double[::1] R_in_sum,
                           double[:, ::1] E,
                           double[:, ::1] O,
                           double[:, ::1] ratio,
                           const double[::1] theta,
                           double[::1] theta_times_ratio,
                           int[::1] idx_list,
                           const double tol,
                           const int max_iter,
                           const double sigma,
                           const int block_size):
                cdef int num_cells = Z_norm.shape[0]
                cdef int num_PCs = Z_norm.shape[1]
                cdef int num_batches = E.shape[0]
                cdef int num_clusters = E.shape[1]
                cdef int i, j, k, iter, num_cells_in_block, pos
                cdef unsigned batch_label
                cdef long state
                cdef double two_over_sigma = 2 / sigma
                cdef double exp_neg_two_over_sigma = exp(-two_over_sigma)
                cdef double norm, objective = -1, old, new
                cdef double[:, ::1] Z_norm_block, R_block
                cdef double past_clustering_objectives[3]
                
                for iter in range(max_iter):
                    # Compute Cluster Centroids
                    matrix_multiply(R, Z_norm, Y_norm, transpose_A=True,
                                    transpose_B=False, alpha=1, beta=0)
                    for i in range(num_clusters):
                        norm = 0
                        for j in range(num_PCs):
                            norm = norm + Y_norm[i, j] * Y_norm[i, j]
                        norm = 1 / sqrt(norm)
                        for j in range(Y_norm.shape[1]):
                            Y_norm[i, j] = Y_norm[i, j] * norm
                    for i in range(num_cells):
                        idx_list[i] = i
                    state = srand(iter)
                    shuffle_array(idx_list, &state)
                    
                    # Update cells blockwise
                    pos = 0
                    Z_norm_block = Z_norm_in
                    R_block = R_in
                    while pos < num_cells:
                        if pos + block_size > num_cells:
                            num_cells_in_block = num_cells - pos
                            Z_norm_block = Z_norm_block[:num_cells_in_block]
                            R_block = R_block[:num_cells_in_block]
                        else:
                            num_cells_in_block = block_size
                        
                        # Remove the cells in this block from E and O
                        R_in_sum[:] = 0
                        for i in range(num_cells_in_block):
                            k = idx_list[pos + i]
                            batch_label = batch_labels[k]
                            for j in range(num_clusters):
                                O[batch_label, j] -= R[k, j]
                                R_in_sum[j] += R[k, j]
                            for j in range(num_PCs):
                                Z_norm_block[i, j] = Z_norm[k, j]
                        for i in range(num_batches):
                            for j in range(num_clusters):
                                E[i, j] -= Pr_b[i] * R_in_sum[j]
                        
                        # Recompute R for the removed cells
                        # Note: the original formula is
                        # exp(-2 / sigma * (1 - Z_norm_block @ Y_norm.T)),
                        # which expands to exp(-2 / sigma) *
                        # exp(2 / sigma * Z_norm_block @ Y_norm.T)). Since
                        # exp(-2 / sigma) is a constant, we fold it into
                        # `ratio`.
                        matrix_multiply(Z_norm_block, Y_norm, R_block,
                                        transpose_A=False, transpose_B=True,
                                        alpha=two_over_sigma, beta=0)
                        for i in range(num_batches):
                            for j in range(num_clusters):
                                ratio[i, j] = exp_neg_two_over_sigma * \
                                    pow((E[i, j] + 1) / (O[i, j] + 1),
                                        theta[i])
                        R_in_sum[:] = 0
                        for i in range(num_cells_in_block):
                            k = idx_list[i + pos]
                            batch_label = batch_labels[k]
                            norm = 0
                            for j in range(num_clusters):
                                R[k, j] = exp(R_block[i, j]) * \
                                          ratio[batch_label, j]
                                norm += R[k, j]
                            norm = 1 / norm
                            for j in range(num_clusters):
                                R[k, j] *= norm
                                R_in_sum[j] += R[k, j]
                                # Add the removed cells back into O
                                O[batch_label, j] += R[k, j]
                        
                        # Add the removed cells back into E
                        for i in range(num_batches):
                            for j in range(num_clusters):
                                E[i, j] += Pr_b[i] * R_in_sum[j]
                        
                        # Move to the next block
                        pos += block_size
                    
                    # Compute the objective and decide whether we've converged
                    matrix_multiply(Z_norm, Y_norm, Z_norm_times_Y_norm,
                                    transpose_A=False, transpose_B=True,
                                    alpha=1, beta=0)
                    objective = compute_objective(
                        Z_norm_times_Y_norm, R, E, O, ratio, theta,
                        theta_times_ratio, sigma, num_cells, num_clusters,
                        num_batches)
                    if iter < 3:
                        past_clustering_objectives[iter] = objective
                    else:
                        old = past_clustering_objectives[0] + \
                            past_clustering_objectives[1] + \
                            past_clustering_objectives[2]
                        new = past_clustering_objectives[1] + \
                            past_clustering_objectives[2] + objective
                        if old - new < tol * abs(old):
                            break
                        else:
                            past_clustering_objectives[0] = \
                                past_clustering_objectives[1]
                            past_clustering_objectives[1] = \
                                past_clustering_objectives[2]
                            past_clustering_objectives[2] = objective
                    PyErr_CheckSignals()
                return objective
            
            def correction(const double[:, ::1] Z,
                           double[:, ::1] Z_hat,
                           const double[:, ::1] R,
                           const double[:, ::1] O,
                           const double ridge_lambda,
                           const unsigned[::1] batch_labels,
                           double[::1] factor,
                           double[:, ::1] P,
                           double[:, ::1] P_t_B_inv,
                           double[:, ::1] inv_mat,
                           double[:, ::1] Phi_t_diag_R_by_X,
                           double[:, ::1] W):
                cdef int num_cells = Z.shape[0]
                cdef int num_PCs = Z.shape[1]
                cdef int num_batches = O.shape[0]
                cdef int num_clusters = O.shape[1]
                cdef int i, j, k
                cdef unsigned batch_label
                cdef double c, c_inv
                
                Z_hat[:] = Z[:]
                
                # Initialize P to the identity matrix
                P[:] = 0
                for i in range(num_batches + 1):
                    P[i, i] = 1
            
                for k in range(num_clusters):
                    # Compute factor, c_inv and P
                    c = 0
                    for i in range(num_batches):
                        factor[i] = 1 / (O[i, k] + ridge_lambda)
                        c += O[i, k] * (1 - factor[i] * O[i, k])
                        P[0, i + 1] = -factor[i] * O[i, k]
                    c_inv = 1 / c
                    
                    # Compute P_t_B_inv
                    P_t_B_inv[:] = 0
                    P_t_B_inv[0, 0] = c_inv
                    for i in range(1, num_batches + 1):
                        P_t_B_inv[i, i] = factor[i - 1]
                        P_t_B_inv[i, 0] = P[0, i] * c_inv
                    
                    # Compute inv_mat
                    matrix_multiply(P_t_B_inv, P, inv_mat, transpose_A=False,
                                    transpose_B=False, alpha=1, beta=0)
                    
                    # Compute Phi_t_diag_R @ X
                    Phi_t_diag_R_by_X[:] = 0
                    for i in range(num_cells):
                        batch_label = batch_labels[i]
                        for j in range(num_PCs):
                            Phi_t_diag_R_by_X[0, j] += Z[i, j] * R[i, k]
                            Phi_t_diag_R_by_X[batch_label + 1, j] += \
                                Z[i, j] * R[i, k]
                            
                    # Compute W
                    matrix_multiply(inv_mat, Phi_t_diag_R_by_X, W,
                                    transpose_A=False, transpose_B=False,
                                    alpha=1, beta=0)
                     
                    # Update Z_hat
                    for i in range(num_cells):
                        batch_label = batch_labels[i]
                        for j in range(num_PCs):
                            Z_hat[i, j] = Z_hat[i, j] - \
                                W[batch_label + 1, j] * R[i, k]
                        
            def normalize_rows(const double[:, ::1] arr, double[:, ::1] out):
                cdef int i, j
                cdef double norm
                for i in range(arr.shape[0]):
                    norm = 0
                    for j in range(arr.shape[1]):
                        norm = norm + arr[i, j] * arr[i, j]
                    norm = 1 / sqrt(norm)
                    for j in range(arr.shape[1]):
                        out[i, j] = arr[i, j] * norm

            def normalize_rows_inplace(double[:, ::1] arr):
                cdef int i, j
                cdef double norm
                for i in range(arr.shape[0]):
                    norm = 0
                    for j in range(arr.shape[1]):
                        norm = norm + arr[i, j] * arr[i, j]
                    norm = 1 / sqrt(norm)
                    for j in range(arr.shape[1]):
                        arr[i, j] = arr[i, j] * norm
        ''')
        initialize = cython_functions['initialize']
        clustering = cython_functions['clustering']
        correction = cython_functions['correction']
        normalize_rows = cython_functions['normalize_rows']
        normalize_rows_inplace = cython_functions['normalize_rows_inplace']
        
        # Get dimensions of everything
        num_cells, num_PCs = Z.shape
        block_size = int(num_cells * block_proportion)
        if num_clusters is None:
            num_clusters = min(100, int(num_cells / 30))
        N_b = bincount(batch_labels, num_bins=batch_labels[-1] + 1)
        num_batches = len(N_b)
        
        # Allocate arrays
        Z_norm = np.empty((num_cells, num_PCs))
        Z_norm_in = np.empty((block_size, num_PCs))
        Z_norm_times_Y_norm = np.empty((num_cells, num_clusters))
        Pr_b = np.empty(num_batches)
        R = np.empty((num_cells, num_clusters))
        R_in = np.empty((block_size, num_clusters))
        R_in_sum = np.zeros(num_clusters)
        E = np.empty((num_batches, num_clusters))
        O = np.zeros((num_batches, num_clusters))
        ratio = np.empty((num_batches, num_clusters))
        theta = np.repeat(theta, num_batches)
        theta_times_ratio = np.empty(num_clusters)
        idx_list = np.empty(num_cells, dtype=np.int32)
        factor = np.empty(num_cells)
        P = np.empty((num_batches + 1, num_batches + 1))
        P_t_B_inv = np.empty((num_batches + 1, num_batches + 1))
        inv_mat = np.empty((num_batches + 1, num_batches + 1))
        Phi_t_diag_R_by_X = np.empty((num_batches + 1, num_PCs))
        W = np.empty((num_batches + 1, num_PCs))
        
        # Run k-means
        normalize_rows(Z, Z_norm)
        kmeans = faiss.Kmeans(num_PCs, num_clusters, seed=seed)
        kmeans.train(Z_norm)
        Y_norm = kmeans.centroids.astype(float)
        normalize_rows_inplace(Y_norm)
        
        # Complete initialization in Cython
        objective = initialize(
            Z_norm=Z_norm, Y_norm=Y_norm,
            Z_norm_times_Y_norm=Z_norm_times_Y_norm, N_b=N_b, Pr_b=Pr_b,
            batch_labels=batch_labels, R=R, R_sum=R_in_sum, E=E, O=O,
            sigma=sigma, ratio=ratio, theta=theta,
            theta_times_ratio=theta_times_ratio, tau=tau)
        
        if verbose:
            print(f'Initialization is complete: objective = {objective:.2f}')
        
        iteration_string = plural('iteration', max_iter_harmony)
        
        for i in range(1, max_iter_harmony + 1):
            prev_objective = objective
            with Timer(f'[iteration {i}] clustering', verbose=verbose):
                objective = clustering(
                    Z_norm=Z_norm, Z_norm_in=Z_norm_in, Y_norm=Y_norm,
                    Z_norm_times_Y_norm=Z_norm_times_Y_norm, Pr_b=Pr_b,
                    batch_labels=batch_labels, R=R, R_in=R_in,
                    R_in_sum=R_in_sum, E=E, O=O, ratio=ratio, theta=theta,
                    theta_times_ratio=theta_times_ratio, idx_list=idx_list,
                    tol=tol_clustering, max_iter=max_iter_clustering,
                    sigma=sigma, block_size=block_size)
            with Timer(f'[iteration {i}] correction', verbose=verbose):
                correction(Z=Z, Z_hat=Z_norm, R=R, O=O,
                           ridge_lambda=ridge_lambda,
                           batch_labels=batch_labels, factor=factor, P=P,
                           P_t_B_inv=P_t_B_inv, inv_mat=inv_mat,
                           Phi_t_diag_R_by_X=Phi_t_diag_R_by_X, W=W)
            
            if verbose:
                print(f'Completed {i} of {max_iter_harmony} '
                      f'{iteration_string}: objective = {objective:.2f}')
            
            if prev_objective - objective < tol_harmony * abs(prev_objective):
                if verbose:
                    print(f'Reached convergence after {i} {iteration_string}')
                break
            
            if i == max_iter_harmony:
                if verbose:
                    print(f'Failed to converge after {i} {iteration_string}; '
                          f'returning the best effort so far')
                break
            
            normalize_rows_inplace(Z_norm)
        
        del batch_labels, Z, Pr_b, theta, R, E, O, Z_norm_in, Y_norm, \
            Z_norm_times_Y_norm, R_in, R_in_sum, ratio, theta_times_ratio, \
            idx_list, factor, P, P_t_B_inv, inv_mat, Phi_t_diag_R_by_X, W
        
        # Store each dataset's Harmony embedding in its obsm
        for dataset_index, (dataset, QC_col, num_cells, end_index) in \
                enumerate(zip(datasets, QC_columns_NumPy,
                              num_cells_per_dataset,
                              num_cells_per_dataset.cumsum())):
            start_index = end_index - num_cells
            dataset_Harmony_embedding = Z_norm[start_index:end_index]
            # If `QC_col` is not None for this dataset, back-project from
            # QCed cells to all cells, filling with NaN
            if QC_col is not None:
                dataset_Harmony_embedding_QCed = dataset_Harmony_embedding
                dataset_Harmony_embedding = np.full(
                    (len(dataset), dataset_Harmony_embedding_QCed.shape[1]),
                    np.nan)
                # noinspection PyUnboundLocalVariable
                dataset_Harmony_embedding[QC_col] = \
                    dataset_Harmony_embedding_QCed
            datasets[dataset_index] = SingleCell(
                X=dataset._X, obs=dataset._obs, var=dataset._var,
                obsm=dataset._obsm | {Harmony_key: dataset_Harmony_embedding},
                varm=self._varm, uns=self._uns)
        return tuple(datasets) if others else datasets[0]
    
    def label_transfer_from(
            self,
            other: SingleCell,
            original_cell_type_column: SingleCellColumn,
            *,
            QC_column: SingleCellColumn | None = 'passed_QC',
            other_QC_column: SingleCellColumn | None = 'passed_QC',
            Harmony_key: str = 'Harmony_PCs',
            cell_type_column: str = 'cell_type',
            confidence_column: str = 'cell_type_confidence',
            next_best_cell_type_column: str = 'next_best_cell_type',
            next_best_confidence_column: str =
                'next_best_cell_type_confidence',
            num_neighbors: int | np.integer = 20,
            num_clusters: int | np.integer | None = None,
            num_probes: int | np.integer | None = None,
            seed: int | np.integer = 0,
            num_threads: int | np.integer | None = 1,
            overwrite: bool = False,
            verbose: bool = True) -> SingleCell:
        """
        Transfer cell-type labels from another dataset to this one, by running
        approximate k-nearest neighbors on the Harmony embeddings from
        `harmonize()`.
                
        For each cell in `self`, the transferred cell-type label is the most
        common cell-type label among the `num_neighbors` cells in `other` with
        the nearest Harmony embeddings. The cell-type confidence is the
        fraction of these neighbors that share this most common cell-type
        label.
        
        Args:
            other: the dataset to transfer cell-type labels from
            original_cell_type_column: a column of `other.obs` containing
                                       cell-type labels. Can be a column name, 
                                       a polars expression, a polars Series, a
                                       1D NumPy array, or a function that takes
                                       in `other` and returns a polars Series
                                       or 1D NumPy array.
            QC_column: an optional Boolean column of `self.obs` indicating
                       which cells passed QC. Can be a column name, a polars 
                       expression, a polars Series, a 1D NumPy array, or a
                       function that takes in `self` and returns a polars 
                       Series or 1D NumPy array. Set to None to include all
                       cells. Cells failing QC will be ignored during the label
                       transfer.
            other_QC_column: an optional Boolean column of `other.obs`
                             indicating which cells passed QC. Can be a column
                             name, a polars expression, a polars Series, a 1D
                             NumPy array, or a function that takes in `other`
                             and returns a polars Series or 1D NumPy array. Set
                             to None to include all cells. Cells failing QC 
                             will be ignored during the label transfer.
            Harmony_key: the key of `self.obsm` and `other.obsm` containing the
                         Harmony embeddings for each dataset
            cell_type_column: the name of a column to be added to `self.obs`
                              indicating each cell's most likely cell type,
                              i.e. the most common cell-type label among the
                              cell's `num_neighbors` nearest neighbors in obs
            confidence_column: the name of a column to be added to `self.obs`
                               indicating each cell's cell-type confidence,
                               i.e. the fraction of the cell's `num_neighbors`
                               nearest neighbors in obs that share the most
                               common cell-type label. If multiple cell types
                               are equally common among the nearest neighbors,
                               tiebreak based on which of them is most common
                               in `other`.
            next_best_cell_type_column: the name of a column to be added to
                                        `self.obs` indicating each cell's
                                        second-most likely cell type, i.e. the
                                        second-most common cell-type label
                                        among the cell's `num_neighbors`
                                        nearest neighbors in obs
            next_best_confidence_column: the name of a column to be
                                                   added to
                                         `self.obs` indicating each cell's
                                         cell-type confidence, i.e. the
                                         fraction of the cell's `num_neighbors`
                                         nearest neighbors in obs that share
                                         the most common cell-type label. If
                                         multiple cell types are equally common
                                         among the nearest neighbors, tiebreak
                                         based on which of them is most common
                                         in `other`.
            num_neighbors: the number of nearest neighbors to use when
                           determining a cell's label. All cell-type
                           confidences will be multiples of
                           `1 / num_neighbors`.
            num_clusters: the number of k-means clusters to use during the
                          nearest-neighbor search. Called `nlist` internally by
                          faiss. Must be between 1 and the number of cells. If
                          None, will be set to
                          `ceil(min(sqrt(num_cells), num_cells / 100))`
                          clusters, i.e. the minimum of the square root of the 
                          number of cells and 1% of the number of cells, 
                          rounding up. The core of the heuristic, 
                          `sqrt(num_cells)`, is on the order of the range
                          recommended by faiss, 4 to 16 times the square root
                          (github.com/facebookresearch/faiss/wiki/
                          Guidelines-to-choose-an-index). However, faiss also
                          recommends using between 39 and 256 data points per
                          centroid when training the k-means clustering used
                          in the k-nearest neighbors search. If there are more
                          than 256, the dataset is automatically subsampled
                          for the k-means step, but if there are fewer than 39,
                          faiss gives a warning. To avoid this warning, we
                          switch to using `num_cells / 100` centroids for small
                          datasets, since 100 is the midpoint of 39 and 256 in
                          log space.
            num_probes: the number of nearest k-means clusters to search for a
                        given cell's nearest neighbors. Called `nprobe`
                        internally by faiss. Must be between 1 and
                        `num_clusters`, and should generally be a small
                        fraction of `num_clusters`. If None, will be set to
                        `min(num_clusters, 10)`.
            seed: the random seed to use when finding nearest neighbors
            num_threads: the number of threads to use for the nearest-neighbors
                         tree construction and search. Set `num_threads=None`
                         to use all available cores (as determined by
                         `os.cpu_count()`).
            overwrite: if True, overwrite `cell_type_column` and/or
                       `confidence_column` if already present in this dataset's
                       obs, instead of raising an error
            verbose: whether to print details of the nearest-neighbor search
        
        Returns:
            `self`, but with two additional columns: `cell_type_column`,
            containing the transferred cell-type labels, and
            `confidence_column`, containing the cell-type confidences.
        """
        with ignore_sigint():
            import faiss
        # Check that `other` is a SingleCell dataset
        check_type(other, 'other', SingleCell, 'a SingleCell dataset')
        # Get `QC_column` from `self` and `other_QC_column` from `other`
        if QC_column is not None:
            QC_column = self._get_column(
                'obs', QC_column, 'QC_column', pl.Boolean,
                allow_missing=QC_column == 'passed_QC')
        if other_QC_column is not None:
            other_QC_column = other._get_column(
                'obs', other_QC_column, 'other_QC_column', pl.Boolean,
                allow_missing=other_QC_column == 'passed_QC')
        # Get the number of cells in `self` and `other`
        num_cells_in_self = len(self._obs) if QC_column is None else \
            QC_column.sum()
        num_cells_in_other = len(other._obs) if other_QC_column is None else \
            other_QC_column.sum()
        # Get `original_cell_type_column` from `other`
        original_cell_type_column = other._get_column(
            'obs', original_cell_type_column, 'original_cell_type_column',
            (pl.Categorical, pl.Enum, pl.String), QC_column=other_QC_column)
        # If `other_QC_column` is not None, filter the cell type labels in
        # `original_cell_type_column` to cells passing QC
        if other_QC_column is not None:
            original_cell_type_column = \
                original_cell_type_column.filter(other_QC_column)
        # Check that `original_cell_type_column` has at least two distinct cell
        # types
        most_common_cell_types = \
            original_cell_type_column.value_counts(sort=True).to_series()
        if len(most_common_cell_types) == 1:
            error_message = (
                f'original_cell_type_column must have at least two distinct '
                f'cell types')
            if other_QC_column is not None:
                error_message += ' after filtering to cells passing QC'
            raise ValueError(error_message)
        # Check that `overwrite` is Boolean
        check_type(overwrite, 'overwrite', bool, 'Boolean')
        # Check that `Harmony_key` is a string and in both `self.obsm` and
        # `other.obsm`
        check_type(Harmony_key, 'Harmony_key', str, 'a string')
        datasets = (self, 'self'), (other, 'other')
        for dataset, dataset_name in datasets:
            if Harmony_key not in dataset._obsm:
                error_message = (
                    f'Harmony_key {Harmony_key!r} is not a column of '
                    f'{dataset_name}.obs; did you forget to run harmonize() '
                    f'before label_transfer_from()? Set overwrite=True to '
                    f'overwrite.')
                raise ValueError(error_message)
        # Check that `cell_type_column`, `confidence_column`,
        # `next_best_cell_type_column` and `next_best_confidence_column` are
        # strings and not already columns of `self.obs`
        for column, column_name in (
                (cell_type_column, 'cell_type_column'),
                (confidence_column, 'confidence_column'),
                (next_best_cell_type_column, 'next_best_cell_type_column'),
                (next_best_confidence_column, 'next_best_confidence_column')):
            check_type(column, column_name, str, 'a string')
            if not overwrite and column in self._obs:
                error_message = (
                    f'{column_name} {column!r} is already a column '
                    f'of obs; did you already run label_transfer_from()? Set '
                    f'overwrite=True to overwrite.')
                raise ValueError(error_message)
        # Check that `num_neighbors` is a positive integer
        check_type(num_neighbors, 'num_neighbors', int, 'a positive integer')
        check_bounds(num_neighbors, 'num_neighbors', 1)
        # Check that `num_clusters` is between 1 and the number of cells, and
        # that `num_probes` is between 1 and the number of clusters. If either
        # is None, set them to their default values.
        if num_clusters is None:
            num_clusters = int(np.ceil(min(np.sqrt(num_cells_in_other),
                                           num_cells_in_other / 100)))
        else:
            check_type(num_clusters, 'num_clusters', int, 'a positive integer')
            if not 1 <= num_clusters <= num_cells_in_other:
                error_message = (
                    f'num_clusters is {num_clusters:,}, but must be ≥ 1 and ≤ '
                    f'the number of cells in other ({num_cells_in_other:,})')
                raise ValueError(error_message)
        if num_probes is None:
            num_probes = min(num_clusters, 10)
        else:
            check_type(num_probes, 'num_probes', int, 'a positive integer')
            if not 1 <= num_probes <= num_clusters:
                error_message = (
                    f'num_probes is {num_probes:,}, but must be ≥ 1 and ≤ '
                    f'num_clusters ({num_clusters:,})')
                raise ValueError(error_message)
        # Check that `seed` is an integer
        check_type(seed, 'seed', int, 'an integer')
        # Check that `num_threads` is a positive integer or None; if None, set
        # to `os.cpu_count()`. Set this as the number of threads for faiss.
        if num_threads is None:
            num_threads = os.cpu_count()
        else:
            check_type(num_threads, 'num_threads', int, 'a positive integer')
            check_bounds(num_threads, 'num_threads', 1)
        faiss.omp_set_num_threads(num_threads)
        # Check that `verbose` is Boolean
        check_type(verbose, 'verbose', bool, 'Boolean')
        # Recode cell types so the most common is 0, the next-most common 1,
        # etc. This has the effect of breaking ties by taking the most common
        # cell type: we pick the first element in case of ties.
        cell_type_to_code = dict(zip(most_common_cell_types, range(
            len(most_common_cell_types))))
        original_cell_type_column = original_cell_type_column\
            .replace_strict(cell_type_to_code, return_dtype=pl.Int32)
        # Get the Harmony embeddings for self and other
        self_Harmony_embeddings = self._obsm[Harmony_key] \
            if QC_column is None else \
            self._obsm[Harmony_key][QC_column.to_numpy()]
        other_Harmony_embeddings = other._obsm[Harmony_key] \
            if other_QC_column is None else \
            other._obsm[Harmony_key][other_QC_column.to_numpy()]
        dim = self_Harmony_embeddings.shape[1]
        if other_Harmony_embeddings.shape[1] != dim:
            error_message = (
                f"the two datasets' Harmony embeddings have different numbers "
                f"of dimensions ({dim:,} vs "
                f"{other_Harmony_embeddings.shape[1]:,}")
            raise ValueError(error_message)
        # Use faiss to get the indices of the num_neighbors nearest
        # neighbors in other for each cell in self
        quantizer = faiss.IndexFlatL2(dim)
        quantizer.verbose = verbose
        index = faiss.IndexIVFFlat(quantizer, dim, num_clusters)
        index.cp.seed = seed
        index.verbose = verbose
        index.cp.verbose = verbose
        # noinspection PyArgumentList
        index.train(other_Harmony_embeddings)
        # noinspection PyArgumentList
        index.add(other_Harmony_embeddings)
        index.nprobe = num_probes
        # noinspection PyArgumentList
        nearest_neighbor_indices = \
            index.search(self_Harmony_embeddings, num_neighbors)[1]
        # Sometimes there aren't enough nearest neighbors for certain cells
        # with `num_probes` probes; if so, double `num_probes` (and threshold
        # to at most `num_clusters`), then re-run nearest-neighbor finding for
        # those cells
        needs_update = nearest_neighbor_indices[:, -1] == -1
        if needs_update.any():
            needs_update_X = self_Harmony_embeddings[needs_update]
            while True:
                num_probes = min(num_probes * 2, num_clusters)
                if verbose:
                    print(f'{len(needs_update_X):,} cells '
                          f'({len(needs_update_X) / len(self._obs):.2f}%) '
                          f'did not have enough neighbors with '
                          f'{index.nprobe:,} probes; re-running '
                          f'nearest-neighbors finding for these cells with '
                          f'{num_probes:,} probes')
                index.nprobe = num_probes
                # noinspection PyArgumentList
                new_indices = index.search(needs_update_X, num_neighbors)[1]
                nearest_neighbor_indices[needs_update] = new_indices
                still_needs_update = new_indices[:, -1] == -1
                if not still_needs_update.any():
                    break
                needs_update[needs_update] = still_needs_update
                needs_update_X = needs_update_X[still_needs_update]
        # Get the cell-type labels of these nearest neighbors (using our
        # integer encoding where the most common cell type is 0, the next-most
        # common 1, etc.)
        nearest_neighbor_cell_types = \
            original_cell_type_column.to_numpy()[nearest_neighbor_indices]
        # Get the two most common cell types for each cell in `self` among its
        # `num_neighbors` nearest neighbors in `other`. Pick the first element
        # in case of ties, which according to our encoding is the most common
        # cell type. Also get the cell-type confidence of these two cell types,
        # i.e. their frequency among the nearest neighbors.
        cell_types = np.empty(num_cells_in_self, dtype=np.uint32)
        confidences = np.empty(num_cells_in_self, dtype=float)
        next_best_cell_types = np.empty(num_cells_in_self, dtype=np.uint32)
        next_best_confidences = np.empty(num_cells_in_self, dtype=float)
        cython_inline(rf'''
            from cython.parallel cimport threadid, prange
            from libc.stdlib cimport free, malloc
            
            def label_transfer(const int[:, ::1] nearest_neighbor_cell_types,
                               const int num_cell_types,
                               unsigned[::1] cell_types,
                               double[::1] confidences,
                               unsigned[::1] next_best_cell_types,
                               double[::1] next_best_confidences,
                               const unsigned num_threads):
                cdef int i, j, start, thread_id, cell_type, count, \
                    most_common_cell_type, second_most_common_cell_type, \
                    max_count, second_max_count
                cdef int num_cells = nearest_neighbor_cell_types.shape[0]
                cdef int num_neighbors = nearest_neighbor_cell_types.shape[1]
                cdef double inv_num_neighbors = 1.0 / num_neighbors
                cdef int[::1] counts
                cdef int* counts_buffer = <int*> malloc(
                    num_cell_types * num_threads * sizeof(int))
                if not counts_buffer:
                    raise MemoryError
                try:
                    counts = \
                        <int[:num_cell_types * num_threads:]> counts_buffer
                    if num_threads == 1:
                        for i in range(num_cells):
                            counts[:] = 0
                            for j in range(num_neighbors):
                                cell_type = nearest_neighbor_cell_types[i, j]
                                counts[cell_type] += 1
                            if counts[0] >= counts[1]:
                                max_count = counts[0]
                                second_max_count = counts[1]
                                most_common_cell_type = 0
                                second_most_common_cell_type = 1
                            else:
                                max_count = counts[1]
                                second_max_count = counts[0]
                                most_common_cell_type = 1
                                second_most_common_cell_type = 0
                            for cell_type in range(2, num_cell_types):
                                count = counts[cell_type]
                                if count > max_count:
                                    second_max_count = max_count
                                    second_most_common_cell_type = \
                                        most_common_cell_type
                                    max_count = count
                                    most_common_cell_type = cell_type
                                elif count > second_max_count:
                                    second_max_count = count
                                    second_most_common_cell_type = cell_type
                            cell_types[i] = most_common_cell_type
                            confidences[i] = max_count * inv_num_neighbors
                            next_best_cell_types[i] = \
                                second_most_common_cell_type
                            next_best_confidences[i] = \
                                second_max_count * inv_num_neighbors
                    else:
                        for i in prange(num_cells, num_threads=num_threads,
                                        nogil=True):
                            thread_id = threadid()
                            start = num_cell_types * thread_id
                            counts[start:start + num_cell_types] = 0
                            for j in range(num_neighbors):
                                cell_type = nearest_neighbor_cell_types[i, j]
                                counts[start + cell_type] += 1
                            if counts[start] >= counts[start + 1]:
                                max_count = counts[start]
                                second_max_count = counts[start + 1]
                                most_common_cell_type = 0
                                second_most_common_cell_type = 1
                            else:
                                max_count = counts[start + 1]
                                second_max_count = counts[start]
                                most_common_cell_type = 1
                                second_most_common_cell_type = 0
                            for cell_type in range(2, num_cell_types):
                                count = counts[start + cell_type]
                                if count > max_count:
                                    second_max_count = max_count
                                    second_most_common_cell_type = \
                                        most_common_cell_type
                                    max_count = count
                                    most_common_cell_type = cell_type
                                elif count > second_max_count:
                                    second_max_count = count
                                    second_most_common_cell_type = cell_type
                            cell_types[i] = most_common_cell_type
                            confidences[i] = max_count * inv_num_neighbors
                            next_best_cell_types[i] = \
                                second_most_common_cell_type
                            next_best_confidences[i] = \
                                second_max_count * inv_num_neighbors
                finally:
                    free(counts_buffer)
            ''')['label_transfer'](
                nearest_neighbor_cell_types=nearest_neighbor_cell_types,
                num_cell_types=len(cell_type_to_code), cell_types=cell_types,
                confidences=confidences,
                next_best_cell_types=next_best_cell_types,
                next_best_confidences=next_best_confidences,
                num_threads=num_threads)
        # Map the cell-type codes back to their labels by constructing a polars
        # Series from the codes, then casting it to an Enum. Also convert
        # cell-type confidences to Series.
        cell_types = pl.Series(cell_type_column, cell_types)\
            .cast(pl.Enum(most_common_cell_types.to_list()))
        confidences = pl.Series(confidence_column, confidences)
        next_best_cell_types = \
            pl.Series(next_best_cell_type_column, next_best_cell_types)\
            .cast(pl.Enum(most_common_cell_types.to_list()))
        next_best_confidences = \
            pl.Series(next_best_confidence_column, next_best_confidences)
        # Add the four columns to obs. If `QC_column` is not None, back-project
        # from QCed cells to all cells, filling with nulls.
        columns = cell_types, confidences, next_best_cell_types, \
            next_best_confidences
        if QC_column is None:
            obs = self._obs.with_columns(columns)
        else:
            # noinspection PyTypeChecker
            expand = lambda series: pl.when(QC_column.name)\
                .then(pl.lit(series).gather(pl.col(QC_column.name).cum_sum()
                                            .cast(pl.Int32) - 1))
            obs = self._obs.with_columns(map(expand, columns))
        # Return a new SingleCell dataset containing the cell-type labels and
        # confidences
        return SingleCell(X=self._X, obs=obs, var=self._var, obsm=self._obsm,
                          varm=self._varm, uns=self._uns)
    
    def get_markers(self,
                    cell_type_column: SingleCellColumn,
                    *,
                    QC_column: SingleCellColumn | None = 'passed_QC',
                    min_detection_rate: int | float | np.integer |
                                        np.floating = 0.25,
                    min_fold_change: int | float | np.integer |
                                     np.floating = 2,
                    all_genes: bool = False,
                    num_threads: int | np.integer | None = 1,
                    verbose: bool = True) -> pl.DataFrame:
        """
        Find "marker genes" that distinguish each cell type from all other cell
        types. This function gives the same result regardless of whether it is
        run before or after normalization.
        
        Marker genes are chosen via an adaptation of the strategy of Fischer
        and Gillis 2021 (ncbi.nlm.nih.gov/pmc/articles/PMC8571500). For a given
        cell type, genes are scored based on a) their "detection rate" in that
        cell type (the fraction of cells of that type that have non-zero count
        for that gene), as well as b) the fold change in detection rate between
        that cell type and every other cell type. Genes must also have a
        detection rate of at least `min_detection_rate` (25% by default) and a
        minimum fold change of at least `min_fold_change` (2-fold by default)
        to be considered as markers.
        
        There is an inherent tradeoff between these two metrics. For instance,
        candidate marker genes with high enough expression to be expressed in
        every cell of a given type (i.e. to have a high detection rate) tend to
        also have at least some expression in other cell types (i.e. a low fold
        change in detection rate).
        
        Thus, marker genes are selected to optimally trade off between these
        two metrics: all genes on the Pareto front of the two metrics (i.e.
        genes for which there is no other gene that does better on both metrics
        simultaneously) are selected as marker genes.
        
        Note that Fischer and Gillis use AUROC versus log2 fold change in
        detection rate, instead of detection rate versus fold change in
        detection rate. However, detection rate is much faster to compute than
        AUROC, and is a very accurate proxy for AUROC: as Figure 1D in their
        paper shows, AUROC is almost perfectly correlated with detection rate
        across marker genes.
        
        Args:
            cell_type_column: a column of obs containing cell-type labels. Can
                              be a column name, a polars expression, a polars
                              Series, a 1D NumPy array, or a function that
                              takes in this SingleCell dataset and returns a
                              polars Series or 1D NumPy array.
            QC_column: an optional Boolean column of obs indicating which cells
                       passed QC. Can be a column name, a polars expression, a
                       polars Series, a 1D NumPy array, or a function that
                       takes in this SingleCell dataset and returns a polars
                       Series or 1D NumPy array. Set to None to include all
                       cells. Cells failing QC will be ignored.
            min_detection_rate: the minimum detection rate required to select
                                a gene as a marker gene; must be positive and
                                less than or equal to 1
            min_fold_change: the minimum fold change in detection rate required
                             to select a gene as a marker gene; must be greater
                             than 1
            all_genes: if True, include all genes in the output, not just
                       marker genes. An additional Boolean column will be
                       included to specify which genes are the marker genes.
                       Note that this option does not change which marker genes
                       are selected, only which information is returned.
            num_threads: the number of threads to use for marker-gene finding.
                         Set `num_threads=None` to use all available cores (as
                         determined by `os.cpu_count()`). For count matrices
                         stored in the usual CSR format, parallelization of the
                         most time-consuming step (calculating detection counts
                         in each cell type) takes place across cell types, so
                         specifying more cores than the number of cell types
                         may not improve peformance.
            verbose: whether to print a warning when X is a csc_array;
                     marker-gene finding may be much slower in this case

        Returns:
            By default, a DataFrame with one row per marker gene, with columns:
            - `'cell_type'`: a cell-type name from `cell_type_column`
            - `'gene'`: a gene symbol from var_names
            - `'detection_rate'`: the gene's detection rate in that cell type
            - `'fold_change'`, the gene's fold change in detection rate
              between that cell type and all other cell types
            If `all_genes=True`, a DataFrame with one row per cell type-gene
            pair, with those four columns plus one other:
            - `'marker'`, a Boolean column listing whether the gene is a marker
              for that cell type
            If `all_genes=False`, marker genes within each cell type will be
            sorted in increasing order of detection rate, and decreasing order
            of fold change.
        
        Note:
            This function may give an incorrect output if the count matrix
            contains explicit zeros (i.e. if `(X.data == 0).any()`): this is
            not checked for, due to speed considerations. In the unlikely event
            that your dataset contains explicit zeros, remove them by running
            `X.eliminate_zeros()` (an in-place operation).
        """
        X = self._X
        # Get the QC column, if not None
        if QC_column is not None:
            QC_column = self._get_column(
                'obs', QC_column, 'QC_column', pl.Boolean,
                allow_missing=QC_column == 'passed_QC')
        # Get the cell-type column
        cell_type_column = \
            self._get_column('obs', cell_type_column, 'cell_type_column',
                             (pl.String, pl.Categorical, pl.Enum),
                             QC_column=QC_column)
        cell_type_column_name = cell_type_column.name
        # Check that `num_threads` is a positive integer or None; if None, set
        # to `os.cpu_count()`
        if num_threads is None:
            num_threads = os.cpu_count()
        else:
            check_type(num_threads, 'num_threads', int, 'a positive integer')
            check_bounds(num_threads, 'num_threads', 1)
        # Check that `min_detection_rate` and `min_fold_change` are numeric and
        # have the correct ranges: 0 < min_detection_rate <= 1,
        # min_fold_change > 1
        check_type(min_detection_rate, 'min_detection_rate',
                   (int, float), 'a positive number less than or equal to 1')
        check_bounds(min_detection_rate, 'min_detection_rate', 0, 1,
                     left_open=True)
        check_type(min_fold_change, 'min_fold_change', (int, float),
                   'a number greater than 1')
        check_bounds(min_fold_change, 'min_fold_change', 1, left_open=True)
        # Check that `all_genes` and `verbose` are Boolean
        check_type(all_genes, 'all_genes', bool, 'Boolean')
        check_type(verbose, 'verbose', bool, 'Boolean')
        # Get the indices corresponding to the cells of each cell type
        # noinspection PyUnboundLocalVariable
        groups = (cell_type_column.to_frame()
                  .with_columns(
                      _SingleCell_group_indices=pl.int_range(pl.len(),
                                                             dtype=pl.Int32))
                  if QC_column is None else
                  pl.DataFrame((cell_type_column, QC_column))
                  .with_columns(
                      _SingleCell_group_indices=pl.int_range(pl.len(),
                                                             dtype=pl.Int32))
                  .filter(QC_column.name))\
            .group_by(cell_type_column_name, maintain_order=True)\
            .agg('_SingleCell_group_indices', _SingleCell_num_cells=pl.len())\
            .sort(cell_type_column_name)
        # Get a cell-type-by-gene matrix of the number of cells of each type
        # with non-zero expression of each gene, i.e. the gene's detection
        # count in that cell type
        num_cell_types = len(groups)
        if num_cell_types == 1:
            error_message = 'cell_type_column only contains one unique value'
            raise ValueError(error_message)
        num_genes = X.shape[1]
        detection_count = np.zeros((num_cell_types, num_genes), dtype=np.int32)
        prange_import = \
            'from cython.parallel cimport prange' if num_threads > 1 else ''
        if isinstance(X, csr_array):
            group_indices = \
                groups['_SingleCell_group_indices'].explode().to_numpy()
            group_ends = \
                groups['_SingleCell_num_cells'].cum_sum().to_numpy()
            cython_inline(f'''
                {prange_import}
                def groupby_getnnz_csr(
                        const {cython_type(X.indices.dtype)}[::1] indices,
                        const {cython_type(X.indptr.dtype)}[::1] indptr,
                        const int[::1] group_indices,
                        const unsigned[::1] group_ends,
                        int[:, ::1] result,
                        const unsigned num_threads):
                    cdef int num_groups = group_ends.shape[0]
                    cdef int group, row, gene
                    cdef unsigned cell
                    # For each group (cell type)...
                    for group in {prange('num_groups', num_threads)}:
                        # For each cell within this group...
                        for cell in range(
                                0 if group == 0 else group_ends[group - 1],
                                group_ends[group]):
                            # Get this cell's row index in the sparse matrix
                            row = group_indices[cell]
                            # For each gene (column) that's non-zero for this
                            # cell...
                            for gene in range(indptr[row], indptr[row + 1]):
                                # Add 1 to the total for this group and gene
                                result[group, indices[gene]] += 1
            ''')['groupby_getnnz_csr'](
                indices=X.indices, indptr=X.indptr,
                group_indices=group_indices, group_ends=group_ends,
                result=detection_count, num_threads=num_threads)
        else:
            if verbose:
                print('Warning: X is a csc_array rather than a csr_array, '
                      'so marker gene finding may be thousands of times '
                      'slower. If you have enough memory, call .tocsr() on '
                      'your SingleCell dataset before finding markers. Set '
                      'verbose=False to silence this warning.')
            group_indices = pl.int_range(X.shape[0], dtype=pl.Int32,
                                         eager=True)\
                .to_frame('_SingleCell_group_indices')\
                .join(groups
                      .select('_SingleCell_group_indices',
                              _SingleCell_index=pl.int_range(pl.len(),
                                                             dtype=pl.Int32))
                      .explode('_SingleCell_group_indices'),
                      on='_SingleCell_group_indices', how='left')\
                ['_SingleCell_index']
            if QC_column is not None:
                group_indices = group_indices.fill_null(-1)
                continue_if_missing = 'if group == -1: continue'
            else:
                continue_if_missing = ''
            group_indices = group_indices.to_numpy()
            cython_inline(f'''
                {prange_import}
                def groupby_getnnz_csc(
                        const {cython_type(X.indices.dtype)}[::1] indices,
                        const {cython_type(X.indptr.dtype)}[::1] indptr,
                        const int[::1] group_indices,
                        int[:, ::1] result,
                        const unsigned num_threads):
                    cdef int row, column, group
                    # For each gene (column of the sparse matrix)...
                    for column in {prange('result.shape[1]', num_threads)}:
                        # For each cell (row) that's non-zero for this gene...
                        for row in range(indptr[column], indptr[column + 1]):
                            # Get the group index for this row (-1 if it failed
                            # QC)
                            group = group_indices[indices[row]]
                            {continue_if_missing}
                            # Add 1 to the total for this group and column
                            result[group, column] += 1
                    ''')['groupby_getnnz_csc'](
                        indices=X.indices, indptr=X.indptr,
                        group_indices=group_indices, result=detection_count,
                        num_threads=num_threads)
        # For each cell type, calculate the detection rate and the fold change
        # of the detection rate. Also, initialize the candidate set of points
        # on the Pareto front to those with detection rate of at least
        # `min_detection_rate` and fold change of at least `min_fold_change`
        total_detection_count = detection_count.sum(axis=0, dtype=np.int32)
        num_cells_per_cell_type = groups['_SingleCell_num_cells'].to_numpy()
        total_num_cells = num_cells_per_cell_type.sum()
        detection_rate = np.empty((num_cell_types, num_genes), dtype=float)
        fold_change = np.empty((num_cell_types, num_genes), dtype=float)
        is_pareto = np.empty((num_cell_types, num_genes), dtype=bool)
        cython_inline(rf'''
            {prange_import}
            def get_detection_rate_and_fold_change(
                    const int[:, ::1] detection_count,
                    const int[::1] total_detection_count,
                    const unsigned[::1] num_cells_per_cell_type,
                    const unsigned total_num_cells,
                    const double min_detection_rate,
                    const double min_fold_change,
                    double[:, ::1] detection_rate,
                    double[:, ::1] fold_change,
                    char[:, ::1] is_pareto,
                    const unsigned num_threads):
                cdef int cell_type, gene
                cdef int count, background_count
                cdef unsigned num_cells, background_num_cells
                cdef double pair_detection_rate, pair_fold_change
                for cell_type in \
                        {prange('detection_count.shape[0]', num_threads)}:
                    num_cells = num_cells_per_cell_type[cell_type]
                    background_num_cells = total_num_cells - num_cells
                    for gene in range(detection_count.shape[1]):
                        count = detection_count[cell_type, gene]
                        pair_detection_rate = <double> count / num_cells
                        background_count = total_detection_count[gene] - count
                        pair_fold_change = pair_detection_rate * \
                            background_num_cells / background_count
                        detection_rate[cell_type, gene] = pair_detection_rate
                        fold_change[cell_type, gene] = pair_fold_change
                        is_pareto[cell_type, gene] = \
                            (pair_detection_rate >= min_detection_rate) & \
                            (pair_fold_change >= min_fold_change)
            ''')['get_detection_rate_and_fold_change'](
                detection_count=detection_count,
                total_detection_count=total_detection_count,
                num_cells_per_cell_type=num_cells_per_cell_type,
                total_num_cells=total_num_cells,
                min_detection_rate=min_detection_rate,
                min_fold_change=min_fold_change, detection_rate=detection_rate,
                fold_change=fold_change, is_pareto=is_pareto,
                num_threads=num_threads)
        # Get the set of genes on the Pareto front of the two metrics for each
        # cell type; these are the marker genes
        cython_inline(rf'''
            {prange_import}
            def pareto_front(double[:, ::1] detection_rate,
                             double[:, ::1] fold_change,
                             char[:, ::1] is_pareto,
                             int num_threads):
                cdef int gene, other_gene, cell_type
                cdef double gene_detection_rate, gene_fold_change
            
                for cell_type in \
                        {prange('detection_rate.shape[0]', num_threads)}:
                    for gene in range(detection_rate.shape[1]):
                        if not is_pareto[cell_type, gene]:
                            continue
                        gene_detection_rate = detection_rate[cell_type, gene]
                        gene_fold_change = fold_change[cell_type, gene]
                        for other_gene in range(detection_rate.shape[1]):
                            if gene == other_gene or \
                                    not is_pareto[cell_type, other_gene]:
                                continue
                            if gene_detection_rate <= \
                                    detection_rate[cell_type, other_gene] and \
                                    gene_fold_change <= \
                                    fold_change[cell_type, other_gene]:
                                is_pareto[cell_type, gene] = 0
                                break
            ''')['pareto_front'](
                detection_rate=detection_rate, fold_change=fold_change,
                is_pareto=is_pareto, num_threads=num_threads)
        # Return a DataFrame of the selected marker genes, or all genes if
        # `all_genes=True`
        cell_types = groups[cell_type_column_name].rename('cell_type')
        genes = self._var[:, 0].rename('gene')
        if all_genes:
            cell_types = pl.select(pl.lit(cell_types).repeat_by(num_genes))\
                .explode('cell_type')\
                .to_series()
            genes = pl.concat([genes] * num_cell_types)
            return pl.DataFrame((
                cell_types, genes,
                pl.Series('marker', is_pareto.ravel()),
                pl.Series('detection_rate', detection_rate.ravel()),
                pl.Series('fold_change', fold_change.ravel())))
        else:
            cell_type_indices, gene_indices = is_pareto.nonzero()
            return pl.DataFrame((
                cell_types[cell_type_indices], genes[gene_indices],
                pl.Series('detection_rate', detection_rate[
                    cell_type_indices, gene_indices].ravel()),
                pl.Series('fold_change', fold_change[
                    cell_type_indices, gene_indices].ravel())))\
                .select(pl.all().sort_by('detection_rate').over('cell_type'))
    
    def plot_markers(self,
                     genes: str | Iterable[str],
                     cell_type_column: SingleCellColumn,
                     filename: str | Path | None = None,
                     *,
                     QC_column: SingleCellColumn | None = 'passed_QC',
                     ax: 'Axes' | None = None,
                     point_size_multiplier: int | float | np.integer |
                                            np.floating = 250,
                     colormap: str | 'Colormap' | dict[Any, Color] = 'RdBu_r',
                     scatter_kwargs: dict[str, Any] | None = None,
                     legend_kwargs: dict[str, Any] | None = None,
                     colorbar_kwargs: dict[str, Any] | None = None,
                     title: str | None = None,
                     title_kwargs: dict[str, Any] | None = None,
                     xlabel: str | None = None,
                     xlabel_kwargs: dict[str, Any] | None = None,
                     ylabel: str | None = None,
                     ylabel_kwargs: dict[str, Any] | None = None,
                     despine: bool = True,
                     savefig_kwargs: dict[str, Any] | None = None,
                     num_threads: int | np.integer | None = 1) -> None:
        """
        Make a dot plot of the detection rate and fold change of a set of genes
        across cell types - the same metrics used to select marker genes in
        `get_markers()`. This function gives the same result regardless of
        whether it is run before or after normalization.
        
        The size of a gene's dot for a cell type represents its
        "detection rate" in that cell type (the fraction of cells of that type
        that have non-zero count for that gene), while its color represents the
        gene's fold change in detection rate between that cell type and every
        other cell type (by default: more positive fold changes are redder,
        while more negative fold changes are bluer).
        
        Args:
            genes: a list of genes to plot: for instance, marker genes found by
                   `get_markers()`, or marker genes from the literature
            cell_type_column: a column of obs containing cell-type labels. Can
                              be a column name, a polars expression, a polars
                              Series, a 1D NumPy array, or a function that
                              takes in this SingleCell dataset and returns a
                              polars Series or 1D NumPy array.
            QC_column: an optional Boolean column of obs indicating which cells
                       passed QC. Can be a column name, a polars expression, a
                       polars Series, a 1D NumPy array, or a function that
                       takes in this SingleCell dataset and returns a polars
                       Series or 1D NumPy array. Set to None to include all
                       cells. Cells failing QC will be ignored.
            num_threads: the number of threads to use when tabulating each
                         gene's detection rate and fold change. Set
                         `num_threads=None` to use all available cores (as
                         determined by `os.cpu_count()`). For count matrices
                         stored in the usual CSR format, parallelization of the
                         most time-consuming step (calculating detection counts
                         in each cell type) takes place across cell types, so
                         specifying more cores than the number of cell types
                         may not improve peformance.
        
        Note:
            This function may give an incorrect output if the count matrix
            contains explicit zeros (i.e. if `(X.data == 0).any()`): this is
            not checked for, due to speed considerations. In the unlikely event
            that your dataset contains explicit zeros, remove them by running
            `X.eliminate_zeros()` (an in-place operation).
        """
        import matplotlib.pyplot as plt
        X = self._X
        # Get `genes` as a polars Series of the same dtype as var_names; make
        # sure all its entries are unique and present in var_names
        genes = to_tuple(genes)
        check_types(genes, 'genes', str, 'strings')
        genes = pl.Series(genes)
        if genes.n_unique() < len(genes):
            error_message = 'genes contains duplicates'
            raise ValueError(error_message)
        var_names = self._var[:, 0]
        if not genes.is_in(var_names).all():
            if not genes.is_in(var_names).any():
                error_message = \
                    'none of the specified genes were found in var_names'
                raise ValueError(error_message)
            else:
                for gene in genes:
                    if gene not in var_names:
                        error_message = (
                            f'one of the specified genes, {gene!r}, was not '
                            f'found in var_names')
                        raise ValueError(error_message)
        if var_names.dtype != pl.String:
            genes = genes.cast(var_names.dtype)
        # Get the QC column, if not None
        if QC_column is not None:
            QC_column = self._get_column(
                'obs', QC_column, 'QC_column', pl.Boolean,
                allow_missing=QC_column == 'passed_QC')
        # Get the cell-type column
        cell_type_column = \
            self._get_column('obs', cell_type_column, 'cell_type_column',
                             (pl.String, pl.Categorical, pl.Enum),
                             QC_column=QC_column)
        cell_type_column_name = cell_type_column.name
        # Check that `num_threads` is a positive integer or None; if None, set
        # to `os.cpu_count()`
        if num_threads is None:
            num_threads = os.cpu_count()
        else:
            check_type(num_threads, 'num_threads', int, 'a positive integer')
            check_bounds(num_threads, 'num_threads', 1)
        # Get the indices corresponding to the cells of each cell type
        # noinspection PyUnboundLocalVariable
        groups = (cell_type_column.to_frame()
                  .with_columns(
                      _SingleCell_group_indices=pl.int_range(pl.len(),
                                                             dtype=pl.Int32))
                  if QC_column is None else
                  pl.DataFrame((cell_type_column, QC_column))
                  .with_columns(
                      _SingleCell_group_indices=pl.int_range(pl.len(),
                                                             dtype=pl.Int32))
                  .filter(QC_column.name))\
            .group_by(cell_type_column_name, maintain_order=True)\
            .agg('_SingleCell_group_indices', _SingleCell_num_cells=pl.len())\
            .sort(cell_type_column_name)
        # Get a cell-type-by-gene matrix of the number of cells of each type
        # with non-zero expression of each gene in `genes`, i.e. the gene's
        # detection count in that cell type
        num_cell_types = len(groups)
        if num_cell_types == 1:
            error_message = 'cell_type_column only contains one unique value'
            raise ValueError(error_message)
        num_genes = len(genes)
        detection_count = np.zeros((num_cell_types, num_genes), dtype=np.int32)
        prange_import = \
            'from cython.parallel cimport prange' if num_threads > 1 else ''
        if isinstance(X, csr_array):
            group_indices = \
                groups['_SingleCell_group_indices'].explode().to_numpy()
            group_ends = \
                groups['_SingleCell_num_cells'].cum_sum().to_numpy()
            # Get an array mapping each gene in `var_names` to its position in
            # `genes` (-1 if missing from `genes`)
            gene_indices = var_names\
                .to_frame()\
                .join(genes
                      .to_frame(var_names.name)
                      .with_columns(index=pl.int_range(pl.len(),
                                                       dtype=pl.Int32)),
                      on=var_names.name, how='left')\
                .select('index')\
                .to_series()
            gene_indices = gene_indices.fill_null(-1).to_numpy()
            cython_inline(f'''
                {prange_import}
                def groupby_getnnz_csr(
                        const {cython_type(X.indices.dtype)}[::1] indices,
                        const {cython_type(X.indptr.dtype)}[::1] indptr,
                        const int[::1] group_indices,
                        const unsigned[::1] group_ends,
                        const int[::1] gene_indices,
                        int[:, ::1] result,
                        const unsigned num_threads):
                    cdef int num_groups = group_ends.shape[0]
                    cdef int group, gene, row, column
                    cdef unsigned cell
                    # For each group (cell type)...
                    for group in {prange('num_groups', num_threads)}:
                        # For each cell within this group...
                        for cell in range(
                                0 if group == 0 else group_ends[group - 1],
                                group_ends[group]):
                            # Get this cell's row index in the sparse matrix
                            row = group_indices[cell]
                            # For each gene (column) that's non-zero for this
                            # cell...
                            for gene in range(indptr[row], indptr[row + 1]):
                                # Get this gene's column index in `result`
                                # (-1 if the gene is not in `genes`)
                                column = gene_indices[indices[gene]]
                                if column == -1: continue
                                # Add 1 to the total for this group and gene
                                result[group, column] += 1
            ''')['groupby_getnnz_csr'](
                indices=X.indices, indptr=X.indptr,
                group_indices=group_indices, group_ends=group_ends,
                gene_indices=gene_indices, result=detection_count,
                num_threads=num_threads)
        else:
            group_indices = pl.int_range(X.shape[0], dtype=pl.Int32,
                                         eager=True)\
                .to_frame('_SingleCell_group_indices')\
                .join(groups
                      .select('_SingleCell_group_indices',
                              _SingleCell_index=pl.int_range(pl.len(),
                                                             dtype=pl.Int32))
                      .explode('_SingleCell_group_indices'),
                      on='_SingleCell_group_indices', how='left')\
                ['_SingleCell_index']
            if QC_column is not None:
                group_indices = group_indices.fill_null(-1)
                continue_if_missing = 'if group == -1: continue'
            else:
                continue_if_missing = ''
            group_indices = group_indices.to_numpy()
            # Get an array mapping each gene in `genes` to its position in
            # `var_names`
            gene_indices = genes\
                .to_frame(var_names.name)\
                .join(var_names
                      .to_frame()
                      .with_columns(index=pl.int_range(pl.len(),
                                                       dtype=pl.Int32)),
                      on=var_names.name, how='left')\
                .select('index')\
                .to_series()
            gene_indices = gene_indices.fill_null(-1).to_numpy()
            cython_inline(f'''
                {prange_import}
                def groupby_getnnz_csc(
                        const {cython_type(X.indices.dtype)}[::1] indices,
                        const {cython_type(X.indptr.dtype)}[::1] indptr,
                        const int[::1] group_indices,
                        const int[::1] gene_indices,
                        int[:, ::1] result,
                        const unsigned num_threads):
                    cdef int gene, cell, group
                    # For each gene...
                    for gene_index in {prange('result.shape[1]', num_threads)}:
                        # Get the index of this gene in the count matrix
                        gene = gene_indices[gene_index]
                        # For each cell (row) that's non-zero for this gene...
                        for cell in range(indptr[gene], indptr[gene + 1]):
                            # Get the group index for this cell (-1 if it
                            # failed QC)
                            group = group_indices[indices[cell]]
                            {continue_if_missing}
                            # Add 1 to the total for this group and gene
                            result[group, gene] += 1
                    ''')['groupby_getnnz_csc'](
                        indices=X.indices, indptr=X.indptr,
                        group_indices=group_indices, gene_indices=gene_indices,
                        result=detection_count, num_threads=num_threads)
        # For each cell type, calculate the detection rate and the fold change
        # of the detection rate. Also, initialize the candidate set of points
        # on the Pareto front to those with detection rate of at least
        # `min_detection_rate` and fold change of at least `min_fold_change`
        total_detection_count = detection_count.sum(axis=0, dtype=np.int32)
        num_cells_per_cell_type = groups['_SingleCell_num_cells'].to_numpy()
        total_num_cells = num_cells_per_cell_type.sum()
        detection_rate = np.empty((num_cell_types, num_genes), dtype=float)
        fold_change = np.empty((num_cell_types, num_genes), dtype=float)
        cython_inline(rf'''
            {prange_import}
            def get_detection_rate_and_fold_change(
                    const int[:, ::1] detection_count,
                    const int[::1] total_detection_count,
                    const unsigned[::1] num_cells_per_cell_type,
                    const unsigned total_num_cells,
                    double[:, ::1] detection_rate,
                    double[:, ::1] fold_change,
                    const unsigned num_threads):
                cdef int cell_type, gene
                cdef int count, background_count
                cdef unsigned num_cells, background_num_cells
                cdef double pair_detection_rate, pair_fold_change
                for cell_type in \
                        {prange('detection_count.shape[0]', num_threads)}:
                    num_cells = num_cells_per_cell_type[cell_type]
                    background_num_cells = total_num_cells - num_cells
                    for gene in range(detection_count.shape[1]):
                        count = detection_count[cell_type, gene]
                        pair_detection_rate = <double> count / num_cells
                        background_count = total_detection_count[gene] - count
                        pair_fold_change = pair_detection_rate * \
                            background_num_cells / background_count
                        detection_rate[cell_type, gene] = pair_detection_rate
                        fold_change[cell_type, gene] = pair_fold_change
            ''')['get_detection_rate_and_fold_change'](
                detection_count=detection_count,
                total_detection_count=total_detection_count,
                num_cells_per_cell_type=num_cells_per_cell_type,
                total_num_cells=total_num_cells, detection_rate=detection_rate,
                fold_change=fold_change, num_threads=num_threads)
        
        # Plot
        cell_types = groups[cell_type_column_name]
        fig = plt.figure(constrained_layout=True)
        ax = plt.gca()
        x, y = np.meshgrid(range(len(genes)), range(len(cell_types)))
        sizes = detection_rate * point_size_multiplier
        
        # Plot the circles
        scatter = ax.scatter(x.ravel(), y.ravel(), s=sizes.ravel(),
                             c=np.log2(fold_change), cmap=colormap,
                             linewidths=0)
        ax.set_aspect('equal')
        
        # Set x and y limits
        padding = 0.5
        ax.set_xlim([-padding, len(genes) - 1 + padding])
        ax.set_ylim([-padding, len(cell_types) - 1 + padding])
        
        # Add x and y ticks and tick labels
        ax.set_yticks(range(len(cell_types)), cell_types)
        ax.set_xticks(range(len(genes)), genes, rotation=30, ha='right',
                      rotation_mode='anchor', position=(0, 0.005))
        if xlabel is not None:
            ax.set_xlabel(xlabel)
        if ylabel is not None:
            ax.set_ylabel(ylabel)
        
        # Add a legend for detection rate; markers should be at intervals of X%
        # (X%, 2X%, 3X%, ...) up to the maximum detection rate (rounded up to
        # the nearest X%).
        max_detection_rate = detection_rate.max()
        interval = 0.2 if max_detection_rate > 0.5 else \
            0.1 if max_detection_rate > 0.1 else 0.05
        max_detection_rate = \
            np.ceil(max_detection_rate / interval) * interval
        legend_point_sizes = \
            np.arange(interval, max_detection_rate + interval, interval)
        legend_elements = [
            plt.Line2D([0], [0], label=f'{100 * size:.0f}%',
                       markersize=np.sqrt(size * point_size_multiplier),
                       marker='o', linestyle='None', markerfacecolor='black',
                       markeredgecolor='None') for size in legend_point_sizes]
        legend = ax.legend(handles=legend_elements, title='Detection rate',
                           bbox_to_anchor=(0.825, 0.975),
                           bbox_transform=fig.transFigure, loc='upper left',
                           frameon=False)
                
        # Add a colorbar for fold change; labels should be at powers of 2
        cbar_ax = plt.gcf().add_axes((0.912, 0.1, 0.04, 0.4))
        cbar = plt.gcf().colorbar(scatter, cax=cbar_ax)
        cbar.ax.set_title('Fold change', size='medium')
        cbar.ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
        cbar.ax.yaxis.set_major_formatter(plt.FuncFormatter(
            lambda x, pos: f'{2 ** x:.2f}'.rstrip('0').rstrip('.')))
        
        # Save
        plt.savefig('foo.pdf')
        plt.close()
        
    def embed(self,
              *,
              QC_column: SingleCellColumn | None = 'passed_QC',
              PC_key: str = 'PCs',
              embedding_key: str | None = 'PaCMAP',
              num_neighbors: int | np.integer = 10,
              num_extra_neighbors: int | np.integer = 50,
              num_clusters: int | np.integer | None = None,
              num_probes: int | np.integer | None = None,
              num_mid_near_pairs: int | np.integer = 5,
              num_further_pairs: int | np.integer = 20,
              num_iterations: int | np.integer |
                              tuple[int | np.integer, int | np.integer,
                                    int | np.integer] | None = None,
              learning_rate: int | float | np.integer | np.floating = 1.0,
              seed: int | np.integer = 0,
              num_threads: int | np.integer | None = 1,
              faster_single_threading: bool = False,
              overwrite: bool = False,
              verbose: bool = True) -> SingleCell:
        """
        Calculate a two-dimensional embedding of this SingleCell dataset
        suitable for plotting with `plot_embedding()`.
        
        Uses PaCMAP (Pairwise Controlled Manifold Approximation;
        github.com/YingfanWang/PaCMAP; arxiv.org/pdf/2012.04456), a faster UMAP
        alternative that also captures global structure better.
        
        This function is intended to be run after `PCA()`; by default, it uses
        `obsm['PCs']` as the input to PaCMAP, and stores the output in
        `obsm['PaCMAP']` as a num_cells × 2 NumPy array.
        
        Args:
            QC_column: an optional Boolean column of obs indicating which cells
                       passed QC. Can be a column name, a polars expression, a
                       polars Series, a 1D NumPy array, or a function that
                       takes in this SingleCell dataset and returns a polars
                       Series or 1D NumPy array. Set to None to include all
                       cells. Cells failing QC will be ignored and have their
                       embeddings set to NaN.
            PC_key: the key of obsm containing a NumPy array of the principal
                    components (or harmonized principal components) to embed
            embedding_key: the key of obsm where a NumPy array of the
                           embeddings will be stored
            num_neighbors: the number of nearest neighbors to consider for
                           local structure preservation
            num_extra_neighbors: the number of extra nearest neighbors (on top
                                 of `num_neighbors`) to search for initially,
                                 before pruning to the `num_neighbors` of these
                                 `num_neighbors + num_extra_neighbors` cells
                                 with the smallest scaled distances. For a pair
                                 of cells `i` and `j`, the scaled distance
                                 between `i` and `j` is its squared Euclidean
                                 distance, divided by `i`'s average Euclidean
                                 distance to its 3rd, 4th, and 5th nearest
                                 neighbors, divided by `j`'s average Euclidean
                                 distance to its 3rd, 4th, and 5th nearest
                                 neighbors. Must be a positive integer or 0.
            num_clusters: the number of k-means clusters to use during the
                          nearest-neighbor search. Called `nlist` internally by
                          faiss. Must be between 1 and the number of cells. If
                          None, will be set to
                          `ceil(min(sqrt(num_cells), num_cells / 100))`
                          clusters, i.e. the minimum of the square root of the 
                          number of cells and 1% of the number of cells, 
                          rounding up. The core of the heuristic, 
                          `sqrt(num_cells)`, is on the order of the range
                          recommended by faiss, 4 to 16 times the square root
                          (github.com/facebookresearch/faiss/wiki/
                          Guidelines-to-choose-an-index). However, faiss also
                          recommends using between 39 and 256 data points per
                          centroid when training the k-means clustering used
                          in the k-nearest neighbors search. If there are more
                          than 256, the dataset is automatically subsampled
                          for the k-means step, but if there are fewer than 39,
                          faiss gives a warning. To avoid this warning, we
                          switch to using `num_cells / 100` centroids for small
                          datasets, since 100 is the midpoint of 39 and 256 in
                          log space.
            num_probes: the number of nearest k-means clusters to search for a
                        given cell's nearest neighbors. Called `nprobe`
                        internally by faiss. Must be between 1 and
                        `num_clusters`, and should generally be a small
                        fraction of `num_clusters`. If None, will be set to
                        `min(num_clusters, 10)`.
            num_mid_near_pairs: the number of mid-near pairs to consider for
                                global structure preservation
            num_further_pairs: the number of further pairs to consider for
                               local and global structure preservation
            num_iterations: the number of iterations/epochs to run PaCMAP for.
                            Can be a length-3 tuple of the number of iterations
                            for each of the 3 stages of PaCMAP, or a single
                            integer of the number of iterations for the third
                            stage (in which case the number of iterations for
                            the first two stages will be set to 100).
            learning_rate: the learning rate of the Adam optimizer for PaCMAP
            seed: the random seed to use for nearest-neighbor finding and
                  PaCMAP
            num_threads: the number of cores to run nearest-neighbor finding
                         and PaCMAP on. Set `num_threads=None` to use all
                         available cores (as determined by `os.cpu_count()`).
            faster_single_threading: if True, use a different order of
                                     operations for single-threaded PaCMAP,
                                     which gives a modest (~15%) boost in
                                     single-threaded performance at the cost of
                                     no longer matching the embedding produced
                                     by the multithreaded version (due to
                                     differences in floating-point round-off
                                     arising from the different order of
                                     operations). Must be False unless
                                     `num_threads=1`.
            overwrite: if True, overwrite `embedding_key` if already present in
                       obsm, instead of raising an error
            verbose: whether to print details of the nearest-neighbors search
                     and PaCMAP construction
        
        Returns:
            A new SingleCell dataset with the PaCMAP embedding added as
            `obsm[embedding_key]`.
        
        Note:
            PaCMAP's original implementation assumes generic input data, so it
            initializes the embedding by standardizing the input data, running
            PCA on it, and taking the first two PCs. Because our input data is
            already PCs (or harmonized PCs), we avoid redundancy by omitting
            this step and initializing the embedding with the first two columns
            of our input data, i.e. the first two PCs.
        """
        with ignore_sigint():
            import faiss
        # Get the QC column, if not None
        if QC_column is not None:
            QC_column = self._get_column(
                'obs', QC_column, 'QC_column', pl.Boolean,
                allow_missing=QC_column == 'passed_QC')
        # Check that `PC_key` is the name of a key in obsm
        check_type(PC_key, 'PC_key', str, 'a string')
        if PC_key not in self._obsm:
            error_message = f'PC_key {PC_key!r} is not a key of obsm'
            if PC_key == 'PCs':
                error_message += \
                    '; did you forget to run PCA() before embed()?'
            raise ValueError(error_message)
        # Get PCs, for QCed cells only if `QC_column` is not None; check that
        # the PCs are float64 and C-contiguous
        PCs = self._obsm[PC_key]
        dtype = PCs.dtype
        if dtype != float:
            error_message = (
                f'obsm[{PC_key!r}].dtype is {dtype!r}, but must be float64')
            raise TypeError(error_message)
        if not PCs.flags['C_CONTIGUOUS']:
            error_message = (
                f'obsm[{PC_key!r}].dtype is not C-contiguous for at least '
                f'one dataset; make it C-contiguous with '
                f'np.ascontiguousarray(dataset.obsm[{PC_key!r}])')
            raise ValueError(error_message)
        if QC_column is not None:
            QCed_NumPy = QC_column.to_numpy()
            PCs = PCs[QCed_NumPy]
        # Check that `overwrite` is Boolean
        check_type(overwrite, 'overwrite', bool, 'Boolean')
        # Check that `embedding_key` is a string and not already a key in obsm
        check_type(embedding_key, 'embedding_key', str, 'a string')
        if not overwrite and embedding_key in self._obsm:
            error_message = (
                f'embedding_key {embedding_key!r} is already a key of obsm; '
                f'did you already run embed()? Set overwrite=True to '
                f'overwrite.')
            raise ValueError(error_message)
        # Check that `num_neighbors`, `num_mid_near_pairs` and
        # `num_further_pairs` are positive integers
        for variable, variable_name in (
                (num_neighbors, 'num_neighbors'),
                (num_mid_near_pairs, 'num_mid_near_pairs'),
                (num_further_pairs, 'num_further_pairs')):
            check_type(variable, variable_name, int, 'a positive integer')
            check_bounds(variable, variable_name, 1)
        # Check that `num_extra_neighbors` is a positive integer or 0
        check_type(num_extra_neighbors, 'num_extra_neighbors', int,
                   'a positive integer or 0')
        check_bounds(num_extra_neighbors, 'num_extra_neighbors', 0)
        # Get the number of total neighbors to search for initially. We want
        # each cell's `num_neighbors + num_extra_neighbors`-nearest neighbors,
        # but we add 1 extra neighbor for the cell itself.
        num_total_neighbors = num_neighbors + num_extra_neighbors + 1
        # Check that `num_clusters` is between 1 and the number of cells, and
        # that `num_probes` is between 1 and the number of clusters. If either
        # is None, set them to their default values.
        num_cells = len(PCs)
        if num_clusters is None:
            num_clusters = int(np.ceil(min(np.sqrt(num_cells), 
                                           num_cells / 100)))
        else:
            check_type(num_clusters, 'num_clusters', int, 'a positive integer')
            if not 1 <= num_clusters <= num_cells:
                error_message = (
                    f'num_clusters is {num_clusters:,}, but must be ≥ 1 and ≤ '
                    f'the number of cells ({num_cells:,})')
                raise ValueError(error_message)
        if num_probes is None:
            num_probes = min(num_clusters, 10)
        else:
            check_type(num_probes, 'num_probes', int, 'a positive integer')
            if not 1 <= num_probes <= num_clusters:
                error_message = (
                    f'num_probes is {num_probes:,}, but must be ≥ 1 and ≤ '
                    f'num_clusters ({num_clusters:,})')
                raise ValueError(error_message)
        # Calculate the number of nearest neighbors, if not specified
        if num_neighbors is None:
            num_neighbors = 10 if num_cells <= 10_000 else \
                int(round(10 + 15 * (np.log10(num_cells) - 4)))
        if num_cells <= num_neighbors:
            error_message = (
                f'the number of cells ({num_cells:,}) must be greater than '
                f'num_neighbors ({num_neighbors:,})')
            raise ValueError(error_message)
        # Check that `num_iterations` is an integer or length-3 tuple of
        # integers, or None
        if num_iterations is not None:
            check_type(num_iterations, 'num_iterations', (int, tuple),
                       'a positive integer or length-3 tuple of positive '
                       'integers')
            if isinstance(num_iterations, tuple):
                if len(num_iterations) != 3:
                    error_message = (
                        f'num_iterations must be a positive integer or '
                        f'length-3 tuple of positive integers, but has length '
                        f'{len(num_iterations):,}')
                    raise ValueError(error_message)
                for step, step_num_iterations in enumerate(num_iterations):
                    check_type(step_num_iterations,
                               f'num_iterations[{step!r}]', int,
                               'a positive integer')
                    check_bounds(step_num_iterations,
                                 f'num_iterations[{step!r}]', 1)
            else:
                check_bounds(num_iterations, 'num_iterations', 1)
                num_iterations = 100, 100, num_iterations
        else:
            num_iterations = 100, 100, 250
        # Check that `learning_rate` is a positive floating-point number
        check_type(learning_rate, 'learning_rate', (int, float),
                   'a positive number')
        check_bounds(learning_rate, 'learning_rate', 0, left_open=True)
        # Check that `seed` is an integer
        check_type(seed, 'seed', int, 'an integer')
        # Check that `num_threads` is a positive integer or None; if None, set
        # to `os.cpu_count()`. Set this as the number of threads for faiss.
        if num_threads is None:
            num_threads = os.cpu_count()
        else:
            check_type(num_threads, 'num_threads', int, 'a positive integer')
            check_bounds(num_threads, 'num_threads', 1)
        faiss.omp_set_num_threads(num_threads)
        # Check that `faster_single_threading` is Boolean, and False unless
        # `num_threads` is 1
        check_type(faster_single_threading, 'faster_single_threading', bool,
                   'Boolean')
        if faster_single_threading and num_threads != 1:
            error_message = \
                'faster_single_threading must be False unless num_threads is 1'
            raise ValueError(error_message)
        # Check that `verbose` is Boolean
        check_type(verbose, 'verbose', bool, 'Boolean')
        # Define Cython functions
        scaled_distance_prange = \
            prange('scaled_distances.shape[0]', num_threads, nogil=False)
        cython_functions = cython_inline(rf'''
        from cython.parallel cimport threadid, prange
        from libc.math cimport sqrt
        from libc.stdlib cimport free, malloc
        from libc.string cimport memcpy
        
        cdef extern from "limits.h":
            cdef int INT_MAX
        
        cdef inline int rand(long* state) noexcept nogil:
            cdef long x = state[0]
            state[0] = x * 6364136223846793005L + 1442695040888963407L
            cdef int s = (x ^ (x >> 18)) >> 27
            cdef int rot = x >> 59
            return (s >> rot) | (s << ((-rot) & 31))
        
        cdef inline long srand(const long seed) noexcept nogil:
            cdef long state = seed + 1442695040888963407L
            rand(&state)
            return state
        
        cdef inline int randint(const int bound, long* state) noexcept nogil:
            cdef int r, threshold = -bound % bound
            while True:
                r = rand(state)
                if r >= threshold:
                    return r % bound

        cdef void quicksort(const double[::1] arr,
                            int[::1] indices,
                            int left,
                            int right) noexcept nogil:
            cdef double pivot_value
            cdef int pivot_index, mid, i, temp
    
            while left < right:
                mid = left + (right - left) // 2
                if arr[indices[mid]] < arr[indices[left]]:
                    temp = indices[left]
                    indices[left] = indices[mid]
                    indices[mid] = temp
                if arr[indices[right]] < arr[indices[left]]:
                    temp = indices[left]
                    indices[left] = indices[right]
                    indices[right] = temp
                if arr[indices[right]] < arr[indices[mid]]:
                    temp = indices[mid]
                    indices[mid] = indices[right]
                    indices[right] = temp
    
                pivot_value = arr[indices[mid]]
                temp = indices[mid]
                indices[mid] = indices[right]
                indices[right] = temp
                pivot_index = left
    
                for i in range(left, right):
                    if arr[indices[i]] < pivot_value:
                        temp = indices[i]
                        indices[i] = indices[pivot_index]
                        indices[pivot_index] = temp
                        pivot_index += 1
    
                temp = indices[right]
                indices[right] = indices[pivot_index]
                indices[pivot_index] = temp
    
                if pivot_index - left < right - pivot_index:
                    quicksort(arr, indices, left, pivot_index - 1)
                    left = pivot_index + 1
                else:
                    quicksort(arr, indices, pivot_index + 1, right)
                    right = pivot_index - 1
    
        cdef inline void argsort(const double[::1] arr,
                                 int[::1] indices) noexcept nogil:
            cdef int i
            for i in range(indices.shape[0]):
                indices[i] = i
            quicksort(arr, indices, 0, indices.shape[0] - 1)
        
        def remove_self_neighbors(long[:, ::1] neighbors,
                                  const unsigned num_threads):
            cdef int i, j
            for i in {prange('neighbors.shape[0]', num_threads)}:
                # If the cell is its own nearest neighbor (almost always), skip
                if neighbors[i, 0] == i:
                    continue
                # Find the position where the cell is listed as its own
                # self-neighbor
                for j in range(1, neighbors.shape[1]):
                    if neighbors[i, j] == i:
                        break
                # Shift all neighbors before it to the right, overwriting it
                while j > 0:
                    neighbors[i, j] = neighbors[i, j - 1]
                    j = j - 1
        
        def get_scaled_distances(const double[:, ::1] X,
                                 const long[:, :] neighbors,
                                 double[:, ::1] scaled_distances,
                                 const unsigned num_threads):
            cdef int i, j, k
            cdef double* sig = \
                <double*> malloc(scaled_distances.shape[0] * sizeof(double))
            if not sig:
                raise MemoryError
            try:
                {"with nogil" if num_threads > 1 else "if True"}:
                    for i in {scaled_distance_prange}:
                        for j in range(scaled_distances.shape[1]):
                            scaled_distances[i, j] = 0
                            for k in range(X.shape[1]):
                                scaled_distances[i, j] = \
                                    scaled_distances[i, j] + \
                                    (X[i, k] - X[j, k]) ** 2
                        sig[i] = (sqrt(scaled_distances[i, 3]) +
                                  sqrt(scaled_distances[i, 4]) +
                                  sqrt(scaled_distances[i, 5])) / 3
                        if sig[i] < 1e-10:
                            sig[i] = 1e-10
                    
                    for i in {scaled_distance_prange}:
                        for j in range(scaled_distances.shape[1]):
                            scaled_distances[i, j] = scaled_distances[i, j] / \
                                sig[i] / sig[neighbors[i, j]]
            finally:
                free(sig)
    
        def get_neighbor_pairs(const double[:, ::1] X,
                               const double[:, ::1] scaled_distances,
                               const long[:, :] neighbors,
                               int[:, ::1] neighbor_pairs,
                               const unsigned num_threads):
            cdef int i, j, thread_id
            cdef int num_neighbors = neighbor_pairs.shape[1]
            cdef int num_total_neighbors = scaled_distances.shape[1]
            cdef int[::1] indices
            cdef int* indices_buffer = <int*> malloc(
                num_total_neighbors * num_threads * sizeof(int))
            if not indices_buffer:
                raise MemoryError
            try:
                indices = \
                    <int[:num_total_neighbors * num_threads:]> indices_buffer
                if num_threads == 1:
                    for i in range(neighbor_pairs.shape[0]):
                        argsort(scaled_distances[i], indices)
                        for j in range(num_neighbors):
                            neighbor_pairs[i, j] = neighbors[i, indices[j]]
                else:
                    for i in prange(X.shape[0], num_threads=num_threads,
                                    nogil=True):
                        thread_id = threadid()
                        argsort(scaled_distances[i],
                                indices[thread_id * num_total_neighbors:
                                        thread_id * num_total_neighbors +
                                        num_total_neighbors])
                        for j in range(num_neighbors):
                            neighbor_pairs[i, j] = neighbors[i, indices[
                                thread_id * num_total_neighbors + j]]
            finally:
                free(indices_buffer)
        
        def sample_mid_near_pairs(const double[:, ::1] X,
                                  int[:, ::1] mid_near_pairs,
                                  const int seed,
                                  const unsigned num_threads):
            cdef int i, j, k, l, thread_id, n = X.shape[0], \
                closest_cell = -1, second_closest_cell = -1
            cdef double squared_distance, smallest, second_smallest
            cdef long state
            cdef int* sampled = <int*> malloc(6 * num_threads * sizeof(int))
            if not sampled:
                raise MemoryError
            try:
                if num_threads == 1:
                    for i in range(n):
                        state = srand(seed + i)
                        for j in range(mid_near_pairs.shape[1]):
                            # Randomly sample 6 cells (which are not the
                            # current cell) and select the 2nd-closest
                            smallest = INT_MAX
                            second_smallest = INT_MAX
                            for k in range(6):
                                while True:
                                    # Sample a random cell...
                                    sampled[k] = randint(n, &state)
                                    # ...that is not this cell...
                                    if sampled[k] == i:
                                        continue
                                    # ...nor a previously sampled cell
                                    for l in range(k):
                                        if sampled[k] == sampled[l]:
                                            break
                                    else:
                                        break
                            for k in range(6):
                                squared_distance = 0
                                for l in range(X.shape[1]):
                                    squared_distance = squared_distance + \
                                        (X[i, l] - X[sampled[k], l]) ** 2
                                if squared_distance < smallest:
                                    second_smallest = smallest
                                    second_closest_cell = closest_cell
                                    smallest = squared_distance
                                    closest_cell = sampled[k]
                                elif squared_distance < second_smallest:
                                    second_smallest = squared_distance
                                    second_closest_cell = sampled[k]
                            mid_near_pairs[i, j] = second_closest_cell
                else:
                    for i in prange(n, num_threads=num_threads, nogil=True):
                        thread_id = threadid()
                        state = srand(seed + i)
                        for j in range(mid_near_pairs.shape[1]):
                            smallest = INT_MAX
                            second_smallest = INT_MAX
                            for k in range(6 * thread_id, 6 * thread_id + 6):
                                while True:
                                    sampled[k] = randint(n, &state)
                                    if sampled[k] == i:
                                        continue
                                    for l in range(6 * thread_id, k):
                                        if sampled[k] == sampled[l]:
                                            break
                                    else:
                                        break
                            for k in range(6 * thread_id, 6 * thread_id + 6):
                                squared_distance = 0
                                for l in range(X.shape[1]):
                                    squared_distance = squared_distance + \
                                        (X[i, l] - X[sampled[k], l]) ** 2
                                if squared_distance < smallest:
                                    second_smallest = smallest
                                    second_closest_cell = closest_cell
                                    smallest = squared_distance
                                    closest_cell = sampled[k]
                                elif squared_distance < second_smallest:
                                    second_smallest = squared_distance
                                    second_closest_cell = sampled[k]
                            mid_near_pairs[i, j] = second_closest_cell
            finally:
                free(sampled)
        
        def sample_further_pairs(const double[:, ::1] X,
                                 const int[:, ::1] neighbor_pairs,
                                 int[:, ::1] further_pairs,
                                 const int seed,
                                 const unsigned num_threads):
            """Sample Further pairs using the given seed."""
            cdef int i, j, k, n = X.shape[0], further_pair_index
            cdef long state
            for i in {prange('n', num_threads)}:
                state = srand(seed + i)
                for j in range(further_pairs.shape[1]):
                    while True:
                        # Sample a random cell...
                        further_pair_index = randint(n, &state)
                        # ...that is not this cell...
                        if further_pair_index == i:
                            continue
                        # ...nor one of its nearest neighbors...
                        for k in range(neighbor_pairs.shape[1]):
                            if further_pair_index == neighbor_pairs[i, k]:
                                break
                        else:
                            # ...nor a previously sampled cell
                            for k in range(j):
                                if further_pair_index == further_pairs[i, k]:
                                    break
                            else:
                                break
                    further_pairs[i, j] = further_pair_index
        
        def reformat_for_parallel(const int[:, ::1] pairs,
                                  int[::1] pair_indices, int[::1] pair_indptr):
            cdef int i, j, k, dest_index
            cdef int* dest_indices
            # Tabulate how often each cell appears in pairs; at a minimum, it
            # will appear `pairs.shape[1]` times (i.e. the number of
            # neighbors), as the `i` in the pair, but it will also appear a
            # variable number of times as the `j` in the pair.
            pair_indptr[0] = 0
            pair_indptr[1:] = pairs.shape[1]
            for i in range(pairs.shape[0]):
                for k in range(pairs.shape[1]):
                    j = pairs[i, k]
                    pair_indptr[j + 1] += 1
            # Take the cumulative sum of the values in `pair_indptr`
            for i in range(2, pair_indptr.shape[0]):
                pair_indptr[i] += pair_indptr[i - 1]
            # Now that we know how many pairs each cell is a part of, do a
            # second pass over `pairs` to populate `pair_indices` with the
            # pairs' indices. Use a temporary buffer, `dest_indices`, to keep
            # track of the index within `pair_indptr` to write each cell's next
            # pair to.
            dest_indices = <int*> malloc(pairs.shape[0] * sizeof(int))
            if not dest_indices:
                raise MemoryError
            try:
                memcpy(dest_indices, &pair_indptr[0],
                       pairs.shape[0] * sizeof(int))
                for i in range(pairs.shape[0]):
                    for k in range(pairs.shape[1]):
                        j = pairs[i, k]
                        pair_indices[dest_indices[i]] = j
                        pair_indices[dest_indices[j]] = i
                        dest_indices[i] += 1
                        dest_indices[j] += 1
            finally:
                free(dest_indices)

        def get_gradients(const double[:, ::1] embedding,
                          const int[:, ::1] neighbor_pairs,
                          const int[:, ::1] mid_near_pairs,
                          const int[:, ::1] further_pairs,
                          const double w_neighbors,
                          const double w_mid_near,
                          double[:, ::1] gradients):
            cdef int i, j, k
            cdef double distance_ij, embedding_ij_0, embedding_ij_1, w
            gradients[:] = 0
            # Nearest-neighbor pairs
            for i in range(neighbor_pairs.shape[0]):
                for k in range(neighbor_pairs.shape[1]):
                    j = neighbor_pairs[i, k]
                    embedding_ij_0 = embedding[i, 0] - embedding[j, 0]
                    embedding_ij_1 = embedding[i, 1] - embedding[j, 1]
                    distance_ij = 1 + embedding_ij_0 ** 2 + embedding_ij_1 ** 2
                    w = w_neighbors * (20 / (10 + distance_ij) ** 2)
                    gradients[i, 0] += w * embedding_ij_0
                    gradients[j, 0] -= w * embedding_ij_0
                    gradients[i, 1] += w * embedding_ij_1
                    gradients[j, 1] -= w * embedding_ij_1
            # Mid-near pairs
            for i in range(mid_near_pairs.shape[0]):
                for k in range(mid_near_pairs.shape[1]):
                    j = mid_near_pairs[i, k]
                    embedding_ij_0 = embedding[i, 0] - embedding[j, 0]
                    embedding_ij_1 = embedding[i, 1] - embedding[j, 1]
                    distance_ij = 1 + embedding_ij_0 ** 2 + embedding_ij_1 ** 2
                    w = w_mid_near * 20000 / (10000 + distance_ij) ** 2
                    gradients[i, 0] += w * embedding_ij_0
                    gradients[j, 0] -= w * embedding_ij_0
                    gradients[i, 1] += w * embedding_ij_1
                    gradients[j, 1] -= w * embedding_ij_1
            # Further pairs
            for i in range(further_pairs.shape[0]):
                for k in range(further_pairs.shape[1]):
                    j = further_pairs[i, k]
                    embedding_ij_0 = embedding[i, 0] - embedding[j, 0]
                    embedding_ij_1 = embedding[i, 1] - embedding[j, 1]
                    distance_ij = 1 + embedding_ij_0 ** 2 + embedding_ij_1 ** 2
                    w = 2 / (1 + distance_ij) ** 2
                    gradients[i, 0] -= w * embedding_ij_0
                    gradients[j, 0] += w * embedding_ij_0
                    gradients[i, 1] -= w * embedding_ij_1
                    gradients[j, 1] += w * embedding_ij_1
        
        def get_gradients_parallel(const double[:, ::1] embedding,
                                   const int[::1] neighbor_pair_indices,
                                   const int[::1] neighbor_pair_indptr,
                                   const int[::1] mid_near_pair_indices,
                                   const int[::1] mid_near_pair_indptr,
                                   const int[::1] further_pair_indices,
                                   const int[::1] further_pair_indptr,
                                   const double w_neighbors,
                                   const double w_mid_near,
                                   double[:, ::1] gradients,
                                   const unsigned num_threads):
            cdef int i, j, k
            cdef double distance_ij, embedding_ij_0, embedding_ij_1, w
            for i in prange(embedding.shape[0], nogil=True,
                            num_threads=num_threads):
                gradients[i, 0] = 0
                gradients[i, 1] = 0
                # Nearest-neighbor pairs
                for k in range(neighbor_pair_indptr[i],
                               neighbor_pair_indptr[i + 1]):
                    j = neighbor_pair_indices[k]
                    embedding_ij_0 = embedding[i, 0] - embedding[j, 0]
                    embedding_ij_1 = embedding[i, 1] - embedding[j, 1]
                    distance_ij = 1 + embedding_ij_0 ** 2 + embedding_ij_1 ** 2
                    w = w_neighbors * (20 / (10 + distance_ij) ** 2)
                    gradients[i, 0] = gradients[i, 0] + w * embedding_ij_0
                    gradients[i, 1] = gradients[i, 1] + w * embedding_ij_1
                # Mid-near pairs
                for k in range(mid_near_pair_indptr[i],
                               mid_near_pair_indptr[i + 1]):
                    j = mid_near_pair_indices[k]
                    embedding_ij_0 = embedding[i, 0] - embedding[j, 0]
                    embedding_ij_1 = embedding[i, 1] - embedding[j, 1]
                    distance_ij = 1 + embedding_ij_0 ** 2 + embedding_ij_1 ** 2
                    w = w_mid_near * 20000 / (10000 + distance_ij) ** 2
                    gradients[i, 0] = gradients[i, 0] + w * embedding_ij_0
                    gradients[i, 1] = gradients[i, 1] + w * embedding_ij_1
                # Further pairs
                for k in range(further_pair_indptr[i],
                               further_pair_indptr[i + 1]):
                    j = further_pair_indices[k]
                    embedding_ij_0 = embedding[i, 0] - embedding[j, 0]
                    embedding_ij_1 = embedding[i, 1] - embedding[j, 1]
                    distance_ij = 1 + embedding_ij_0 ** 2 + embedding_ij_1 ** 2
                    w = 2 / (1 + distance_ij) ** 2
                    gradients[i, 0] = gradients[i, 0] - w * embedding_ij_0
                    gradients[i, 1] = gradients[i, 1] - w * embedding_ij_1
        
        def update_embedding_adam(double[:, ::1] embedding,
                                  const double[:, ::1] gradients,
                                  double[:, ::1] momentum,
                                  double[:, ::1] velocity,
                                  const double beta1,
                                  const double beta2,
                                  double learning_rate,
                                  const int iteration,
                                  const unsigned num_threads):
            cdef int i
            learning_rate = \
                learning_rate * sqrt(1 - beta2 ** (iteration + 1)) / \
                (1 - beta1 ** (iteration + 1))
            for i in {prange('embedding.shape[0]', num_threads)}:
                momentum[i, 0] += \
                    (1 - beta1) * (gradients[i, 0] - momentum[i, 0])
                velocity[i, 0] += \
                    (1 - beta2) * (gradients[i, 0] ** 2 - velocity[i, 0])
                embedding[i, 0] -= learning_rate * momentum[i, 0] / \
                                   (sqrt(velocity[i, 0]) + 1e-7)
                momentum[i, 1] += \
                    (1 - beta1) * (gradients[i, 1] - momentum[i, 1])
                velocity[i, 1] += \
                    (1 - beta2) * (gradients[i, 1] ** 2 - velocity[i, 1])
                embedding[i, 1] -= learning_rate * momentum[i, 1] / \
                                   (sqrt(velocity[i, 1]) + 1e-7)
            ''')
        remove_self_neighbors = cython_functions['remove_self_neighbors']
        get_scaled_distances = cython_functions['get_scaled_distances']
        get_neighbor_pairs = cython_functions['get_neighbor_pairs']
        sample_mid_near_pairs = cython_functions['sample_mid_near_pairs']
        sample_further_pairs = cython_functions['sample_further_pairs']
        update_embedding_adam = cython_functions['update_embedding_adam']
        # Calculate each cell's `num_neighbors + num_extra_neighbors`-nearest
        # neighbors with faiss. (`num_total_neighbors` is
        # `num_neighbors + num_extra_neighbors + 1`, where the `+ 1` is for the
        # cell itself. We exclude the cell itself below.)
        with Timer('Calculating nearest neighbors', verbose=verbose):
            dim = PCs.shape[1]
            quantizer = faiss.IndexFlatL2(dim)
            quantizer.verbose = verbose
            index = faiss.IndexIVFFlat(quantizer, dim, num_clusters)
            index.cp.seed = seed
            index.verbose = verbose
            index.cp.verbose = verbose
            # noinspection PyArgumentList
            index.train(PCs)
            # noinspection PyArgumentList
            index.add(PCs)
            index.nprobe = num_probes
            # noinspection PyArgumentList
            nearest_neighbor_indices = \
                index.search(PCs, num_total_neighbors)[1]
            # Sometimes there aren't enough nearest neighbors for certain cells
            # with `num_probes` probes; if so, double `num_probes` (and
            # threshold to at most `num_clusters`), then re-run
            # nearest-neighbor finding for those cells
            needs_update = nearest_neighbor_indices[:, -1] == -1
            if needs_update.any():
                needs_update_X = PCs[needs_update]
                while True:
                    num_probes = min(num_probes * 2, num_clusters)
                    if verbose:
                        print(f'{len(needs_update_X):,} cells '
                              f'({len(needs_update_X) / len(self._obs):.2f}%)'
                              f' did not have enough neighbors with '
                              f'{index.nprobe:,} probes; re-running '
                              f'nearest-neighbors finding for these cells '
                              f'with {num_probes:,} probes')
                    index.nprobe = num_probes
                    # noinspection PyArgumentList
                    new_indices = \
                        index.search(needs_update_X, num_total_neighbors)[1]
                    nearest_neighbor_indices[needs_update] = new_indices
                    still_needs_update = new_indices[:, -1] == -1
                    if not still_needs_update.any():
                        break
                    needs_update[needs_update] = still_needs_update
                    needs_update_X = needs_update_X[still_needs_update]
        if verbose:
            percent = \
                (nearest_neighbor_indices[:, 0] == range(num_cells)).mean()
            print(f'{100 * percent:.3f}% of cells are correctly detected '
                  f'as their own nearest neighbors (a measure of the '
                  f'quality of the k-nearest neighbors search)')
        # Remove self-neighbors from each cell's list of nearest neighbors.
        # These are almost always in the 0th column, but occasionally later due
        # to the inaccuracy of the nearest-neighbors search. This leaves us
        # with `num_neighbors + num_extra_neighbors` nearest neighbors.
        remove_self_neighbors(nearest_neighbor_indices, num_threads)
        nearest_neighbor_indices = nearest_neighbor_indices[:, 1:]
        # Get scaled distances between each cell and its nearest neighbors
        scaled_distances = np.empty_like(nearest_neighbor_indices, dtype=float)
        get_scaled_distances(PCs, nearest_neighbor_indices, scaled_distances,
                             num_threads)
        # Select the `num_neighbors` of the `num_total_neighbors`
        # nearest-neighbor pairs with the lowest scaled distances
        neighbor_pairs = np.empty((num_cells, num_neighbors), dtype=np.int32)
        get_neighbor_pairs(PCs, scaled_distances, nearest_neighbor_indices,
                           neighbor_pairs, num_threads)
        del scaled_distances, nearest_neighbor_indices
        # Sample mid-near pairs
        mid_near_pairs = np.empty((num_cells, num_mid_near_pairs),
                                  dtype=np.int32)
        sample_mid_near_pairs(PCs, mid_near_pairs, seed, num_threads)
        # Sample further pairs
        further_pairs = np.empty((num_cells, num_further_pairs),
                                 dtype=np.int32)
        sample_further_pairs(PCs, neighbor_pairs, further_pairs,
                             seed + mid_near_pairs.size, num_threads)
        # If multithreaded, or single-threaded with
        # `faster_single_threading=False`, reformat the three lists of pairs to
        # allow deterministic parallelism. Specifically, transform pairs of
        # cell indices from the original format of a 2D array `pairs` where
        # `pairs[i]` contains all js for which (i, j) is a pair, to a pair of
        # 1D arrays `pair_indices` and `pair_indptr` forming a sparse matrix,
        # where `pair_indices[pair_indptr[i]:pair_indptr[i + 1]]` contains all
        # js for which (i, j) is a pair or (j, i) is a pair. `pair_indices`
        # must have length `2 * pairs.size`, since each pair will appear twice,
        # once for (i, j) and once for (j, i). `pair_indptr` must have length
        # equal to the number of cells plus one, just like for scipy sparse
        # matrices.
        if faster_single_threading:
            get_gradients = cython_functions['get_gradients']
        else:
            reformat_for_parallel = cython_functions['reformat_for_parallel']
            
            neighbor_pair_indices = np.empty(2 * neighbor_pairs.size,
                                              dtype=np.int32)
            neighbor_pair_indptr = np.empty(num_cells + 1, dtype=np.int32)
            reformat_for_parallel(neighbor_pairs, neighbor_pair_indices,
                                  neighbor_pair_indptr)
            del neighbor_pairs
            mid_near_pair_indices = \
                np.empty(2 * mid_near_pairs.size, dtype=np.int32)
            mid_near_pair_indptr = np.empty(num_cells + 1, dtype=np.int32)
            reformat_for_parallel(mid_near_pairs, mid_near_pair_indices,
                                  mid_near_pair_indptr)
            del mid_near_pairs
            further_pair_indices = \
                np.empty(2 * further_pairs.size, dtype=np.int32)
            further_pair_indptr = np.empty(num_cells + 1, dtype=np.int32)
            reformat_for_parallel(further_pairs, further_pair_indices,
                                  further_pair_indptr)
            del further_pairs
            get_gradients = cython_functions['get_gradients_parallel']
        # Initialize the embedding, gradients, and other optimizer parameters
        embedding = 0.01 * PCs[:, :2]
        gradients = np.zeros_like(embedding, dtype=float)
        momentum = np.zeros_like(embedding, dtype=float)
        velocity = np.zeros_like(embedding, dtype=float)
        w_mid_near_init = 1000
        beta1 = 0.9
        beta2 = 0.999
        # Optimize the embedding
        for iteration in range(sum(num_iterations)):
            num_phase_1_iterations, num_phase_2_iterations = num_iterations[:2]
            if iteration < num_phase_1_iterations:
                w_mid_near = \
                    (1 - iteration / num_phase_1_iterations) * \
                    w_mid_near_init + iteration / num_phase_1_iterations * 3
                w_neighbors = 2
            elif iteration < num_phase_1_iterations + num_phase_2_iterations:
                w_mid_near = 3
                w_neighbors = 3
            else:
                w_mid_near = 0
                w_neighbors = 1
            # Calculate gradients
            if faster_single_threading:
                # noinspection PyUnboundLocalVariable
                get_gradients(embedding, neighbor_pairs, mid_near_pairs,
                              further_pairs, w_neighbors, w_mid_near,
                              gradients)
            else:
                # noinspection PyUnboundLocalVariable
                get_gradients(embedding, neighbor_pair_indices,
                              neighbor_pair_indptr, mid_near_pair_indices,
                              mid_near_pair_indptr, further_pair_indices,
                              further_pair_indptr, w_neighbors, w_mid_near,
                              gradients, num_threads)
            # Update the embedding based on the gradients, via the Adam
            # optimizer
            update_embedding_adam(embedding, gradients, momentum, velocity,
                                  beta1, beta2, learning_rate, iteration,
                                  num_threads)
        # If `QC_column` is not None, back-project from QCed cells to all
        # cells, filling with NaN
        if QC_column is not None:
            embedding_QCed = embedding
            embedding = np.full((len(self), embedding_QCed.shape[1]), np.nan)
            # noinspection PyUnboundLocalVariable
            embedding[QCed_NumPy] = embedding_QCed
        # noinspection PyTypeChecker
        return SingleCell(X=self._X, obs=self._obs, var=self._var,
                          obsm=self._obsm | {embedding_key: embedding},
                          varm=self._varm, uns=self._uns)
    
    # noinspection PyUnresolvedReferences
    def plot_embedding(
            self,
            color_column: SingleCellColumn | None,
            filename: str | Path | None = None,
            *,
            cells_to_plot_column: SingleCellColumn | None = 'passed_QC',
            embedding_key: str = 'PaCMAP',
            ax: 'Axes' | None = None,
            point_size: int | float | np.integer | np.floating | str |
                        None = None,
            sort_by_frequency: bool = False,
            palette: str | 'Colormap' | dict[Any, Color] = None,
            palette_kwargs: dict[str, Any] | None = None,
            default_color: Color = 'lightgray',
            scatter_kwargs: dict[str, Any] | None = None,
            label: bool = False,
            label_kwargs: dict[str, Any] | None = None,
            legend: bool = True,
            legend_kwargs: dict[str, Any] | None = None,
            colorbar: bool = True,
            colorbar_kwargs: dict[str, Any] | None = None,
            title: str | None = None,
            title_kwargs: dict[str, Any] | None = None,
            xlabel: str | None = 'Component 1',
            xlabel_kwargs: dict[str, Any] | None = None,
            ylabel: str | None = 'Component 2',
            ylabel_kwargs: dict[str, Any] | None = None,
            xlim: tuple[int | float | np.integer | np.floating,
                        int | float | np.integer | np.floating] | None = None,
            ylim: tuple[int | float | np.integer | np.floating,
                        int | float | np.integer | np.floating] | None = None,
            despine: bool = True,
            savefig_kwargs: dict[str, Any] | None = None) -> None:
        """
        Plot an embedding created by `embed()`, using Matplotlib.
        
        Requires the colorspacious package. Install via:
        mamba install -y colorspacious
        
        Args:
            color_column: an optional column of obs indicating how to color
                          each cell in the plot. Can be a column name, a polars
                          expression, a polars Series, a 1D NumPy array, or a
                          function that takes in this SingleCell dataset and
                          returns a polars Series or 1D NumPy array. Can be
                          discrete (e.g. cell-type labels), specified as a
                          String/Categorical/Enum column, or quantitative (e.g.
                          the number of UMIs per cell), specified as an
                          integer/floating-point column. Missing (null) cells
                          will be plotted with the color `default_color`. Set
                          to None to use `default_color` for all cells.
            filename: the file to save to. If None, generate the plot but do
                      not save it, which allows it to be shown interactively or
                      modified further (e.g. by adding a title or axis labels)
                      before saving.
            cells_to_plot_column: an optional Boolean column of obs indicating
                                  which cells to plot. Can be a column name, a
                                  polars expression, a polars Series, a 1D
                                  NumPy array, or a function that takes in this
                                  SingleCell dataset and returns a polars
                                  Series or 1D NumPy array. Set to None to plot
                                  all cells passing QC.
            embedding_key: the key of obsm containing a NumPy array of the
                           embedding to plot
            ax: the Matplotlib axes to save the plot onto; if None, create a
                new figure with Matpotlib's constrained layout and plot onto it
            point_size: the size of the points for each cell; defaults to
                        30,000 divided by the number of cells, one quarter of
                        scanpy's default. Can be a single number, or the name
                        of a column of obs to make each point a different size.
            sort_by_frequency: if True, assign colors and sort the legend in
                               order of decreasing frequency; if False (the
                               default), use natural sorted order
                               (en.wikipedia.org/wiki/Natural_sort_order).
                               Cannot be True unless `palette` is None and
                               `color_column` is discrete; if `palette` is
                               not None, the plot order is determined by the
                               order of the keys in `palette`.
            palette: a string or Colormap object indicating the Matplotlib
                     colormap to use; or, if `color_column` is discrete, a
                     dictionary mapping values in `color_column` to Matplotlib
                     colors (cells with values of `color_column` that are not
                     in the dictionary will be plotted in the color
                     `default_color`). Defaults to `plt.rcParams['image.cmap']`
                     (`'viridis'` by default) if `color_column` is continous,
                     or the colors from `generate_palette()` if `color_column`
                     is discrete (with colors assigned in decreasing order of
                     frequency). Cannot be specified if `color_column` is None.
            palette_kwargs: a dictionary of keyword arguments to be passed to
                            `generate_palette()`. Can only be specified when
                            `color_column` is discrete and `palette` is None.
            default_color: the default color to plot cells in when
                           `color_column` is None, or when certain cells have
                           missing (null) values for `color_column`, or when
                           `palette` is a dictionary and some cells have values
                           of `color_column` that are not in the dictionary
            scatter_kwargs: a dictionary of keyword arguments to be passed to
                            `ax.scatter()`, such as:
                            - `rasterized`: whether to convert the scatter plot
                              points to a raster (bitmap) image when saving to
                              a vector format like PDF. Defaults to True,
                              instead of the Matplotlib default of False.
                            - `marker`: the shape to use for plotting each cell
                            - `norm`, `vmin`, and `vmax`: control how the
                              numbers in `color_column` are converted to
                              colors, if `color_column` is numeric
                            - `alpha`: the transparency of each point
                            - `linewidths` and `edgecolors`: the width and
                              color of the borders around each marker. These
                              are absent by default (`linewidths=0`), unlike
                              Matplotlib's default. Both arguments can be
                              either single values or sequences.
                            - `zorder`: the order in which the cells are
                              plotted, with higher values appearing on top of
                              lower ones.
                            Specifying `s`, `c`/`color`, or `cmap` will raise
                            an error, since these arguments conflict with the
                            `point_size`, `color_column`, and `colormap`
                            arguments, respectively.
            label: whether to label cells with each distinct value of
                   `color_column`. Labels will be placed at the median x and y
                   position of the points with that color. Can only be True
                   when `color_column` is discrete. When set to True, you may
                   also want to set `legend=False` to avoid redundancy.
            label_kwargs: a dictionary of keyword arguments to be passed to
                          `ax.text()` when adding labels to control the text
                          properties, such as:
                           - `color` and `size` to modify the text color/size
                           - `verticalalignment` and `horizontalalignment` to
                             control vertical and horizontal alignment. By
                             default, unlike Matplotlib's default, these are
                             both set to `'center'`.
                           - `path_effects` to set properties for the border
                             around the text. By default, set to
                             `matplotlib.patheffects.withStroke(
                                  linewidth=3, foreground='white', alpha=0.75)`
                             instead of Matplotlib's default of None, to put a
                             semi-transparent white border around the labels
                             for better contrast.
                          Can only be specified when `label=True`.
            legend: whether to add a legend for each value in `color_column`.
                    Ignored unless `color_column` is discrete.
            legend_kwargs: a dictionary of keyword arguments to be passed to
                           `ax.legend()` to modify the legend, such as:
                           - `loc`, `bbox_to_anchor`, and `bbox_transform` to
                             set its location. By default, `loc` is set to
                             `center left` and `bbox_to_anchor` to `(1, 0.5)`
                             to put the legend to the right of the plot,
                             anchored at the middle.
                           - `ncols` to set its number of columns. By
                             default, set to
                             `obs[color_column].n_unique() // 16 + 1` to have
                             at most 16 items per column.
                           - `prop`, `fontsize`, and `labelcolor` to set its
                             font properties
                           - `facecolor` and `framealpha` to set its background
                             color and transparency
                           - `frameon=True` or `edgecolor` to add or color its
                             border (`frameon` is False by default, unlike
                             Matplotlib's default of True)
                           - `title` to add a legend title
                           Can only be specified when `color_column` is
                           discrete and `legend=True`.
            colorbar: whether to add a colorbar. Ignored unless `color_column`
                      is quantitative.
            colorbar_kwargs: a dictionary of keyword arguments to be passed to
                             `ax.colorbar()`, such as:
                             - `location`: `'left'`, `'right'`, `'top'`, or
                               `'bottom'`
                             - `orientation`: `'vertical'` or `'horizontal'`
                             - `fraction`: the fraction of the axes to
                               allocate to the colorbar (default: 0.15)
                             - `shrink`: the fraction to multiply the size of
                               the colorbar by (default: 0.5, unlike
                               Matplotlib's default of 1)
                             - `aspect`: the ratio of the colorbar's long to
                               short dimensions (default: 20)
                             - `pad`: the fraction of the axes between the
                               colorbar and the rest of the figure (default:
                               0.05 if vertical, 0.15 if horizontal)
                             Can only be specified when `color_column` is
                             quantitative and `colorbar=True`.
            title: the title of the plot, or None to not add a title
            title_kwargs: a dictionary of keyword arguments to be passed to
                          ax.title() to modify the title; see `label_kwargs`
                          for examples
            xlabel: the x-axis label, or None to not label the x-axis
            xlabel_kwargs: a dictionary of keyword arguments to be passed to
                           ax.set_xlabel() to modify the x-axis label
            ylabel: the y-axis label, or None to not label the y-axis
            ylabel_kwargs: a dictionary of keyword arguments to be passed to
                           ay.set_ylabel() to modify the y-axis label
            xlim: a length-2 tuple of the left and right x-axis limits, or None
                  to set the limits based on the data
            ylim: a length-2 tuple of the bottom and top y-axis limits, or None
                  to set the limits based on the data
            despine: whether to remove the top and right spines (borders of the
                     plot area) from the plot
            savefig_kwargs: a dictionary of keyword arguments to be passed to
                            `plt.savefig()`, such as:
                            - `dpi`: defaults to 300 instead of Matplotlib's
                              default of 150
                            - `bbox_inches`: the bounding box of the portion of
                              the figure to save; defaults to `'tight'` (crop
                              out any blank borders) instead of Matplotlib's
                              default of None (save the entire figure)
                            - `pad_inches`: the number of inches of padding to
                              add on each of the four sides of the figure when
                              saving. Defaults to `'layout'` (use the padding
                              from the constrained layout engine, when `ax` is
                              not None) instead of Matplotlib's default of 0.1.
                            - `transparent`: whether to save with a transparent
                              background; defaults to True if saving to a PDF
                              (i.e. when `filename` ends with `'.pdf'`) and
                              False otherwise, instead of Matplotlib's default
                              of always being False.
                            Can only be specified when `filename` is not None.
        """
        import matplotlib.pyplot as plt
        # Get `cells_to_plot_column`, if not None
        if cells_to_plot_column is not None:
            cells_to_plot_column = self._get_column(
                'obs', cells_to_plot_column, 'cells_to_plot_column',
                pl.Boolean,
                custom_error='cells_to_plot_column {} is not a column of obs; '
                             'set cells_to_plot_column=None to include all '
                             'cells')  # TODO fix
        # If `color_column` is not None, check that it either discrete
        # (Categorical, Enum, or String) or quantitative (integer or
        # floating-point). If discrete, require at least two distinct values.
        if color_column is not None:
            color_column = self._get_column(
                'obs', color_column, 'color_column',
                (pl.Categorical, pl.Enum, pl.String, 'integer',
                 'floating-point'), allow_null=True,
                QC_column=cells_to_plot_column)
            unique_color_column = color_column.unique().drop_nulls()
            dtype = color_column.dtype
            discrete = dtype in (pl.Categorical, pl.Enum, pl.String)
            if discrete and len(unique_color_column) == 1:
                error_message = (
                    f'color_column {color_column!r} must have at least two '
                    f'distinct values when its data type is '
                    f'{dtype.base_type()!r}')
                raise ValueError(error_message)
        # If `filename` is not None, check that it is a string or pathlib.Path
        # and that its base directory exists; if `filename` is None, make sure
        # `savefig_kwargs` is also None
        if filename is not None:
            check_type(filename, 'filename', (str, Path),
                       'a string or pathlib.Path')
            directory = os.path.dirname(filename)
            if directory and not os.path.isdir(directory):
                error_message = (
                    f'{filename} refers to a file in the directory '
                    f'{directory!r}, but this directory does not exist')
                raise NotADirectoryError(error_message)
            filename = str(filename)
        elif savefig_kwargs is not None:
            error_message = 'savefig_kwargs must be None when filename is None'
            raise ValueError(error_message)
        # Check that `embedding_key` is the name of a key in obsm
        check_type(embedding_key, 'embedding_key', str, 'a string')
        if embedding_key not in self._obsm:
            error_message = (
                f'embedding_key {embedding_key!r} is not a key of obsm; '
                f'did you forget to run embed() before plot_embedding()?')
            raise ValueError(error_message)
        # Check that the embedding `embedding_key` references is 2D.
        embedding = self._obsm[embedding_key]
        if embedding.shape[1] != 2:
            error_message = (
                f'the embedding at obsm[{embedding_key!r}] is '
                f'{embedding.shape[1]:,}-dimensional, but must be '
                f'2-dimensional to be plotted')
            raise ValueError(error_message)
        # If `cells_to_plot_column` is not None, subset to these cells
        if cells_to_plot_column is not None:
            embedding = embedding[cells_to_plot_column.to_numpy()]
            if color_column is not None:
                color_column = color_column.filter(cells_to_plot_column)
                unique_color_column = color_column.unique().drop_nulls()
        # Check that the embedding does not contain NaNs
        if np.isnan(embedding).any():
            error_message = \
                f'the embedding at obsm[{embedding_key!r}] contains NaNs; '
            if cells_to_plot_column is None and QC_column is not None:
                error_message += (
                    'did you forget to set QC_column to None in embed(), to '
                    'match the fact that you set cells_to_plot_column to '
                    'None in plot_embedding()?')
            else:
                error_message += (
                    'does your cells_to_plot_column contain cells that were '
                    'excluded by the QC_column used in embed()?')
            raise ValueError(error_message)
        # For each of the kwargs arguments, if the argument is not None, check
        # that it is a dictionary and that all its keys are strings.
        for kwargs, kwargs_name in (
                (palette_kwargs, 'palette_kwargs'),
                (scatter_kwargs, 'scatter_kwargs'),
                (label_kwargs, 'label_kwargs'),
                (legend_kwargs, 'legend_kwargs'),
                (colorbar_kwargs, 'colorbar_kwargs'),
                (title_kwargs, 'title_kwargs'),
                (xlabel_kwargs, 'xlabel_kwargs'),
                (ylabel_kwargs, 'ylabel_kwargs'),
                (savefig_kwargs, 'savefig_kwargs')):
            if kwargs is not None:
                check_type(kwargs, kwargs_name, dict, 'a dictionary')
                for key in kwargs:
                    if not isinstance(key, str):
                        error_message = (
                            f'all keys of {kwargs_name} must be strings, but '
                            f'it contains a key of type '
                            f'{type(key).__name__!r}')
                        raise TypeError(error_message)
        # If point_size is None, default to 30,000 / num_cells; otherwise,
        # require it to be a positive number or the name of a numeric column of
        # obs with all-positive numbers
        num_cells = \
            len(self) if cells_to_plot_column is None else len(embedding)
        if point_size is None:
            # noinspection PyUnboundLocalVariable
            point_size = 30_000 / num_cells
        else:
            check_type(point_size, 'point_size', (int, float, str),
                       'a positive number or string')
            if isinstance(point_size, (int, float)):
                check_bounds(point_size, 'point_size', 0, left_open=True)
            else:
                if point_size not in self._obs:
                    error_message = \
                        f'point_size {point_size!r} is not a column of obs'
                    raise ValueError(error_message)
                point_size = self._obs[point_size]
                if not (point_size.dtype.is_integer() or
                        point_size.dtype.is_float()):
                    error_message = (
                        f'the point_size column, obs[{point_size!r}], must '
                        f'have an integer or floating-point data type, but '
                        f'has data type {point_size.dtype.base_type()!r}')
                    raise TypeError(error_message)
                if point_size.min() <= 0:
                    error_message = (
                        f'the point_size column, obs[{point_size!r}], does '
                        f'not have all-positive elements')
                    raise ValueError(error_message)
        # If `sort_by_frequency=True`, ensure `palette` is None and
        # `color_column` is discrete
        check_type(sort_by_frequency, 'sort_by_frequency', bool, 'Boolean')
        if sort_by_frequency:
            if palette is not None:
                error_message = (
                    f'sort_by_frequency must be False when palette is '
                    f'specified')
                raise ValueError(error_message)
            # noinspection PyUnboundLocalVariable
            if color_column is None or not discrete:
                error_message = (
                    f'sort_by_frequency must be False when color_column is '
                    f'{"None" if color_column is None else "continuous"}')
                raise ValueError(error_message)
        # If `palette` is not None, check that it is a string, Colormap
        # object, or dictionary where all keys are in `color_column` and all
        # values are valid Matplotlib colors (and normalize these to hex
        # codes). If None and `color_column` is discrete, assign colors via
        # `generate_palette()`, in natural sort order, or decreasing order of
        # frequency if `sort_by_frequency=True`. Also make sure `palette` and
        # `palette_kwargs` are None when `color_column` is None.
        if palette is not None:
            if color_column is None:
                error_message = \
                    'palette must be None when color_column is None'
                raise ValueError(error_message)
            if palette_kwargs is not None:
                error_message = \
                    'palette_kwargs must be None when palette is specified'
                raise ValueError(error_message)
            check_type(palette, 'palette',
                       (str, plt.matplotlib.colors.Colormap, dict),
                       'a string, matplotlib Colormap object, or dictionary')
            if isinstance(palette, str):
                palette = plt.colormaps[palette]
            elif isinstance(palette, dict):
                # noinspection PyUnboundLocalVariable
                if not discrete:
                    error_message = (
                        'palette cannot be a dictionary when color_column is '
                        'continuous')
                    raise ValueError(error_message)
                for key, value in palette.items():
                    if not isinstance(key, str):
                        error_message = (
                            f'all keys of palette must be strings, but it '
                            f'contains a key of type {type(key).__name__!r}')
                        raise TypeError(error_message)
                    # noinspection PyUnboundLocalVariable
                    if key not in unique_color_column:
                        error_message = (
                            f'palette is a dictionary containing the key '
                            f'{key!r}, which is not one of the values in '
                            f'obs[{color_column!r}]')
                        raise ValueError(error_message)
                    if not plt.matplotlib.colors.is_color_like(value):
                        error_message = \
                            f'palette[{key!r}] is not a valid Matplotlib color'
                        raise ValueError(error_message)
                    palette[key] = plt.matplotlib.colors.to_hex(value)
        else:
            if color_column is not None and discrete:
                if palette_kwargs is None:
                    palette_kwargs = {}
                # noinspection PyUnboundLocalVariable
                color_order = color_column\
                    .value_counts(sort=True)\
                    .to_series()\
                    .drop_nulls() if sort_by_frequency else \
                        sorted(unique_color_column,
                               key=lambda color_label: [
                                   int(text) if text.isdigit() else
                                   text.lower() for text in
                                   re.split('([0-9]+)', color_label)])
                palette = generate_palette(len(color_order), **palette_kwargs)
                palette = dict(zip(color_order, palette))
            elif palette_kwargs is not None:
                error_message = (
                    f'palette_kwargs must be None when color_column is '
                    f'{"None" if color_column is None else "continuous"}')
                raise ValueError(error_message)
        # Check that `default_color` is a valid Matplotlib color, and convert
        # it to hex
        if not plt.matplotlib.colors.is_color_like(default_color):
            error_message = 'default_color is not a valid Matplotlib color'
            raise ValueError(error_message)
        default_color = plt.matplotlib.colors.to_hex(default_color)
        # Override the defaults for certain keys of `scatter_kwargs`
        default_scatter_kwargs = dict(rasterized=True, linewidths=0)
        scatter_kwargs = default_scatter_kwargs | scatter_kwargs \
            if scatter_kwargs is not None else default_scatter_kwargs
        # Check that `scatter_kwargs` does not contain the `s`, `c`/`color`, or
        # `cmap` keys
        if 's' in scatter_kwargs:
            error_message = (
                "'s' cannot be specified as a key in scatter_kwargs; specify "
                "the point_size argument instead")
            raise ValueError(error_message)
        for key in 'c', 'color', 'cmap':
            if key in scatter_kwargs:
                error_message = (
                    f'{key!r} cannot be specified as a key in '
                    f'scatter_kwargs; specify the color_column, palette, '
                    f'palette_kwargs, and/or default_color arguments instead')
                raise ValueError(error_message)
        # If `label=True`, require `color_column` to be discrete.
        # If `label=False`, `label_kwargs` must be None.
        check_type(label, 'label', bool, 'Boolean')
        if label:
            if color_column is None:
                error_message = 'color_column cannot be None when label=True'
                raise ValueError(error_message)
            if not discrete:
                error_message = \
                    'color_column cannot be continuous when label=True'
                raise ValueError(error_message)
        elif label_kwargs is not None:
            error_message = 'label_kwargs must be None when label=False'
            raise ValueError(error_message)
        # Only add a legend if `legend=True` and `color_column` is discrete.
        # If not adding a legend, `legend_kwargs` must be None.
        check_type(legend, 'legend', bool, 'Boolean')
        add_legend = legend and color_column is not None and discrete
        if not add_legend and legend_kwargs is not None:
            error_message = (
                f'legend_kwargs must be None when color_column is '
                f'{"None" if color_column is None else "continuous"}')
            raise ValueError(error_message)
        # Only add a colorbar if `colorbar=True` and `color_column` is
        # continuous. If not adding a colorbar, `colorbar_kwargs` must be None.
        check_type(colorbar, 'colorbar', bool, 'Boolean')
        add_colorbar = colorbar and color_column is not None and not discrete
        if not add_colorbar and colorbar_kwargs is not None:
            error_message = (
                f'colorbar_kwargs must be None when color_column is '
                f'{"None" if color_column is None else "discrete"}')
            raise ValueError(error_message)
        # `title` must be a string or None; if None, `title_kwargs` must be as
        # well; ditto for `xlabel` and `ylabel`
        for arg, arg_name, arg_kwargs in (
                (title, 'title', title_kwargs),
                (xlabel, 'xlabel', xlabel_kwargs),
                (ylabel, 'ylabel', ylabel_kwargs)):
            if arg is not None:
                check_type(arg, arg_name, str, 'a string')
            elif arg_kwargs is not None:
                error_message = \
                    f'{arg_name}_kwargs must be None when {arg_name} is None'
                raise ValueError(error_message)
        # `xlim` and `ylim` must be length-2 tuples or None, with the first
        # element less than the second
        for arg, arg_name in (xlim, 'xlim'), (ylim, 'ylim'):
            if arg is not None:
                check_type(arg, arg_name, tuple, 'a length-2 tuple')
                if len(arg) != 2:
                    error_message = (
                        f'{arg_name} must be a length-2 tuple, but has length '
                        f'{len(arg):,}')
                    raise ValueError(error_message)
                if arg[0] >= arg[1]:
                    error_message = \
                        f'{arg_name}[0] must be less than {arg_name}[1]'
                    raise ValueError(error_message)
        # If `color_column` is None, plot all cells in `default_color`. If
        # `palette` is a dictionary, generate an explicit list of colors to
        # plot each cell in. If `palette` is a Colormap, just pass it as the
        # cmap` argument. If `palette` is missing and `color_column` is
        # continuous, set it to `plt.rcParams['image.cmap']` ('viridis' by
        # default)
        if color_column is None:
            c = default_color
            cmap = None
        elif isinstance(palette, dict):
            # Note: `replace_strict(..., default=default_color)` fills both
            # missing values and values missing from `palette` with
            # `default_color`
            c = color_column\
                .replace_strict(palette, default=default_color,
                                return_dtype=pl.String)\
                .to_numpy()
            cmap = None
        else:
            # Need to `copy()` because `set_bad()` is in-place
            c = color_column.to_numpy()
            if palette is not None:
                cmap = palette.copy()
                # noinspection PyUnresolvedReferences
                cmap.set_bad(default_color)
            else:  # `color_column` is continuous
                cmap = plt.rcParams['image.cmap']
        # If `ax` is None, create a new figure with `constrained_layout=True`;
        # otherwise, check that it is a Matplotlib axis
        make_new_figure = ax is None
        try:
            if make_new_figure:
                plt.figure(constrained_layout=True)
                ax = plt.gca()
            else:
                check_type(ax, 'ax', plt.Axes, 'a Matplotlib axis')
            # Make a scatter plot of the embedding with equal x-y aspect ratios
            scatter = ax.scatter(embedding[:, 0], embedding[:, 1],
                                 s=point_size, c=c, cmap=cmap,
                                 **scatter_kwargs)
            ax.set_aspect('equal')
            # Add the title, axis labels and axis limits
            if title is not None:
                if title_kwargs is None:
                    ax.set_title(title)
                else:
                    ax.set_title(title, **title_kwargs)
            if xlabel is not None:
                if xlabel_kwargs is None:
                    ax.set_xlabel(xlabel)
                else:
                    ax.set_xlabel(xlabel, **xlabel_kwargs)
            if ylabel is not None:
                if ylabel_kwargs is None:
                    ax.set_ylabel(ylabel)
                else:
                    ax.set_ylabel(ylabel, **ylabel_kwargs)
            if xlim is not None:
                ax.set_xlim(*xlim)
            if ylim is not None:
                ax.set_ylim(*ylim)
            # Add the legend; override the defaults for certain values of
            # `legend_kwargs`
            if add_legend:
                default_legend_kwargs = dict(
                    loc='center left', bbox_to_anchor=(1, 0.5), frameon=False,
                    ncols=len(unique_color_column) // 16 + 1)
                legend_kwargs = default_legend_kwargs | legend_kwargs \
                    if legend_kwargs is not None else default_legend_kwargs
                if isinstance(palette, dict):
                    for color_label, color in palette.items():
                        ax.scatter([], [], c=color, label=color_label,
                                   **scatter_kwargs)
                    plt.legend(**legend_kwargs)
                else:
                    plt.legend(*scatter.legend_elements(), **legend_kwargs)
            # Add the colorbar; override the defaults for certain keys of
            # `colorbar_kwargs`
            if add_colorbar:
                default_colorbar_kwargs = dict(shrink=0.5)
                colorbar_kwargs = default_colorbar_kwargs | colorbar_kwargs \
                    if colorbar_kwargs is not None else default_colorbar_kwargs
                plt.colorbar(scatter, ax=ax, **colorbar_kwargs)
            # Label cells; override the defaults for certain keys of
            # `label_kwargs`
            if label:
                from matplotlib.patheffects import withStroke
                if label_kwargs is None:
                    label_kwargs = {}
                # noinspection PyUnresolvedReferences
                label_kwargs |= dict(
                    horizontalalignment=label_kwargs.pop(
                        'horizontalalignment',
                        label_kwargs.pop('ha', 'center')),
                    verticalalignment=label_kwargs.pop(
                        'verticalalignment',
                        label_kwargs.pop('va', 'center')),
                    path_effects=[withStroke(linewidth=3, foreground='white',
                                             alpha=0.75)])
                for color_label in unique_color_column:
                    ax.text(*np.median(embedding[color_column == color_label],
                                       axis=0), color_label, **label_kwargs)
            # Despine, if specified
            if despine:
                spines = ax.spines
                spines['top'].set_visible(False)
                spines['right'].set_visible(False)
            # Save; override the defaults for certain keys of `savefig_kwargs`
            if filename is not None:
                default_savefig_kwargs = \
                    dict(dpi=300, bbox_inches='tight', pad_inches='layout',
                         transparent=filename is not None and
                                     filename.endswith('.pdf'))
                savefig_kwargs = default_savefig_kwargs | savefig_kwargs \
                    if savefig_kwargs is not None else default_savefig_kwargs
                plt.savefig(filename, **savefig_kwargs)
                if make_new_figure:
                    plt.close()
        except:
            # If we made a new figure, make sure to close it if there's an
            # exception (but not if there was no error and `filename` is None,
            # in case the user wants to modify it further before saving)
            if make_new_figure:
                plt.close()
            raise
    
    def aggregate_var(self, ID_column, path_to_gtf):
        
        agg_matrix = pl.concat([pl.DataFrame(self.X.toarray().T), self.var[ID_column].to_frame()], how = "horizontal")\
            .group_by(ID_column)\
            .agg(pl.sum("*"))
        
        gtf = pl.read_csv(
            path_to_gtf, 
            separator = "\t", comment_prefix = "##",
            new_columns=[
                "chromosome",
                "source",
                "feature",
                "start",
                "end",
                "score",
                "strand",
                "frame",
                "attribute"]).filter(pl.col("feature") == "gene").with_columns(
                associated_gene = pl.col("attribute").str.extract(r'gene_name "([^;]*)";'),
                gene_type = pl.col("attribute").str.extract(r'gene_type "([^;]*)";'),
            ).drop("attribute")
        
        agg_matrix = agg_matrix\
            .cast({"associated_gene": pl.String})\
            .join(gtf["associated_gene", "gene_type"], how = "left", on = "associated_gene")

        gene_count_matrix = SingleCell(
            X = csr_array(agg_matrix.select(pl.selectors.starts_with("column_")).to_numpy().T),
            obs = self.obs,
            var = agg_matrix.drop(pl.selectors.starts_with("column_"))
        )          
        return gene_count_matrix

class Pseudobulk:
    """
    A pseudobulked single-cell dataset resulting from calling `pseudobulk()`
    on a SingleCell dataset.
    
    Has slots for:
    - X: a dict of NumPy arrays of counts per cell and gene for each cell type
    - obs: a dict of polars DataFrames of sample metadata for each cell type
    - var: a dict of polars DataFrames of gene metadata for each cell type
    as well as `obs_names` and `var_names`, aliases for a dict of obs[:, 0] and
    var[:, 0] for each cell type, and `cell_types`, a tuple of cell types.
    
    Supports iteration:
    - `for cell_type in pseudobulk:` yields the cell type names, as does
      `for cell_type in pseudobulk.keys():`
    - `for X, obs, var in pseudobulk.values():` yields the X, obs and var for
       each cell type
    - `for cell_type, (X, obs, var) in pseudobulk.items():` yields both the
      name and the X, obs and var for each cell type
    - `for X in pseudobulk.iter_X():` yields just the X for each cell type
    - `for X in pseudobulk.iter_obs():` yields just the obs for each cell type
    - `for X in pseudobulk.iter_var():` yields just the var for each cell type
    """
    def __init__(self,
                 X: dict[str, np.ndarray[2, np.integer | np.floating]] |
                    str | Path,
                 obs: dict[str, pl.DataFrame] = None,
                 var: dict[str, pl.DataFrame] = None) -> None:
        """
        Load a saved Pseudobulk dataset, or create one from an in-memory count
        matrix + metadata for each cell type. The latter functionality is
        mainly for internal use; most users will create new pseudobulk datasets
        by calling `pseudobulk()` on a SingleCell dataset.
        
        Args:
            X: a {cell type: NumPy array} dict of counts or log CPMs, or a
               directory to load a saved Pseudobulk dataset from (see save())
            obs: a {cell type: polars DataFrame} dict of metadata per sample.
                 The first column must be String, Enum or Categorical.
            var: a {cell type: polars DataFrame} dict of metadata per gene.
                 The first column must be String, Enum or Categorical.
        """
        if isinstance(X, dict):
            if obs is None:
                error_message = (
                    'obs is None, but since X is a dictionary, obs must also '
                    'be a dictionary')
                raise TypeError(error_message)
            if var is None:
                error_message = (
                    'var is None, but since X is a dictionary, var must also '
                    'be a dictionary')
                raise TypeError(error_message)
            if not X:
                error_message = 'X is an empty dictionary'
                raise ValueError(error_message)
            if X.keys() != obs.keys():
                error_message = \
                    'X and obs must have the same keys (cell types)'
                raise ValueError(error_message)
            if X.keys() != var.keys():
                error_message = \
                    'X and var must have the same keys (cell types)'
                raise ValueError(error_message)
            for cell_type in X:
                if not isinstance(cell_type, str):
                    error_message = (
                        f'all keys of X (cell types) must be strings, but X '
                        f'contains a key of type {type(cell_type).__name__!r}')
                    raise TypeError(error_message)
                check_type(X[cell_type], f'X[{cell_type!r}]', np.ndarray,
                           'a NumPy array')
                if X[cell_type].ndim != 2:
                    error_message = (
                        f'X[{cell_type!r}] is a {X[cell_type].ndim:,}-'
                        f'dimensional NumPy array, but must be 2-dimensional')
                    raise ValueError(error_message)
                check_type(obs[cell_type], f'obs[{cell_type!r}]', pl.DataFrame,
                           'a polars DataFrame')
                check_type(var[cell_type], f'var[{cell_type!r}]', pl.DataFrame,
                           'a polars DataFrame')
            self._X = X
            self._obs = obs
            self._var = var
        elif isinstance(X, (str, Path)):
            X = str(X)
            if not os.path.exists(X):
                error_message = f'Pseudobulk directory {X!r} does not exist'
                raise FileNotFoundError(error_message)
            cell_types = [line.rstrip('\n') for line in
                          open(f'{X}/cell_types.txt')]
            self._X = {cell_type: np.load(
                os.path.join(X, f'{cell_type.replace("/", "-")}.X.npy'))
                for cell_type in cell_types}
            self._obs = {cell_type: pl.read_parquet(
                os.path.join(X, f'{cell_type.replace("/", "-")}.obs.parquet'))
                for cell_type in cell_types}
            self._var = {cell_type: pl.read_parquet(
                os.path.join(X, f'{cell_type.replace("/", "-")}.var.parquet'))
                for cell_type in cell_types}
        else:
            error_message = (
                f'X must be a dictionary of NumPy arrays or a directory '
                f'containing a saved Pseudobulk dataset, but has type '
                f'{type(X).__name__!r}')
            raise ValueError(error_message)
        for cell_type in self:
            if len(self._obs[cell_type]) == 0:
                error_message = \
                    f'len(obs[{cell_type!r}]) is 0: no samples remain'
                raise ValueError(error_message)
            if len(self._var[cell_type]) == 0:
                error_message = \
                    f'len(var[{cell_type!r}]) is 0: no genes remain'
                raise ValueError(error_message)
            if len(self._obs[cell_type]) != len(self._X[cell_type]):
                error_message = (
                    f'len(obs[{cell_type!r}]) is '
                    f'{len(self._obs[cell_type]):,}, but '
                    f'len(X[{cell_type!r}]) is {len(X[cell_type]):,}')
                raise ValueError(error_message)
            if len(self._var[cell_type]) != self._X[cell_type].shape[1]:
                error_message = (
                    f'len(var[{cell_type!r}]) is '
                    f'{len(self._var[cell_type]):,}, but '
                    f'X[{cell_type!r}].shape[1] is '
                    f'{self._X[cell_type].shape[1]:,}')
                raise ValueError(error_message)
            if self._obs[cell_type][:, 0].dtype not in \
                    (pl.String, pl.Categorical, pl.Enum):
                error_message = (
                    f'the first column of obs[{cell_type!r}] '
                    f'({self._obs[cell_type].columns[0]!r}) must be String, '
                    f'Categorical or Enum, but has data type '
                    f'{self._obs[cell_type][:, 0].dtype.base_type()!r}')
                raise ValueError(error_message)
            if self._var[cell_type][:, 0].dtype not in \
                    (pl.String, pl.Categorical, pl.Enum):
                error_message = (
                    f'the first column of var[{cell_type!r}] '
                    f'({self._var[cell_type].columns[0]!r}) must be String, '
                    f'Categorical or Enum, but has data type '
                    f'{self._var[cell_type][:, 0].dtype.base_type()!r}')
                raise ValueError(error_message)
    
    @staticmethod
    def _setter_check(new: dict[str, np.ndarray[2, np.integer | np.floating]
                                     | pl.DataFrame],
                      old: dict[str, np.ndarray[2, np.integer | np.floating]
                                     | pl.DataFrame],
                      name: str) -> None:
        """
        When setting X, obs or var, raise an error if the new value is not a
        dictionary, the new keys (cell types) differ from the old ones, or
        the new values differ in length (or shape, in the case of X) from the
        old ones. For obs and var, also check that the first column is String,
        Categorical or Enum.
        
        Args:
            new: the new X, obs or var
            old: the old X, obs or var
            name: the name of the field: 'X', 'obs' or 'var'
        """
        if not isinstance(new, dict):
            error_message = (
                f'new {name} must be a dictionary, but has type '
                f'{type(new).__name__!r}')
            raise TypeError(error_message)
        if new.keys() != old.keys():
            error_message = (
                f'new {name} has different cell types (keys) from the old '
                f'{name}')
            raise ValueError(error_message)
        if name == 'X':
            for cell_type in new:
                check_type(new[cell_type], f'X[{cell_type!r}]', np.ndarray,
                           'a NumPy array')
                new_shape = new[cell_type].shape
                old_shape = old[cell_type].shape
                if new_shape != old_shape:
                    error_message = (
                        f'new X is {new_shape.shape[0]:,} × '
                        f'{new_shape.shape[1]:,}, but old X is '
                        f'{old_shape.shape[0]:,} × {old_shape.shape[1]:,}')
                    raise ValueError(error_message)
        else:
            for cell_type in new:
                check_type(new[cell_type], f'{name}[{cell_type!r}]',
                           pl.DataFrame, 'a polars DataFrame')
                if new[cell_type][:, 0].dtype not in (
                        pl.String, pl.Categorical, pl.Enum):
                    error_message = (
                        f'the first column of {name}[{cell_type!r}] '
                        f'({new[cell_type].columns[0]!r}) must be String, '
                        f'Categorical or Enum, but the first column of the '
                        f'new {name} has data type '
                        f'{new[cell_type][:, 0].dtype.base_type()!r}')
                    raise ValueError(error_message)
                if len(new) != len(old):
                    error_message = (
                        f'new {name} has length {len(new):,}, but old {name} '
                        f'has length {len(old):,}')
                    raise ValueError(error_message)
    
    @property
    def X(self) -> dict[str, np.ndarray[2, np.integer | np.floating]]:
        return self._X
    
    @X.setter
    def X(self, X: dict[str, np.ndarray[2, np.integer | np.floating]]) -> \
            None:
        self._setter_check(X, self._X, 'X')
        self._X = X
    
    @property
    def obs(self) -> dict[str, pl.DataFrame]:
        return self._obs
    
    @obs.setter
    def obs(self, obs: dict[str, pl.DataFrame]) -> None:
        self._setter_check(obs, self._obs, 'obs')
        self._obs = obs

    @property
    def var(self) -> dict[str, pl.DataFrame]:
        return self._var
    
    @var.setter
    def var(self, var: dict[str, pl.DataFrame]) -> None:
        self._setter_check(var, self._var, 'var')
        self._var = var

    @property
    def obs_names(self) -> dict[str, pl.Series]:
        return {cell_type: obs[:, 0] for cell_type, obs in self._obs.items()}
    
    @property
    def var_names(self) -> dict[str, pl.Series]:
        return {cell_type: var[:, 0] for cell_type, var in self._var.items()}
    
    def set_obs_names(self, column: str) -> Pseudobulk:
        """
        Sets a column as the new first column of obs, i.e. the obs_names.
        
        Args:
            column: the column name in obs; must have String, Categorical, or
                    Enum data type

        Returns:
            A new Pseudobulk dataset with `column` as the first column of each
            cell type's obs. If `column` is already the first column for every
            cell type, return this dataset unchanged.
        """
        check_type(column, 'column', str, 'a string')
        if all(column == cell_type_obs.columns[0]
               for cell_type_obs in self._obs.values()):
            return self
        obs = {}
        for cell_type, cell_type_obs in self._obs.items():
            if column not in cell_type_obs:
                error_message = \
                    f'{column!r} is not a column of obs[{cell_type!r}]'
                raise ValueError(error_message)
            check_dtype(cell_type_obs, f'obs[{column!r}]',
                        (pl.String, pl.Categorical, pl.Enum))
            obs[cell_type] = cell_type_obs.select(column, pl.exclude(column))
        return Pseudobulk(X=self._X, obs=obs, var=self._var)
    
    def set_var_names(self, column: str) -> Pseudobulk:
        """
        Sets a column as the new first column of var, i.e. the var_names.
        
        Args:
            column: the column name in var; must have String, Categorical, or
                    Enum data type

        Returns:
            A new Pseudobulk dataset with `column` as the first column of each
            cell type's var. If `column` is already the first column for every
            cell type, return this dataset unchanged.
        """
        check_type(column, 'column', str, 'a string')
        if all(column == cell_type_var.columns[0]
               for cell_type_var in self._var.values()):
            return self
        var = {}
        for cell_type, cell_type_var in self._var.items():
            if column not in cell_type_var:
                error_message = \
                    f'{column!r} is not a column of var[{cell_type!r}]'
                raise ValueError(error_message)
            check_dtype(cell_type_var, f'var[{column!r}]',
                        (pl.String, pl.Categorical, pl.Enum))
            var[cell_type] = cell_type_var.select(column, pl.exclude(column))
        return Pseudobulk(X=self._X, obs=self._obs, var=var)

    def keys(self) -> KeysView[str]:
        """
        Get a KeysView (like you would get from `dict.keys()`) of this
        Pseudobulk dataset's cell types. `for cell_type in pb.keys()` is
        equivalent to `for cell_type in pb`.
        
        Returns:
            A KeysView of the cell types.
        """
        return self._X.keys()
    
    def values(self) -> ValuesView[tuple[np.ndarray[2, np.integer |
                                                       np.floating],
                                       pl.DataFrame, pl.DataFrame]]:
        """
        Get a ValuesView (like you would get from `dict.values()`) of
        `(X, obs, var)` tuples for each cell type in this Pseudobulk dataset.
        
        Returns:
            A ValuesView of `(X, obs, var)` tuples for each cell type.
        """
        return {cell_type: (self._X[cell_type], self._obs[cell_type],
                            self._var[cell_type])
                for cell_type in self._X}.values()
    
    def items(self) -> ItemsView[str, tuple[np.ndarray[2, np.integer |
                                                          np.floating],
                                            pl.DataFrame, pl.DataFrame]]:
        """
        Get an ItemsView (like you would get from `dict.items()`) of
        `(cell_type, (X, obs, var))` tuples for each cell type in this
        Pseudobulk dataset.
        
        Yields:
            An ItemsView of `(cell_type, (X, obs, var))` tuples for each cell
            type.
        """
        return {cell_type: (self._X[cell_type], self._obs[cell_type],
                            self._var[cell_type])
                for cell_type in self._X}.items()
    
    def iter_X(self) -> Iterable[np.ndarray[2, np.integer | np.floating]]:
        """
        Iterate over each cell type's X.
        
        Yields:
            X for each cell type.
        """
        for cell_type in self:
            yield self._X[cell_type]
    
    def iter_obs(self) -> Iterable[pl.DataFrame]:
        """
        Iterate over each cell type's obs.
        
        Yields:
            obs for each cell type.
        """
        for cell_type in self:
            yield self._obs[cell_type]
    
    def iter_var(self) -> Iterable[pl.DataFrame]:
        """
        Iterate over each cell type's var.
        
        Yields:
            var for each cell type.
        """
        for cell_type in self:
            yield self._var[cell_type]
    
    def map_X(self, function: Callable[[np.ndarray[2, np.integer |
                                                      np.floating], ...],
                                        np.ndarray[2, np.integer |
                                                      np.floating]],
              *args: Any, **kwargs: Any) -> Pseudobulk:
        """
        Apply a function to each cell type's X.
        
        Args:
            function: the function to apply
            *args: the positional arguments to the function
            **kwargs: the keyword arguments to the function

        Returns:
            A new Pseudobulk dataset where the function has been applied to
            each cell type's X.
        """
        return Pseudobulk(X={cell_type: function(self._X[cell_type], *args,
                                                 **kwargs)
                             for cell_type in self},
                          obs=self._obs, var=self._var)
    
    def map_obs(self, function: Callable[[pl.DataFrame, ...], pl.DataFrame],
                *args: Any, **kwargs: Any) -> Pseudobulk:
        """
        Apply a function to each cell type's obs.
        
        Args:
            function: the function to apply
            *args: the positional arguments to the function
            **kwargs: the keyword arguments to the function

        Returns:
            A new Pseudobulk dataset where the function has been applied to
            each cell type's obs.
        """
        return Pseudobulk(X=self._X,
                          obs={cell_type: function(self._obs[cell_type], *args,
                                                   **kwargs)
                               for cell_type in self},
                          var=self._var)
    
    def map_var(self, function: Callable[[pl.DataFrame, ...], pl.DataFrame],
                *args: Any, **kwargs: Any) -> Pseudobulk:
        """
        Apply a function to each cell type's var.
        
        Args:
            function: the function to apply
            *args: the positional arguments to the function
            **kwargs: the keyword arguments to the function

        Returns:
            A new Pseudobulk dataset where the function has been applied to
            each cell type's var.
        """
        return Pseudobulk(X=self._X, obs=self._obs,
                          var={cell_type: function(self._var[cell_type], *args,
                                                   **kwargs)
                               for cell_type in self})
        
    def __eq__(self, other: Pseudobulk) -> bool:
        """
        Test for equality with another Pseudobulk dataset.
        
        Args:
            other: the other Pseudobulk dataset to test for equality with

        Returns:
            Whether the two Pseudobulk datasets are identical.
        """
        if not isinstance(other, Pseudobulk):
            error_message = (
                f'the left-hand operand of `==` is a Pseudobulk dataset, but '
                f'the right-hand operand has type {type(other).__name__!r}')
            raise TypeError(error_message)
        # noinspection PyUnresolvedReferences
        return tuple(self.keys()) == tuple(other.keys()) and \
            all(obs.equals(other_obs) for obs, other_obs in
                zip(self._obs.values(), other._obs.values())) and \
            all(var.equals(other_var) for var, other_var in
                zip(self._var.values(), other._var.values())) and \
            all(np.array_equal(X, other_X, equal_nan=True) for X, other_X in
                zip(self._X.values(), other._X.values()))
    
    def __or__(self, other: Pseudobulk) -> Pseudobulk:
        """
        Combine the cell types of this Pseudobulk dataset with another. The
        two datasets must have non-overlapping cell types.
        
        Args:
            other: the other Pseudobulk dataset to combine with this one

        Returns:
            A Pseudobulk dataset with each of the cell types in the first
            Pseudobulk dataset, followed by each of the cell types in the
            second.
        """
        if not isinstance(other, Pseudobulk):
            error_message = (
                f'the left-hand operand of `|` is a Pseudobulk dataset, but '
                f'the right-hand operand has type {type(other).__name__!r}')
            raise TypeError(error_message)
        if self.keys() & other.keys():
            error_message = (
                'the left- and right-hand operands of `|` are Pseudobulk '
                'datasets that share some cell types')
            raise ValueError(error_message)
        return Pseudobulk(X=self._X | other._X, obs=self._obs | other._obs,
                          var=self._var | other._var)
    
    def __contains__(self, cell_type: str) -> bool:
        """
        Check if this Pseudobulk dataset contains the specified cell type.
        
        Args:
            cell_type: the cell type

        Returns:
            Whether the cell type is present in the Pseudobulk dataset.
        """
        check_type(cell_type, 'cell_type', str, 'a string')
        return cell_type in self._X
    
    @staticmethod
    def _getitem_error(item: Indexer) -> None:
        """
        Raise an error if the indexer is invalid.
        
        Args:
            item: the indexer
        """
        types = tuple(type(elem).__name__ for elem in to_tuple(item))
        if len(types) == 1:
            types = types[0]
        error_message = (
            f'Pseudobulk indices must be a cell-type string, a length-1 tuple '
            f'of (cell_type,), a length-2 tuple of (cell_type, samples), or a '
            f'length-3 tuple of (cell_type, samples, genes). Samples and '
            f'genes must each be a string or integer; a slice of strings or '
            f'integers; or a list, NumPy array, or polars Series of strings, '
            f'integers, or Booleans. You indexed with: {types}.')
        raise ValueError(error_message)
    
    @staticmethod
    def _getitem_by_string(df: pl.DataFrame, string: str) -> int:
        """
        Get the index where df[:, 0] == string, raising an error if no rows or
        multiple rows match.
        
        Args:
            df: a DataFrame (obs or var)
            string: the string to find the index of in the first column of df

        Returns:
            The integer index of the string within the first column of df.
        """
        first_column = df.columns[0]
        try:
            return df\
                .select(pl.int_range(pl.len(), dtype=pl.Int32)
                        .alias('__Pseudobulk_getitem'), first_column)\
                .row(by_predicate=pl.col(first_column) == string)\
                [0]
        except pl.exceptions.NoRowsReturnedError:
            raise KeyError(string)
    
    @staticmethod
    def _getitem_process(item: Indexer, index: int, df: pl.DataFrame) -> \
            list[int] | slice | pl.Series:
        """
        Process an element of an item passed to __getitem__().
        
        Args:
            item: the item
            index: the index of the element to process
            df: the DataFrame (obs or var) to process the element with respect
                to

        Returns:
            A new indexer indicating the rows/columns to index.
        """
        subitem = item[index]
        if is_integer(subitem):
            return [subitem]
        elif isinstance(subitem, str):
            return [Pseudobulk._getitem_by_string(df, subitem)]
        elif isinstance(subitem, slice):
            start = subitem.start
            stop = subitem.stop
            step = subitem.step
            if isinstance(start, str):
                start = Pseudobulk._getitem_by_string(df, start)
            elif start is not None and not is_integer(start):
                Pseudobulk._getitem_error(item)
            if isinstance(stop, str):
                stop = Pseudobulk._getitem_by_string(df, stop)
            elif stop is not None and not is_integer(stop):
                Pseudobulk._getitem_error(item)
            if step is not None and not is_integer(step):
                Pseudobulk._getitem_error(item)
            return slice(start, stop, step)
        elif isinstance(subitem, (list, np.ndarray, pl.Series)):
            if not isinstance(subitem, pl.Series):
                subitem = pl.Series(subitem)
            if subitem.is_null().any():
                error_message = (
                    'indexer contains null entries; this may happen when '
                    'indexing with a mix of strings and integers, due to a '
                    'bug in Polars: github.com/pola-rs/polars/issues/11156')
                raise ValueError(error_message)
            if subitem.dtype == pl.String or subitem.dtype == \
                    pl.Categorical or subitem.dtype == pl.Enum:
                indices = subitem\
                    .to_frame(df.columns[0])\
                    .join(df.with_columns(_Pseudobulk_index=pl.int_range(
                              pl.len(), dtype=pl.Int32)),
                          on=df.columns[0], how='left')\
                    ['_Pseudobulk_index']
                if indices.null_count():
                    error_message = subitem.filter(indices.is_null())[0]
                    raise KeyError(error_message)
                return indices
            elif subitem.dtype.is_integer() or subitem.dtype == pl.Boolean:
                return subitem
            else:
                Pseudobulk._getitem_error(item)
        else:
            Pseudobulk._getitem_error(item)
    
    def __getitem__(self, item: Indexer | tuple[str, Indexer, Indexer]) -> \
            Pseudobulk:
        """
        Subset to specific cell type(s), sample(s), and/or gene(s).
        
        Index with a tuple of `(cell_types, samples, genes)`. If `samples` and
        `genes` are integers, arrays/lists/slices of integers, or arrays/lists
        of Booleans, the result will be a Pseudobulk dataset subset to
        `X[samples, genes]`, `obs[samples]`, and `var[genes]` for each of the
        cell types in `cell_types`. However, `samples` and/or `genes` can
        instead be strings (or arrays or slices of strings), in which case they
        refer to the first column of obs and/or var, respectively.
        
        Examples:
        - Subset to one cell type:
          pseudobulk['Astro']
        - Subset to multiple cell types:
          pseudobulk[['Astro', 'Micro']]
        - Subset to one cell type and sample, for all genes:
          pb['Astro', 'H19.30.002']
          pb['Astro', 2]
        - Subset to one gene, for all cell types and samples:
          pb[:, :, 'APOE']
          pb[:, :, 13196]
        - Subset to one cell type, sample and gene:
          pb['Astro', 'H18.30.002', 'APOE']
          pb['Astro', 2, 13196]
        - Subset to one cell type and a range of samples and genes:
          pb['Astro', 'H18.30.002':'H19.33.004', 'APOE':'TREM2']
          pb['Astro', 'H18.30.002':'H19.33.004', 13196:34268]
        - Subset to one a cell type and specific samples and genes:
          pb['Astro', ['H18.30.002', 'H19.33.004']]
          pb['Astro', :, pl.Series(['APOE', 'TREM2'])]
          pb['Astro', ('H18.30.002', 'H19.33.004'),
             np.array(['APOE', 'TREM2'])]
        
        Args:
            item: the item to index with

        Returns:
            A new Pseudobulk dataset subset to the specified cell types,
            samples, and/or genes.
        """
        if isinstance(item, tuple):
            if not 1 <= len(item) <= 3:
                self._getitem_error(item)
            cell_types = to_tuple(item[0])
        elif isinstance(item, list):
            cell_types = to_tuple(item)
        elif isinstance(item, str):
            cell_types = item,
        else:
            self._getitem_error(item)
        # noinspection PyUnboundLocalVariable
        for cell_type in cell_types:
            if cell_type not in self:
                if isinstance(cell_type, str):
                    error_message = (
                        f'tried to select {cell_type!r}, which is not a cell '
                        f'type in this Pseudobulk')
                    raise ValueError(error_message)
                else:
                    error_message = (
                        f'tried to select a non-existent cell type of type '
                        f'{type(cell_type).__name__!r}')
                    raise TypeError(error_message)
        if not isinstance(item, tuple) or len(item) == 1:
            return Pseudobulk(X={cell_type: self._X[cell_type]
                                 for cell_type in cell_types},
                              obs={cell_type: self._obs[cell_type]
                                   for cell_type in cell_types},
                              var={cell_type: self._var[cell_type]
                                   for cell_type in cell_types})
        X, obs, var = {}, {}, {}
        for cell_type in cell_types:
            rows = self._getitem_process(item, 1, self._obs[cell_type])
            if isinstance(rows, pl.Series):
                obs[cell_type] = self._obs[cell_type].filter(rows) \
                    if rows.dtype == pl.Boolean else self._obs[cell_type][rows]
                rows = rows.to_numpy()
            else:
                obs[cell_type] = self._obs[cell_type][rows]
            if len(item) == 2:
                X[cell_type] = self._X[cell_type][rows]
                var[cell_type] = self._var[cell_type]
            else:
                columns = self._getitem_process(item, 2, self._var[cell_type])
                if isinstance(columns, pl.Series):
                    var[cell_type] = self._var[cell_type].filter(columns) \
                        if columns.dtype == pl.Boolean \
                        else self._var[cell_type][columns]
                    columns = columns.to_numpy()
                else:
                    var[cell_type] = self._var[cell_type][columns]
                X[cell_type] = self._X[cell_type][rows, columns] \
                    if isinstance(rows, slice) or \
                       isinstance(columns, slice) else \
                    self._X[cell_type][np.ix_(rows, columns)]
        return Pseudobulk(X=X, obs=obs, var=var)
    
    def sample(self, cell_type: str, sample: str) -> np.ndarray[1, Any]:
        """
        Get the row of X[cell_type] corresponding to a single sample, based on
        the sample's name in obs_names.
        
        Args:
            cell_type: the cell type to retrieve the row of X from
            sample: the name of the sample in obs_names
        
        Returns:
            The corresponding row of X[cell_type], as a dense 1D NumPy array
            with zeros included.
        """
        row_index = Pseudobulk._getitem_by_string(self._obs[cell_type], sample)
        return self._X[cell_type][[row_index]].toarray().squeeze()
    
    def gene(self, cell_type: str, gene: str) -> np.ndarray[1, Any]:
        """
        Get the column of X[cell_type] corresponding to a single gene, based on
        the gene's name in var_names.
        
        Args:
            cell_type: the cell type to retrieve the row of X from
            gene: the name of the gene in var_names
        
        Returns:
            The corresponding column of X[cell_type], as a dense 1D NumPy array
            with zeros included.
        """
        column_index = \
            Pseudobulk._getitem_by_string(self._var[cell_type], gene)
        return self._X[cell_type][:, [column_index]].toarray().squeeze()
    
    def __iter__(self) -> Iterable[str]:
        """
        Iterate over the cell types of this Pseudobulk dataset.
        `for cell_type in pb` is equivalent to `for cell_type in pb.keys()`.
        
        Returns:
            An iterator over the cell types.
        """
        return iter(self._X)
    
    def __len__(self) -> dict[str, int]:
        """
        Get the number of samples in each cell type of this Pseudobulk dataset.
        
        Returns:
            A dictionary mapping each cell type to its number of samples.
        """
        return {cell_type: len(X_cell_type)
                for cell_type, X_cell_type in self._X.items()}
    
    def __repr__(self) -> str:
        """
        Get a string representation of this Pseudobulk dataset.
        
        Returns:
            A string summarizing the dataset.
        """
        min_num_samples = min(len(obs) for obs in self._obs.values())
        max_num_samples = max(len(obs) for obs in self._obs.values())
        min_num_genes = min(len(var) for var in self._var.values())
        max_num_genes = max(len(var) for var in self._var.values())
        samples_string = \
            f'{min_num_samples:,} {plural("sample", max_num_samples)}' \
            if min_num_samples == max_num_samples else \
            f'{min_num_samples:,}-{max_num_samples:,} samples'
        genes_string = \
            f'{min_num_genes:,} {plural("gene", max_num_genes)}' \
            if min_num_genes == max_num_genes else \
            f'{min_num_genes:,}-{max_num_genes:,} genes'
        return f'Pseudobulk dataset with {len(self._X):,} cell ' \
               f'{"types, each" if len(self._X) > 1 else "type,"} with ' \
               f'{samples_string} (obs) and {genes_string} (var)\n' + \
            fill(f'    Cell types: {", ".join(self._X)}',
                 width=os.get_terminal_size().columns,
                 subsequent_indent=' ' * 17)
    
    @property
    def shape(self) -> dict[str, tuple[int, int]]:
        """
        Get the shape of each cell type in this Pseudobulk dataset.
        
        Returns:
            A dictionary mapping each cell type to a length-2 tuple where the
            first element is the number of samples, and the second is the
            number of genes.
        """
        return {cell_type: X_cell_type.shape
                for cell_type, X_cell_type in self._X.items()}
    
    def save(self, directory: str | Path, overwrite: bool = False) -> None:
        """
        Saves a Pseudobulk dataset to `directory` (which must not exist unless
        `overwrite=True`, and will be created) with three files per cell type:
        the X at f'{cell_type}.X.npy', the obs at f'{cell_type}.obs.parquet',
        and the var at f'{cell_type}.var.parquet'. Also saves a text file,
        cell_types.txt, containing the cell types.
        
        Args:
            directory: the directory to save the Pseudobulk dataset to
            overwrite: if False, raises an error if the directory exists; if
                       True, overwrites files inside it as necessary
        """
        check_type(directory, 'directory', (str, Path),
                   'a string or pathlib.Path')
        directory = str(directory)
        if not overwrite and os.path.exists(directory):
            error_message = (
                f'directory {directory!r} already exists; set overwrite=True '
                f'to overwrite')
            raise FileExistsError(error_message)
        os.makedirs(directory, exist_ok=overwrite)
        with open(os.path.join(directory, 'cell_types.txt'), 'w') as f:
            print('\n'.join(self._X), file=f)
        for cell_type in self._X:
            escaped_cell_type = cell_type.replace('/', '-')
            np.save(os.path.join(directory, f'{escaped_cell_type}.X.npy'),
                    self._X[cell_type])
            self._obs[cell_type].write_parquet(
                os.path.join(directory, f'{escaped_cell_type}.obs.parquet'))
            self._var[cell_type].write_parquet(
                os.path.join(directory, f'{escaped_cell_type}.var.parquet'))
    
    def copy(self, deep: bool = False) -> Pseudobulk:
        """
        Make a deep (if deep=True) or shallow copy of this Pseudobulk dataset.
        
        Returns:
            A copy of the Pseudobulk dataset. Since polars DataFrames are
            immutable, obs[cell_type] and var[cell_type] will always point to
            the same underlying data as the original for all cell types. The
            only difference when deep=True is that X[cell_type] will point to a
            fresh copy of the data, rather than the same data. Watch out: when
            deep=False, any modifications to X[cell_type] will modify both
            copies!
        """
        check_type(deep, 'deep', bool, 'Boolean')
        return Pseudobulk(X={cell_type: cell_type_X.copy()
                             for cell_type, cell_type_X in self._X.items()}
                            if deep else self._X, obs=self._obs, var=self._var)

    def concat_obs(self,
                   datasets: Pseudobulk,
                   *more_datasets: Pseudobulk,
                   flexible: bool = False) -> Pseudobulk:
        """
        Concatenate the samples of multiple Pseudobulk datasets. All datasets
        must have the same cell types.
        
        By default, all datasets must have the same var. They must also have
        the same columns in obs, with the same data types.
        
        Conversely, if `flexible=True`, subset to genes present in all datasets
        (according to the first column of var, i.e. `var_names`) before
        concatenating. Subset to columns of var that are identical in all
        datasets after this subsetting. Also, subset to columns of obs that are
        present in all datasets, and have the same data types. All datasets'
        `obs_names` must have the same name and dtype, and similarly for
        `var_names`.
        
        The one exception to the obs "same data type" rule: if a column is Enum
        in some datasets and Categorical in others, or Enum in all datasets but
        with different categories in each dataset, that column will be retained
        as an Enum column (with the union of the categories) in the
        concatenated obs.
        
        Args:
            datasets: one or more Pseudobulk datasets to concatenate with this
                      one
            *more_datasets: additional Pseudobulk datasets to concatenate with
                            this one, specified as positional arguments
            flexible: whether to subset to genes and columns of obs and var
                      common to all datasets before concatenating, rather than
                      raising an error on any mismatches
        
        Returns:
            The concatenated Pseudobulk dataset.
        """
        # Check inputs
        if isinstance(datasets, Pseudobulk):
            datasets = datasets,
        datasets = (self,) + datasets + more_datasets
        if len(datasets) == 1:
            error_message = \
                'need at least one other Pseudobulk dataset to concatenate'
            raise ValueError(error_message)
        check_types(datasets[1:], 'datasets', Pseudobulk,
                    'Pseudobulk datasets')
        check_type(flexible, 'flexible', bool, 'Boolean')
        # Check that cell types match across all datasets
        if not all(self.keys() == dataset.keys() for dataset in datasets[1:]):
            error_message = \
                'not all Pseudobulk datasets have the same cell types'
            raise ValueError(error_message)
        # Perform either flexible or non-flexible concatenation
        X = {}
        obs = {}
        var = {}
        for cell_type in self:
            if flexible:
                # Check that `obs_names` and `var_names` have the same name and
                # data type for each cell type across all datasets
                obs_names_name = self._obs[cell_type][:, 0].name
                if not all(dataset._obs[cell_type][:, 0] == obs_names_name
                           for dataset in datasets[1:]):
                    error_message = (
                        f'[{cell_type!r}] not all Pseudobulk datasets have '
                        f'the same name for the first column of obs (the '
                        f'obs_names column)')
                    raise ValueError(error_message)
                var_names_name = self._var[cell_type][:, 0].name
                if not all(dataset._var[cell_type][:, 0].name == var_names_name
                           for dataset in datasets[1:]):
                    error_message = (
                        f'[{cell_type!r}] not all Pseudobulk datasets have '
                        f'the same name for the first column of var (the '
                        f'var_names column)')
                    raise ValueError(error_message)
                obs_names_dtype = self._obs[cell_type][:, 0].dtype
                if not all(dataset._obs[cell_type][:, 0].dtype ==
                           obs_names_dtype for dataset in datasets[1:]):
                    error_message = (
                        f'[{cell_type!r}] not all Pseudobulk datasets have '
                        f'the same data type for the first column of obs (the '
                        f'obs_names column)')
                    raise TypeError(error_message)
                var_names_dtype = self._var[cell_type][:, 0].dtype
                if not all(dataset._var[cell_type][:, 0].dtype ==
                           var_names_dtype for dataset in datasets[1:]):
                    error_message = (
                        f'[{cell_type!r}] not all Pseudobulk datasets have '
                        f'the same data type for the first column of var (the '
                        f'var_names column)')
                    raise TypeError(error_message)
                # Subset to genes in common across all datasets
                genes_in_common = self._var[cell_type][:, 0]\
                    .filter(self._var[cell_type][:, 0]
                            .is_in(pl.concat([dataset._var[cell_type][:, 0]
                                              for dataset in datasets[1:]])))
                if len(genes_in_common) == 0:
                    error_message = (
                        f'[{cell_type!r}] no genes are shared across all '
                        f'Pseudobulk datasets')
                    raise ValueError(error_message)
                cell_type_X = []
                cell_type_var = []
                for dataset in datasets:
                    gene_indices = dataset._getitem_process(
                        genes_in_common, 1, dataset._var[cell_type])
                    cell_type_X.append(
                        dataset._X[cell_type][:, gene_indices.to_numpy()])
                    cell_type_var.append(dataset._var[cell_type][gene_indices])
                # Subset to columns of var that are identical in all datasets
                # after this subsetting
                var_columns_in_common = [
                    column.name for column in cell_type_var[0][:, 1:]
                    if all(column.name in dataset_cell_type_var and
                           dataset_cell_type_var[column.name].equals(column)
                           for dataset_cell_type_var in cell_type_var[1:])]
                cell_type_var = cell_type_var[0]
                cell_type_var = cell_type_var.select(cell_type_var.columns[0],
                                                     var_columns_in_common)
                # Subset to columns of obs that are present in all datasets,
                # and have the same data types. Also include columns of obs
                # that are Enum in some datasets and Categorical in others, or
                # Enum in all datasets but with different categories in each
                # dataset; cast these to Categorical.
                obs_mismatched_categoricals = {
                    column for column, dtype in self._obs[cell_type][:, 1:]
                    .select(pl.col(pl.Categorical, pl.Enum)).schema.items()
                    if all(column in dataset._obs[cell_type] and
                           dataset._obs[cell_type][column].dtype in
                           (pl.Categorical, pl.Enum)
                           for dataset in datasets[1:]) and
                       not all(dataset._obs[cell_type][column].dtype == dtype
                               for dataset in datasets[1:])}
                obs_columns_in_common = [
                    column
                    for column, dtype in islice(
                        self._obs[cell_type].schema.items(), 1, None)
                    if column in obs_mismatched_categoricals or
                       all(column in dataset[cell_type]._obs and
                           dataset._obs[cell_type][column].dtype == dtype
                           for dataset in datasets[1:])]
                cast_dict = {column: pl.Enum(
                    pl.concat([dataset._obs[cell_type][column]
                              .cat.get_categories() for dataset in datasets])
                    .unique(maintain_order=True))
                    for column in obs_mismatched_categoricals}
                cell_type_obs = [
                    dataset._obs[cell_type]
                    .cast(cast_dict)
                    .select(obs_columns_in_common) for dataset in datasets]
            else:  # non-flexible
                # Check that all var are identical
                cell_type_var = self._var[cell_type]
                for dataset in datasets[1:]:
                    if not dataset._var[cell_type].equals(cell_type_var):
                        error_message = (
                            f'[{cell_type!r}] all Pseudobulk datasets must '
                            f'have the same var, unless flexible=True')
                        raise ValueError(error_message)
                # Check that all obs have the same columns and data types
                schema = self._obs[cell_type].schema
                for dataset in datasets[1:]:
                    if dataset._obs[cell_type].schema != schema:
                        error_message = (
                            f'[{cell_type!r}] all Pseudobulk datasets must '
                            f'have the same columns in obs, with the same '
                            f'data types, unless flexible=True')
                        raise ValueError(error_message)
                cell_type_X = [dataset._X[cell_type] for dataset in datasets]
                cell_type_obs = [dataset._obs[cell_type]
                                 for dataset in datasets]
            # Concatenate
            X[cell_type] = np.vstack(cell_type_X)
            obs[cell_type] = pl.concat(cell_type_obs)
            var[cell_type] = cell_type_var
        return Pseudobulk(X=X, obs=obs, var=var)

    def concat_var(self,
                   datasets: Pseudobulk,
                   *more_datasets: Pseudobulk,
                   flexible: bool = False) -> Pseudobulk:
        """
        Concatenate the genes of multiple Pseudobulk datasets. All datasets
        must have the same cell types.
        
        By default, all datasets must have the same obs. They must also have
        the same columns in var, with the same data types.
        
        Conversely, if `flexible=True`, subset to cells present in all
        datasets (according to the first column of obs, i.e. `obs_names`)
        before concatenating. Subset to columns of obs that are identical in
        all datasets after this subsetting. Also, subset to columns of var that
        are present in all datasets, and have the same data types. All
        datasets' `obs_names` must have the same name and dtype, and similarly
        for `var_names`.
        
        The one exception to the var "same data type" rule: if a column is Enum
        in some datasets and Categorical in others, or Enum in all datasets but
        with different categories in each dataset, that column will be retained
        as an Enum column (with the union of the categories) in the
        concatenated var.
        
        Args:
            datasets: one or more Pseudobulk datasets to concatenate with this
                      one
            *more_datasets: additional Pseudobulk datasets to concatenate with
                            this one, specified as positional arguments
            flexible: whether to subset to cells and columns of obs and var
                      common to all datasets before concatenating, rather than
                      raising an error on any mismatches
        
        Returns:
            The concatenated Pseudobulk dataset.
        """
        # Check inputs
        if isinstance(datasets, Pseudobulk):
            datasets = datasets,
        datasets = (self,) + datasets + more_datasets
        if len(datasets) == 1:
            error_message = \
                'need at least one other Pseudobulk dataset to concatenate'
            raise ValueError(error_message)
        check_types(datasets[1:], 'datasets', Pseudobulk,
                    'Pseudobulk datasets')
        check_type(flexible, 'flexible', bool, 'Boolean')
        # Check that cell types match across all datasets
        if not all(self.keys() == dataset.keys() for dataset in datasets[1:]):
            error_message = \
                'not all Pseudobulk datasets have the same cell types'
            raise ValueError(error_message)
        # Perform either flexible or non-flexible concatenation
        X = {}
        obs = {}
        var = {}
        for cell_type in self:
            if flexible:
                # Check that `var_names` and `obs_names` have the same name and
                # data type for each cell type across all datasets
                var_names_name = self._var[cell_type][:, 0].name
                if not all(dataset._var[cell_type][:, 0] == var_names_name
                           for dataset in datasets[1:]):
                    error_message = (
                        f'[{cell_type!r}] not all Pseudobulk datasets have '
                        f'the same name for the first column of var (the '
                        f'var_names column)')
                    raise ValueError(error_message)
                obs_names_name = self._obs[cell_type][:, 0].name
                if not all(dataset._obs[cell_type][:, 0].name == obs_names_name
                           for dataset in datasets[1:]):
                    error_message = (
                        f'[{cell_type!r}] not all Pseudobulk datasets have '
                        f'the same name for the first column of obs (the '
                        f'obs_names column)')
                    raise ValueError(error_message)
                var_names_dtype = self._var[cell_type][:, 0].dtype
                if not all(dataset._var[cell_type][:, 0].dtype ==
                           var_names_dtype for dataset in datasets[1:]):
                    error_message = (
                        f'[{cell_type!r}] not all Pseudobulk datasets have '
                        f'the same data type for the first column of var (the '
                        f'var_names column)')
                    raise TypeError(error_message)
                obs_names_dtype = self._obs[cell_type][:, 0].dtype
                if not all(dataset._obs[cell_type][:, 0].dtype ==
                           obs_names_dtype for dataset in datasets[1:]):
                    error_message = (
                        f'[{cell_type!r}] not all Pseudobulk datasets have '
                        f'the same data type for the first column of obs (the '
                        f'obs_names column)')
                    raise TypeError(error_message)
                # Subset to genes in common across all datasets
                genes_in_common = self._obs[cell_type][:, 0]\
                    .filter(self._obs[cell_type][:, 0]
                            .is_in(pl.concat([dataset._obs[cell_type][:, 0]
                                              for dataset in datasets[1:]])))
                if len(genes_in_common) == 0:
                    error_message = (
                        f'[{cell_type!r}] no genes are shared across all '
                        f'Pseudobulk datasets')
                    raise ValueError(error_message)
                cell_type_X = []
                cell_type_obs = []
                for dataset in datasets:
                    gene_indices = dataset._getitem_process(
                        genes_in_common, 1, dataset._obs[cell_type])
                    cell_type_X.append(
                        dataset._X[cell_type][:, gene_indices.to_numpy()])
                    cell_type_obs.append(dataset._obs[cell_type][gene_indices])
                # Subset to columns of obs that are identical in all datasets
                # after this subsetting
                obs_columns_in_common = [
                    column.name for column in cell_type_obs[0][:, 1:]
                    if all(column.name in dataset_cell_type_obs and
                           dataset_cell_type_obs[column.name].equals(column)
                           for dataset_cell_type_obs in cell_type_obs[1:])]
                cell_type_obs = cell_type_obs[0]
                cell_type_obs = cell_type_obs.select(cell_type_obs.columns[0],
                                                     obs_columns_in_common)
                # Subset to columns of var that are present in all datasets,
                # and have the same data types. Also include columns of var
                # that are Enum in some datasets and Categorical in others, or
                # Enum in all datasets but with different categories in each
                # dataset; cast these to Categorical.
                var_mismatched_categoricals = {
                    column for column, dtype in self._var[cell_type][:, 1:]
                    .select(pl.col(pl.Categorical, pl.Enum)).schema.items()
                    if all(column in dataset._var[cell_type] and
                           dataset._var[cell_type][column].dtype in
                           (pl.Categorical, pl.Enum)
                           for dataset in datasets[1:]) and
                       not all(dataset._var[cell_type][column].dtype == dtype
                               for dataset in datasets[1:])}
                var_columns_in_common = [
                    column
                    for column, dtype in islice(
                        self._var[cell_type].schema.items(), 1, None)
                    if column in var_mismatched_categoricals or
                       all(column in dataset[cell_type]._var and
                           dataset._var[cell_type][column].dtype == dtype
                           for dataset in datasets[1:])]
                cast_dict = {column: pl.Enum(
                    pl.concat([dataset._var[cell_type][column]
                              .cat.get_categories() for dataset in datasets])
                    .unique(maintain_order=True))
                    for column in var_mismatched_categoricals}
                cell_type_var = [
                    dataset._var[cell_type]
                    .cast(cast_dict)
                    .select(var_columns_in_common) for dataset in datasets]
            else:  # non-flexible
                # Check that all obs are identical
                cell_type_obs = self._obs[cell_type]
                for dataset in datasets[1:]:
                    if not dataset._obs[cell_type].equals(cell_type_obs):
                        error_message = (
                            f'[{cell_type!r}] all Pseudobulk datasets must '
                            f'have the same obs, unless flexible=True')
                        raise ValueError(error_message)
                # Check that all var have the same columns and data types
                schema = self._var[cell_type].schema
                for dataset in datasets[1:]:
                    if dataset._var[cell_type].schema != schema:
                        error_message = (
                            f'[{cell_type!r}] all Pseudobulk datasets must '
                            f'have the same columns in var, with the same '
                            f'data types, unless flexible=True')
                        raise ValueError(error_message)
                cell_type_X = [dataset._X[cell_type] for dataset in datasets]
                cell_type_var = [dataset._var[cell_type]
                                 for dataset in datasets]
            # Concatenate
            X[cell_type] = np.hstack(cell_type_X)
            var[cell_type] = pl.concat(cell_type_var)
            obs[cell_type] = cell_type_obs
        return Pseudobulk(X=X, obs=obs, var=var)
    
    def _get_columns(self,
                     obs_or_var_name: Literal['obs', 'var'],       
                     columns: PseudobulkColumn | None |
                              Sequence[PseudobulkColumn | None],
                     variable_name: str,
                     dtypes: pl.datatypes.classes.DataTypeClass | str |
                             tuple[pl.datatypes.classes.DataTypeClass | str,
                                   ...],
                     custom_error: str | None = None,
                     allow_None: bool = True,
                     allow_null: bool = False,
                     cell_types: Sequence[str] | None = None) -> \
            dict[str, pl.Series | None]:
        """
        Get a column of the same length as obs or var for each cell type.
        
        Args:
            obs_or_var_name: the name of the DataFrame the column is with
                             respect to, i.e. `'obs'` or `'var'`
            columns: a string naming a column of each cell type's obs/var, a
                     polars expression that evaluates to a single column when 
                     applied to each cell type's obs/var, a polars Series or 
                     NumPy array of the same length as each cell type's 
                     obs/var, or a function that takes in two arguments, `self`
                     and a cell type, and returns a polars Series or NumPy 
                     array of the same length as obs/var. Or, a sequence of
                     any combination of these for each cell type. May also be
                     None (or a Sequence containing None) if `allow_None=True`.
            variable_name: the name of the variable corresponding to `columns`
            dtypes: the required dtype(s) of the column
            custom_error: a custom error message for when (an element of)
                          `columns` is a string and is not found in obs/var;
                          use `{}` as a placeholder for the name of the column
            allow_None: whether to allow `columns` or its elements to be None
            allow_null: whether to allow `columns` to contain null values
            cell_types: a list of cell types; if None, use all cell types. If
                        specified and `column` is a Sequence, `column` and
                        `cell_types` should have the same length.
        
        Returns:
            A dictionary mapping each cell type to a polars Series of the same
            length as the cell type's obs/var. Or, if `columns` is None (or if
            some elements are None), a dict where some or all values are None.
        """
        obs_or_var = self._obs if obs_or_var_name == 'obs' else self._var
        if cell_types is None:
            cell_types = self._X
        if columns is None:
            if not allow_None:
                error_message = f'{variable_name} is None'
                raise TypeError(error_message)
            return {cell_type: None for cell_type in cell_types}
        columns_dict = {}
        if isinstance(columns, str):
            for cell_type in cell_types:
                if columns not in obs_or_var[cell_type]:
                    error_message = (
                        f'{columns!r} is not a column of '
                        f'{obs_or_var_name}[{cell_type!r}]'
                        if custom_error is None else
                        custom_error.format(f'{columns!r}'))
                    raise ValueError(error_message)
                columns_dict[cell_type] = obs_or_var[cell_type][columns]
        elif isinstance(columns, pl.Expr):
            for cell_type in cell_types:
                columns_dict[cell_type] = obs_or_var[cell_type].select(columns)
                if columns_dict[cell_type].width > 1:
                    error_message = (
                        f'{variable_name} is a polars expression that expands '
                        f'to {columns_dict[cell_type].width:,} columns rather '
                        f'than 1 for cell type {cell_type!r}')
                    raise ValueError(error_message)
                columns_dict[cell_type] = columns_dict[cell_type].to_series()
        elif isinstance(columns, pl.Series):
            for cell_type in cell_types:
                if len(columns) != len(obs_or_var[cell_type]):
                    error_message = (
                        f'{variable_name} is a polars Series of length '
                        f'{len(columns):,}, which differs from the length of '
                        f'{obs_or_var_name}[{cell_type!r}] '
                        f'({len(obs_or_var[cell_type]):,})')
                    raise ValueError(error_message)
                columns_dict[cell_type] = columns
        elif isinstance(columns, np.ndarray):
            for cell_type in cell_types:
                if len(columns) != len(obs_or_var[cell_type]):
                    error_message = (
                        f'{variable_name} is a NumPy array of length '
                        f'{len(columns):,}, which differs from the length of '
                        f'{obs_or_var_name}[{cell_type!r}] '
                        f'({len(obs_or_var[cell_type]):,})')
                    raise ValueError(error_message)
                columns_dict[cell_type] = pl.Series(variable_name, columns)
        elif callable(columns):
            function = columns
            for cell_type in cell_types:
                columns = function(self, cell_type)
                if isinstance(columns, np.ndarray):
                    if columns.ndim != 1:
                        error_message = (
                            f'{variable_name} is a function that returns a '
                            f'{columns.ndim:,}D NumPy array, but must return '
                            f'a polars Series or 1D NumPy array')
                        raise ValueError(error_message)
                    columns = pl.Series(variable_name, columns)
                elif not isinstance(columns, pl.Series):
                    error_message = (
                        f'{variable_name} is a function that returns a '
                        f'variable of type {type(columns).__name__}, but must '
                        f'return a polars Series or 1D NumPy array')
                    raise TypeError(error_message)
                if len(columns) != len(obs_or_var[cell_type]):
                    error_message = (
                        f'{variable_name} is a function that returns a column '
                        f'of length {len(columns):,} for cell type '
                        f'{cell_type!r}, which differs from the length of '
                        f'{obs_or_var_name}[{cell_type!r}] '
                        f'({len(obs_or_var[cell_type]):,})')
                    raise ValueError(error_message)
                columns_dict[cell_type] = columns
        elif isinstance(columns, Sequence):
            if len(columns) != len(cell_types):
                error_message = (
                    f'{variable_name} is a sequence of length '
                    f'{len(columns):,}, which differs from the number of cell '
                    f'types ({len(cell_types):,})')
                raise ValueError(error_message)
            if not allow_None and any(column is None for column in columns):
                error_message = \
                    f'{variable_name} contains an element that is None'
                raise TypeError(error_message)
            for index, (column, cell_type) in \
                    enumerate(zip(columns, cell_types)):
                if isinstance(column, str):
                    if column not in obs_or_var[cell_type]:
                        error_message = (
                            f'{column!r} is not a column of '
                            f'{obs_or_var_name}[{cell_type!r}]'
                            if custom_error is None else
                            custom_error.format(f'{column!r}'))
                        raise ValueError(error_message)
                    columns_dict[cell_type] = obs_or_var[cell_type][column]
                elif isinstance(column, pl.Expr):
                    columns_dict[cell_type] = \
                        obs_or_var[cell_type].select(column)
                    if columns[cell_type].width > 1:
                        error_message = (
                            f'{variable_name}[{index}] is a polars expression '
                            f'that expands to {columns[cell_type].width:,} '
                            f'columns rather than 1 for cell type '
                            f'{cell_type!r}')
                        raise ValueError(error_message)
                    columns_dict[cell_type] = \
                        columns_dict[cell_type].to_series()
                elif isinstance(column, pl.Series):
                    if len(column) != len(obs_or_var[cell_type]):
                        error_message = (
                            f'{variable_name}[{index}] is a polars Series of '
                            f'length {len(column):,}, which differs from the '
                            f'length of {obs_or_var_name}[{cell_type!r}] '
                            f'({len(obs_or_var[cell_type]):,})')
                        raise ValueError(error_message)
                    columns_dict[cell_type] = column
                elif isinstance(column, np.ndarray):
                    if len(column) != len(obs_or_var[cell_type]):
                        error_message = (
                            f'{variable_name}[{index}] is a NumPy array of '
                            f'length {len(column):,}, which differs from the '
                            f'length of {obs_or_var_name}[{cell_type!r}] '
                            f'({len(obs_or_var[cell_type]):,})')
                        raise ValueError(error_message)
                    columns_dict[cell_type] = pl.Series(variable_name, column)
                elif callable(columns):
                    column = column(self, cell_type)
                    if isinstance(column, np.ndarray):
                        if column.ndim != 1:
                            error_message = (
                                f'{variable_name}[{index}] is a function that '
                                f'returns a {column.ndim:,}D NumPy array, but '
                                f'must return a polars Series or 1D NumPy '
                                f'array')
                            raise ValueError(error_message)
                        column = pl.Series(variable_name, column)
                    elif not isinstance(column, pl.Series):
                        error_message = (
                            f'{variable_name}[{index}] is a function that '
                            f'returns a variable of type '
                            f'{type(column).__name__}, but must return a '
                            f'polars Series or 1D NumPy array')
                        raise TypeError(error_message)
                    if len(column) != len(obs_or_var[cell_type]):
                        error_message = (
                            f'{variable_name}[{index}] is a function that '
                            f'returns a column of length {len(column):,} for '
                            f'cell type {cell_type!r}, which differs from the '
                            f'length of {obs_or_var_name}[{cell_type!r}] '
                            f'({len(obs_or_var[cell_type]):,})')
                        raise ValueError(error_message)
                    columns_dict[cell_type] = column
                else:
                    error_message = (
                        f'{variable_name}[{index}] must be a string column '
                        f'name, a polars expression or Series, a 1D NumPy '
                        f'array, or a function that returns any of these when '
                        f'applied to this Pseudobulk dataset and a given cell '
                        f'type, but has type {type(column).__name__!r}')
                    raise TypeError(error_message)
        else:
            error_message = (
                f'{variable_name} must be a string column name, a polars '
                f'expression or Series, a 1D NumPy array, or a function that '
                f'returns any of these when applied to this Pseudobulk '
                f'dataset and a given cell type, but has type '
                f'{type(columns).__name__!r}')
            raise TypeError(error_message)
        # Check dtypes
        if not isinstance(dtypes, tuple):
            dtypes = dtypes,
        for cell_type, column in columns_dict.items():
            base_type = column.dtype.base_type()
            for expected_type in dtypes:
                if base_type == expected_type or expected_type == 'integer' \
                        and base_type in pl.INTEGER_DTYPES or \
                        expected_type == 'floating-point' and \
                        base_type in pl.FLOAT_DTYPES:
                    break
            else:
                if len(dtypes) == 1:
                    dtypes = str(dtypes[0])
                elif len(dtypes) == 2:
                    dtypes = ' or '.join(map(str, dtypes))
                else:
                    dtypes = ', '.join(map(str, dtypes[:-1])) + ', or ' + \
                             str(dtypes[-1])
                error_message = (
                    f'{variable_name} must be {dtypes}, but has data type '
                    f'{base_type!r} for cell type {cell_type!r}')
                raise TypeError(error_message)
        # Check nulls, if `allow_null=False`
        if not allow_null:
            for cell_type, column in columns_dict.items():
                null_count = column.null_count()
                if null_count > 0:
                    error_message = (
                        f'{variable_name} contains {null_count:,} '
                        f'{plural("null value", null_count)} for cell type '
                        f'{cell_type!r}, but must not contain any')
                    raise ValueError(error_message)
        return columns_dict
    
    def filter_obs(self,
                   *predicates: str | pl.Expr | pl.Series |
                                Iterable[str | pl.Expr | pl.Series] | bool |
                                list[bool] | np.ndarray[1, np.bool_],
                   **constraints: Any) -> Pseudobulk:
        """
        Equivalent to `df.filter()` from polars, but applied to both obs and X
        for each cell type.
        
        Args:
            *predicates: one or more column names, expressions that evaluate to
                         Boolean Series, Boolean Series, lists of Booleans,
                         and/or 1D Boolean NumPy arrays
            **constraints: column filters: `name=value` filters to samples
                           where the column named `name` has the value `value`
        
        Returns:
            A new Pseudobulk dataset filtered to samples passing all the
            Boolean filters in `predicates` and `constraints`.
        """
        X = {}
        obs = {}
        for cell_type in self:
            obs[cell_type] = self._obs[cell_type]\
                .with_columns(__Pseudobulk_index=pl.int_range(pl.len(),
                                                              dtype=pl.Int32))\
                .filter(*predicates, **constraints)
            X[cell_type] = self._X[cell_type][
                obs[cell_type]['__Pseudobulk_index'].to_numpy()]
            obs[cell_type] = obs[cell_type].drop('__Pseudobulk_index')
        return Pseudobulk(X=X, obs=obs, var=self._var)
    
    def filter_var(self,
                   *predicates: pl.Expr | pl.Series | str |
                                Iterable[pl.Expr | pl.Series | str] | bool |
                                list[bool] | np.ndarray[1, np.bool_],
                   **constraints: Any) -> Pseudobulk:
        """
        Equivalent to `df.filter()` from polars, but applied to both var and X
        for each cell type.
        
        Args:
            *predicates: one or more column names, expressions that evaluate to
                         Boolean Series, Boolean Series, lists of Booleans,
                         and/or 1D Boolean NumPy arrays
            **constraints: column filters: `name=value` filters to genes
                           where the column named `name` has the value `value`
        
        Returns:
            A new Pseudobulk dataset filtered to genes passing all the
            Boolean filters in `predicates` and `constraints`.
        """
        X = {}
        var = {}
        for cell_type in self:
            var[cell_type] = self._var[cell_type]\
                .with_columns(__Pseudobulk_index=pl.int_range(pl.len(),
                                                              dtype=pl.Int32))\
                .filter(*predicates, **constraints)
            X[cell_type] = self._X[cell_type][
                :, var[cell_type]['__Pseudobulk_index'].to_numpy()]
            var[cell_type] = var[cell_type].drop('__Pseudobulk_index')
        return Pseudobulk(X=X, obs=self._obs, var=var)
    
    def select_obs(self,
                   *exprs: Scalar | pl.Expr | pl.Series |
                           Iterable[Scalar | pl.Expr | pl.Series],
                   **named_exprs: Scalar | pl.Expr | pl.Series) -> Pseudobulk:
        """
        Equivalent to `df.select()` from polars, but applied to each cell
        type's obs. obs_names will be automatically included as the first
        column, if not included explicitly.
        
        Args:
            *exprs: column(s) to select, specified as positional arguments.
                    Accepts expression input. Strings are parsed as column
                    names, other non-expression inputs are parsed as literals.
            **named_exprs: additional columns to select, specified as keyword
                           arguments. The columns will be renamed to the
                           keyword used.
        
        Returns:
            A new Pseudobulk dataset with
            obs[cell_type]=obs[cell_type].select(*exprs, **named_exprs) for all
            cell types in obs, and obs_names as the first column unless already
            included explicitly.
        """
        obs = {}
        for cell_type, cell_type_obs in self._obs.items():
            new_cell_type_obs = cell_type_obs.select(*exprs, **named_exprs)
            if cell_type_obs.columns[0] not in new_cell_type_obs:
                new_cell_type_obs = \
                    new_cell_type_obs.select(cell_type_obs[:, 0], pl.all())
            obs[cell_type] = new_cell_type_obs
        return Pseudobulk(X=self._X, obs=obs, var=self._var)
    
    def select_var(self,
                   *exprs: Scalar | pl.Expr | pl.Series |
                           Iterable[Scalar | pl.Expr | pl.Series],
                   **named_exprs: Scalar | pl.Expr | pl.Series) -> Pseudobulk:
        """
        Equivalent to `df.select()` from polars, but applied to each cell
        type's var. var_names will be automatically included as the first
        column, if not included explicitly.
        
        Args:
            *exprs: column(s) to select, specified as positional arguments.
                    Accepts expression input. Strings are parsed as column
                    names, other non-expression inputs are parsed as literals.
            **named_exprs: additional columns to select, specified as keyword
                           arguments. The columns will be renamed to the
                           keyword used.
        
        Returns:
            A new Pseudobulk dataset with
            var[cell_type]=var[cell_type].select(*exprs, **named_exprs) for all
            cell types in var, and var_names as the first column unless already
            included explicitly.
        """
        var = {}
        for cell_type, cell_type_var in self._var.items():
            new_cell_type_var = cell_type_var.select(*exprs, **named_exprs)
            if cell_type_var.columns[0] not in new_cell_type_var:
                new_cell_type_var = \
                    new_cell_type_var.select(cell_type_var[:, 0], pl.all())
            var[cell_type] = new_cell_type_var
        return Pseudobulk(X=self._X, obs=self._obs, var=var)
    
    def select_cell_types(self, cell_types: str, *more_cell_types: str) -> \
            Pseudobulk:
        """
        Create a new Pseudobulk dataset subset to the cell type(s) in
        `cell_types` and `more_cell_types`.
        
        Args:
            cell_types: cell type(s) to select
            *more_cell_types: additional cell types to select, specified as
                              positional arguments
        
        Returns:
            A new Pseudobulk dataset subset to the specified cell type(s).
        """
        cell_types = to_tuple(cell_types) + more_cell_types
        check_types(cell_types, 'cell_types', str, 'strings')
        for cell_type in cell_types:
            if cell_type not in self:
                error_message = (
                    f'tried to select {cell_type!r}, which is not a cell type '
                    f'in this Pseudobulk')
                raise ValueError(error_message)
        return Pseudobulk(X={cell_type: self._X[cell_type]
                             for cell_type in cell_types},
                          obs={cell_type: self._obs[cell_type]
                               for cell_type in cell_types},
                          var={cell_type: self._var[cell_type]
                               for cell_type in cell_types})
    
    def with_columns_obs(self,
                         *exprs: Scalar | pl.Expr | pl.Series |
                                 Iterable[Scalar | pl.Expr | pl.Series],
                         **named_exprs: Scalar | pl.Expr | pl.Series) -> \
            Pseudobulk:
        """
        Equivalent to `df.with_columns()` from polars, but applied to each cell
        type's obs.
        
        Args:
            *exprs: column(s) to add, specified as positional arguments.
                    Accepts expression input. Strings are parsed as column
                    names, other non-expression inputs are parsed as literals.
            **named_exprs: additional columns to add, specified as keyword
                           arguments. The columns will be renamed to the
                           keyword used.
        
        Returns:
            A new Pseudobulk dataset with
            obs[cell_type]=obs[cell_type].with_columns(*exprs, **named_exprs)
            for all cell types in obs.
        """
        return Pseudobulk(X=self._X, obs={
            cell_type: obs.with_columns(*exprs, **named_exprs)
            for cell_type, obs in self._obs.items()}, var=self._var)
    
    def with_columns_var(self,
                         *exprs: Scalar | pl.Expr | pl.Series |
                                 Iterable[Scalar | pl.Expr | pl.Series],
                         **named_exprs: Scalar | pl.Expr | pl.Series) -> \
            Pseudobulk:
        """
        Equivalent to `df.with_columns()` from polars, but applied to each cell
        type's var.
        
        Args:
            *exprs: column(s) to add, specified as positional arguments.
                    Accepts expression input. Strings are parsed as column
                    names, other non-expression inputs are parsed as literals.
            **named_exprs: additional columns to add, specified as keyword
                           arguments. The columns will be renamed to the
                           keyword used.
        
        Returns:
            A new Pseudobulk dataset with
            var[cell_type]=var[cell_type].with_columns(*exprs, **named_exprs)
            for all cell types in var.
        """
        return Pseudobulk(X=self._X, obs=self._obs, var={
            cell_type: var.with_columns(*exprs, **named_exprs)
            for cell_type, var in self._var.items()})

    def drop_obs(self,
                 columns: pl.type_aliases.ColumnNameOrSelector |
                          Iterable[pl.type_aliases.ColumnNameOrSelector],
                 *more_columns: pl.type_aliases.ColumnNameOrSelector) -> \
            Pseudobulk:
        """
        Create a new Pseudobulk dataset with `columns` and `more_columns`
        removed from obs.
        
        Args:
            columns: columns(s) to drop
            *more_columns: additional columns to drop, specified as
                           positional arguments
        
        Returns:
            A new Pseudobulk dataset with the column(s) removed.
        """
        columns = to_tuple(columns) + more_columns
        return Pseudobulk(X=self._X,
                          obs={cell_type: self._obs[cell_type].drop(columns)
                               for cell_type in self}, var=self._var)

    def drop_var(self,
                 columns: pl.type_aliases.ColumnNameOrSelector |
                          Iterable[pl.type_aliases.ColumnNameOrSelector],
                 *more_columns: pl.type_aliases.ColumnNameOrSelector) -> \
            Pseudobulk:
        """
        Create a new Pseudobulk dataset with `columns` and `more_columns`
        removed from var.
        
        Args:
            columns: columns(s) to drop
            *more_columns: additional columns to drop, specified as
                           positional arguments
        
        Returns:
            A new Pseudobulk dataset with the column(s) removed.
        """
        columns = to_tuple(columns) + more_columns
        return Pseudobulk(X=self._X, obs=self._obs,
                          var={cell_type: self._var[cell_type].drop(columns)
                               for cell_type in self})
    
    def drop_cell_types(self, cell_types: str, *more_cell_types: str) -> \
            Pseudobulk:
        """
        Create a new Pseudobulk dataset with `cell_types` and `more_cell_types`
        removed. Raises an error if all cell types would be dropped.
        
        Args:
            cell_types: cell type(s) to drop
            *more_cell_types: additional cell types to drop, specified as
                              positional arguments
        
        Returns:
            A new Pseudobulk dataset with the cell type(s) removed.
        """
        cell_types = set(to_tuple(cell_types)) | set(more_cell_types)
        check_types(cell_types, 'cell_types', str, 'strings')
        # noinspection PyTypeChecker
        original_cell_types = set(self)
        if not cell_types < original_cell_types:
            if cell_types == original_cell_types:
                error_message = 'all cell types would be dropped'
                raise ValueError(error_message)
            for cell_type in cell_types:
                if cell_type not in original_cell_types:
                    error_message = (
                        f'tried to drop {cell_type!r}, which is not a cell '
                        f'type in this Pseudobulk')
                    raise ValueError(error_message)
        new_cell_types = \
            [cell_type for cell_type in self if cell_type not in cell_types]
        return Pseudobulk(X={cell_type: self._X[cell_type]
                             for cell_type in new_cell_types},
                          obs={cell_type: self._obs[cell_type]
                               for cell_type in new_cell_types},
                          var={cell_type: self._var[cell_type]
                               for cell_type in new_cell_types})
    
    def rename_obs(self, mapping: dict[str, str] | Callable[[str], str]) -> \
            Pseudobulk:
        """
        Create a new Pseudobulk dataset with column(s) of obs renamed for each
        cell type.
        
        Args:
            mapping: the renaming to apply, either as a dictionary with the old
                     names as keys and the new names as values, or a function
                     that takes an old name and returns a new name
        
        Returns:
            A new Pseudobulk dataset with the column(s) of obs renamed.
        """
        return Pseudobulk(X=self._X, obs={
            cell_type: self._obs[cell_type].rename(mapping)
            for cell_type in self}, var=self._var)
    
    def rename_var(self, mapping: dict[str, str] | Callable[[str], str]) -> \
            Pseudobulk:
        """
        Create a new Pseudobulk dataset with column(s) of var renamed for each
        cell type.
        
        Args:
            mapping: the renaming to apply, either as a dictionary with the old
                     names as keys and the new names as values, or a function
                     that takes an old name and returns a new name
        
        Returns:
            A new Pseudobulk dataset with the column(s) of var renamed.
        """
        return Pseudobulk(X=self._X, obs=self._obs, var={
            cell_type: self._var[cell_type].rename(mapping)
            for cell_type in self})
    
    def rename_cell_types(self,
                          mapping: dict[str, str] | Callable[[str], str]) -> \
            Pseudobulk:
        """
        Create a new Pseudobulk dataset with cell type(s) renamed.
        
        Args:
            mapping: the renaming to apply, either as a dictionary with the old
                     cell type names as keys and the new names as values, or a
                     function that takes an old name and returns a new name
        
        Returns:
            A new Pseudobulk dataset with the cell type(s) renamed.
        """
        if isinstance(mapping, dict):
            new_cell_types = [mapping.get(cell_type, cell_type)
                              for cell_type in self._X]
        elif callable(mapping):
            new_cell_types = [mapping(cell_type) for cell_type in self._X]
        else:
            raise TypeError(f'mapping must be a dictionary or function, but '
                            f'has type {type(mapping).__name__!r}')
        return Pseudobulk(X={new_cell_type: X
                             for new_cell_type, X in
                             zip(new_cell_types, self._X.values())},
                          obs={new_cell_type: obs
                               for new_cell_type, obs in
                               zip(new_cell_types, self._obs.values())},
                          var={new_cell_type: var
                               for new_cell_type, var in
                               zip(new_cell_types, self._var.values())})
    
    def cast_X(self, dtype: np._typing.DTypeLike) -> Pseudobulk:
        """
        Cast each cell type's X to the specified data type.
        
        Args:
            dtype: a NumPy data type

        Returns:
            A new Pseudobulk dataset with each cell type's X cast to the
            specified data type.
        """
        return Pseudobulk(X={cell_type: self._X[cell_type].astype(dtype)
                             for cell_type in self},
                          obs=self._obs, var=self._var)
    
    def cast_obs(self,
                 dtypes: Mapping[pl.type_aliases.ColumnNameOrSelector |
                                 pl.type_aliases.PolarsDataType,
                                 pl.type_aliases.PolarsDataType] |
                         pl.type_aliases.PolarsDataType,
                 *,
                 strict: bool = True) -> Pseudobulk:
        """
        Cast column(s) of each cell type's obs to the specified data type(s).
        
        Args:
            dtypes: a mapping of column names (or selectors) to data types, or
                    a single data type to which all columns will be cast
            strict: whether to raise an error if a cast could not be done (for
                    instance, due to numerical overflow)

        Returns:
            A new Pseudobulk dataset with column(s) of each cell type's obs
            cast to the specified data type(s).
        """
        return Pseudobulk(X=self._X,
                          obs={cell_type: self._obs[cell_type].cast(
                              dtypes, strict=strict) for cell_type in self},
                          var=self._var)
    
    def cast_var(self,
                 dtypes: Mapping[pl.type_aliases.ColumnNameOrSelector |
                                 pl.type_aliases.PolarsDataType,
                                 pl.type_aliases.PolarsDataType] |
                         pl.type_aliases.PolarsDataType,
                 *,
                 strict: bool = True) -> Pseudobulk:
        """
        Cast column(s) of each cell type's var to the specified data type(s).
        
        Args:
            dtypes: a mapping of column names (or selectors) to data types, or
                    a single data type to which all columns will be cast
            strict: whether to raise an error if a cast could not be done (for
                    instance, due to numerical overflow)

        Returns:
            A new Pseudobulk dataset with column(s) of each cell type's var
            cast to the specified data type(s).
        """
        return Pseudobulk(X=self._X,
                          obs=self._obs,
                          var={cell_type: self._var[cell_type].cast(
                              dtypes, strict=strict) for cell_type in self})
    
    def join_obs(self,
                 other: pl.DataFrame,
                 on: str | pl.Expr | Sequence[str | pl.Expr] | None = None,
                 *,
                 left_on: str | pl.Expr | Sequence[str | pl.Expr] |
                          None = None,
                 right_on: str | pl.Expr | Sequence[str | pl.Expr] |
                           None = None,
                 suffix: str = '_right',
                 validate: Literal['m:m', 'm:1', '1:m', '1:1'] = 'm:m',
                 join_nulls: bool = False,
                 coalesce: bool = True) -> Pseudobulk:
        """
        Left join each cell type's obs with another DataFrame.
        
        Args:
            other: a polars DataFrame to join each cell type's obs with
            on: the name(s) of the join column(s) in both DataFrames
            left_on: the name(s) of the join column(s) in obs
            right_on: the name(s) of the join column(s) in `other`
            suffix: a suffix to append to columns with a duplicate name
            validate: checks whether the join is of the specified type. Can be:
                      - 'm:m' (many-to-many): the default, no checks performed.
                      - '1:1' (one-to-one): check that none of the values in
                        the join column(s) appear more than once in obs or more
                        than once in `other`.
                      - '1:m' (one-to-many): check that none of the values in
                        the join column(s) appear more than once in obs.
                      - 'm:1' (many-to-one): check that none of the values in
                        the join column(s) appear more than once in `other`.
            join_nulls: whether to include null as a valid value to join on.
                        By default, null values will never produce matches.
            coalesce: if True, coalesce each of the pairs of join columns
                      (the columns in `on` or `left_on`/`right_on`) from obs
                      and `other` into a single column, filling missing values
                      from one with the corresponding values from the other.
                      If False, include both as separate columns, adding
                      `suffix` to the join columns from `other`.
        
        Returns:
            A new Pseudobulk dataset with the columns from `other` joined to
            each cell type's obs.
        
        Note:
            If a column of `on`, `left_on` or `right_on` is Enum in obs and
            Categorical in `other` (or vice versa), or Enum in both but with
            different categories in each, that pair of columns will be
            automatically cast to a common Enum data type (with the union of
            the categories) before joining.
        """
        # noinspection PyTypeChecker
        check_type(other, 'other', pl.DataFrame, 'a polars DataFrame')
        if on is None:
            if left_on is None and right_on is None:
                error_message = (
                    f"either 'on' or both of 'left_on' and 'right_on' must be "
                    f"specified")
                raise ValueError(error_message)
            elif left_on is None:
                error_message = \
                    'right_on is specified, so left_on must be specified'
                raise ValueError(error_message)
            elif right_on is None:
                error_message = \
                    'left_on is specified, so right_on must be specified'
                raise ValueError(error_message)
        else:
            if left_on is not None:
                error_message = "'on' is specified, so 'left_on' must be None"
                raise ValueError(error_message)
            if right_on is not None:
                error_message = "'on' is specified, so 'right_on' must be None"
                raise ValueError(error_message)
        obs = {}
        for cell_type in self:
            left = self._obs[cell_type]
            right = other
            if on is None:
                left_columns = left.select(left_on)
                right_columns = right.select(right_on)
            else:
                left_columns = left.select(on)
                right_columns = right.select(on)
            left_cast_dict = {}
            right_cast_dict = {}
            for left_column, right_column in zip(left_columns, right_columns):
                left_dtype = left_column.dtype
                right_dtype = right_column.dtype
                if left_dtype == right_dtype:
                    continue
                if (left_dtype == pl.Enum or left_dtype == pl.Categorical) \
                        and (right_dtype == pl.Enum or
                             right_dtype == pl.Categorical):
                    common_dtype = \
                        pl.Enum(pl.concat([left_column.cat.get_categories(),
                                           right_column.cat.get_categories()])
                                .unique(maintain_order=True))
                    left_cast_dict[left_column.name] = common_dtype
                    right_cast_dict[right_column.name] = common_dtype
                else:
                    error_message = (
                        f'obs[{cell_type!r}][{left_column.name!r}] has data '
                        f'type {left_dtype.base_type()!r}, but '
                        f'other[{cell_type!r}][{right_column.name!r}] has '
                        f'data type {right_dtype.base_type()!r}')
                    raise TypeError(error_message)
            if left_cast_dict is not None:
                left = left.cast(left_cast_dict)
                right = right.cast(right_cast_dict)
            obs[cell_type] = \
                left.join(right, on=on, how='left', left_on=left_on,
                          right_on=right_on, suffix=suffix, validate=validate,
                          join_nulls=join_nulls, coalesce=coalesce)
            if len(obs[cell_type]) > len(self._obs[cell_type]):
                other_on = to_tuple(right_on if right_on is not None else on)
                assert other.select(other_on).is_duplicated().any()
                duplicate_column = other_on[0] if len(other_on) == 1 else \
                    next(column for column in other_on
                         if other[column].is_duplicated().any())
                error_message = (
                    f'other[{duplicate_column!r}] contains duplicate values, '
                    f'so it must be deduplicated before being joined on')
                raise ValueError(error_message)
        return Pseudobulk(X=self._X, obs=obs, var=self._var)
    
    def join_var(self,
                 other: pl.DataFrame,
                 on: str | pl.Expr | Sequence[str | pl.Expr] | None = None,
                 *,
                 left_on: str | pl.Expr | Sequence[str | pl.Expr] |
                          None = None,
                 right_on: str | pl.Expr | Sequence[str | pl.Expr] |
                           None = None,
                 suffix: str = '_right',
                 validate: Literal['m:m', 'm:1', '1:m', '1:1'] = 'm:m',
                 join_nulls: bool = False,
                 coalesce: bool = True) -> Pseudobulk:
        """
        Join each cell type's var with another DataFrame.
        
        Args:
            other: a polars DataFrame to join each cell type's var with
            on: the name(s) of the join column(s) in both DataFrames
            left_on: the name(s) of the join column(s) in var
            right_on: the name(s) of the join column(s) in `other`
            suffix: a suffix to append to columns with a duplicate name
            validate: checks whether the join is of the specified type. Can be:
                      - 'm:m' (many-to-many): the default, no checks performed.
                      - '1:1' (one-to-one): check that none of the values in
                        the join column(s) appear more than once in var or more
                        than once in `other`.
                      - '1:m' (one-to-many): check that none of the values in
                        the join column(s) appear more than once in var.
                      - 'm:1' (many-to-one): check that none of the values in
                        the join column(s) appear more than once in `other`.
            join_nulls: whether to include null as a valid value to join on.
                        By default, null values will never produce matches.
            coalesce: if True, coalesce each of the pairs of join columns
                      (the columns in `on` or `left_on`/`right_on`) from obs
                      and `other` into a single column, filling missing values
                      from one with the corresponding values from the other.
                      If False, include both as separate columns, adding
                      `suffix` to the join columns from `other`.
        
        Returns:
            A new Pseudobulk dataset with the columns from `other` joined to
            each cell type's var.
        
        Note:
            If a column of `on`, `left_on` or `right_on` is Enum in obs and
            Categorical in `other` (or vice versa), or Enum in both but with
            different categories in each, that pair of columns will be
            automatically cast to a common Enum data type (with the union of
            the categories) before joining.
        """
        check_type(other, 'other', pl.DataFrame, 'a polars DataFrame')
        if on is None:
            if left_on is None and right_on is None:
                error_message = (
                    "either 'on' or both of 'left_on' and 'right_on' must be "
                    "specified")
                raise ValueError(error_message)
            elif left_on is None:
                error_message = \
                    'right_on is specified, so left_on must be specified'
                raise ValueError(error_message)
            elif right_on is None:
                error_message = \
                    'left_on is specified, so right_on must be specified'
                raise ValueError(error_message)
        else:
            if left_on is not None:
                error_message = "'on' is specified, so 'left_on' must be None"
                raise ValueError(error_message)
            if right_on is not None:
                error_message = "'on' is specified, so 'right_on' must be None"
                raise ValueError(error_message)
        var = {}
        for cell_type in self:
            left = self._var[cell_type]
            right = other
            if on is None:
                left_columns = left.select(left_on)
                right_columns = right.select(right_on)
            else:
                left_columns = left.select(on)
                right_columns = right.select(on)
            left_cast_dict = {}
            right_cast_dict = {}
            for left_column, right_column in zip(left_columns, right_columns):
                left_dtype = left_column.dtype
                right_dtype = right_column.dtype
                if left_dtype == right_dtype:
                    continue
                if (left_dtype == pl.Enum or left_dtype == pl.Categorical) \
                        and (right_dtype == pl.Enum or
                             right_dtype == pl.Categorical):
                    common_dtype = \
                        pl.Enum(pl.concat([left_column.cat.get_categories(),
                                           right_column.cat.get_categories()])
                                .unique(maintain_order=True))
                    left_cast_dict[left_column.name] = common_dtype
                    right_cast_dict[right_column.name] = common_dtype
                else:
                    error_message = (
                        f'var[{cell_type!r}][{left_column.name!r}] has data '
                        f'type {left_dtype.base_type()!r}, but '
                        f'other[{cell_type!r}][{right_column.name!r}] has '
                        f'data type {right_dtype.base_type()!r}')
                    raise TypeError(error_message)
            if left_cast_dict is not None:
                left = left.cast(left_cast_dict)
                right = right.cast(right_cast_dict)
            var[cell_type] = \
                left.join(right, on=on, how='left', left_on=left_on,
                          right_on=right_on, suffix=suffix, validate=validate,
                          join_nulls=join_nulls, coalesce=coalesce)
            if len(var[cell_type]) > len(self._var[cell_type]):
                other_on = to_tuple(right_on if right_on is not None else on)
                assert other.select(other_on).is_duplicated().any()
                duplicate_column = other_on[0] if len(other_on) == 1 else \
                    next(column for column in other_on
                         if other[column].is_duplicated().any())
                error_message = (
                    f'other[{duplicate_column!r}] contains duplicate values, '
                    f'so it must be deduplicated before being joined on')
                raise ValueError(error_message)
        return Pseudobulk(X=self._X, obs=self._obs, var=var)
    
    def peek_obs(self, cell_type: str | None = None, row: int = 0) -> None:
        """
        Print a row of obs (the first row, by default) for a cell type (the
        first cell type, by default) with each column on its own line.
        
        Args:
            cell_type: the cell type to print the row for, or None to use the
                       first cell type
            row: the index of the row to print
        """
        if cell_type is None:
            cell_type = next(iter(self._obs))
        else:
            check_type(cell_type, 'cell_type', str, 'a string')
        check_type(row, 'row', int, 'an integer')
        with pl.Config(tbl_rows=-1):
            print(self._obs[cell_type][row].unpivot(variable_name='column'))
    
    def peek_var(self, cell_type: str | None = None, row: int = 0) -> None:
        """
        Print a row of var (the first row, by default) with each column on its
        own line.
        
        Args:
            cell_type: the cell type to print the row for, or None to use the
                       first cell type
            row: the index of the row to print
        """
        if cell_type is None:
            cell_type = next(iter(self._var))
        else:
            check_type(cell_type, 'cell_type', str, 'a string')
        check_type(row, 'row', int, 'an integer')
        with pl.Config(tbl_rows=-1):
            print(self._var[cell_type][row].unpivot(variable_name='column'))
    
    def subsample_obs(self,
                      n: int | np.integer | None = None,
                      *,
                      fraction: int | float | np.integer | np.floating |
                                None = None,
                      by_column: PseudobulkColumn | None |
                                 Sequence[PseudobulkColumn | None] = None,
                      subsample_column: str | None = None,
                      seed: int | np.integer = 0,
                      overwrite: bool = False) -> Pseudobulk:
        """
        Subsample a specific number or fraction of samples.
        
        Args:
            n: the number of samples to return; mutually exclusive with
               `fraction`
            fraction: the fraction of samples to return; mutually exclusive
                      with `n`
            by_column: an optional String, Categorical, Enum, or integer column
                       of obs to subsample by. Can be None, a column name, a
                       polars expression, a polars Series, a 1D NumPy array, or
                       a function that takes in this Pseudobulk dataset and a
                       cell type and returns a polars Series or 1D NumPy array.
                       Can also be a sequence of any combination of these for
                       each cell type. Specifying `by_column` ensures that the
                       same fraction of cells with each value of `by_column`
                       are subsampled. When combined with `n`, to make sure the
                       total number of samples is exactly `n`, some of the
                       smallest groups may be oversampled by one element, or
                       some of the largest groups can be undersampled by one
                       element. Can contain null entries: the corresponding
                       samples will not be included in the result.
            subsample_column: an optional name of a Boolean column to add to
                              obs indicating the subsampled genes; if None,
                              subset to these genes instead
            seed: the random seed to use when subsampling
            overwrite: if True, overwrite `subsample_column` if already present
                       in obs, instead of raising an error. Must be False when
                       `subsample_column` is None.
        
        Returns:
            A new Pseudobulk dataset subset to the subsampled cells, or if
            `subsample_column` is not None, the full dataset with
            `subsample_column` added to obs.
        """
        if n is not None:
            check_type(n, 'n', int, 'a positive integer')
            check_bounds(n, 'n', 1)
        elif fraction is not None:
            check_type(fraction, 'fraction', float,
                       'a floating-point number between 0 and 1')
            check_bounds(fraction, 'fraction', 0, 1, left_open=True,
                         right_open=True)
        else:
            error_message = 'one of n and fraction must be specified'
            raise ValueError(error_message)
        if n is not None and fraction is not None:
            error_message = 'only one of n and fraction must be specified'
            raise ValueError(error_message)
        by_column = self._get_columns(
            'obs', by_column, 'by_column',
            (pl.String, pl.Categorical, pl.Enum, 'integer'), allow_null=True)
        check_type(overwrite, 'overwrite', bool, 'Boolean')
        if subsample_column is not None:
            check_type(subsample_column, 'subsample_column', str, 'a string')
            if not overwrite:
                for cell_type, obs in self._obs.items():
                    if subsample_column in obs:
                        error_message = (
                            f'subsample_column {subsample_column!r} is '
                            f'already a column of obs[{cell_type!r}]')
                        raise ValueError(error_message)
        elif overwrite:
            error_message = (
                'overwrite must be False when subsample_column is None; did '
                'you already run subsample_obs()? Set overwrite=True to '
                'overwrite.')
            raise ValueError(error_message)
        check_type(seed, 'seed', int, 'an integer')
        by = lambda expr, cell_type: \
            expr if by_column[cell_type] is None else \
            expr.over(by_column[cell_type])
        if by_column is not None and n is not None:
            # Reassign n to be a vector of sample sizes per group, broadcast to
            # the length of obs. The total sample size should exactly match the
            # original n; if necessary, oversample the smallest groups or
            # undersample the largest groups to make this happen.
            cell_type_n = {}
            for cell_type, cell_type_by_column in by_column.items():
                if cell_type_by_column is None:
                    cell_type_n[cell_type] = n
                else:
                    by_frame = cell_type_by_column.to_frame()
                    by_name = cell_type_by_column.name
                    group_counts = by_frame\
                        .group_by(by_name)\
                        .agg(pl.len(), n=(n * pl.len() / len(by_column))
                                         .round().cast(pl.Int32))\
                        .drop_nulls(by_name)
                    diff = n - group_counts['n'].sum()
                    if diff != 0:
                        group_counts = group_counts\
                            .sort('len', descending=diff < 0)\
                            .with_columns(n=pl.col.n +
                                            pl.int_range(pl.len(),
                                                         dtype=pl.Int32)
                                            .lt(abs(diff)).cast(pl.Int32) *
                                            pl.lit(diff).sign())
                    cell_type_n[cell_type] = \
                        group_counts.join(by_frame, on=by_name)['n']
        # noinspection PyUnboundLocalVariable,PyUnresolvedReferences
        expressions = {
            cell_type: pl.int_range(pl.len(), dtype=pl.Int32)
                       .shuffle(seed=seed)
                       .pipe(by, cell_type=cell_type)
                       .lt((cell_type_n[cell_type] if by_column is not None
                            else n) if fraction is None else
                           fraction * pl.len().pipe(by, cell_type=cell_type))
                       for cell_type in self}
        if subsample_column is None:
            X = {}
            obs = {}
            for cell_type in self:
                obs[cell_type] = self._obs[cell_type]\
                    .with_columns(__Pseudobulk_index=pl.int_range(
                        pl.len(), dtype=pl.Int32))\
                    .filter(expressions[cell_type])
                X[cell_type] = self._X[cell_type][
                    obs[cell_type]['__Pseudobulk_index'].to_numpy()]
                obs[cell_type] = obs[cell_type].drop('__Pseudobulk_index')
            return Pseudobulk(X=X, obs=obs, var=self._var)
        else:
            return Pseudobulk(X=self._X, obs={
                cell_type: obs.with_columns(expressions[cell_type]
                                            .alias(subsample_column))
                for cell_type, obs in self._obs.items()}, var=self._var)
    
    def subsample_var(self,
                      n: int | np.integer | None = None,
                      *,
                      fraction: int | float | np.integer | np.floating |
                                None = None,
                      by_column: PseudobulkColumn | None |
                                 Sequence[PseudobulkColumn | None] = None,
                      subsample_column: str | None = None,
                      seed: int | np.integer = 0,
                      overwrite: bool = False) -> Pseudobulk:
        """
        Subsample a specific number or fraction of genes.
        
        Args:
            n: the number of genes to return; mutually exclusive with
               `fraction`
            fraction: the fraction of genes to return; mutually exclusive with
                      `n`
            by_column: an optional String, Categorical, Enum, or integer column
                       of var to subsample by. Can be None, a column name, a
                       polars expression, a polars Series, a 1D NumPy array, or
                       a function that takes in this Pseudobulk dataset and a
                       cell type and returns a polars Series or 1D NumPy array.
                       Can also be a sequence of any combination of these for
                       each cell type. Specifying `by_column` ensures that the
                       same fraction of genes with each value of `by_column`
                       are subsampled. When combined with `n`, to make sure the
                       total number of samples is exactly `n`, some of the
                       smallest groups may be oversampled by one element, or
                       some of the largest groups may be undersampled by one
                       element. Can contain null entries: the corresponding
                       genes will not be included in the result.
            subsample_column: an optional name of a Boolean column to add to
                              var indicating the subsampled genes; if None,
                              subset to these genes instead
            seed: the random seed to use when subsampling
            overwrite: if True, overwrite `subsample_column` if already present
                       in var, instead of raising an error. Must be False when
                       `subsample_column` is None.

        Returns:
            A new Pseudobulk dataset subset to the subsampled genes, or if
            `subsample_column` is not None, the full dataset with
            `subsample_column` added to var.
        """
        if n is not None:
            check_type(n, 'n', int, 'a positive integer')
            check_bounds(n, 'n', 1)
        elif fraction is not None:
            check_type(fraction, 'fraction', float,
                       'a floating-point number between 0 and 1')
            check_bounds(fraction, 'fraction', 0, 1, left_open=True,
                         right_open=True)
        else:
            error_message = 'one of n and fraction must be specified'
            raise ValueError(error_message)
        if n is not None and fraction is not None:
            error_message = 'only one of n and fraction must be specified'
            raise ValueError(error_message)
        by_column = self._get_columns(
            'var', by_column, 'by_column',
            (pl.String, pl.Categorical, pl.Enum, 'integer'), allow_null=True)
        check_type(overwrite, 'overwrite', bool, 'Boolean')
        if subsample_column is not None:
            check_type(subsample_column, 'subsample_column', str, 'a string')
            if not overwrite:
                for cell_type, var in self._var.items():
                    if subsample_column in var:
                        error_message = (
                            f'subsample_column {subsample_column!r} is '
                            f'already a column of var[{cell_type!r}]')
                        raise ValueError(error_message)
        elif overwrite:
            error_message = (
                'overwrite must be False when subsample_column is None; did '
                'you already run subsample_var()? Set overwrite=True to '
                'overwrite.')
            raise ValueError(error_message)
        check_type(seed, 'seed', int, 'an integer')
        by = lambda expr, cell_type: \
            expr if by_column[cell_type] is None else \
            expr.over(by_column[cell_type])
        if by_column is not None and n is not None:
            # Reassign n to be a vector of sample sizes per group, broadcast to
            # the length of var. The total sample size should exactly match the
            # original n; if necessary, oversample the smallest groups or
            # undersample the largest groups to make this happen.
            cell_type_n = {}
            for cell_type, cell_type_by_column in by_column.items():
                if cell_type_by_column is None:
                    cell_type_n[cell_type] = n
                else:
                    by_frame = cell_type_by_column.to_frame()
                    by_name = cell_type_by_column.name
                    group_counts = by_frame\
                        .group_by(by_name)\
                        .agg(pl.len(), n=(n * pl.len() / len(by_column))
                                         .round().cast(pl.Int32))\
                        .drop_nulls(by_name)
                    diff = n - group_counts['n'].sum()
                    if diff != 0:
                        group_counts = group_counts\
                            .sort('len', descending=diff < 0)\
                            .with_columns(n=pl.col.n +
                                            pl.int_range(pl.len(),
                                                       dtype=pl.Int32)
                                            .lt(abs(diff)).cast(pl.Int32) *
                                            pl.lit(diff).sign())
                    cell_type_n[cell_type] = \
                        group_counts.join(by_frame, on=by_name)['n']
        # noinspection PyUnboundLocalVariable,PyUnresolvedReferences
        expressions = {
            cell_type: pl.int_range(pl.len(), dtype=pl.Int32)
                       .shuffle(seed=seed)
                       .pipe(by, cell_type=cell_type)
                       .lt((cell_type_n[cell_type] if by_column is not None
                            else n) if fraction is None else
                           fraction * pl.len().pipe(by, cell_type=cell_type))
                       for cell_type in self}
        if subsample_column is None:
            X = {}
            var = {}
            for cell_type in self:
                var[cell_type] = self._var[cell_type]\
                    .with_columns(__Pseudobulk_index=pl.int_range(
                        pl.len(), dtype=pl.Int32))\
                    .filter(expressions[cell_type])
                X[cell_type] = self._X[cell_type][
                    :, var[cell_type]['__Pseudobulk_index'].to_numpy()]
                var[cell_type] = var[cell_type].drop('__Pseudobulk_index')
            return Pseudobulk(X=X, obs=self._obs, var=var)
        else:
            return Pseudobulk(X=self._X, obs=self._obs, var={
                cell_type: var.with_columns(expressions[cell_type]
                                            .alias(subsample_column))
                for cell_type, var in self._var.items()})
    
    def pipe(self,
             function: Callable[[Pseudobulk, ...], Any],
             *args: Any,
             **kwargs: Any) -> Any:
        """
        Apply a function to a Pseudobulk dataset.
        
        Args:
            function: the function to apply
            *args: the positional arguments to the function
            **kwargs: the keyword arguments to the function

        Returns:
            function(self, *args, **kwargs)
        """
        return function(self, *args, **kwargs)
    
    def pipe_X(self,
               function: Callable[[np.ndarray[2, np.integer | np.floating],
                                   ...],
                                  np.ndarray[2, np.integer | np.floating]],
               *args: Any,
               **kwargs: Any) -> Pseudobulk:
        """
        Apply a function to a Pseudobulk dataset's X.
        
        Args:
            function: the function to apply to X
            *args: the positional arguments to the function
            **kwargs: the keyword arguments to the function

        Returns:
            A new Pseudobulk dataset where the function has been applied to X.
        """
        return Pseudobulk(X={cell_type: function(self._X, *args, **kwargs)
                             for cell_type in self},
                          obs=self._obs, var=self._var)
    
    def pipe_obs(self,
                 function: Callable[[pl.DataFrame, ...], pl.DataFrame],
                 *args: Any,
                 **kwargs: Any) -> Pseudobulk:
        """
        Apply a function to a Pseudobulk dataset's obs for each cell type.
        
        Args:
            function: the function to apply to each cell type's obs
            *args: the positional arguments to the function
            **kwargs: the keyword arguments to the function

        Returns:
            A new Pseudobulk dataset where the function has been applied to
            each cell type's obs.
        """
        return Pseudobulk(X=self._X, obs={
            cell_type: function(self._obs[cell_type], *args, **kwargs)
            for cell_type in self}, var=self._var)
    
    def pipe_var(self,
                 function: Callable[[pl.DataFrame, ...], pl.DataFrame],
                 *args: Any,
                 **kwargs: Any) -> Pseudobulk:
        """
        Apply a function to a Pseudobulk dataset's var for each cell type.
        
        Args:
            function: the function to apply to each cell type's var
            *args: the positional arguments to the function
            **kwargs: the keyword arguments to the function

        Returns:
            A new Pseudobulk dataset where the function has been applied to
            each cell type's var.
        """
        return Pseudobulk(X=self._X, obs=self._obs, var={
            cell_type: function(self._var[cell_type], *args, **kwargs)
            for cell_type in self})
    
    @staticmethod
    def _too_few_samples(obs: pl.DataFrame,
                         case_control_mask: pl.Series,
                         min_samples: int | np.integer,
                         cell_type: str,
                         verbose: bool) -> bool:
        """
        Skip cell types with fewer than `min_samples` cases or `min_samples`
        controls, or with fewer than `min_samples` total samples if
        non-case-control.
        
        Args:
            obs: the cell type's obs, after filtering
            case_control_mask: the case-control mask, after filtering
            min_samples: filter to cell types with at least this many cases and
                         this many controls, or with at least this many total
                         samples if `case_control_mask` is None
            cell_type: the name of the cell type
            verbose: whether to explain why the cell type is being skipped, if
                     it is

        Returns:
            Whether this cell type has too few samples and should be skipped.
        """
        if case_control_mask is not None:
            num_cases = case_control_mask.sum()
            if num_cases < min_samples:
                if verbose:
                    print(f'[{cell_type}] Skipping this cell type '
                          f'because it has only {num_cases:,} '
                          f'{plural("case", num_cases)} after filtering, '
                          f'which is fewer than min_samples '
                          f'({min_samples:,})')
                return True
            num_controls = len(case_control_mask) - num_cases - \
                           case_control_mask.null_count()
            if num_controls < min_samples:
                if verbose:
                    print(f'[{cell_type}] Skipping this cell type '
                          f'because it has only {num_controls:,} '
                          f'{plural("control", num_controls)} after '
                          f'filtering, which is fewer than min_samples '
                          f'{min_samples:,})')
                return True
        else:
            num_samples = len(obs)
            if num_samples < min_samples:
                if verbose:
                    print(f'[{cell_type}] Skipping this cell type because '
                          f'it has only {num_samples:,} '
                          f'{plural("sample", num_samples)} after '
                          f'filtering, which is fewer than min_samples '
                          f'({min_samples:,})')
                return True
        return False
    
    def qc(self,
           case_control_column: PseudobulkColumn | None |
                                Sequence[PseudobulkColumn | None],
           *,
           custom_filter: PseudobulkColumn | None |
                          Sequence[PseudobulkColumn | None] = None,
           min_samples: int | np.integer = 2,
           min_cells: int | np.integer | None = 10,
           max_standard_deviations: int | float | np.integer | np.floating |
                                    None = 3,
           min_nonzero_fraction: int | float | np.integer | np.floating |
                                 None = 0.8,
           error_if_negative_counts: bool = True,
           allow_float: bool = False,
           verbose: bool = True) -> Pseudobulk:
        """
        Subsets each cell type to samples passing quality control (QC).
        This is different from `SingleCell.qc()`, which (for memory efficiency)
        just adds a Boolean column to obs of which cells passed QC.
        
        Filters, in order, to:
        - samples that have at least `min_cells` cells of that type (default:
          10), and pass the `custom_filter` if specified
        - samples where the number of genes with 0 counts is at most
          `max_standard_deviations` standard deviations above the mean
          (default: 3)
        - genes with at least 1 count in `100 * min_nonzero_fraction`% of
          controls AND `100 * min_nonzero_fraction`% of cases (default: 80%),
          or if `case_control_column` is None, at least one count in
          `100 * min_nonzero_fraction`% of samples
        
        Also filters to cell types with at least `min_samples` cases and
        `min_samples` controls after applying the above sample-level filters,
        or `min_samples` total samples if `case_control_column` is None.
        
        Args:
            case_control_column: an optional column of obs with case-control
                                 labels; set to None for non-case-control data.
                                 Can be None, a column name, a polars
                                 expression, a polars Series, a 1D NumPy array,
                                 or a function that takes in this Pseudobulk
                                 dataset and a cell type and returns a polars
                                 Series or 1D NumPy array. Can also be a
                                 sequence of any combination of these for each
                                 cell type. Not used when
                                 `min_nonzero_fraction` is None. Must be
                                 Boolean, integer, floating-point, or Enum with
                                 cases = 1/True and controls = 0/False. Can
                                 contain null entries: the corresponding
                                 samples will fail QC.
            custom_filter: an optional Boolean column of obs containing a
                           filter to apply on top of the other QC filters; True
                           elements will be kept. Can be None, a column name, a
                           polars expression, a polars Series, a 1D NumPy
                           array, or a function that takes in this Pseudobulk
                           dataset and a cell type and returns a polars Series
                           or 1D NumPy array. Can also be a sequence of any
                           combination of these for each cell type.
            min_samples: filter to cell types with at least this many cases and
                         this many controls, or with at least this many total
                         samples if `case_control_column` is None
            min_cells: if not None, filter to samples with ≥ this many cells of
                       each cell type
            max_standard_deviations: if not None, filter to samples where the
                                     number of genes with 0 counts is at most
                                     this many standard deviations above the
                                     mean
            min_nonzero_fraction: if not None, filter to genes with at least
                                  one count in this fraction of controls AND
                                  this fraction of cases (or if
                                  `case_control_column` is None, at least one
                                  count in this fraction of samples)
            error_if_negative_counts: if True, raise an error if any counts are
                                      negative
            allow_float: if False, raise an error if `X.dtype` is
                         floating-point (suggesting the user may not be using
                         the raw counts); if True, disable this sanity check
            verbose: whether to print how many samples and genes were filtered
                     out at each step of the QC process
        
        Returns:
            A new Pseudobulk dataset with each cell type's X, obs and var
            subset to samples and genes passing QC.
        """
        # Check inputs
        case_control_column = self._get_columns(
            'obs', case_control_column, 'case_control_column',
            (pl.Boolean, 'integer', 'floating-point', pl.Enum),
            allow_null=True)
        custom_filter = self._get_columns(
            'obs', custom_filter, 'custom_filter', pl.Boolean)
        check_type(min_samples, 'min_samples', int,
                   'an integer greater than or equal to 2')
        check_bounds(min_samples, 'min_samples', 2)
        if min_cells is not None:
            check_type(min_cells, 'min_cells', int, 'a positive integer')
            check_bounds(min_cells, 'min_cells', 1)
        if max_standard_deviations is not None:
            check_type(max_standard_deviations, 'max_standard_deviations',
                       (int, float), 'a positive number')
            check_bounds(max_standard_deviations, 'max_standard_deviations', 0,
                         left_open=True)
        if min_nonzero_fraction is not None:
            check_type(min_nonzero_fraction, 'min_nonzero_fraction',
                       (int, float), 'a number between 0 and 1, inclusive')
            check_bounds(min_nonzero_fraction, 'min_nonzero_fraction', 0, 1)
        check_type(error_if_negative_counts, 'error_if_negative_counts', bool,
                   'Boolean')
        check_type(allow_float, 'allow_float', bool, 'Boolean')
        check_type(verbose, 'verbose', bool, 'Boolean')
        # If `error_if_negative_counts=True`, raise an error if X has any
        # negative values
        if error_if_negative_counts:
            for cell_type in self:
                if self._X[cell_type].ravel().min() < 0:
                    error_message = f'X[{cell_type!r}] has negative counts'
                    raise ValueError(error_message)
        # If `allow_float=False`, raise an error if `X` is floating-point
        if not allow_float:
            for cell_type in self:
                dtype = self._X[cell_type].dtype
                if np.issubdtype(dtype, np.floating):
                    error_message = (
                        f"qc() requires raw counts but X[{cell_type!r}].dtype "
                        f"is {dtype!r}, a floating-point data type; if you "
                        f"are sure that all values are raw integer counts, "
                        f"i.e. that (X[{cell_type!r}].data == "
                        f"X[{cell_type!r}].data.astype(int)).all(), then set "
                        f"allow_float=True (or just cast X to an integer data "
                        f"type).")
                    raise TypeError(error_message)
        if verbose:
            print()
        X_qced, obs_qced, var_qced = {}, {}, {}
        for cell_type, (X, obs, var) in self.items():
            if verbose:
                print(f'[{cell_type}] Starting with {len(obs):,} samples and '
                      f'{len(var):,} genes.')
            # Get the case-control mask
            if case_control_column is not None and \
                    case_control_column[cell_type] is not None:
                case_control_mask = case_control_column[cell_type]
                if case_control_mask.dtype != pl.Boolean:
                    if case_control_mask.dtype == pl.Enum:
                        categories = case_control_mask.cat.get_categories()
                        if len(categories) != 2:
                            error_message = (
                                f'case_control_column is an Enum column '
                                f'with {len(categories):,} categor'
                                f'{"y" if len(categories) == 1 else "ies"} '
                                f'for cell type {cell_type!r}, but must have '
                                f'2 (cases = 1, controls = 0)')
                            raise ValueError(error_message)
                        case_control_mask = case_control_mask.to_physical()
                    else:
                        unique_labels = case_control_mask.unique().drop_nulls()
                        num_unique_labels = len(unique_labels)
                        if num_unique_labels != 2:
                            plural_string = \
                                plural('unique value', num_unique_labels)
                            error_message = (
                                f'case_control_column is a numeric column '
                                f'with {num_unique_labels:,} '
                                f'{plural_string} for cell type '
                                f'{cell_type!r}, but must have 2 '
                                f'(cases = 1, controls = 0)')
                            raise ValueError(error_message)
                        if not unique_labels.sort().equals(pl.Series([0, 1])):
                            error_message = (
                                f'case_control_column is a numeric column '
                                f'with 2 unique values for cell type '
                                f'{cell_type!r}, {unique_labels[0]} and '
                                f'{unique_labels[1]}, but must have '
                                f'cases = 1 and controls = 0')
                            raise ValueError(error_message)
                    case_control_mask = case_control_mask.cast(pl.Boolean)
            else:
                case_control_mask = None
            # Check if we have enough samples for this cell type
            if Pseudobulk._too_few_samples(obs, case_control_mask, min_samples,
                                           cell_type, verbose):
                continue
            # Get the custom filter
            if custom_filter is not None:
                sample_mask = custom_filter[cell_type]
                if sample_mask is not None:
                    if verbose:
                        print(f'[{cell_type}] {sample_mask.sum():,} samples '
                              f'remain after applying the custom filter.')
            else:
                sample_mask = None
            # Get a mask of samples with at least `min_cells` cells of this
            # cell type
            if min_cells is not None:
                if verbose:
                    print(f'[{cell_type}] Filtering to samples with at least '
                          f'{min_cells} {cell_type} cells...')
                if sample_mask is None:
                    sample_mask = obs['num_cells'] >= min_cells
                else:
                    sample_mask &= obs['num_cells'] >= min_cells
                if verbose:
                    print(f'[{cell_type}] {sample_mask.sum():,} samples '
                          f'remain after filtering to samples with at least '
                          f'{min_cells} {cell_type} cells.')
            # Apply the custom filter and/or `min_cells` filters, if either
            # were specified
            if sample_mask is not None:
                # noinspection PyUnresolvedReferences
                X = X[sample_mask.to_numpy()]
                obs = obs.filter(sample_mask)
                if case_control_mask is not None:
                    case_control_mask = case_control_mask.filter(sample_mask)
                # Check if we still have enough samples for this cell type,
                # after applying these two filters
                if Pseudobulk._too_few_samples(obs, case_control_mask,
                                               min_samples, cell_type,
                                               verbose):
                    continue
            # Filter to samples where the number of genes with 0 counts is less
            # than `max_standard_deviations` standard deviations above the mean
            if max_standard_deviations is not None:
                if verbose:
                    print(f'[{cell_type}] Filtering to samples where the '
                          f'number of genes with 0 counts is '
                          f'<{max_standard_deviations} standard deviations '
                          f'above the mean...')
                num_zero_counts = (X == 0).sum(axis=1, dtype=np.int32)
                sample_mask = num_zero_counts < num_zero_counts.mean() + \
                              max_standard_deviations * num_zero_counts.std()
                X = X[sample_mask]
                sample_mask = pl.Series(sample_mask)
                obs = obs.filter(sample_mask)
                if case_control_mask is not None:
                    case_control_mask = case_control_mask.filter(sample_mask)
                if verbose:
                    print(f'[{cell_type}] {len(obs):,} samples remain after '
                          f'filtering to samples where the number of genes '
                          f'with 0 counts is <{max_standard_deviations} '
                          f'standard deviations above the mean.')
                # Check if we have still enough samples for this cell type,
                # after applying this filter
                if Pseudobulk._too_few_samples(obs, case_control_mask,
                                               min_samples, cell_type,
                                               verbose):
                    continue
            # Filter to genes with at least 1 count in
            # `100 * min_nonzero_fraction`% of controls AND
            # `100 * min_nonzero_fraction`% of cases, or if
            # `case_control_column` is None (for this cell type), at least
            # one count in `100 * min_nonzero_fraction`% of samples
            if min_nonzero_fraction is not None:
                if case_control_mask is not None:
                    if verbose:
                        print(f'[{cell_type}] Filtering to genes with at '
                              f'least one count in '
                              f'{100 * min_nonzero_fraction}% of cases and '
                              f'{100 * min_nonzero_fraction}% of controls...')
                    gene_mask = \
                        (np.quantile(X[case_control_mask.fill_null(False)
                                                    .to_numpy()],
                                     1 - min_nonzero_fraction, axis=0) > 0) & \
                        (np.quantile(X[(~case_control_mask).fill_null(False)
                                       .to_numpy()],
                                     1 - min_nonzero_fraction, axis=0) > 0)
                    X = X[:, gene_mask]
                    var = var.filter(gene_mask)
                    if verbose:
                        print(f'[{cell_type}] {len(var):,} genes remain '
                              f'after filtering to genes with at least one '
                              f'count in {100 * min_nonzero_fraction}% of '
                              f'cases and {100 * min_nonzero_fraction}% of '
                              f'controls.')
                else:
                    if verbose:
                        print(f'[{cell_type}] Filtering to genes with at '
                              f'least one count in '
                              f'{100 * min_nonzero_fraction}% of samples...')
                    gene_mask = np.quantile(X, 1 - min_nonzero_fraction,
                                            axis=0) > 0
                    X = X[:, gene_mask]
                    var = var.filter(gene_mask)
                    if verbose:
                        print(f'[{cell_type}] {len(var):,} genes remain '
                              f'after filtering to genes with at least one '
                              f'count in {100 * min_nonzero_fraction}% of '
                              f'samples.')
            X_qced[cell_type] = X
            obs_qced[cell_type] = obs
            var_qced[cell_type] = var
            if verbose:
                print()
        return Pseudobulk(X=X_qced, obs=obs_qced, var=var_qced)
    
    @staticmethod
    def _calc_norm_factors(X: np.ndarray[2, np.integer | np.floating],
                           *,
                           logratio_trim: int | float | np.integer |
                                          np.floating = 0.3,
                           sum_trim: int | float | np.integer |
                                     np.floating = 0.05,
                           do_weighting: bool = True,
                           A_cutoff: int | float | np.integer |
                                     np.floating = -1e10) -> \
            np.ndarray[2, np.floating]:
        """
        A drop-in replacement for edgeR's calcNormFactors with method='TMM'.
        
        Results were verified to match edgeR to within floating-point error.
        
        Does not support the lib.size and refColumn arguments to
        calcNormFactors; these are both assumed to be NULL (the default) and
        will always be calculated internally.
        
        Args:
            X: a matrix of raw (read) counts
            logratio_trim: the amount of trim to use on log-ratios ("M"
                           values); must be between 0 and 1
            sum_trim: the amount of trim to use on the combined absolute levels
                      ("A" values); must be between 0 and 1
            do_weighting: whether to compute (asymptotic binomial precision)
                          weights
            A_cutoff: the cutoff on "A" values to use before trimming
        
        Returns:
            A 1D NumPy array with the norm factors for each column of X.
        """
        # Check inputs
        check_type(logratio_trim, 'logratio_trim', float,
                   'a floating-point number')
        check_bounds(logratio_trim, 'logratio_trim', 0, 1, left_open=True,
                     right_open=True)
        check_type(sum_trim, 'sum_trim', float, 'a floating-point number')
        check_bounds(sum_trim, 'sum_trim', 0, 1, left_open=True,
                     right_open=True)
        check_type(do_weighting, 'do_weighting', bool, 'Boolean')
        check_type(A_cutoff, 'A_cutoff', float, 'a floating-point number')
        
        # Degenerate cases
        if X.shape[0] == 0 or X.shape[1] == 1:
            return np.ones(X.shape[1])
        
        # Remove all-zero rows
        any_non_zero = (X != 0).any(axis=1)
        if not any_non_zero.all():
            X = X[any_non_zero]
        
        # Calculate library sizes
        lib_size = X.sum(axis=0)
        
        # Determine which column is the reference column
        f75 = np.quantile(X, 0.75, axis=0) / lib_size
        if f75.min() == 0:
            import warnings
            warning_message = 'one or more quantiles are zero'
            warnings.warn(warning_message)
        ref_column = np.argmax(np.sqrt(X).sum(axis=0)) \
            if np.median(f75) < 1e-20 else \
            np.argmin(np.abs(f75 - f75.mean()))
        
        with np.errstate(divide='ignore', invalid='ignore'):
            # Calculate the log ratio of expression accounting for library size
            normalized_X = X / lib_size
            logR_ = np.log2(normalized_X / normalized_X[:, [ref_column]])
            
            # Calculate absolute expression
            log_normalized_X = np.log2(normalized_X)
            absE_ = 0.5 * (log_normalized_X +
                           log_normalized_X[:, [ref_column]])
            
            # Calculate estimated asymptotic variance
            if do_weighting:
                sum_of_reciprocals = 1 / X + 1 / lib_size
                v_ = sum_of_reciprocals + sum_of_reciprocals[:, [ref_column]]
        
        # Remove infinite values, cutoff based on A
        finite_ = (logR_ != -np.inf) & (absE_ > A_cutoff)
        
        # Calculate the normalization factors
        factors = np.empty(X.shape[1])
        for i in range(X.shape[1]):
            finite = finite_[:, i]
            logR = logR_[finite, i]
            absE = absE_[finite, i]
            if np.abs(logR).max() < 1e-6:
                factors[i] = 1
                continue
            n = len(logR)
            loL = int(n * logratio_trim)
            hiL = n - loL
            loS = int(n * sum_trim)
            hiS = n - loS
            logR_rank = rankdata(logR)
            absE_rank = rankdata(absE)
            keep = (logR_rank >= loL + 1) & (logR_rank <= hiL) & \
                   (absE_rank >= loS + 1) & (absE_rank <= hiS)
            factors[i] = 2 ** np.average(logR[keep],
                                         weights=1 / v_[finite, i][keep]
                                                 if do_weighting else None)
        
        # Results will be missing if the two libraries share no features with
        # positive counts; in this case, set to unity
        np.nan_to_num(factors, copy=False, nan=1)
        
        # Factors should multiply to one
        factors /= np.exp(np.log(factors).mean())
        return factors
    
    def CPM(self) -> Pseudobulk:
        """
        Calculate counts per million for each cell type.

        Returns:
            A new Pseudobulk dataset containing the CPMs.
        """
        CPMs = {}
        for cell_type in self:
            X = self._X[cell_type]
            library_size = X.sum(axis=1) * self._calc_norm_factors(X.T)
            CPMs[cell_type] = X / library_size[:, None] * 1e6
        return Pseudobulk(X=CPMs, obs=self._obs, var=self._var)
    
    def log_CPM(self,
                *,
                prior_count: int | float | np.integer |
                             np.floating = 2) -> Pseudobulk:
        """
        Calculate log counts per million for each cell type.
        
        Do NOT run this before DE(), since DE() already runs it internally.
        
        Based on the R translation of edgeR's C++ cpm() code at
        bioinformatics.stackexchange.com/a/4990.
        
        Results were verified to match edgeR to within floating-point error.
        
        Args:
            prior_count: the pseudocount to add before log-transforming. In the
                         current version of edgeR, prior.count is now 2 instead
                         of the old value of 0.5: code.bioconductor.org/browse/
                         edgeR/blob/RELEASE_3_18/R/cpm.R
        
        Returns:
            A new Pseudobulk dataset containing the log(CPMs).
        """
        check_type(prior_count, 'prior_count', (int, float),
                   'a positive number')
        check_bounds(prior_count, 'prior_count', 0, left_open=True)
        log_CPMs = {}
        for cell_type in self:
            X = self._X[cell_type]
            library_size = X.sum(axis=1) * self._calc_norm_factors(X.T)
            pseudocount = prior_count * library_size / library_size.mean()
            library_size += 2 * pseudocount
            log_CPMs[cell_type] = np.log2(X + pseudocount[:, None]) - \
                np.log2(library_size[:, None]) + np.log2(1e6)
        return Pseudobulk(X=log_CPMs, obs=self._obs, var=self._var)

    def regress_out_obs(self,
                        covariate_columns: Sequence[
                            PseudobulkColumn |
                            Sequence[PseudobulkColumn | None]],
                        *,
                        error_if_int: bool = True) -> Pseudobulk:
        """
        Regress out covariates from obs. Must be run after log_CPM().

        Args:
            covariate_columns: a sequence of columns of obs to regress out.
                               Each element of the sequence can be a column
                               name, a polars expression, a polars Series, a
                               1D NumPy array, or a function that takes in this
                               Pseudobulk dataset and a cell type and returns a
                               polars Series or 1D NumPy array. Each element of
                               the sequence may also itself be a sequence of
                               any combination of these for each cell type, or
                               None to not include that covariate for that cell
                               type.
            error_if_int: if True, raise an error if `X.dtype` is integer
                          (indicating the user may not have run log_CPM() yet)

        Returns:
            A new Pseudobulk dataset with covariates regressed out.
        """
        # Check inputs
        if covariate_columns is None:
            error_message = 'covariate_columns is None'
            raise TypeError(error_message)
        for index, column in enumerate(covariate_columns):
            if column is None:
                error_message = f'covariate_columns[{index}] is None'
                raise TypeError(error_message)
        covariate_columns = [
            self._get_columns('obs', column, f'covariate_columns[{index}]',
                              ('integer', 'floating-point', pl.Categorical,
                               pl.Enum))
            for index, column in enumerate(covariate_columns)]
        check_type(error_if_int, 'error_if_int', bool, 'Boolean')
        # For each cell type...
        residuals = {}
        for cell_type, (X, obs, var) in self.items():
            # If error_if_int=True, raise an error if X has an integer dtype
            if error_if_int and np.issubdtype(X.dtype, np.integer):
                error_message = (
                    f'X[{cell_type!r}].dtype is {str(X.dtype)!r}, an integer '
                    f'data type; did you forget to run log_CPM() before '
                    f'regress_out()?')
                raise ValueError(error_message)
            # Get the covariates for this cell type
            covariates = pl.DataFrame([column[cell_type]
                                       for column in covariate_columns])
            # Convert the covariates to a Numpy array, and add an intercept
            covariates = covariates.to_numpy()
            if not np.issubdtype(covariates.dtype, np.number):
                error_message = (
                    f'obs[{cell_type!r}].select(covariate_columns) must be '
                    f'convertible to a numeric NumPy array, but converted to '
                    f'an array of data type {str(covariates.dtype)!r}')
                raise TypeError(error_message)
            covariates = np.column_stack(
                (np.ones(len(covariates), covariates.dtype), covariates))
            # Regress out the covariates; silence warnings with rcond=None
            beta = np.linalg.lstsq(covariates, X, rcond=None)[0]
            # Calculate the residuals
            residuals[cell_type] = X - covariates @ beta
        # Return a new Pseudobulk datasets with the residuals
        return Pseudobulk(X=residuals, obs=self._obs, var=self._var)
    
    def regress_out_var(self,
                        covariate_columns: Sequence[
                            PseudobulkColumn |
                            Sequence[PseudobulkColumn | None]],
                        *,
                        error_if_int: bool = True) -> Pseudobulk:
        """
        Regress out covariates from var. Must be run after log_CPM().

        Args:
            covariate_columns: a sequence of columns of var to regress out.
                               Each element of the sequence can be a column
                               name, a polars expression, a polars Series, a
                               1D NumPy array, or a function that takes in this
                               Pseudobulk dataset and a cell type and returns a
                               polars Series or 1D NumPy array. Each element of
                               the sequence may also itself be a sequence of
                               any combination of these for each cell type, or
                               None to not include that covariate for that cell
                               type.
            error_if_int: if True, raise an error if `X.dtype` is integer
                          (indicating the user may not have run log_CPM() yet)

        Returns:
            A new Pseudobulk dataset with covariates regressed out.
        """
        # Check inputs
        if covariate_columns is None:
            error_message = 'covariate_columns is None'
            raise TypeError(error_message)
        for index, column in enumerate(covariate_columns):
            if column is None:
                error_message = f'covariate_columns[{index}] is None'
                raise TypeError(error_message)
        covariate_columns = [
            self._get_columns('var', column, f'covariate_columns[{index}]',
                              ('integer', 'floating-point', pl.Categorical,
                               pl.Enum))
            for index, column in enumerate(covariate_columns)]
        check_type(error_if_int, 'error_if_int', bool, 'Boolean')
        # For each cell type...
        residuals = {}
        for cell_type, (X, obs, var) in self.items():
            # If error_if_int=True, raise an error if X has an integer dtype
            if error_if_int and np.issubdtype(X.dtype, np.integer):
                error_message = (
                    f'X[{cell_type!r}].dtype is {str(X.dtype)!r}, an integer '
                    f'data type; did you forget to run log_CPM() before '
                    f'regress_out()?')
                raise ValueError(error_message)
            # Get the covariates for this cell type
            covariates = pl.DataFrame([column[cell_type]
                                       for column in covariate_columns])
            # Convert the covariates to a Numpy array, and add an intercept
            covariates = covariates.to_numpy()
            if not np.issubdtype(covariates.dtype, np.number):
                error_message = (
                    f'var[{cell_type!r}].select(covariate_columns) must be '
                    f'convertible to a numeric NumPy array, but converted to '
                    f'an array of data type {str(covariates.dtype)!r}')
                raise TypeError(error_message)
            covariates = np.column_stack(
                (np.ones(len(covariates), covariates.dtype), covariates))
            # Regress out the covariates; silence warnings with rcond=None
            beta = np.linalg.lstsq(covariates, X.T, rcond=None)[0]
            # Calculate the residuals
            residuals[cell_type] = (X.T - covariates @ beta).T
        # Return a new Pseudobulk datasets with the residuals
        return Pseudobulk(X=residuals, obs=self._obs, var=self._var)
    
    # A slightly reformatted version of the voomByGroup source code from
    # github.com/YOU-k/voomByGroup/blob/main/voomByGroup.R, which is available
    # under the MIT license. Copyright (c) 2023 Yue You.
    _voomByGroup_source_code = r'''
    voomByGroup <- function (counts, group = NULL, design = NULL,
                             lib.size = NULL, dynamic = NULL,
                             normalize.method = "none", span = 0.5,
                             save.plot = FALSE, print = TRUE, plot = c("none",
                             "all", "separate", "combine"),
                             col.lines = NULL, pos.legend = c("inside",
                             "outside", "none"), fix.y.axis = FALSE, ...) {
      out <- list()
      if (is(counts, "DGEList")) {
        out$genes <- counts$genes
        out$targets <- counts$samples
        if(is.null(group))
          group <- counts$samples$group
        if (is.null(lib.size))
          lib.size <- with(counts$samples, lib.size * norm.factors)
        counts <- counts$counts
      }
      else {
        isExpressionSet <-
          suppressPackageStartupMessages(is(counts, "ExpressionSet"))
        if (isExpressionSet) {
          if (length(Biobase::fData(counts)))
            out$genes <- Biobase::fData(counts)
          if (length(Biobase::pData(counts)))
            out$targets <- Biobase::pData(counts)
          counts <- Biobase::exprs(counts)
        }
        else {
          counts <- as.matrix(counts)
        }
      }
      if (nrow(counts) < 2L)
        stop("Need at least two genes to fit a mean-variance trend")
      # Library size
      if(is.null(lib.size))
        lib.size <- colSums(counts)
      # Group
      if(is.null(group))
        group <- rep("Group1", ncol(counts))
      group <- as.factor(group)
      intgroup <- as.integer(group)
      levgroup <- levels(group)
      ngroups <- length(levgroup)
      # Design matrix
      if (is.null(design)) {
        design <- matrix(1L, ncol(counts), 1)
        rownames(design) <- colnames(counts)
        colnames(design) <- "GrandMean"
      }
      # Dynamic
      if (is.null(dynamic)) {
        dynamic <- rep(FALSE, ngroups)
      }
      # voom by group
      if(print)
        cat("Group:\n")
      E <- w <- counts
      xy <- line <- as.list(rep(NA, ngroups))
      names(xy) <- names(line) <- levgroup
      for (lev in 1L:ngroups) {
        if(print)
          cat(lev, levgroup[lev], "\n")
        i <- intgroup == lev
        countsi <- counts[, i]
        libsizei <- lib.size[i]
        designi <- design[i, , drop = FALSE]
        QR <- qr(designi)
        if(QR$rank<ncol(designi))
          designi <- designi[,QR$pivot[1L:QR$rank], drop = FALSE]
        if(ncol(designi)==ncol(countsi))
          designi <- matrix(1L, ncol(countsi), 1)
        voomi <- voom(counts = countsi, design = designi, lib.size = libsizei,
                      normalize.method = normalize.method, span = span,
                      plot = FALSE, save.plot = TRUE, ...)
        E[, i] <- voomi$E
        w[, i] <- voomi$weights
        xy[[lev]] <- voomi$voom.xy
        line[[lev]] <- voomi$voom.line
      }
      #voom overall
      if (TRUE %in% dynamic){
        voom_all <- voom(counts = counts, design = design, lib.size = lib.size,
                         normalize.method = normalize.method, span = span,
                         plot = FALSE, save.plot = TRUE, ...)
        E_all <- voom_all$E
        w_all <- voom_all$weights
        xy_all <- voom_all$voom.xy
        line_all <- voom_all$voom.line
        dge <- DGEList(counts)
        disp <- estimateCommonDisp(dge)
        disp_all <- disp$common
      }
      # Plot, can be "both", "none", "separate", or "combine"
      plot <- plot[1]
      if(plot!="none"){
        disp.group <- c()
        for (lev in levgroup) {
          dge.sub <- DGEList(counts[,group == lev])
          disp <- estimateCommonDisp(dge.sub)
          disp.group[lev] <- disp$common
        }
        if(plot %in% c("all", "separate")){
          if (fix.y.axis == TRUE) {
            yrange <- sapply(levgroup, function(lev){
              c(min(xy[[lev]]$y), max(xy[[lev]]$y))
            }, simplify = TRUE)
            yrange <- c(min(yrange[1,]) - 0.1, max(yrange[2,]) + 0.1)
          }
          for (lev in 1L:ngroups) {
            if (fix.y.axis == TRUE){
              plot(xy[[lev]], xlab = "log2( count size + 0.5 )",
                   ylab = "Sqrt( standard deviation )", pch = 16, cex = 0.25,
                   ylim = yrange)
            } else {
              plot(xy[[lev]], xlab = "log2( count size + 0.5 )",
                   ylab = "Sqrt( standard deviation )", pch = 16, cex = 0.25)
            }
            title(paste("voom: Mean-variance trend,", levgroup[lev]))
            lines(line[[lev]], col = "red")
            legend("topleft", bty="n", paste("BCV:",
              round(sqrt(disp.group[lev]), 3)), text.col="red")
          }
        }
        
        if(plot %in% c("all", "combine")){
          if(is.null(col.lines))
            col.lines <- 1L:ngroups
          if(length(col.lines)<ngroups)
            col.lines <- rep(col.lines, ngroups)
          xrange <- unlist(lapply(line, `[[`, "x"))
          xrange <- c(min(xrange)-0.3, max(xrange)+0.3)
          yrange <- unlist(lapply(line, `[[`, "y"))
          yrange <- c(min(yrange)-0.1, max(yrange)+0.3)
          plot(1L,1L, type="n", ylim=yrange, xlim=xrange,
               xlab = "log2( count size + 0.5 )",
               ylab = "Sqrt( standard deviation )")
          title("voom: Mean-variance trend")
          if (TRUE %in% dynamic){
            for (dy in which(dynamic)){
              line[[dy]] <- line_all
              disp.group[dy] <- disp_all
              levgroup[dy] <- paste0(levgroup[dy]," (all)")
            }
          }
          for (lev in 1L:ngroups)
            lines(line[[lev]], col=col.lines[lev], lwd=2)
          pos.legend <- pos.legend[1]
          disp.order <- order(disp.group, decreasing = TRUE)
          text.legend <-
            paste(levgroup, ", BCV: ", round(sqrt(disp.group), 3), sep="")
          if(pos.legend %in% c("inside", "outside")){
            if(pos.legend=="outside"){
              plot(1,1, type="n", yaxt="n", xaxt="n", ylab="", xlab="",
                   frame.plot=FALSE)
              legend("topleft", text.col=col.lines[disp.order],
                     text.legend[disp.order], bty="n")
            } else {
              legend("topright", text.col=col.lines[disp.order],
                     text.legend[disp.order], bty="n")
            }
          }
        }
      }
      # Output
      if (TRUE %in% dynamic){
        E[,intgroup %in% which(dynamic)] <-
          E_all[,intgroup %in% which(dynamic)]
        w[,intgroup %in% which(dynamic)] <-
          w_all[,intgroup %in% which(dynamic)]
      }
      out$E <- E
      out$weights <- w
      out$design <- design
      if(save.plot){
        out$voom.line <- line
        out$voom.xy <- xy
      }
      new("EList", out)
    }
    '''
    
    def DE(self,
           label_column: PseudobulkColumn | Sequence[PseudobulkColumn],
           covariate_columns: Sequence[PseudobulkColumn | None |
                                       Sequence[PseudobulkColumn | None]] |
                              None,
           *,
           case_control: bool = True,
           cell_types: str | Iterable[str] | None = None,
           excluded_cell_types: str | Iterable[str] | None = None,
           library_size_as_covariate: bool = True,
           num_cells_as_covariate: bool = True,
           return_voom_info: bool = True,
           allow_float: bool = False,
           verbose: bool = True) -> DE:
        """
        Perform differential expression (DE) on a Pseudobulk dataset with
        limma-voom. Uses voomByGroup when case_control=True, which is better
        than regular voom for case-control DE.
        
        Loosely based on the `de_pseudobulk()` function from
        github.com/tluquez/utils/blob/main/utils.R, which is itself based on
        github.com/neurorestore/Libra/blob/main/R/pseudobulk_de.R.

        Args:
            label_column: the column of obs to calculate DE with respect to.
                          Can be a column name, a polars expression, a polars
                          Series, a 1D NumPy array, or a function that takes in
                          this Pseudobulk dataset and a cell type and returns a
                          polars Series or 1D NumPy array. Can also be a
                          sequence of any combination of these for each cell
                          type. If `case_control=True`, must be Boolean,
                          integer, floating-point, or Enum with cases = 1/True
                          and controls = 0/False. If `case_control=False`, must
                          be integer or floating-point.
            covariate_columns: an optional sequence of columns of obs to use
                               as covariates, or None to not include
                               covariates. Each element of the sequence can be
                               a column name, a polars expression, a polars
                               Series, a 1D NumPy array, or a function that
                               takes in this Pseudobulk dataset and a cell type
                               and returns a polars Series or 1D NumPy array.
                               Each element of the sequence may also itself be
                               a sequence of any combination of these for each
                               cell type, or None to not include that covariate
                               for that cell type.
            case_control: whether the analysis is case-control or with respect
                          to a quantitative variable.
                          If True, uses voomByGroup instead of regular voom,
                          and uses `label_column` as the `group` argument to
                          calcNormFactors().
            cell_types: one or more cell types to test for differential
                        expression; if None, test all cell types. Mutually
                        exclusive with `excluded_cell_types`.
            excluded_cell_types: cell types to exclude when testing for
                                 differential expression; mutually exclusive
                                 with `cell_types`
            library_size_as_covariate: whether to include the log2 of the
                                       library size, calculated according to
                                       the method of edgeR's calcNormFactors(),
                                       as anadditional covariate
            num_cells_as_covariate: whether to include the log2 of the
                                    `'num_cells'` column of obs, i.e. the
                                    number of cells that went into each
                                    sample's pseudobulk in each cell type, as
                                    an additional covariate
            return_voom_info: whether to include the voom weights and voom plot
                              data in the returned DE object; set to False for
                              reduced runtime if you do not need to use the
                              voom weights or generate voom plots
            allow_float: if False, raise an error if `X.dtype` is
                         floating-point (suggesting the user may not be using
                         the raw counts, e.g. due to accidentally having run
                         log_CPM() already); if True, disable this sanity check
            verbose: whether to print out details of the DE estimation

        Returns:
            A DE object with a `table` attribute containing a polars DataFrame
            of the DE results. If `return_voom_info=True`, also includes a
            `voom_weights` attribute containing a {cell_type: DataFrame}
            dictionary of voom weights, and a `voom_plot_data` attribute
            containing a {cell_type: DataFrame} dictionary of info necessary to
            construct a voom plot with `DE.plot_voom()`.
        """
        # Import required Python and R packages
        from ryp import r, to_py, to_r
        r('suppressPackageStartupMessages(library(edgeR))')
        # Source voomByGroup code
        if case_control:
            r(self._voomByGroup_source_code)
        # Check inputs
        label_column = self._get_columns(
            'obs', label_column, 'label_column',
            (pl.Boolean, 'integer', 'floating-point', pl.Enum)
            if case_control else ('integer', 'floating-point'),
            allow_None=False)
        if covariate_columns is not None:
            covariate_columns = [
                self._get_columns('obs', column, f'covariate_columns[{index}]',
                                  ('integer', 'floating-point', pl.Categorical,
                                   pl.Enum))
                for index, column in enumerate(covariate_columns)]
        check_type(case_control, 'case_control', bool, 'Boolean')
        if cell_types is not None:
            if excluded_cell_types is not None:
                error_message = (
                    'cell_types and excluded_cell_types cannot both be '
                    'specified')
                raise ValueError(error_message)
            is_string = isinstance(cell_types, str)
            cell_types = to_tuple(cell_types)
            if len(cell_types) == 0:
                error_message = 'cell_types is empty'
                raise ValueError(error_message)
            check_types(cell_types, 'cell_types', str, 'strings')
            for cell_type in cell_types:
                if cell_type not in self._X:
                    if is_string:
                        error_message = (
                            f'cell_types is {cell_type!r}, which is not a '
                            f'cell type in this Pseudobulk dataset')
                        raise ValueError(error_message)
                    else:
                        error_message = (
                            f'one of the elements of cell_types, '
                            f'{cell_type!r}, is not a cell type in this '
                            f'Pseudobulk dataset')
                        raise ValueError(error_message)
        elif excluded_cell_types is not None:
            excluded_cell_types = to_tuple(excluded_cell_types)
            check_types(excluded_cell_types, 'cell_types', str, 'strings')
            for cell_type in excluded_cell_types:
                if cell_type not in self._X:
                    if excluded_cell_types:
                        error_message = (
                            f'excluded_cell_types is {cell_type!r}, which is '
                            f'not a cell type in this Pseudobulk dataset')
                        raise ValueError(error_message)
                    else:
                        error_message = (
                            f'one of the elements of excluded_cell_types, '
                            f'{cell_type!r}, is not a cell type in this '
                            f'Pseudobulk dataset')
                        raise ValueError(error_message)
            cell_types = [cell_type for cell_type in self._X
                          if cell_type not in excluded_cell_types]
            if len(cell_types) == 0:
                error_message = \
                    'all cell types were excluded by excluded_cell_types'
                raise ValueError(error_message)
        else:
            cell_types = self._X
        check_type(library_size_as_covariate, 'library_size_as_covariate',
                   bool, 'Boolean')
        check_type(num_cells_as_covariate, 'num_cells_as_covariate', bool,
                   'Boolean')
        if num_cells_as_covariate:
            for cell_type in self:
                if 'num_cells' not in self._obs[cell_type]:
                    error_message = (
                        f"num_cells_as_covariate is True, but 'num_cells' is "
                        f"not a column of obs[{cell_type!r}]")
                    raise ValueError(error_message)
        check_type(allow_float, 'allow_float', bool, 'Boolean')
        check_type(verbose, 'verbose', bool, 'Boolean')
        # Compute DE for each cell type
        DE_results = {}
        if return_voom_info:
            voom_weights = {}
            voom_plot_data = {}
        for cell_type in cell_types:
            X = self._X[cell_type]
            obs = self._obs[cell_type]
            var = self._var[cell_type]
            with Timer(f'[{cell_type}] Calculating DE', verbose=verbose):
                # If `allow_float=False`, raise an error if `X` is
                # floating-point
                if not allow_float and np.issubdtype(X.dtype, np.floating):
                    error_message = (
                        f"DE() requires raw counts but X[{cell_type!r}].dtype "
                        f"is {str(X.dtype)!r}, a floating-point data type. If "
                        f"you are sure that all values are integers, i.e. "
                        f"that X[{cell_type!r}].data == X[{cell_type!r}].data"
                        f".astype(int)).all(), then set allow_float=True (or "
                        f"just cast X to an integer data type). "
                        f"Alternatively, did you accidentally run log_CPM() "
                        f"before DE()?")
                    raise TypeError(error_message)
                # Get the DE labels and covariates for this cell type
                DE_labels = label_column[cell_type]
                covariates = pl.DataFrame([column[cell_type]
                                           for column in covariate_columns]) \
                    if covariate_columns is not None else pl.DataFrame()
                if num_cells_as_covariate:
                    covariates = \
                        covariates.with_columns(obs['num_cells'].log(2))
                if case_control:
                    # If `case_control=True`, check that the DE labels have
                    # only two unique values
                    if DE_labels.dtype != pl.Boolean:
                        if DE_labels.dtype == pl.Enum:
                            categories = DE_labels.cat.get_categories()
                            if len(categories) != 2:
                                suffix = 'y' if len(categories) == 1 else 'ies'
                                error_message = (
                                    f'[{cell_type}] label_column is an Enum '
                                    f'column with {len(categories):,} '
                                    f'categor{suffix}, but must have 2 (cases '
                                    f'and controls)')
                                raise ValueError(error_message)
                            DE_labels = DE_labels.to_physical()
                        else:
                            unique_labels = DE_labels.unique()
                            num_unique_labels = len(unique_labels)
                            if num_unique_labels != 2:
                                plural_string = \
                                    plural('unique value', num_unique_labels)
                                error_message = (
                                    f'[{cell_type}] label_column is a numeric '
                                    f'column with {num_unique_labels:,} '
                                    f'{plural_string}, but must have 2 '
                                    f'(cases = 1, controls = 0) unless '
                                    f'case_control=False')
                                raise ValueError(error_message)
                            if not unique_labels.sort()\
                                    .equals(pl.Series([0, 1])):
                                error_message = (
                                    f'[{cell_type}] label_column is a numeric '
                                    f'column with 2 unique values, '
                                    f'{unique_labels[0]} and '
                                    f'{unique_labels[1]}, but must have '
                                    f'cases = 1 and controls = 0')
                                raise ValueError(error_message)
                # Get the design matrix
                if verbose:
                    print(f'[{cell_type}] Generating design matrix...')
                design_matrix = \
                    obs.select(pl.lit(1).alias('intercept'), DE_labels)
                if covariates.width:
                    design_matrix = pl.concat([
                        design_matrix,
                        covariates.to_dummies(covariates.select(
                            pl.col(pl.Categorical, pl.Enum)).columns,
                            drop_first=True)],
                        how='horizontal')
                # Estimate library sizes
                if verbose:
                    print(f'[{cell_type}] Estimating library sizes...')
                library_size = X.sum(axis=1) * self._calc_norm_factors(X.T)
                if library_size.min() == 0:
                    error_message = f'[{cell_type}] some library sizes are 0'
                    raise ValueError(error_message)
                if library_size_as_covariate:
                    design_matrix = design_matrix\
                        .with_columns(library_size=np.log2(library_size))
                # Check that the design matrix has more rows than columns, and
                # is full-rank
                if verbose:
                    print(f'[{cell_type}] Sanity-checking the design matrix')
                if design_matrix.width >= design_matrix.height:
                    error_message = (
                        f'[{cell_type}] the design matrix must have more rows '
                        f'(samples) than columns (one plus the number of '
                        f'covariates), but has {design_matrix.height:,} rows '
                        f'and {design_matrix.width:,} columns; either reduce '
                        f'the number of covariates or exclude this cell type '
                        f'with e.g. excluded_cell_types={cell_type!r}')
                    raise ValueError(error_message)
                if np.linalg.matrix_rank(design_matrix.to_numpy()) < \
                        design_matrix.width:
                    error_message = (
                        f'[{cell_type}] the design matrix is not full-rank; '
                        f'some of your covariates are linear combinations of '
                        f'other covariates')
                    raise ValueError(error_message)
                try:
                    # Convert the expression matrix, design matrix, and library
                    # sizes to R
                    if verbose:
                        print(f'[{cell_type}] Converting the expression '
                              f'matrix, design matrix and library sizes to '
                              f'R...')
                    to_r(X.T, '.Pseudobulk.X.T', rownames=var[:, 0],
                         colnames=obs['ID'])
                    to_r(design_matrix, '.Pseudobulk.design.matrix',
                         rownames=obs['ID'])
                    to_r(library_size, '.Pseudobulk.library.size',
                         rownames=obs['ID'])
                    # Run voom
                    to_r(return_voom_info, 'save.plot')
                    if case_control:
                        if verbose:
                            print(f'[{cell_type}] Running voomByGroup...')
                        to_r(DE_labels, '.Pseudobulk.DE.labels',
                             rownames=obs['ID'])
                        r('.Pseudobulk.voom.result = voomByGroup('
                          '.Pseudobulk.X.T, .Pseudobulk.DE.labels, '
                          '.Pseudobulk.design.matrix, '
                          '.Pseudobulk.library.size, save.plot=save.plot, '
                          'print=FALSE)')
                    else:
                        if verbose:
                            print(f'[{cell_type}] Running voom...')
                        r('.Pseudobulk.voom.result = voom(.Pseudobulk.X.T, '
                          '.Pseudobulk.design.matrix, '
                          '.Pseudobulk.library.size, save.plot=save.plot)')
                    if return_voom_info:
                        # noinspection PyUnboundLocalVariable
                        voom_weights[cell_type] = \
                            to_py('.Pseudobulk.voom.result$weights',
                                  index='gene')
                        # noinspection PyUnboundLocalVariable
                        voom_plot_data[cell_type] = pl.DataFrame({
                            f'{prop}_{dim}_{case}': to_py(
                                f'.Pseudobulk.voom.result$voom.{prop}$'
                                f'`{case_label}`${dim}', format='numpy')
                            for prop in ('xy', 'line') for dim in ('x', 'y')
                            for case, case_label in zip(
                                (False, True), ('FALSE', 'TRUE')
                                if DE_labels.dtype == pl.Boolean else (0, 1))}
                            if case_control else {
                            f'{prop}_{dim}': to_py(
                                f'.Pseudobulk.voom.result$voom.{prop}${dim}',
                                format='numpy')
                            for prop in ('xy', 'line') for dim in ('x', 'y')})
                    # Run limma
                    if verbose:
                        print(f'[{cell_type}] Running lmFit...')
                    r('.Pseudobulk.lmFit.result = lmFit('
                      '.Pseudobulk.voom.result, .Pseudobulk.design.matrix)')
                    if verbose:
                        print(f'[{cell_type}] Running eBayes...')
                    r('.Pseudobulk.eBayes.result = eBayes('
                      '.Pseudobulk.lmFit.result, trend=FALSE, robust=FALSE)')
                    # Get results table
                    if verbose:
                        print(f'[{cell_type}] Running topTable...')
                    to_r(DE_labels.name, '.Pseudobulk.coef')
                    r('.Pseudobulk.topTable.result = topTable('
                      '.Pseudobulk.eBayes.result, coef=.Pseudobulk.coef, '
                      'number=Inf, adjust.method="none", sort.by="P", '
                      'confint=TRUE)')
                    if verbose:
                        print(f'[{cell_type}] Collating results...')
                    DE_results[cell_type] = \
                        to_py('.Pseudobulk.topTable.result', index='gene')\
                        .select('gene',
                                logFC=pl.col.logFC,
                                SE=to_py('.Pseudobulk.eBayes.result$s2.post')
                                   .sqrt() *
                                   to_py('.Pseudobulk.eBayes.result$stdev.'
                                         'unscaled[,1]', index=False),
                                LCI=pl.col('CI.L'),
                                UCI=pl.col('CI.R'),
                                AveExpr=pl.col.AveExpr,
                                P=pl.col('P.Value'),
                                Bonferroni=bonferroni(pl.col('P.Value')),
                                FDR=fdr(pl.col('P.Value')))
                finally:
                    r('rm(list = Filter(exists, c(".Pseudobulk.X.T", '
                      '".Pseudobulk.DE.labels", ".Pseudobulk.design.matrix", '
                      '".Pseudobulk.library.size", ".Pseudobulk.voom.result", '
                      '".Pseudobulk.lmFit.result", '
                      '".Pseudobulk.eBayes.result", ".Pseudobulk.coef", '
                      '".Pseudobulk.topTable.result")))')
        # Concatenate across cell types
        table = pl.concat([
            cell_type_DE_results
            .select(pl.lit(cell_type).alias('cell_type'), pl.all())
            for cell_type, cell_type_DE_results in DE_results.items()])
        if return_voom_info:
            return DE(table, case_control, voom_weights, voom_plot_data)
        else:
            return DE(table, case_control)


class DE:
    """
    Differential expression results returned by Pseudobulk.DE().
    """
    
    def __init__(self,
                 table: pl.DataFrame,
                 case_control: bool | None = None,
                 voom_weights: dict[str, pl.DataFrame] | None = None,
                 voom_plot_data: dict[str, pl.DataFrame] | None = None) -> \
            None:
        """
        Initialize the DE object.
        
        Args:
            table: a polars DataFrame containing the DE results, with columns:
                   - cell_type: the cell type in which DE was tested
                   - gene: the gene for which DE was tested
                   - logFC: the log2 fold change of the gene, i.e. its effect
                            size
                   - SE: the standard error of the effect size
                   - LCI: the lower 95% confidence interval of the effect size
                   - UCI: the upper 95% confidence interval of the effect size
                   - AveExpr: the gene's average expression in this cell type,
                              in log CPM
                   - P: the DE p-value
                   - Bonferroni: the Bonferroni-corrected DE p-value
                   - FDR: the FDR q-value for the DE
                   Or, a directory containing a DE object saved with `save()`.
            case_control: whether the analysis is case-control or with respect
                          to a quantitative variable. Must be specified unless
                          `table` is a directory.
            voom_weights: an optional {cell_type: DataFrame} dictionary of voom
                         weights, where rows are genes and columns are samples.
                         The first column of each cell type's DataFrame,
                         'gene', contains the gene names.
            voom_plot_data: an optional {cell_type: DataFrame} dictionary of
                            info necessary to construct a voom plot with
                            `DE.plot_voom()`
        """
        if isinstance(table, pl.DataFrame):
            check_type(case_control, 'case_control', bool, 'Boolean')
            if voom_weights is not None:
                if voom_plot_data is None:
                    error_message = (
                        'voom_plot_data must be specified when voom_weights '
                        'is specified')
                    raise ValueError(error_message)
                check_type(voom_weights, 'voom_weights', dict, 'a dictionary')
                if voom_weights.keys() != voom_plot_data.keys():
                    error_message = (
                        'voom_weights and voom_plot_data must have matching '
                        'keys (cell types)')
                    raise ValueError(error_message)
                for key in voom_weights:
                    if not isinstance(key, str):
                        error_message = (
                            f'all keys of voom_weights and voom_plot_data '
                            f'must be strings (cell types), but they contain '
                            f'a key of type {type(key).__name__!r}')
                        raise TypeError(error_message)
            if voom_plot_data is not None:
                if voom_weights is None:
                    error_message = (
                        'voom_weights must be specified when voom_plot_data '
                        'is specified')
                    raise ValueError(error_message)
                check_type(voom_plot_data, 'voom_plot_data', dict,
                           'a dictionary')
        elif isinstance(table, (str, Path)):
            table = str(table)
            if not os.path.exists(table):
                error_message = f'DE object directory {table!r} does not exist'
                raise FileNotFoundError(error_message)
            cell_types = [line.rstrip('\n') for line in
                          open(f'{table}/cell_types.txt')]
            voom_weights = {cell_type: pl.read_parquet(
                os.path.join(table, f'{cell_type.replace("/", "-")}.'
                                    f'voom_weights.parquet'))
                for cell_type in cell_types}
            voom_plot_data = {cell_type: pl.read_parquet(
                os.path.join(table, f'{cell_type.replace("/", "-")}.'
                                    f'voom_plot_data.parquet'))
                for cell_type in cell_types}
            # noinspection PyUnresolvedReferences
            case_control = \
                next(iter(voom_plot_data.values())).columns[0].count('_') == 2
            table = pl.read_parquet(os.path.join(table, 'table.parquet'))
        else:
            error_message = (
                f'table must be a polars DataFrame or a directory (string or '
                f'pathlib.Path) containing a saved DE object, but has type '
                f'{type(table).__name__!r}')
            raise TypeError(error_message)
        self.table = table
        self.case_control = case_control
        self.voom_weights = voom_weights
        self.voom_plot_data = voom_plot_data
    
    def __repr__(self) -> str:
        """
        Get a string representation of this DE object.
        
        Returns:
            A string summarizing the object.
        """
        num_cell_types = self.table['cell_type'].n_unique()
        descr = (
            f'DE object with {len(self.table):,} '
            f'{"entries" if len(self.table) != 1 else "entry"} across '
            f'{num_cell_types:,} {plural("cell type", num_cell_types)}')
        return descr
    
    def __eq__(self, other: DE) -> bool:
        """
        Test for equality with another DE object.
        
        Args:
            other: the other DE object to test for equality with

        Returns:
            Whether the two DE objects are identical.
        """
        if not isinstance(other, DE):
            error_message = (
                f'the left-hand operand of `==` is a DE object, but '
                f'the right-hand operand has type {type(other).__name__!r}')
            raise TypeError(error_message)
        return self.table.equals(other.table) and \
            self.case_control == other.case_control and \
            (other.voom_weights is None if self.voom_weights is None else
             self.voom_weights.keys() == other.voom_weights.keys() and
             all(self.voom_weights[cell_type].equals(
                     other.voom_weights[cell_type]) and
                 self.voom_plot_data[cell_type].equals(
                     other.voom_plot_data[cell_type])
                 for cell_type in self.voom_weights))
    
    def save(self, directory: str | Path, overwrite: bool = False) -> None:
        """
        Save a DE object to `directory` (which must not exist unless
        `overwrite=True`, and will be created) with the table at table.parquet.
        
        If the DE object contains voom info (i.e. was created with
        `return_voom_info=True` in `Pseudobulk.DE()`, the default), also saves
        each cell type's voom weights and voom plot data to
        f'{cell_type}_voom_weights.parquet' and
        f'{cell_type}_voom_plot_data.parquet', as well as a text file,
        cell_types.txt, containing the cell types.
        
        Args:
            directory: the directory to save the DE object to
            overwrite: if False, raises an error if the directory exists; if
                       True, overwrites files inside it as necessary
        """
        check_type(directory, 'directory', (str, Path),
                   'a string or pathlib.Path')
        directory = str(directory)
        if not overwrite and os.path.exists(directory):
            error_message = (
                f'directory {directory!r} already exists; set overwrite=True '
                f'to overwrite')
            raise FileExistsError(error_message)
        os.makedirs(directory, exist_ok=overwrite)
        self.table.write_parquet(os.path.join(directory, 'table.parquet'))
        if self.voom_weights is not None:
            with open(os.path.join(directory, 'cell_types.txt'), 'w') as f:
                print('\n'.join(self.voom_weights), file=f)
            for cell_type in self.voom_weights:
                escaped_cell_type = cell_type.replace('/', '-')
                self.voom_weights[cell_type].write_parquet(
                    os.path.join(directory, f'{escaped_cell_type}.'
                                            f'voom_weights.parquet'))
                self.voom_plot_data[cell_type].write_parquet(
                    os.path.join(directory, f'{escaped_cell_type}.'
                                            f'voom_plot_data.parquet'))
    
    def get_hits(self,
                 significance_column: str = 'FDR',
                 threshold: int | float | np.integer | np.floating = 0.05,
                 num_top_hits: int | np.integer | None = None) -> pl.DataFrame:
        """
        Get all (or the top) differentially expressed genes.
        
        Args:
            significance_column: the name of a Boolean column of self.table to
                                 determine significance from
            threshold: the significance threshold corresponding to
                       significance_column
            num_top_hits: the number of top hits to report for each cell type;
                          if None, report all hits

        Returns:
            The `table` attribute of this DE object, subset to (top) DE hits.
        """
        check_type(significance_column, 'significance_column', str, 'a string')
        if significance_column not in self.table:
            error_message = 'significance_column is not a column of self.table'
            raise ValueError(error_message)
        check_dtype(self.table[significance_column],
                    f'self.table[{significance_column!r}]', 'floating-point')
        check_type(threshold, 'threshold', (int, float),
                   'a number > 0 and ≤ 1')
        check_bounds(threshold, 'threshold', 0, 1, left_open=True)
        if num_top_hits is not None:
            check_type(num_top_hits, 'num_top_hits', int, 'a positive integer')
            check_bounds(num_top_hits, 'num_top_hits', 1)
        return self.table\
            .filter(pl.col(significance_column) < threshold)\
            .pipe(lambda df: df.group_by('cell_type', maintain_order=True)
                  .head(num_top_hits) if num_top_hits is not None else df)
    
    def get_num_hits(self,
                     significance_column: str = 'FDR',
                     threshold: int | float | np.integer |
                                np.floating = 0.05) -> pl.DataFrame:
        """
        Get the number of differentially expressed genes in each cell type.
        
        Args:
            significance_column: the name of a Boolean column of self.table to
                                 determine significance from
            threshold: the significance threshold corresponding to
                       significance_column

        Returns:
            A DataFrame with one row per cell type and two columns:
            'cell_type' and 'num_hits'.
        """
        check_type(significance_column, 'significance_column', str, 'a string')
        if significance_column not in self.table:
            error_message = 'significance_column is not a column of self.table'
            raise ValueError(error_message)
        check_dtype(self.table[significance_column],
                    f'self.table[{significance_column!r}]', 'floating-point')
        check_type(threshold, 'threshold', (int, float),
                   'a number > 0 and ≤ 1')
        check_bounds(threshold, 'threshold', 0, 1, left_open=True)
        return self.table\
            .filter(pl.col(significance_column) < threshold)\
            .group_by('cell_type', maintain_order=True)\
            .agg(num_hits=pl.len())\
            .sort('cell_type')
    
    # noinspection PyUnresolvedReferences
    def plot_voom(self,
                  cell_type: str,
                  filename: str | Path | None = None,
                  *,
                  ax: 'Axes' | None = None,
                  point_color: Color = '#666666',
                  case_point_color: Color = '#ff6666',
                  point_size: int | float | np.integer | np.floating = 1,
                  case_point_size: int | float | np.integer | np.floating = 1,
                  line_color: Color = '#000000',
                  case_line_color: Color = '#ff0000',
                  line_width: int | float | np.integer | np.floating = 1.5,
                  case_line_width: int | float | np.integer |
                                   np.floating = 1.5,
                  scatter_kwargs: dict[str, Any] | None = None,
                  case_scatter_kwargs: dict[str, Any] | None = None,
                  plot_kwargs: dict[str, Any] | None = None,
                  case_plot_kwargs: dict[str, Any] | None = None,
                  legend_labels: list[str] |
                                 tuple[str, str] = ('Controls', 'Cases'),
                  legend_kwargs: dict[str, Any] | None = None,
                  xlabel: str = 'Average log2(count + 0.5)',
                  xlabel_kwargs: dict[str, Any] | None = None,
                  ylabel: str = 'sqrt(standard deviation)',
                  ylabel_kwargs: dict[str, Any] | None = None,
                  title: bool | str | dict[str, str] |
                         Callable[[str], str] = False,
                  title_kwargs: dict[str, Any] | None = None,
                  despine: bool = True,
                  savefig_kwargs: dict[str, Any] | None = None) -> None:
        """
        Generate a voom plot for a cell type that differential expression was
        calculated for.
        
        Voom plots consist of a scatter plot with one point per gene. They
        visualize how the mean expression of each gene across samples (x)
        relates to its variation across samples (y). The plot also includes a
        LOESS (also called LOWESS) fit, a type of non-linear curve fit, of the
        mean-variance (x-y) trend.
        
        Specifically, the x position of a gene's point is the average, across
        samples, of the base-2 logarithm of the gene's count in each sample
        (plus a pseudocount of 0.5): in other words, mean(log2(count + 0.5)).
        The y position is the square root of the standard deviation, across
        samples, of the gene's log counts per million after regressing out,
        across samples, the differential expression design matrix.
        
        For case-control differential expression (`case_control=True` in
        `Pseudobulk.DE()`), voom is run separately for cases and controls
        ("voomByGroup"), and so the voom plot will show a separate LOESS
        trendline for each of the two groups, with the points and trendlines
        for the two groups shown in different colors.
        
        Args:
            cell_type: the cell type to generate the voom plot for
            filename: the file to save to. If None, generate the plot but do
                      not save it, which allows it to be shown interactively or
                      modified further (e.g. by adding a title or axis labels)
                      before saving.
            ax: the Matplotlib axes to save the plot onto; if None, create a
                new figure with Matpotlib's constrained layout and plot onto it
            point_color: the color of the points in the voom plot; if
                         case-control, only points for controls will be plotted
                         in this color
            case_point_color: the color of the points for cases; ignored for
                              non-case-control differential expression
            point_size: the size of the points in the voom plot; if
                        case-control, only the control points will be plotted
                        with this size
            case_point_size: the size of the points for cases; ignored for
                             non-case-control differential expression
            line_color: the color of the LOESS trendline in the voom plot; if
                        case-control, only the control trendline will be
                        plotted in this color
            case_line_color: the color of the LOESS trendline for cases;
                             ignored for non-case-control differential
                             expression
            line_width: the width of the LOESS trendline in the voom plot; if
                        case-control, only the control trendline will be
                        plotted with this width
            case_line_width: the width of the LOESS trendline for cases;
                             ignored for non-case-control differential
                             expression
            scatter_kwargs: a dictionary of keyword arguments to be passed to
                            `ax.scatter()`, such as:
                            - `rasterized`: whether to convert the scatter plot
                              points to a raster (bitmap) image when saving to
                              a vector format like PDF. Defaults to True,
                              instead of the Matplotlib default of False.
                            - `marker`: the shape to use for plotting each cell
                            - `norm`, `vmin`, and `vmax`: control how the
                              numbers in `color_column` are converted to
                              colors, if `color_column` is numeric
                            - `alpha`: the transparency of each point
                            - `linewidths` and `edgecolors`: the width and
                              color of the borders around each marker. These
                              are absent by default (`linewidths=0`), unlike
                              Matplotlib's default. Both arguments can be
                              either single values or sequences.
                            - `zorder`: the order in which the cells are
                              plotted, with higher values appearing on top of
                              lower ones.
                            Specifying `s` or `c`/`color` will raise an error,
                            since these arguments conflict with the
                            `point_size` and `point_color` arguments,
                            respectively.
                            If case-control and `case_scatter_kwargs` is not
                            None, these settings only apply to control points.
            case_scatter_kwargs: a dictionary of keyword arguments to be passed
                                 to `ax.scatter()` for case points. Like for
                                 `scatter_kwargs`, `rasterized=True` is the
                                 default, and specifying `s` or `c`/`color`
                                 will raise an error. If None and
                                 `scatter_kwargs` is not None, the settings in
                                 `scatter_kwargs` apply to all points. Can only
                                 be specified for case-control differential
                                 expression.
            plot_kwargs: a dictionary of keyword arguments to be passed to
                         `ax.plot()` when plotting the trendlines, such as
                         `linestyle` for dashed trendlines. Specifying
                         `color`/`c` or `linewidth` will raise an error, since
                         these arguments conflict with the `line_color` and
                         `line_width` arguments, respectively.
            case_plot_kwargs: a dictionary of keyword arguments to be passed to
                              `ax.plot()` when plotting the case trendlines.
                              Specifying `color`/`c` or `linewidth` will raise
                              an error, like for `plot_kwargs`. If None and
                              `plot_kwargs` is not None, the settings in
                              `plot_kwargs` apply to all points. Can only be
                              specified for case-control differential
                              expression.
            legend_labels: a two-element tuple or list of labels for controls
                           and cases (in that order) in the legend, or None to
                           not include a legend. Ignored for non-case-control
                           differential expression.
            legend_kwargs: a dictionary of keyword arguments to be passed to
                           `ax.legend()` to modify the legend, such as:
                           - `loc`, `bbox_to_anchor`, and `bbox_transform` to
                             set its location.
                           - `prop`, `fontsize`, and `labelcolor` to set its
                             font properties
                           - `facecolor` and `framealpha` to set its background
                             color and transparency
                           - `frameon=True` or `edgecolor` to add or color
                             its border (`frameon` is False by default,
                             unlike Matplotlib's default of True)
                           - `title` to add a legend title
                           Can only be specified for case-control differential
                           expression.
            xlabel: the x-axis label, or None to not label the x-axis
            xlabel_kwargs: a dictionary of keyword arguments to be passed to
                           `ax.set_xlabel()` to control the text properties,
                           such as `color` and `size` to modify the text
                           color/size
            ylabel: the y-axis label, or None to not label the y-axis
            ylabel_kwargs: a dictionary of keyword arguments to be passed to
                           `ax.set_ylabel()` to control the text properties,
                           such as `color` and `size` to modify the text
                           color/size
            title: what to use as the title. If False, do not include a
                   title. If True, use the cell type as a title. If a string,
                   use the string as the title for every cell type. If a
                   dictionary, use `title[cell_type]` as the title; every cell
                   type must be present in the dictionary. If a function, use
                   `title(cell_type)` as the title.
            title_kwargs: a dictionary of keyword arguments to be passed to
                          `ax.set_title()` to control the text properties, such
                          as `color` and `size` to modify the text color/size.
                          Cannot be specified when `title=False`.
            despine: whether to remove the top and right spines (borders of the
                     plot area) from the voom plot
            savefig_kwargs: a dictionary of keyword arguments to be passed to
                            `plt.savefig()`, such as:
                            - `dpi`: defaults to 300 instead of Matplotlib's
                              default of 150
                            - `bbox_inches`: the bounding box of the portion of
                              the figure to save; defaults to 'tight' (crop out
                              any blank borders) instead of Matplotlib's
                              default of None (save the entire figure)
                            - `pad_inches`: the number of inches of padding to
                              add on each of the four sides of the figure when
                              saving. Defaults to 'layout' (use the padding
                              from the constrained layout engine) instead of
                              Matplotlib's default of 0.1.
                            - `transparent`: whether to save with a transparent
                              background; defaults to True if saving to a PDF
                              (i.e. when `PNG=False`) and False if saving to
                              a PNG, instead of Matplotlib's default of always
                              being False.
        """
        import matplotlib.pyplot as plt
        # Check that this DE object contains `voom_plot_data`
        if self.voom_plot_data is None:
            error_message = (
                'this DE object does not contain the voom_plot_data '
                'attribute, which is necessary to generate voom plots; re-run '
                'Pseudobulk.DE() with return_voom_info=True to include this '
                'attribute')
            raise AttributeError(error_message)
        # Check that `cell_type` is a cell type in this DE object
        check_type(cell_type, 'cell_type', str, 'a string')
        if cell_type not in self.voom_plot_data:
            error_message = \
                f'cell_type {cell_type!r} is not a cell type in this DE object'
            raise ValueError(error_message)
        # If `filename` is not None, check that it is a string or pathlib.Path
        # and that its base directory exists; if `filename` is None, make sure
        # `savefig_kwargs` is also None
        if filename is not None:
            check_type(filename, 'filename', (str, Path),
                       'a string or pathlib.Path')
            directory = os.path.dirname(filename)
            if directory and not os.path.isdir(directory):
                error_message = (
                    f'{filename} refers to a file in the directory '
                    f'{directory!r}, but this directory does not exist')
                raise NotADirectoryError(error_message)
            filename = str(filename)
        elif savefig_kwargs is not None:
            error_message = 'savefig_kwargs must be None when filename is None'
            raise ValueError(error_message)
        # Check that each of the colors are valid Matplotlib colors, and
        # convert them to hex
        for color, color_name in ((point_color, 'point_color'),
                                  (line_color, 'line_color'),
                                  (case_point_color, 'case_point_color'),
                                  (case_line_color, 'case_line_color')):
            if not plt.matplotlib.colors.is_color_like(color):
                error_message = f'{color_name} is not a valid Matplotlib color'
                raise ValueError(error_message)
        point_color = plt.matplotlib.colors.to_hex(point_color)
        line_color = plt.matplotlib.colors.to_hex(line_color)
        case_point_color = plt.matplotlib.colors.to_hex(case_point_color)
        case_line_color = plt.matplotlib.colors.to_hex(case_line_color)
        # Check that point sizes are positive numbers
        check_type(point_size, 'point_size', (int, float), 'a positive number')
        check_bounds(point_size, 'point_size', 0, left_open=True)
        check_type(case_point_size, 'case_point_size', (int, float),
                   'a positive number')
        check_bounds(case_point_size, 'case_point_size', 0, left_open=True)
        # For each of the kwargs arguments, if the argument is not None, check
        # that it is a dictionary and that all its keys are strings.
        for kwargs, kwargs_name in (
                (scatter_kwargs, 'scatter_kwargs'),
                (case_scatter_kwargs, 'case_scatter_kwargs'),
                (plot_kwargs, 'plot_kwargs'),
                (case_plot_kwargs, 'case_plot_kwargs'),
                (legend_kwargs, 'legend_kwargs'),
                (xlabel_kwargs, 'xlabel_kwargs'),
                (ylabel_kwargs, 'ylabel_kwargs'),
                (title_kwargs, 'title_kwargs')):
            if kwargs is not None:
                check_type(kwargs, kwargs_name, dict, 'a dictionary')
                for key in kwargs:
                    if not isinstance(key, str):
                        error_message = (
                            f'all keys of {kwargs_name} must be strings, but '
                            f'it contains a key of type '
                            f'{type(key).__name__!r}')
                        raise TypeError(error_message)
        # Check that `case_scatter_kwargs` and `case_plot_kwargs` are None for
        # non-case-control differential expression. If None, use the settings
        # from `scatter_kwargs` and `plot_kwargs`, respectively. Also set
        # `plot_kwargs` to {} if it is None.
        if case_scatter_kwargs is None:
            case_scatter_kwargs = scatter_kwargs
        elif not self.case_control:
            error_message = (
                'case_scatter_kwargs can only be specified for case-control '
                'differential expression')
            raise ValueError(error_message)
        if plot_kwargs is None:
            plot_kwargs = {}
        if case_plot_kwargs is None:
            case_plot_kwargs = plot_kwargs
        elif not self.case_control:
            error_message = (
                'case_plot_kwargs can only be specified for case-control '
                'differential expression')
            raise ValueError(error_message)
        # Override the defaults for certain keys of `scatter_kwargs` and
        # `case_scatter_kwargs`
        default_scatter_kwargs = dict(rasterized=True, linewidths=0)
        scatter_kwargs = default_scatter_kwargs | scatter_kwargs \
            if scatter_kwargs is not None else default_scatter_kwargs
        case_scatter_kwargs = default_scatter_kwargs | case_scatter_kwargs \
            if case_scatter_kwargs is not None else default_scatter_kwargs
        # Check that `scatter_kwargs` and `case_scatter_kwargs` do not contain
        # the `s` or `c`/`color` keys, and that `plot_kwargs` and
        # `case_plot_kwargs` do not contain the `c`/`color` or `linewidth` keys
        for plot, kwargs_set in enumerate(
                (((scatter_kwargs, 'scatter_kwargs'),
                  (case_scatter_kwargs, 'case_scatter_kwargs')),
                 ((plot_kwargs, 'plot_kwargs'),
                  (case_plot_kwargs, 'case_plot_kwargs')))):
            for kwargs, kwargs_name in kwargs_set:
                if kwargs is None:
                    continue
                for bad_key, alternate_argument in (
                        ('linewidth', 'line_width') if plot else
                        ('s', 'point_size'),
                        ('c', 'point_color' if plot else 'line_color'),
                        ('color', 'point_color' if plot else 'line_color')):
                    if bad_key in kwargs:
                        error_message = (
                            f'{bad_key!r} cannot be specified as a key in '
                            f'{kwargs_name}; specify the {alternate_argument} '
                            f'argument instead')
                        raise ValueError(error_message)
        # Check that `legend_labels` is a two-element tuple or list of strings.
        # Only add a legend if case-control and `legend_labels` is not None.
        check_type(legend_labels, 'legend_labels', (tuple, list),
                   'a length-2 tuple or list of strings')
        if len(legend_labels) != 2:
            error_message = (
                f'legend_labels must have a length of 2, but has a length of '
                f'{len(legend_labels):,}')
            raise ValueError(error_message)
        check_type(legend_labels[0], 'legend_labels[0]', str, 'a string')
        check_type(legend_labels[1], 'legend_labels[1]', str, 'a string')
        add_legend = self.case_control and legend_labels is not None
        # Override the defaults for certain values of `legend_kwargs`; check
        # that it is None for non-case-control differential expression
        default_legend_kwargs = dict(frameon=False)
        if legend_kwargs is not None:
            if not self.case_control:
                error_message = (
                    'legend_kwargs can only be specified for case-control '
                    'differential expression')
                raise ValueError(error_message)
            legend_kwargs = default_legend_kwargs | legend_kwargs
        else:
            legend_kwargs = default_legend_kwargs
        # Check that `xlabel` and `ylabel` are strings, or None
        if xlabel is not None:
            check_type(xlabel, 'xlabel', str, 'a string')
        if ylabel is not None:
            check_type(ylabel, 'ylabel', str, 'a string')
        # Check that `title` is Boolean, a string, a dictionary where all keys
        # are cell types and every cell type is present, or a function
        check_type(title, 'title', (bool, str, dict, Callable),
                   'Boolean, a string, a dictionary, or a function')
        if isinstance(title, dict):
            if len(title) != len(self.voom_plot_data) or \
                    set(title) != set(self.voom_plot_data):
                error_message = (
                    'when title is a dictionary, all its keys must be cell '
                    'types, and every cell type must be present')
                raise ValueError(error_message)
        # Check that `title_kwargs` is None when `title=False`
        if title is False and title_kwargs is not None:
            error_message = 'title_kwargs cannot be specified when title=False'
            raise ValueError(error_message)
        # Override the defaults for certain values of `savefig_kwargs`
        default_savefig_kwargs = \
            dict(dpi=300, bbox_inches='tight', pad_inches='layout',
                 transparent=filename is not None and
                             filename.endswith('.pdf'))
        savefig_kwargs = default_savefig_kwargs | savefig_kwargs \
            if savefig_kwargs is not None else default_savefig_kwargs
        # Get the data to generate the voom plot from (scatter-plot points and
        # LOESS trendlines)
        voom_plot_data = self.voom_plot_data[cell_type]
        # If `ax` is None, create a new figure with `constrained_layout=True`;
        # otherwise, check that it is a Matplotlib axis
        make_new_figure = ax is None
        try:
            if make_new_figure:
                plt.figure(constrained_layout=True)
                ax = plt.gca()
            else:
                check_type(ax, 'ax', plt.Axes, 'a Matplotlib axis')
            if self.case_control:
                if add_legend:
                    legend_patches = []
                for case in False, True:
                    # Plot the scatter plot for cases and controls
                    ax.scatter(voom_plot_data[f'xy_x_{case}'],
                               voom_plot_data[f'xy_y_{case}'],
                               s=case_point_size if case else point_size,
                               c=case_point_color if case else point_color,
                               **(case_scatter_kwargs if case else
                                  scatter_kwargs))
                    # Plot the LOESS trendline for cases and controls
                    ax.plot(voom_plot_data[f'line_x_{case}'],
                            voom_plot_data[f'line_y_{case}'],
                            c=case_line_color if case else line_color,
                            linewidth=case_line_width if case else line_width,
                            **(case_plot_kwargs if case else plot_kwargs))
                    # Create a rectangle for cases and controls for the legend,
                    # where the border matches the color of the trendline and
                    # the fill matches the color of the scatter plot points
                    if add_legend:
                        # noinspection PyUnboundLocalVariable
                        # noinspection PyUnresolvedReferences
                        legend_patches.append(
                            plt.matplotlib.patches.Patch(
                                facecolor=case_point_color if case else
                                point_color,
                                edgecolor=case_line_color if case else
                                line_color,
                                linewidth=case_line_width if case else
                                line_width,
                                label=legend_labels[case]))
                # Add the legend; override the defaults for certain values of
                # `legend_kwargs`
                if add_legend:
                    ax.legend(handles=legend_patches, **legend_kwargs)
            else:
                # Plot the scatter plot
                ax.scatter(voom_plot_data['xy_x'], voom_plot_data['xy_y'],
                           s=point_size, c=point_color, **scatter_kwargs)
                # Plot the LOESS trendline
                ax.plot(voom_plot_data['line_x'], voom_plot_data['line_y'],
                         c=line_color, linewidth=line_width, **plot_kwargs)
            # Add the title, axis labels and axis limits
            if xlabel_kwargs is None:
                xlabel_kwargs = {}
            if ylabel_kwargs is None:
                ylabel_kwargs = {}
            ax.set_xlabel(xlabel, **xlabel_kwargs)
            ax.set_ylabel(ylabel, **ylabel_kwargs)
            if title is not False:
                if title_kwargs is None:
                    title_kwargs = {}
                # noinspection PyCallingNonCallable
                ax.set_title(title[cell_type] if isinstance(title, dict)
                             else title if isinstance(title, str) else
                             title(cell_type) if isinstance(title, Callable)
                             else cell_type, **title_kwargs)
            # Despine, if specified
            if despine:
                spines = ax.spines
                spines['top'].set_visible(False)
                spines['right'].set_visible(False)
            # Save; override the defaults for certain keys of `savefig_kwargs`
            if filename is not None:
                default_savefig_kwargs = \
                    dict(dpi=300, bbox_inches='tight', pad_inches='layout',
                         transparent=filename is not None and
                                     filename.endswith('.pdf'))
                savefig_kwargs = default_savefig_kwargs | savefig_kwargs \
                    if savefig_kwargs is not None else default_savefig_kwargs
                plt.savefig(filename, **savefig_kwargs)
                if make_new_figure:
                    plt.close()
        except:
            # If we made a new figure, make sure to close it if there's an
            # exception (but not if there was no error and `filename` is None,
            # in case the user wants to modify it further before saving)
            if make_new_figure:
                plt.close()
            raise
