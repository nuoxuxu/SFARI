from __future__ import annotations
import os
import polars as pl
import signal
import subprocess
import sys
from contextlib import contextmanager
from functools import cache, reduce, wraps
from typing import Any, Iterable
pl.enable_string_cache()

###############################################################################
# [1] General utilities
###############################################################################


@contextmanager
def Timer(message=None, verbose=True):
    """
    Use "with Timer(message):" to time the code inside the with block. Based on
    preshing.com/20110924/timing-your-code-using-pythons-with-statement
    
    Args:
        message: a message to print when starting the with block (with "..."
                 after) and ending the with block (with the time after)
        verbose: if False, disables the Timer. This is useful to conditionally
                 run the Timer based on the value of a boolean variable.
    """
    if verbose:
        from timeit import default_timer
        if message is not None:
            print(f'{message}...')
        start = default_timer()
        aborted = False
        try:
            yield
        except Exception as e:
            aborted = True
            raise e
        finally:
            end = default_timer()
            duration = end - start
            
            days = int(duration // 86400)
            hours = int((duration % 86400) // 3600)
            minutes = int((duration % 3600) // 60)
            seconds = int(duration % 60)
            milliseconds = int((duration * 1000) % 1000)
            microseconds = int((duration * 1000000) % 1000)
            nanoseconds = int((duration * 1000000000) % 1000)
            
            time_parts = []
            if days > 0:
                time_parts.append(f'{days} {plural("day", days)}')
            if hours > 0:
                time_parts.append(f'{hours}h')
            if minutes > 0:
                time_parts.append(f'{minutes}m')
            if seconds > 0:
                time_parts.append(f'{seconds}s')
            if milliseconds > 0:
                time_parts.append(f'{milliseconds}ms')
            if microseconds > 0:
                time_parts.append(f'{microseconds}µs')
            if nanoseconds > 0:
                time_parts.append(f'{nanoseconds}ns')
            
            time_str = \
                ' '.join(time_parts[:2]) if time_parts else 'less than 1ns'
            
            print(f'{message if message is not None else "Command"} '
                  f'{"aborted after" if aborted else "took"} '
                  f'{time_str}')
    else:
        yield  # no-op


@contextmanager
def cd(path, *, create_if_missing=False):
    """
    Use "with cd(path):" to temporarly change directory inside the with block
    
    Args:
        path: the directory to change to temporarily
        create_if_missing: whether to create the directory if missing
    """
    if create_if_missing:
        os.makedirs(path, exist_ok=True)
    original_dir = os.getcwd()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(original_dir)


@contextmanager
def SuppressMessages():
    """
    Use "with SuppressMessages:" to suppress stdout inside the with block
    """
    try:
        sys.stdout = open(os.devnull, "w")
        yield
    finally:
        sys.stdout = sys.__stdout__


def raise_error_if_on_compute_node(message=None):
    """
    Raises an error if the user is on a compute node
    
    Args:
        message: Message to print when user is not on a login node; if None,
                 prints a default message
    """
    import re
    import socket
    if re.search('(nia|nc|nl|ng)[0-9]', socket.gethostname()):
        # nia matches the Niagara compute nodes; nc/nl/ng match the Narval ones
        import inspect
        calling_function = inspect.currentframe().f_back.f_code.co_name
        raise RuntimeError(message if message is not None else
                           f'{calling_function}() needs internet access! Run '
                           f'once on the login node to download required '
                           f'files, then re-run here on this compute node')


def check_cluster(cluster):
    """
    Check the cluster the user is on.
    
    Returns:
        The cluster: "narval" or "niagara". Raises an error if $CLUSTER is not
        set to one of those two.
    """
    if cluster is None:
        raise RuntimeError('The environment variable $CLUSTER is not set; it '
                           'must be set to "narval" or "niagara"')
    if cluster != 'narval' and cluster != 'niagara':
        raise RuntimeError(f"The environment variable $CLUSTER is set to "
                           f"{cluster!r}, but must be set to 'narval' or "
                           f"'niagara'")


@cache
def get_base_data_directory():
    """
    Get the location of the "base" data directory, where data will be stored.
    
    Returns:
        The path of the base data directory
    """
    cluster = os.environ.get('CLUSTER')
    return '/home/wainberg/projects/def-wainberg' if cluster == 'narval' else \
        '/scratch/w/wainberg/wainberg' if cluster == 'niagara' else '.'


def print_large_objects(N=20):
    """
    Print the N largest objects in the user's Python session.
    
    Args:
        N: the number of largest objects to print
    """
    # There are 2 main options for accurately calculating Python object sizes:
    # - objsize.get_deep_size from pypistats.org/packages/objsize
    # - pympler.asizeof.asizeof from pypistats.org/packages/pympler
    from pympler.asizeof import asizeof
    mem = pl.DataFrame({'Variable': list(globals()),
                        'Size': [asizeof(eval(key)) for key in globals()]})\
        .with_columns(pl.all().sort_by('Size', descending=True))
    with pl.Config(tbl_hide_dataframe_shape=True):
        print(mem.head(N))


def escape(x):
    """
    Escapes a string or polars Series or expression by removing leading and
    trailing whitespace and converting internal groups of 1+ whitespace
    characters to a single underscore.
    
    Args:
        x: the string, Series of expression to escape

    Returns:
        The escaped version of x.
    """
    if isinstance(x, str):
        import re
        return re.sub(r'\W+', '_', x).strip('_')
    else:
        return x.str.replace(r'\W+', '_').str.strip_chars('_')


def plural(string, count):
    """
    Adds an s to the end of string, unless count is 1 or -1.
    
    Args:
        string: a string
        count: a count

    Returns:
        string, with an s at the end fo count is 1 or -1
    """
    return string if abs(count) == 1 else f'{string}s'


def is_integer(variable):
    """
    Check if variable is an integer. Count integer NumPy "generics" (scalars)
    like np.int32(0) as integers.
    
    Args:
        variable: the variable to be checked

    Returns:
        Whether variable is an integer.
    """
    if isinstance(variable, int):
        return True
    import numpy as np
    return isinstance(variable, np.generic) and \
        np.issubdtype(variable.dtype, np.integer)


def check_type(variable: Any, variable_name: str,
               expected_types: type | tuple[type, ...],
               expected_type_name: str):
    """
    Raise a TypeError if `variable` is not of the expected type.
    
    Args:
        variable: the variable to be checked
        variable_name: the name of the variable, used in the error message
        expected_types: the expected type or types (specifying int, float, or
                        bool also implicitly includes their NumPy equivalents)
        expected_type_name: the name of the expected type, used in the error
                            message (e.g. 'a polars DataFrame')
    """
    if isinstance(variable, expected_types):
        return
    if not isinstance(expected_types, tuple):
        expected_types = expected_types,
    for t in expected_types:
        if t in (int, float, bool):
            import numpy as np
            if isinstance(variable, np.generic):
                if t == int and np.issubdtype(variable.dtype, np.integer) or \
                    t == float and \
                        np.issubdtype(variable.dtype, np.floating) or \
                        variable.dtype == np.bool_:
                    return
    error_message = (
        f'{variable_name} must be {expected_type_name}, but has type '
        f'{type(variable).__name__!r}')
    raise TypeError(error_message)


def check_types(variable: Iterable[Any], variable_name: str,
                expected_types: type | tuple[type, ...],
                expected_type_name: str):
    """
    Raise a TypeError if not all elements of `variable` are strings.
    
    Unlike `check_type`, specifying int, float, or bool does not implicitly
    include their NumPy equivalents.
    
    Args:
        variable: the variable to be checked
        variable_name: the name of the variable, used in the error message
        expected_types: the expected type or types
        expected_type_name: the name of the expected type, used in the error
                            message (e.g. 'polars DataFrames')
    """
    for element in variable:
        if not isinstance(element, expected_types):
            error_message = (
                f'all elements of {variable_name} must be '
                f'{expected_type_name}, but it contains an element of type '
                f'{type(element).__name__!r}')
            raise TypeError(error_message)


def check_dtype(series: pl.Series, series_name: str,
                expected_dtypes: pl.datatypes.classes.DataTypeClass | str |
                                 tuple[pl.datatypes.classes.DataTypeClass | str,
                                      ...]):
    """
    Raise a TypeError if series is not of the expected polars dtype.
    
    Args:
        series: the polars Series to be checked
        series_name: the name of the variable, used in the error message
        expected_dtypes: the expected dtype or dtypes. Specify the string
                        'integer' to include all integer dtypes, and
                        'floating-point' to include all floating-point dtypes.
    """
    base_type = series.dtype.base_type()
    if not isinstance(expected_dtypes, tuple):
        expected_dtypes = expected_dtypes,
    for expected_type in expected_dtypes:
        if base_type == expected_type or expected_type == 'integer' and \
                base_type in pl.INTEGER_DTYPES or \
                expected_type == 'floating-point' and \
                base_type in pl.FLOAT_DTYPES:
            return
    if len(expected_dtypes) == 1:
        expected_dtypes = str(expected_dtypes[0])
    elif len(expected_dtypes) == 2:
        expected_dtypes = ' or '.join(map(str, expected_dtypes))
    else:
        expected_dtypes = ', '.join(map(str, expected_dtypes[:-1])) + \
                          ', or ' + str(expected_dtypes[-1])
    error_message = (
        f'{series_name} must be {expected_dtypes}, but has data type '
        f'{base_type!r}')
    raise TypeError(error_message)


def check_bounds(variable, variable_name, lower_bound=None, upper_bound=None,
                 *, left_open=False, right_open=False):
    """
    Check whether variable is between lower bound and upper bound, inclusive.
    
    Args:
        variable: the variable to be checked
        variable_name: the name of the variable, used in the error message
        lower_bound: the smallest allowed value for variable, or None to have
                     no lower bound
        upper_bound: the largest allowed value for variable, or None to have no
                     upper bound
        left_open: if True, require variable to be strictly greater than
                   lower_bound, rather than >= lower_bound; has no effect if
                   lower_bound is None
        right_open: if True, require variable to be strictly less than
                    upper_bound, rather than <= upper_bound; has no effect if
                    upper_bound is None
    """
    if lower_bound is None and upper_bound is None:
        error_message = 'lower_bound and upper_bound cannot both be None'
        raise ValueError(error_message)
    if lower_bound is not None and (variable <= lower_bound if left_open
                                    else variable < lower_bound) or \
            upper_bound is not None and (variable >= upper_bound if right_open
                                         else variable > upper_bound):
        error_message = f'{variable_name} is {variable:,}, but must be'
        if lower_bound is not None:
            error_message += f' {">" if left_open else "≥"} {lower_bound:,}'
            if upper_bound is not None:
                error_message += ' and'
        if upper_bound is not None:
            error_message += f' {"<" if right_open else "≤"} {upper_bound:,}'
        raise ValueError(error_message)


def check_R_variable_name(
        R_variable_name: str,
        variable_name: str,
        R_keywords: set[str] = {
            'if', 'else', 'repeat', 'while', 'function', 'for', 'in', 'next',
            'break', 'TRUE', 'FALSE', 'NULL', 'Inf', 'NaN', 'NA',
            'NA_integer_', 'NA_real_', 'NA_complex_', 'NA_character_',
            '...'}) -> None:
    """
    Raise an error if R_variable_name is not a valid variable name in R.
    
    Args:
        R_variable_name: the R variable name to be checked
        variable_name: the name of the Python variable the R variable name
                       `R_variable_name` is stored in
        R_keywords: the set of R reserved keywords to check against
    """
    import re
    if not R_variable_name:
        error_message = f'{variable_name} is an empty string'
        raise ValueError(error_message)
    if R_variable_name[0] == '.':
        if len(R_variable_name) > 1 and R_variable_name[1].isdigit():
            error_message = (
                f'{variable_name} {R_variable_name!r} starts with a period '
                f'followed by a digit, which is not a valid R variable name')
            raise ValueError(error_message)
    elif not R_variable_name[0].isidentifier():
        error_message = (
            f'{variable_name} {R_variable_name!r} must start with a letter, '
            f'number, period or underscore')
        raise ValueError(error_message)
    if not re.fullmatch(r'[\w.]*', R_variable_name[1:]):
        invalid_characters = \
            sorted(set(re.findall(r'[^\w.]',
                                  ''.join(dict.fromkeys(R_variable_name)))))
        if len(invalid_characters) == 1:
            description = f"the character '{invalid_characters[0]}'"
        else:
            description = f"the characters " + ", ".join(
                f"'{character}'" for character in invalid_characters) + \
                f" and '{invalid_characters[-1]}'"
        error_message = (
            f'{variable_name} {R_variable_name!r} contains {description}, but '
            f'must contain only letters, numbers, periods and underscores')
        raise ValueError(error_message)
    if R_variable_name in R_keywords or (R_variable_name.startswith('..') and
                                         R_variable_name[2:].isdigit()):
        error_message = (
            f'{variable_name} {R_variable_name!r} is a reserved keyword in R, '
            f'and cannot be used as a variable name')
        raise ValueError(error_message)


class ProcessPool(object):
    """
    Like multiprocessing.Pool but 1) child processes ignore KeyboardInterrupts
    and 2) apply_async() is called submit() and takes actual *args and **kwargs
    rather than a list of args and a dictionary of kwargs (like
    ProcessPoolExecutor.submit() from the concurrent.futures module)

    Attributes:
        pool (multiprocessing.Pool): the underlying multiprocessing.Pool object
    """
    def __init__(self, max_concurrent, start_method='forkserver'):
        """
        Sets up the multiprocessing.Pool object underlying this ProcessPool.
        
        Args:
            max_concurrent: The number of worker processes to be spawned by the
                            Pool object, i.e. the maximum number of concurrent
                            processes that will be allowed to run at once.
            start_method: How worker processes will be started. Possible
                          values are 'fork', 'spawn', 'forkserver'. For
                          details, see docs.python.org/3/library/
                          multiprocessing.html#contexts-and-start-methods.
        """
        import multiprocessing
        self.pool = multiprocessing.get_context(start_method)\
            .Pool(max_concurrent, initializer=self.ignore_keyboard_interrupts)
    
    @staticmethod
    def ignore_keyboard_interrupts():
        """
        When provided as an initializer to a multiprocessing.Pool object,
        this function tells the Pool to ignore KeyboardInterrupts (i.e. SIGINT)
        so that pressing Ctrl + C doesn't kill all your background processes.
        """
        signal.signal(signal.SIGINT, signal.SIG_IGN)
    
    def submit(self, func, *args, **kwargs):
        """
        Submits a function to the process pool. A thin wrapper over
        Pool.apply_async() that allows the user to pass actual *args and
        **kwargs rather than a list of args and a dictionary of kwargs.
        
        Args:
            func: the function to be submitted
            *args: positional arguments to be passed to function
            **kwargs: keyword arguments to be passed to function

        Returns:
            The multiprocessing.pool.AsyncResult object returned by
            Pool.apply_async().
        """
        return self.pool.apply_async(func, args, kwargs)


@cache
def get_process_pool(max_concurrent, start_method='forkserver'):
    """
    Creates a ProcessPool for a given value of max_concurrent and start_method,
    which due to the @cache decorator will persist across multiple calls to
    functions like run_background() or run_function_background(), so long as
    they call this function with the same max_concurrent and start_method.
    
    Args:
        max_concurrent: The number of worker processes to be spawned by the
                        ProcessPool, i.e. the maximum number of concurrent
                        processes that will be allowed to run at once.
        start_method: How worker processes will be started. Possible
                      values are 'fork', 'spawn', 'forkserver'. For
                      details, see docs.python.org/3/library/
                      multiprocessing.html#contexts-and-start-methods.

    Returns:
        A ProcessPool for the given values of max_concurrent and start_method.
        A new ProcessPool will be created the first time this function is run
        for a given value of max_concurrent and start_method, which will then
        be cached for all subsequent calls with the same max_concurrent and
        start_method.
    """
    return ProcessPool(max_concurrent, start_method=start_method)


@cache
def cython_inline(code, boundscheck=False, cdivision=True,
                  initializedcheck=False, wraparound=False,
                  debug_symbols=False, include_dirs=None, libraries=None,
                  verbose=False, **other_cython_settings):
    """
    A drop-in replacement for cython.inline that supports cimports. It turns on
    the major Cython optimizations (boundscheck=False, cdivision=True,
    initializedcheck=False, wraparound=False), sets quiet=True for quiet
    compilation, and sets language_level=3 for full Python 3 compatibility.
    
    Args:
        code: a string of Cython code to compile
        boundscheck: whether to perform array bounds checking when indexing;
                     always affects array/memoryview indexing, but also affects
                     list, tuple, and string indexing when wraparound=False
        cdivision: whether to use C-style rather than Python-style division and
                   remainder operations; disabling leads to a ~35% speed
                   penalty for these operations
        initializedcheck: whether to check whether memoryviews and C++ classes
                          are initialized before using them
        wraparound: whether to support Python-style negative indexing
        debug_symbols: whether to add debug symbols, so you can run tools like
                       gdb or valgrind; slows down your code
        include_dirs: an optional tuple of include directories of libraries to
                      link against; `np.get_include()` will always be included
        libraries: an optional tuple of libraries to link against,
                   e.g. ('hdf5',)
        verbose: if True, print Cython's compilation logs
        **other_cython_settings: other Cython settings, which will be written
                                 into the source code as #cython compiler
                                 directives
    
    Returns:
        The {function_name: function} dictionary of compiled functions that
        would be returned by cython.inline().
    """
    from hashlib import md5
    from inspect import getmembers
    from textwrap import dedent
    # ~ is read-only on Niagara compute nodes, so build in CYTHON_CACHE_DIR in
    # scratch instead
    cython_cache_dir = os.path.abspath(os.environ.get(
        'CYTHON_CACHE_DIR', os.path.expanduser('~/.cython')))
    os.makedirs(cython_cache_dir, exist_ok=True)
    # Remove extra levels of indentation from the code string (since it's
    # usually defined inside a function, so there's at least one extra level of
    # indentation that would cause a syntax error if not removed) and remove
    # any leading newlines if present (since when users define code strings,
    # they usually use triple-quoted strings, and the code usually doesn't
    # start until the line after the three opening quotes, leading to a single
    # leading newline)
    settings = dict(language_level=3, boundscheck=boundscheck,
                    cdivision=cdivision, initializedcheck=initializedcheck,
                    wraparound=wraparound, **{'warn.undeclared': True})
    settings.update(other_cython_settings)
    code = ''.join(f'#cython: {setting_name}={setting}\n'
                   for setting_name, setting in settings.items()) + \
           dedent(code)
    # Make a short alphabetic module name by taking the code string's MD5 hash
    # and converting the hexadegimal digits to letters (0 -> a, 1 -> b, ...,
    # 9 --> j, a --> k, ..., f --> p)
    module_name = ''.join(chr(ord(c) + (49 if c <= '9' else 10))
                          for c in md5(code.encode('utf-8')).hexdigest())
    code_file = os.path.join(cython_cache_dir, f'{module_name}.pyx')
    # Try to import the module; build it if it does not exist
    sys.path.append(cython_cache_dir)
    try:
        module = __import__(module_name)
    except ModuleNotFoundError:
        # Create the code file
        with open(code_file, 'w') as f:
            print(code, file=f)
        # Write a build script to a temp file based on the module name
        build_file = os.path.join(cython_cache_dir, f'{module_name}_build.py')
        with open(build_file, 'w') as f:
            import numpy as np
            if include_dirs is not None:
                include_dirs = (np.get_include(),) + include_dirs
            else:
                include_dirs = np.get_include(),
            include_dirs = \
                '[' + ', '.join(f'{include_dir!r}'
                                for include_dir in include_dirs) + ']'
            if libraries is not None:
                libraries = \
                    f'[' + ', '.join(f'{library!r}'
                                     for library in libraries) + ']'
            print(dedent(f'''
                from setuptools import Extension, setup
                from Cython.Build import cythonize
                setup(name='{module_name}', ext_modules=cythonize([
                    Extension('{module_name}', ['{code_file}'],
                              language='c++',
                              include_dirs={include_dirs},
                              libraries={libraries},
                              extra_compile_args=['-Ofast', '-march=native',
                                                  '-funroll-loops',
                                                  '-fopenmp', '-Werror'],
                              extra_link_args=['-Ofast', '-fopenmp'])],
                    build_dir='{cython_cache_dir}'))'''), file=f)
        # Build the code (note: `sys.executable` is the location of Python)
        try:
            run(f'cd {cython_cache_dir} && '
                f'{"CFLAGS=-g " if debug_symbols else ""}'
                f'{sys.executable} {build_file} build_ext --inplace'
                f'{"" if verbose else " > /dev/null"}')
        except subprocess.CalledProcessError:
            try:
                os.unlink(code_file)
            except OSError:
                pass
            raise
        # Remove the temp file
        os.unlink(build_file)
        # Try again
        module = __import__(module_name)
    finally:
        sys.path = sys.path[:-1]
    # Create a dict of all the Cython functions defined in the module
    function_dict = {function_name: function
                     for function_name, function in getmembers(module)
                     if repr(function).startswith('<cyfunction')}
    # Return the dict of Cython functions
    return function_dict

    
def prange(range_body, num_threads, *, nogil=True):
    """
    Generates Cython code for either range() or prange() (depending on whether
    num_threads > 1), pasting the code in range_body inside the (p)range().

    Args:
        range_body: a string of Cython code that will be pasted inside the
                    range() or prange()
        num_threads: the number of threads; set `num_threads=None` to use all
                     available cores (as determined by `os.cpu_count()`)
        nogil: whether to release the Global Interpreter Lock (GIL); if already
               released (e.g. via a `with nogil` block), set to False

    Returns:
        Cython code for the prange() or range().
    """
    if num_threads is None:
        num_threads = os.cpu_count()
    return f'prange({range_body}, {"nogil=True, " if nogil else ""}' \
           f'num_threads=num_threads)' \
        if num_threads > 1 else f'range({range_body})'


@cache
def cython_type(dtype):
    """
    Converts a NumPy dtype or string representation of a dtype to its
    corresponding Cython type. Raises a TypeError if dtype isn't recognized or
    is not a Cython type

    Args:
        dtype: An NumPy dtype object (e.g. np.float32) or string representation
               of a dtype (e.g. 'float32').

    Returns:
        str: Corresponding Cython type as a string.
    """
    import numpy as np
    cython_types = {
        'i1': 'char', 'u1': 'unsigned char', 'i2': 'short',
        'u2': 'unsigned short', 'i4': 'int', 'u4': 'unsigned int',
        'i8': 'long', 'u8': 'unsigned long', 'i16': 'long long',
        'u16': 'unsigned long long', 'f2': 'float16', 'f4': 'float',
        'f8': 'double', 'f16': 'long double', 'c8': 'complex float',
        'c16': 'complex double', 'c32': 'complex long double', 'b1': 'bint'}
    try:
        dtype_string = np.dtype(dtype).str[1:]
    except TypeError:
        raise TypeError(f'{dtype!r} is not a valid NumPy dtype')
    try:
        cython_type = cython_types[dtype_string]
    except KeyError:
        raise TypeError(f'{dtype!r} does not correspond to any Cython type')
    return cython_type


def debug(turn_on=True, *, third_party=False):
    """
    Turns on "debug mode", or turns it off if turn_on=False.
    
    In debug mode, whenever you get an error inside a function, local variables
    from inside the function are automatically copied into the global
    namespace, instead of just being discarded.

    Of course, the error may happen many layers of functions deep, so do this
    for every stack frame (nested function call)! Go from the outermost stack
    frame to the innermost, so that variables in inner stack frames overwrite
    variables with the same name from outer stack frames.
    
    However, do not include variables from third-party library code (i.e. code
    files in your miniforge/mambaforge directory), unless
    include_library_variables=True. utils.py is not considered library code!
    
    Implementation details:
    - Uses sys.modules['__main__'].__dict__ (or get_ipython().user_global_ns
      for IPython) instead of globals(): globals() is a module-level variable
      and each module has its own globals(), so we'd only be modifying
      utils.py's globals() and not the REPL's globals()!
    - Imported modules are in f_globals but not f_locals, but functions' local
      variables are in f_locals but not f_globals, so include both.
    
    Args:
        turn_on: whether to turn on (if True) or turn off (if False) debug mode
        third_party: if True, copies variables from third-party library code,
                     not just those from your code
    """
    def add_variables_to_globals(traceback, global_namespace):
        while traceback is not None:
            if third_party or ('miniforge3' not in
                    traceback.tb_frame.f_code.co_filename and 'mambaforge'
                    not in traceback.tb_frame.f_code.co_filename):
                global_namespace.update(traceback.tb_frame.f_globals)
                global_namespace.update(traceback.tb_frame.f_locals)
            traceback = traceback.tb_next
        # Inside polars, the module object "pl" is sometimes reassigned to the
        # module "polars._reexport"; reset it here
        if third_party and 'pl' in global_namespace:
            import polars as pl
            global_namespace['pl'] = pl
    
    try:
        # noinspection PyUnresolvedReferences
        ipython = get_ipython()
    except NameError:
        # Not IPython
        if turn_on:
            def excepthook(exception_class, exception, traceback):
                global_namespace = sys.modules['__main__'].__dict__
                add_variables_to_globals(traceback, global_namespace)
                sys.__excepthook__(exception_class, exception, traceback)
            sys.excepthook = excepthook
        else:
            sys.excepthook = sys.__excepthook__  # reset
    else:
        # IPython
        if turn_on:
            def excepthook(shell, etype, evalue, tb, tb_offset=None):
                global_namespace = ipython.user_global_ns
                add_variables_to_globals(tb, global_namespace)
                shell.showtraceback((etype, evalue, tb), tb_offset=tb_offset)
            ipython.set_custom_exc((Exception,), excepthook)
        else:
            ipython.set_custom_exc((), None)  # reset


def to_tuple(value):
    """
    Cast Iterables (except str/bytes) to tuple, but box non-Iterables (and
    str/bytes) in a length-1 tuple.
    
    Args:
        value: a value

    Returns:
        value as a tuple
    """
    from collections.abc import Iterable
    return tuple(value) if isinstance(value, Iterable) and not \
        isinstance(value, (str, bytes)) else (value,)


def reload(module):
    """
    A drop-in replacement for `importlib.reload()` that also updates existing
    objects' methods (though not class or instance variables).
    
    Specifically, after reloading the given module, it updates all type objects
    in the caller's global namespace that have the same name as a class defined
    in the reloaded module. The update modifies the existing type objects'
    methods to match those of the newly reloaded classes. It then updates the
    `__class__` attribute of all objects defined anywhere (according to
    `gc.get_objects()`) so that `isinstance()` checks don't break.
    
    This function is mainly meant to be called interactively, in which case it
    will update all objects defined in your current Python session.

    Args:
        module: the module to reload

    Returns:
        The newly reloaded module.
    """
    import gc
    import importlib
    from types import FunctionType
    
    def deleted(key):
        error_message = \
            f'method {key!r} no longer exists after running reload()'
        raise AttributeError(error_message)
    
    # Sometimes the reload doesn't work on the first try (due to caching?) so
    # try twice
    for _ in range(2):
        # Perform the standard reload
        module = importlib.reload(module)
        module_name = module.__name__
        module_dict = module.__dict__
        
        # Get the calling frame's globals
        # noinspection PyUnresolvedReferences
        calling_frame_globals = sys._getframe(1).f_globals
        
        # Construct the set of classes to update
        types_to_update = {
            cls for cls in map(type, calling_frame_globals.values())
            if cls.__module__ == module_name and cls.__name__ in module_dict}
        
        # Update each class
        for cls in types_to_update:
            class_dict = cls.__dict__
            class_name = cls.__name__
            new_class = module_dict[class_name]
            new_class_dict = new_class.__dict__
            
            # Add/update methods
            for key, value in new_class_dict.items():
                if isinstance(value, FunctionType):
                    setattr(cls, key, value)
            
            # Remove methods that no longer exist
            # (note: it's not possible to fully remove the method because of
            # the way Python caches method lookups for performance)
            for key in class_dict.keys() - new_class_dict.keys():
                if isinstance(class_dict[key], FunctionType):
                    setattr(cls, key, (
                        lambda key: lambda *args, **kwargs: deleted(key))(key))
            
            # Update __class__ attribute of all instances of the class to point
            # to the new class, so that `isinstance()` checks don't break
            for obj in gc.get_objects():
                if type(obj).__name__ == class_name:
                    obj.__class__ = new_class
    
    return module

###############################################################################
# [2] Polars
###############################################################################


def save_npy(df, filename):
    """
    Saves df to filename in NumPy's .npy binary format.
    
    Args:
        df: a polars DataFrame; all columns must have the same numeric dtype
        filename: a filename to save to. df's data will be saved to
                  f'{filename.removesuffix(".npy")}.npy', and its columns to
                  f'{filename.removesuffix(".npy")}.columns'.
    """
    if df.is_empty():
        raise ValueError('df is empty!')
    dtypes = set(df.dtypes)
    if len(dtypes) > 1 or dtypes.pop() not in pl.NUMERIC_DTYPES:
        raise ValueError('All columns of df must have the same numeric dtype '
                         'to save with save_npy()')
    import numpy as np
    prefix = filename.removesuffix('.npy')
    np.save(f'{prefix}.npy', df.to_numpy())
    pl.DataFrame(df.columns).write_csv(f'{prefix}.columns',
                                       include_header=False)


def load_npy(filename):
    """
    Loads a polars DataFrame saved via save_npy().
    
    Args:
        filename: a filename to load from. The DataFrame's data will be loaded
                  from f'{filename.removesuffix(".npy")}.npy', and its columns
                  from f'{filename.removesuffix(".npy")}.columns'.
    Returns:
        The polars DataFrame.
    """
    import numpy as np
    prefix = filename.removesuffix('.npy')
    return pl.from_numpy(np.load(f'{prefix}.npy'), pl.read_csv(
        f'{prefix}.columns', has_header=False).to_series().to_list())


def print_df(df, num_rows=-1, num_columns=-1):
    """
    Prints the entirety of a polars DataFrame without truncating.
    
    Args:
        df: the DataFrame to print
        num_rows: the number of rows to print (-1 to print all rows)
        num_columns: the number of columns to print (-1 to print all columns)
    """
    with pl.Config(tbl_rows=num_rows, tbl_cols=num_columns):
        print(df)


def print_row(df, row_number=0):
    """
    Prints a row of a polars DataFrame with each column's value on its own line
    alongside its column header. Similar to df.glimpse() but less cluttered and
    only showing one row.
    
    Args:
        df: the DataFrame to print the row of
        row_number: which row to print (by default, the first)
    """
    print_df(df[row_number].unpivot(variable_name='column'))


def filter_columns(df, predicates, *more_predicates):
    """
    Selects columns from a polars DataFrame where all the boolean expressions
    in predicates evaluate to True, like filter() but for columns instead of
    rows. Use it in method chains, e.g. df.pipe(filter_columns,
    pl.all().n_unique() > 1). See github.com/pola-rs/polars/issues/11254.
    
    Args:
        df: a polars DataFrame
        predicates: the boolean expressions to filter on
        *more_predicates: additional boolean expressions, specified as
                          positional arguments

    Returns:
        df, filtered to the columns where all the boolean expressions in
        predicates evaluate to True.
    """
    predicates = to_tuple(predicates) + more_predicates
    boolean_expression = reduce(lambda a, b: a & b, predicates)
    return df.pipe(lambda df: df.select(df.select(boolean_expression)
                                        .unpivot()
                                        .filter(pl.col.value)
                                        ['variable']
                                    .to_list()))


def map_df(df, map_col, other_df, key_col, value_col, *,
           retain_missing=False):
    """
    Maps df[map_col] based on the mapping other_df[key_col] -> 
    other_df[value_col].

    In other words, for each element of df[map_col], check if it's in 
    other_df[key_col], and if so, replace it with the corresponding entry of 
    other_df[value_col].

    Equivalent to df.with_columns(pl.col(map_col).replace_strict(dict(zip(
    other_df[key_col], other_df[value_col])), default=pl.first() if
    retain_missing else None)).
    
    Implementation detail: uses a join, but prefixes other_df's key_col and
    value_col with "__MAP_DF_" to handle the possibility that map_col might
    have the same name as key_col or value_col.
    
    Args:
        df: a polars DataFrame
        map_col: a column in df
        other_df: another polars DataFrame
        key_col: a column in other_df with the mapping keys; all values must be
                 unique, although this is not checked, for speed
        value_col: a column in other_df with the mapping values
        retain_missing: if False, sets elements of map_col that don't appear in
                        key_col to null; if True, leaves them unchanged

    Returns:
        df with map_col transformed so that each of its values that are in
        key_col are transformed to the corresponding value in value_col.
    """
    if key_col == value_col:
        raise ValueError(f'Both key_col and value_col are set to the column '
                         f'name "{key_col}"')
    if isinstance(df, pl.LazyFrame):
        other_df = other_df.lazy()
    if isinstance(other_df, pl.LazyFrame):
        df = df.lazy()
    prefix = '__MAP_DF_'
    df = df.join(other_df.select(pl.col(key_col, value_col)
                                 .name.prefix(prefix)),
                 left_on=map_col, right_on=prefix + key_col, how='left')
    if retain_missing:
        df = df.with_columns(pl.col(prefix + value_col)
                             .fill_null(pl.col(map_col)))
    df = df.with_columns(pl.col(prefix + value_col).alias(map_col))\
        .drop(prefix + value_col)
    return df


def pivot_longer(df, column_groups):
    """
    Transforms df from wide to long format by pivoting multiple lists of
    columns (all of which must have the same length) at the same time, i.e.
    stacking them on top of each other while repeating all the remaining
    columns that aren't listed. You can specify the same column multiple times.
    Specify None in place of a column name to insert nulls instead of a column.
    
    Example:
    >>> df = pl.DataFrame({'A': [1, 2], 'B1': [3, 4], 'B2': [5, 6],
    ...                    'C1': [7, 8], 'C2': [9, 10], 'D': [10, 11]})
    >>> df.pipe(pivot_longer, {'B': ['B1', 'B2'], 'C': ['C1', 'C2'],
    ...                        'D': ['D', None]})
     A  B  C   D
     1  3  7   10
     2  4  8   11
     1  5  9   null
     2  6  10  null
    shape: (4, 4)
    
    A is repeated twice because it's not listed, B1 is stacked on top of B2, C1
    on top of C2, and D on top of an all-nulls column.
    
    Args:
        df: the DataFrame to pivot
        column_groups: a dictionary mapping new column names to groups of old
                       column names that will be stacked to form the new
                       column; all lists must be the same length

    Returns:
        the pivoted DataFrame; columns not listed in column_groups will appear
        first, followed by the columns in column_groups in the order they are
        listed there
    """
    group_lengths = set(map(len, column_groups.values()))
    if len(group_lengths) != 1:
        if len(column_groups) == 0:
            raise ValueError('column_groups is empty')
        else:
            raise ValueError(f'column_groups has values() with a mix of '
                             f'lengths: '
                             f'{", ".join(map(str, sorted(group_lengths)))}')
    # noinspection PyTypeChecker
    listed_columns = set.union(*map(set, column_groups.values()))
    unlisted_columns = [column for column in df.columns
                        if column not in listed_columns]
    return pl.concat([df.select(unlisted_columns, **{
        group_name: group[index] for group_name, group in
        column_groups.items()}) for index in range(group_lengths.pop())])


def polars_numpy_autoconvert(use_columns_from=None):
    """
    Convert polars DataFrames or Series to NumPy arrays before calling func,
    then convert returned NumPy arrays back to polars DataFrames or Series.

    Requires all DataFrames to have the same column names, unless
    use_columns_from is not None, in which case the column names are taken from
    the column(s) specified there.
    
    Args:
        use_columns_from: if a string, the column names applied to the 
                          function's outputs will be taken from its input
                          argument with that name. If a list of strings, the
                          column names will be taken from the first argument in
                          the list that is not None.

    Returns:
        A decorator that can be applied to any function to turn on
        autoconversion of its arguments.
    """
    def decorator(func):
        @wraps(func)
        def polars_to_numpy_wrapper(*args, **kwargs):
            import numpy as np
            columns = None
            # When use_columns_from is not None, get column names from the
            # column(s) it specifies
            if use_columns_from is not None:
                from inspect import signature
                arguments = tuple(signature(func).bind(
                    *args, **kwargs).arguments)
                for argument_name in ([use_columns_from]
                                      if isinstance(use_columns_from, str) else
                                      use_columns_from):
                    try:
                        argument_index = arguments.index(argument_name)
                    except ValueError:
                        continue
                    argument = args[argument_index] \
                        if argument_index < len(args) else \
                        tuple(kwargs.items())[argument_index - len(args)][1]
                    if argument is not None:
                        break
                else:
                    use_columns_joined = \
                        ', '.join(f"{col}" for col in use_columns_from)
                    raise ValueError(f'All of {use_columns_joined} are None!')
                if not isinstance(argument, (np.ndarray, pl.Series,
                                             pl.DataFrame)):
                    raise ValueError(f'The "{argument_name}" argument must '
                                     f'be a NumPy array or polars Series or '
                                     f'DataFrame!')
                if isinstance(argument, pl.DataFrame):
                    columns = argument.columns
            # Convert polars DataFrames to NumPy arrays
            args = list(args)
            for arg_index, arg in enumerate(args):
                if isinstance(arg, pl.DataFrame):
                    if use_columns_from is None:
                        if columns is None:
                            columns = arg.columns
                        else:
                            if columns != arg.columns:
                                raise ValueError('Two arguments are polars '
                                                 'DataFrames with '
                                                 'non-matching columns!')
                    args[arg_index] = arg.to_numpy()
            for kwarg_name, kwarg in kwargs.items():
                if isinstance(kwarg, pl.DataFrame):
                    if use_columns_from is None:
                        if columns is None:
                            columns = kwarg.columns
                        else:
                            if columns != kwarg.columns:
                                raise ValueError('Two arguments are polars '
                                                 'DataFrames with '
                                                 'non-matching columns!')
                    kwargs[kwarg_name] = kwarg.to_numpy()
            # Run the function
            return_values = func(*args, **kwargs)
            # Convert polars DataFrames to NumPy arrays
            # noinspection PyUnboundLocalVariable
            prefix = 'All input DataFrames have' \
                if use_columns_from is None else \
                f'The "{argument_name}" argument has'
            if isinstance(return_values, tuple) and columns is not None:
                return_values = list(return_values)
                for return_index, return_value in enumerate(return_values):
                    if isinstance(return_value, np.ndarray):
                        if return_value.shape[1] != len(columns):
                            error_message = (
                                f'{prefix} {len(columns):,} columns, but a '
                                f'returned NumPy array has '
                                f'{return_value.shape[1]:,} columns!')
                            raise ValueError(error_message)
                        return_values[return_index] = \
                            pl.from_numpy(return_value, columns)
                return tuple(return_values)
            elif columns is not None:
                if isinstance(return_values, np.ndarray):
                    if return_values.shape[1] != len(columns):
                        raise ValueError(f'{prefix} {len(columns):,} columns, '
                                         f'but the returned NumPy array has '
                                         f'{return_values.shape[1]:,} '
                                         f'columns!')
                    return_values = pl.from_numpy(return_values, columns)
            return return_values
        return polars_to_numpy_wrapper
    return decorator


###############################################################################
# [3] Interfacing with bash
###############################################################################


def thread_string(num_threads):
    """
    Generates a string for setting relevant environment variables to limit the
    number of threads for a called process.
    
    Args:
        num_threads: The number of threads.

    Returns:
        A string for setting the relevant environment variables to num_threads.
    """
    return f'export MKL_NUM_THREADS={num_threads}; ' \
           f'export OMP_NUM_THREADS={num_threads}; ' \
           f'export OPENBLAS_NUM_THREADS={num_threads}; ' \
           f'export NUMEXPR_MAX_THREADS={num_threads}; ' \
           if num_threads is not None else ''


def run(cmd, *, log_file=None, unbuffered=False, pipefail=True,
        num_threads=None, **kwargs):
    """
    Runs a bash code segment interactively.
    
    Args:
        cmd: the command to be run
        log_file: a filename to log stdout/stderr to, in addition to printing
        unbuffered: set to True for unbuffered I/O (currently broken when cmd
                    contains multiple commands)
        pipefail: set to False when piping commands to head to avoid errors
        num_threads: set to a positive integer to limit how many threads cmd
                     uses, or to None to not limit the number of threads
        **kwargs: passed on to subprocess.run()

    Returns:
        The CompletedProcess object returned by subprocess.run().
    """
    run_kwargs = dict(check=True, shell=True, executable='/bin/bash')
    run_kwargs.update(**kwargs)
    return subprocess.run(
        f'{thread_string(num_threads)}'
        f'set -eu{"o pipefail" if pipefail else ""}; '
        f'{"stdbuf -i0 -o0 -e0 " if unbuffered else ""}{cmd}'
        f'{f" 2>&1 | tee {log_file}" if log_file is not None else ""}',
        **run_kwargs)


def run_background(cmd, *, log_file=None, unbuffered=False, pipefail=True,
                   num_threads=1, max_concurrent=os.cpu_count(), **kwargs):
    """
    Runs a bash code segment in the background, on the same node. Uses Python's
    multiprocessing module to ensure no more than max_concurrent processes
    run at the same time, to avoid overloading the node.
    
    Background processes will stop when (or shortly after) Python exits; if
    this is problematic, call `.get()` on the return value of this function to
    wait for the process to finish and retrieve its result.
    
    Args:
        cmd: same as run()
        log_file: same as run()
        unbuffered: same as run()
        pipefail: same as run()
        num_threads: same as run(), but the default is now 1, not None
        max_concurrent: if not None, at most max_concurrent processes will run
                        at the same time across all calls to run_background()
        **kwargs: same as run()
    
    Returns:
        An object representing the started background process; specifically,
        the multiprocessing.pool.AsyncResult object returned by
        multiprocessing.Pool.apply_async(). Call `.get()` on this object to
        wait for the process to finish and retrieve its result, as a
        `subprocess.CompletedProcess` object (the same object returned by
        `run()`).
    """
    run_kwargs = dict(shell=True, executable='/bin/bash')
    run_kwargs.update(**kwargs)
    pool = get_process_pool(max_concurrent)
    return pool.submit(
        subprocess.run,
        f'{thread_string(num_threads)}'
        f'set -eu{"o pipefail" if pipefail else ""}; '
        f'{"stdbuf -i0 -o0 -e0 " if unbuffered else ""}{cmd} '
        f'&> {log_file if log_file is not None else "/dev/null"}',
        **run_kwargs)


def run_slurm(cmd, *, job_name='job', log_file=None,
              CPUs=80 if os.environ.get('CLUSTER') == 'niagara' else 1, GPUs=0,
              days=None, hours=None, memory=None,
              partition='compute' if os.environ.get('CLUSTER') == 'niagara'
                         else None,
              verbose=False):
    """
    Runs a bash code segment on another node via slurm. Unlike with
    run_background(), slurm jobs will not stop when Python exits.
    
    Args:
        cmd: same as run()
        job_name: the slurm job name (e.g. for viewing with squeue); cannot be
                  'interactive', since that name is reserved for interactive
                  jobs created with the `n` command
        log_file: same as run()
        CPUs: the number of CPUs to request; 1 by default on Narval but can be
              up to 64; 80 on Niagara since an entire node must be requested
        GPUs: the number of GPUs to request; 0 by default, but can be up to 4
              on Narval; must be 0 on Niagara
        days: the number of days to request; must be an integer between 1 and 7
              on Narval and 1 on Niagara, or None; mutually exclusive with
              hours
        hours: the number of hours to request; must be an integer between 1 and
               168, i.e. 7 days, on Narval, or between 1 and 24 (if partition
               == "compute") or exactly 1 (if partition == "debug") on Niagara;
               mutually exclusive with days. If neither days nor hours are
               specified, default to 1 day (or 1 hour for the Niagara debug
               partition)
        memory: the total amount of memory to request, in GiB (gibibytes, i.e.
                1024 * 1024 * 1024 bytes); defaults to 4 GiB per CPU requested
                on Narval (up to a maximum of 249 GiB, since the most common
                node has 255,000 MiB ~= 249.02 GiB); must be between 1 and 4000
                or None on Narval and None on Niagara
        partition: the Niagara partition to request (e.g. compute, debug);
                   must be None on Narval
        verbose: if True, print the amount of resources requested before
                 requesting the node
    """
    if days is not None and hours is not None:
        raise ValueError('days and hours are mutually exclusive; at least one '
                         'must be None')
    if job_name == 'interactive':
        raise ValueError("job_name cannot be 'interactive', since that name "
                         "is reserved for interactive jobs created with the "
                         "'n' command")
    cluster = os.environ.get('CLUSTER')
    check_cluster(cluster)
    if cluster == 'narval':
        if CPUs not in range(1, 65):
            raise ValueError(f'CPUs was set to "{CPUs}", but must be an '
                             f'integer between 1 and 64 on Narval')
        if GPUs not in range(0, 5):
            raise ValueError(f'GPUs was set to "{GPUs}", but must be an '
                             f'integer between 0 and 4 on Narval')
        if days is not None and days not in range(1, 8):
            raise ValueError(f'days was set to "{days}", but must be an '
                             f'integer between 1 and 7, or None, on Narval')
        if hours is not None and hours not in range(1, 168):
            raise ValueError(f'hours was set to "{hours}", but must be an '
                             f'integer between 1 and 168, or None, on Narval')
        if memory is None:
            memory = min(4 * CPUs, 249)
        elif memory not in range(1, 4001):
            raise ValueError(f'memory was set to "{memory}", but must be an '
                             f'integer between 1 and 4000 on Narval')
        if partition is not None:
            raise ValueError('partition must be None on Narval')
    else:  # cluster == 'niagara'
        if CPUs != 80:
            raise ValueError('CPUs must be set to 80 (the default) on Niagara')
        if GPUs != 0:
            raise ValueError('GPUs must be set to 0 (the default) on Niagara')
        if days is not None and days != 1:
            raise ValueError(f'days was set to "{days}", but must be 1 or '
                             f'None on Niagara')
        if hours is not None and hours not in range(1, 25):
            raise ValueError(f'days was set to "{days}", but must be an '
                             f'integer between 1 and 24, or None, on Niagara')
        if memory is not None:
            raise ValueError('memory must be None on Niagara')
        if partition != 'compute' and partition != 'debug':
            raise ValueError('partition must be "compute" or "debug" on '
                             'Niagara')
        if partition == 'debug' and (days == 1 or
                                     hours is not None and hours > 1):
            raise ValueError('For partition="debug", do not set days or '
                             'hours, or set hours=1')
    runtime = f'{days}-00:00:00' if days is not None else \
        f'{hours}:00:00' if hours is not None else '1-00:00:00' \
            if partition != 'debug' else '1:00:00'
    if verbose:
        memory_description = \
            f' and {memory} GiB memory' if memory is not None else ''
        partition_description = \
            f' on the {partition} partition' if partition is not None else ''
        print(f'Requesting {CPUs} {plural("CPU", CPUs)}, {GPUs} '
              f'{plural("GPU", GPUs)} {memory_description} for {runtime}'
              f'{partition_description}...')
    job_name = job_name.replace(' ', '_')
    from tempfile import NamedTemporaryFile
    try:
        with NamedTemporaryFile('w', dir=os.environ.get('SCRATCH', '.'),
                                suffix='.sh', delete=False) as temp_file:
            partition_settings = f'#SBATCH -p {partition}\n' \
                if partition is not None else ''
            account_settings = '#SBATCH --account=def-wainberg\n' \
                if cluster == 'narval' else ''
            memory_settings = f'#SBATCH --mem {memory}G\n' \
                if memory is not None else ""
            print(
                f'#!/bin/bash\n'
                f'{partition_settings}'
                f'{account_settings}'
                f'#SBATCH -N 1\n'
                f'{f"#SBATCH -G {GPUs}" if GPUs > 0 else ""}'
                f'#SBATCH -n {CPUs}\n'
                f'{memory_settings}'
                f'#SBATCH -t {runtime}\n'
                f'#SBATCH -J {job_name}\n'
                f'{f"#SBATCH -o {log_file}" if log_file is not None else ""}\n'
                f'set -euo pipefail; {cmd}\n',
                file=temp_file)
        sbatch = '.sbatch' if cluster == 'niagara' else 'sbatch'
        # noinspection PyUnresolvedReferences
        sbatch_message = run(f'{sbatch} {temp_file.name}',
                             stdout=subprocess.PIPE)\
            .stdout.decode().rstrip('\n')
        print(f'{sbatch_message} ("{job_name}")')
    finally:
        try:
            os.unlink(temp_file.name)
        except NameError:
            pass


def run_function_background(function, *, log_file=None, unbuffered=False,
                            pipefail=True, num_threads=1,
                            max_concurrent=os.cpu_count(),
                            function_kwargs={}, **kwargs):
    """
    Runs a standalone Python function in the background, similarly to
    run_background().
    
    Background processes will stop when (or shortly after) Python exits; if
    this is problematic, call `.get()` on the return value of this function to
    wait for the process to finish and retrieve its result.
    
    Args:
        function: the function; must import everything it needs and write its
                  output to disk
        log_file: same as run_background()
        unbuffered: same as run_background()
        pipefail: same as run_background()
        num_threads: same as run_background()
        max_concurrent: same as run_background()
        function_kwargs: passed to function as keyword arguments
        **kwargs: same as run_background()
    
    Returns:
        An object representing the started background process; specifically,
        the multiprocessing.pool.AsyncResult object returned by
        multiprocessing.Pool.apply_async(). Call `.get()` on this object to
        wait for the process to finish and retrieve its result, as a
        `subprocess.CompletedProcess` object (the same object returned by
        `run()`).
    """
    from dill.source import getname, getsource
    return run_background(
        f'python -u << EOF &> '
        f'{log_file if log_file is not None else "/dev/null"}\n'
        f'{getsource(function)}\n{getname(function)}'
        f'({", ".join(f"{k}={repr(v)}" for k, v in function_kwargs.items())})'
        f'\nEOF',
        # have to put log_file after the EOF, so set it to None here
        log_file=None, unbuffered=unbuffered, pipefail=pipefail,
        num_threads=num_threads, max_concurrent=max_concurrent, **kwargs)


def run_function_slurm(function, *, job_name='job', log_file=None,
                       CPUs=1 if os.environ.get('CLUSTER') == 'narval' else 80,
                       GPUs=0, days=None, hours=None, memory=None,
                       partition=None if os.environ.get('CLUSTER') ==
                                 'narval' else 'compute',
                       verbose=False, function_kwargs={}):
    """
    Runs a standalone Python function on another node via slurm.
    
    Args:
        function: same as run_function()
        job_name: same as run_slurm()
        log_file: same as run_slurm()
        CPUs: same as run_slurm()
        GPUs: same as run_slurm()
        days: same as run_slurm()
        hours: same as run_slurm()
        memory: same as run_slurm()
        partition: same as run_slurm()
        verbose: same as run_slurm()
        function_kwargs: same as run_function()
    """
    from dill.source import getname, getsource
    run_slurm(
        f'python -u << EOF\n{getsource(function)}\n{getname(function)}'
        f'({", ".join(f"{k}={repr(v)}" for k, v in function_kwargs.items())})'
        f'\nEOF',
        job_name=job_name, log_file=log_file, CPUs=CPUs, GPUs=GPUs, days=days,
        hours=hours, memory=memory, partition=partition, verbose=verbose)


def read_csv_from_command(cmd, *, run_kwargs={}, **kwargs):
    """
    Read a columnar file with polars from the output of a bash command.
    
    Args:
        cmd: the bash command
        run_kwargs: keyword arguments to utils.run()
        **kwargs: keyword arguments to pl.read_csv()

    Returns:
        A polars DataFrame with the contents of the columnar file produced by
        the bash command.
    """
    from io import BytesIO
    # noinspection PyTypeChecker,PyArgumentList
    return pl.read_csv(BytesIO(run(cmd, stdout=subprocess.PIPE, **run_kwargs)
                               .stdout), **kwargs)

        
def read_csv_delim_whitespace(whitespace_delimited_file, **kwargs):
    """
    Reads a whitespace-delimited file with polars, mimicking the behavior of
    delim_whitespace=True in pandas.read_csv().
    
    Explanation of the sed command:
    s/^[[:space:]]\\+// removes leading whitespace
    s/[[:space:]]\\+$// removes trailing whitespace
    s/[[:space:]]\\+/\t/g replaces internal whitespace with tabs
    
    Args:
        whitespace_delimited_file: the whitespace-delimited file
        **kwargs: keyword arguments to pl.read_csv()
    
    Returns:
        A polars DataFrame with the contents of the whitespace-delimted file.
    """
    if not os.path.exists(whitespace_delimited_file):
        error_message = \
            f'No such file or directory: {whitespace_delimited_file}'
        raise FileNotFoundError(error_message)
    if whitespace_delimited_file.endswith('.gz'):
        whitespace_delimited_file = f'<(zcat {whitespace_delimited_file})'
    return read_csv_from_command(
        f"sed 's/^[[:space:]]\\+//; s/[[:space:]]\\+$//; "
        f"s/[[:space:]]\\+/\t/g' {whitespace_delimited_file}",
        separator='\t', **kwargs)


###############################################################################
# [4] Plotting
###############################################################################


def use_font(font_name):
    """
    Tells matplotlib to use the specified font name. A font file containing the
    named font must be available in the user's ~/.fonts directory.
    
    Args:
        font_name: the name of the font to use
    """
    from matplotlib import rcParams
    from matplotlib.font_manager import fontManager, findSystemFonts
    for font_file in findSystemFonts(os.path.expanduser('~/.fonts')):
        fontManager.addfont(font_file)
    # Make sure the font is there
    fontManager.findfont(font_name, fallback_to_default=False)
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = [font_name]


def savefig(filename, *, despine=False, **kwargs):
    """
    A drop-in replacement for plt.savefig(). Saves a matplotlib figure to
    filename, but additionally:
    - removes the right and top spines for a cleaner-looking plot (if
      despine=True)
    - saves with dpi=450 for higher-resolution plots
    - saves with bbox_inches='tight' and pad_inches=0 to avoid extra whitespace
      around plots
    - saves PDF files with transparent=True so transparency info isn't
      discarded
    - closes the plot with plt.close() to avoid memory leaks
    
    Args:
        filename: the filename to save the matplotlib figure to
        despine: whether to remove the right and top spines using seaborn's
                 despine() function
        **kwargs:
    """
    import matplotlib.pyplot as plt
    if despine:
        spines = plt.gca().spines
        spines['top'].set_visible(False)
        spines['right'].set_visible(False)
    all_kwargs = dict(dpi=450, bbox_inches='tight', pad_inches=0,
                      transparent=filename.endswith('pdf'))
    all_kwargs.update(kwargs)
    plt.savefig(filename, **all_kwargs)
    plt.close()


def mantissa_and_exponent(x, base=10):
    """
    Given a floating-point number, returns its mantissa and exponent in the
    specified base (e.g. base 10, base 2).
    
    Args:
        x: the floating-point number
        base: the base

    Returns:
        A tuple of mantissa and exponent.
    """
    import math
    exponent = math.floor(math.log(x, base))
    mantissa = x / base ** exponent
    return mantissa, exponent


def scientific_notation(number, num_decimals, *, threshold=None, latex=False):
    """
    Converts number to a string, using scientific notation with num_decimals
    decimal places if it is below threshold in absolute value.
    
    Args:
        number: the number to convert to scientific notation
        num_decimals: the number of decimal places to format the number with
        threshold: the threshold below which to use scientific notation, versus
                   just displaying the number normally
        latex: whether to use Latex notation, compatible with matplotlib

    Returns:
        The number formatted as a string in scientific notation.
    """
    if (threshold is not None and abs(number) >= threshold) or number == 0:
        return f'{number:.{num_decimals}f}'
    else:
        mantissa, exponent = mantissa_and_exponent(number)
        if f'{mantissa:.{num_decimals}e}' == '1e+01':
            mantissa = 1
            exponent += 1
        if latex:
            return f'{mantissa:.{num_decimals}f} × ' \
                rf'10$^\mathregular{{{exponent}}}$'
        else:
            return f'{mantissa:.{num_decimals}f} × 10^{exponent}'


def grouped_barplot_with_errorbars(data, lower_CIs, upper_CIs, **kwargs):
    """
    Makes a grouped barplot with error bars with pandas. See pandas.pydata.org/
    pandas-docs/stable/user_guide/visualization.html#visualization-errorbars.
    
    Args:
        data: a polars DataFrame of data points: one row per group, one column
              per bar in each group, wich each column corresponding to a
              specific color: rows = groups, columns = bars-in-groups = colors
        lower_CIs: a DataFrame of lower confidence intervals, same size as data
        upper_CIs: a DataFrame of upper confidence intervals, same size as data
        **kwargs: keyword arguments passed on to pd.DataFrame.plot.bar
    
    Consider other visualizations, like stackoverflow.com/a/58041227/1397061.
    """
    import numpy as np
    if data.shape != lower_CIs.shape or data.shape != upper_CIs.shape:
        raise ValueError('data, lower_CIs and upper_CIs must all be the same '
                         'shape!')
    data = data.to_pandas()
    lower_CIs = lower_CIs.to_pandas()
    upper_CIs = upper_CIs.to_pandas()
    data.plot.bar(yerr=np.stack([data - lower_CIs, upper_CIs - data],
                                axis=1).T, **kwargs)


def generate_palette(num_colors, *, lightness_range=(100 / 3, 200 / 3),
                     chroma_range=(50, 100), hue_range=None,
                     first_color='#008cb9', stride=5):
    """
    Generate a maximally perceptually distinct color palette.
    Great for when you want dozens of colors on the same plot!
    
    The first color in the palette is `first_color`. The second color is the
    color that's most perceptually distinct from `first_color`, i.e. has the
    largest distance from it in the perceptually uniform CAM02-UCS color space.
    The third color is the color that has the largest distance from either of
    the first two colors, i.e. the color that maximizes the minimum distance
    to any of the colors currently in the palette. And so on.
    
    Requires the colorspacious package. Install via:
    mamba install -y colorspacious
    
    An optimized version of github.com/taketwo/glasbey that only generates R,
    G, and B values of (0, 5, 10, ..., 255) instead of (0, 1, 2, ..., 255).
    You can change this stride (by default 5) with the `stride` paramter.
    
    For all i < j, generate_palette(i) == generate_palette(j)[:i], so if you 
    need one plot with 10 colors and another with 5, just run this function 
    once for 10 and then do [:5] on the returned palette.
    
    Useful discussion: github.com/holoviz/colorcet/issues/11
    Parameter examples: colorcet.holoviz.org/assets/images/named.png
    
    Args:
        num_colors: the number of colors to include in the palette
        lightness_range: a two-element tuple with the lightness range of colors
                         to generate, or None to take the full range: (0, 100)
        chroma_range: a two-element tuple with the chroma range of colors to
                      generate, or None to take the full range: (0, 100). Grays
                      have low chroma, and vivid colors have high chroma.
        hue_range: a two-element tuple with the hue range of colors to
                   generate, or None to take the full range: (0, 360). Red is
                   at 0°, green at 120°, and blue at 240°. Because it wraps
                   around, the first element of the tuple can be greater than
                   the second, unlike for `lightness_range` and `chroma_range`.
        first_color: a hex code (the `#` symbol followed by 6 hex digits) for
                     the first color of the palette.
        stride: as an optimization, consider only RGB colors where R, G, and B
                are all multiples of this value. Must be a small divisor of 255:
                1, 3, 5, 15, or 17. Set to 1 for the best possible solution, at 
                orders of magnitude more computational cost.
    
    Returns:
        A list of hex codes like #A06B72, with `first_color` as the first hex
        code. This palette can be passed to seaborn's color_palette() function,
        among other possible uses.
    """
    import numpy as np
    from colorspacious import cspace_convert
    from matplotlib.colors import to_hex, to_rgb
    # Check ranges
    for argument, argument_name, max_value in (
            (lightness_range, 'lightness_range', 100),
            (chroma_range, 'chroma_range', 100),
            (hue_range, 'hue_range', 360)):
        if argument is not None:
            check_type(argument, argument_name, tuple, 'a two-element tuple')
            if len(argument) != 2:
                error_message = (
                    f'{argument_name} must be a two-element tuple, but has '
                    f'{len(argument):,} elements')
                raise ValueError(error_message)
            for i in range(2):
                check_type(argument[i], f'{argument_name}[i]', (int, float), 
                           f'a number between 0 and {max_value}, inclusive')
            if argument[0] < 0:
                error_message = f'{argument_name}[0] must be ≥ 0'
                raise ValueError(error_message)
            if argument[1] > max_value:
                error_message = f'{argument_name}[1] must be ≤ {max_value}'
                raise ValueError(error_message)
            if argument is not hue_range and argument[0] > argument[1]:
                error_message = \
                    f'{argument_name}[0] must be ≤ {argument_name}[1]'
                raise ValueError(error_message)
    # Check `first_color`, and convert to the perceptually uniform CAM02-UCS
    # color space
    check_type(first_color, 'first_color', str,
               'a string containing a hex code')
    if not first_color:
        error_message = 'first_color is an empty string'
        raise ValueError(error_message)
    if first_color[0] != '#':
        error_message = 'first_color must start with "#"'
        raise ValueError(error_message)
    first_color = to_rgb(first_color)
    first_color = cspace_convert(first_color, 'sRGB1', 'CAM02-UCS')
    if lightness_range is not None or chroma_range is not None or \
            hue_range is not None:
        lightness, chroma, hue = \
            cspace_convert(first_color, 'CAM02-UCS', 'JCh')
        if lightness_range is not None and \
                not lightness_range[0] <= lightness <= lightness_range[1]:
            error_message = (
                f'first_color has a lightness of {lightness}, outside the '
                f'specified lightness_range of {lightness_range}')
            raise ValueError(error_message)
        if chroma_range is not None and \
                not chroma_range[0] <= chroma <= chroma_range[1]:
            error_message = (
                f'first_color has a chroma of {chroma}, outside the specified '
                f'chroma_range of {chroma_range}')
            raise ValueError(error_message)
        if hue_range is not None and not (hue_range[0] <= hue <= hue_range[1]
                                          if hue_range[0] <= hue_range[1] else
                                          hue_range[0] <= hue or 
                                          hue <= hue_range[1]):
            error_message = (
                f'first_color has a hue of {hue}, outside the specified '
                f'hue_range of {hue_range}')
            raise ValueError(error_message)
    # Check `stride`
    check_type(stride, 'stride', int, 'one of the integers 1, 3, 5, 15, or 17')
    if stride not in (1, 3, 5, 15, 17):
        error_message = 'stride must be 1, 3, 5, 15, or 17'
        raise ValueError(error_message)
    # Generate a lookup table with all possible RGB colors where R, G and B are
    # multiples of 5, encoded in CAM02-UCS space. Table rows correspond to 
    # individual RGB colors; columns correspond to J', a', and b' components.
    rgb = np.arange(0, 256, stride)
    colors = np.empty([len(rgb)] * 3 + [3])
    colors[..., 0] = rgb[:, None, None]
    colors[..., 1] = rgb[None, :, None]
    colors[..., 2] = rgb[None, None, :]
    colors = colors.reshape(-1, 3)
    colors = cspace_convert(colors, 'sRGB255', 'CAM02-UCS')
    # Remove colors outside the specified lightness, chroma and/or hue ranges
    if lightness_range is not None or chroma_range is not None or \
            hue_range is not None:
        jch = cspace_convert(colors, 'CAM02-UCS', 'JCh')
        mask = np.ones(len(colors), dtype=bool)
        if lightness_range is not None:
            mask &= (jch[:, 0] >= lightness_range[0]) & \
                    (jch[:, 0] <= lightness_range[1])
        if chroma_range is not None:
            mask &= (jch[:, 1] >= chroma_range[0]) & \
                    (jch[:, 1] <= chroma_range[1])
        if hue_range is not None:
            if hue_range[0] <= hue_range[1]:
                mask &= (jch[:, 2] >= hue_range[0]) & \
                        (jch[:, 2] <= hue_range[1])
            else:
                mask &= (jch[:, 2] >= hue_range[0]) | \
                        (jch[:, 2] <= hue_range[1])
        colors = colors[mask]
    # Initialize the palette to `first_color`, then iteratively add the color
    # that's farthest away from all other colors (i.e. with the maximum min
    # distance to any color already in the palette)
    palette = [first_color]
    distances = np.full(len(colors), np.inf)
    while len(palette) < num_colors:
        # Update palette-colors distances to account for the color just added
        distance_to_newest_color = \
            np.linalg.norm((colors - palette[-1]), axis=1)
        np.minimum(distances, distance_to_newest_color, distances)
        # Add the color with the new maximum distance
        palette.append(colors[distances.argmax()])
    # Convert the generated palette to sRGB1 format
    palette = cspace_convert(palette, 'CAM02-UCS', 'sRGB1')
    # Clip palette to [0, 1], in case some colors are slightly out-of-range
    palette = palette.clip(0, 1)
    # Convert RGB to hex
    palette = np.apply_along_axis(to_hex, 1, palette)
    return palette


def ordered_legend(order, *, ax=None, **kwargs):
    """
    A drop-in replacement for ax.legend() that lets you specify an explicit
    ordering of the items. You can run it even if a legend has already been
    created, which will reorder the items after the fact.
    
    For instance, ordered_legend(['Top', 'Middle', 'Bottom']) makes the first
    legend entry 'Top', the second 'Middle', and the third 'Bottom'.
    
    Adapted from stackoverflow.com/a/35926913/1397061.
    
    Args:
        order: the ordering of the items
        ax: an axis; if None, use plt.gca()
        **kwargs: keyword arguments to be passed to ax.legend()
    """
    import matplotlib.pyplot as plt
    import numpy as np
    if ax is None:
        ax = plt.gca()
    ax.legend(**kwargs)
    handles, labels = ax.get_legend_handles_labels()
    keys = dict(zip(order, range(len(order))))
    labels, handles = zip(*sorted(zip(labels, handles),
                                  key=lambda t: keys.get(t[0], np.inf)))
    ax.legend(handles, labels, **kwargs)


def manhattan_plot(sumstats, *, clumped_variants=None, genome_build=None,
                   SNP_col='SNP', chrom_col='CHROM', bp_col='BP',
                   ref_col='REF', alt_col='ALT', p_col='P',
                   max_p=1e-3, max_p_to_label=5e-8,
                   p_thresholds={5e-8: ('#D62728', 'dashed')},
                   text_padding=3, text_padding_overrides={}, 
                   text_size='small', text_size_overrides={},
                   text_overrides={},
                   gene_annotation_dir=f'{get_base_data_directory()}/'
                                       f'gene-annotations'):
    """
    Make a Manhattan plot of the variants in sumstats; highlight lead variants.
    If clumped_variants is not None, highlight the LD clumps as well.
    If genome_build is not None, label each lead variant with its nearest
    coding gene(s) from that genome build.
    After running, call savefig() to save the plot, or customize further first.
    Style inspired by ncbi.nlm.nih.gov/pmc/articles/PMC6481311/figure/F1.
    
    Args:
        sumstats: the summary statistics to plot
        clumped_variants: a DataFrame of clumped variants from ld_clump() to
                          highlight lead variants on the Manhattan plot, or
                          None to skip
        genome_build: if not None, a genome build used to label each lead
                      variant with its nearest coding gene(s)
        SNP_col: the name of the variant ID column in sumstats
        chrom_col: the name of the chromosome column in sumstats
        bp_col: the name of the base-pair position column in sumstats
        ref_col: the name of the reference allele column in sumstats
        alt_col: the name of the alternate allele column in sumstats
        p_col: the name of the p-value column in sumstats
        max_p: defines the bottom y limit of the Manhattan plot; variants with
               p-values >= this value will not be plotted
        max_p_to_label: variants with p-values >= this value will not be
                        labeled; only has an effect if genome_build is not None
        p_thresholds: a dictionary of p-value thresholds to draw horizontal
                      lines at; values are tuples of (line color, line style)
        text_padding: the amount of vertical and horizontal padding between
                      lead variants and labels in points, either as a single 
                      number or as a pair of (horizontal, vertical) numbers; 
                      only used if genome_build is not None
        text_padding_overrides: A dictionary of {gene_name: padding} used to
                                override text_padding for specific gene labels.
                                gene_name must be a string with a single gene
                                name or a comma-separated sequence of gene
                                names. padding must be a single number or a
                                pair of numbers, with the same semantics as the
                                text_padding argument. Only used if
                                genome_build is not None.
        text_size: the text size for the gene labels, in any format recognized 
                   by Matplotlib; only used if genome_build is not None
        text_size_overrides: A dictionary of {gene_name: size} used to
                             override text_size for specific gene labels.
                             gene_name must be a string with a single gene name
                             or a comma-separated sequence of gene names. Only
                             used if genome_build is not None.
        text_overrides: A dictionary of {gene_name: label} used to override
                        specific gene labels.
        gene_annotation_dir: the directory where the coding gene locations
                             returned by get_coding_genes() will be cached.
                             Must be run on the login node to generate this
                             cache, if it doesn't exist. Because generating the
                             cache can take a long time, you will probably want
                             to leave this argument at its default value.
    """
    import matplotlib.pyplot as plt
    import numpy as np
    if sumstats.is_empty():
        raise ValueError(f'sumstats is empty!')
    if SNP_col not in sumstats:
        raise ValueError(f'{SNP_col!r} not in sumstats; specify SNP_col')
    if chrom_col not in sumstats:
        raise ValueError(f'{chrom_col!r} not in sumstats; specify chrom_col')
    if bp_col not in sumstats:
        raise ValueError(f'{bp_col!r} not in sumstats; specify bp_col')
    if p_col not in sumstats:
        raise ValueError(f'{p_col!r} not in sumstats; specify p_col')
    if genome_build is not None:
        check_valid_genome_build(genome_build=genome_build)
    # Subset sumstats to p < max_p; standardize chromosome names; alternate
    # each chromosome in a different color (dark blue then light blue); offset
    # each chromosome by the cumulative number of base pairs from the start of
    # chr1 to the start of that chromosome
    sumstats = sumstats\
        .filter(pl.col(p_col) < max_p)\
        .with_columns(pl.col(chrom_col)
                      .pipe(standardize_chromosomes, omit_chr_prefix=True))\
        .with_columns(color=pl.col(chrom_col).cast(pl.Categorical)
                      .to_physical().mod(2)
                      .replace_strict({0: '#0A6FA5', 1: '#008FCD'}))
    # Get the total number of bps to the start and end of each chromosome
    cumulative_bp = sumstats\
        .group_by(chrom_col, maintain_order=True)\
        .agg(pl.max(bp_col))\
        .with_columns(end=pl.col(bp_col).cum_sum())\
        .with_columns(start=pl.col.end.shift().fill_null(0))\
        .drop(bp_col)
    # Get x and y coordinates to plot
    sumstats = sumstats\
        .join(cumulative_bp, on=chrom_col, how='left')\
        .with_columns(x=pl.col.start + pl.col(bp_col),
                      y=-pl.col(p_col).log10())
    # Plot horizontal lines indicating significance thresholds
    for p_threshold, (color, linestyle) in p_thresholds.items():
        plt.axhline(y=-np.log10(p_threshold), color=color, linestyle=linestyle,
                    zorder=-1)
    # Overplot three times: first non-clumped variants, then clumped
    # variants, then lead variants. Rasterize non-clumped and (optionally)
    # clumped variants for quick plotting and rendering.
    w, h = plt.gcf().get_size_inches()
    plt.gcf().set_size_inches(1.2 * w, h)
    if clumped_variants is None:
        plt.scatter(sumstats['x'], sumstats['y'], c=sumstats['color'], s=1,
                    rasterized=True)
    else:
        non_clump_sumstats = sumstats.filter(
            ~pl.col(SNP_col).is_in(clumped_variants[SNP_col]))
        clump_sumstats = sumstats.filter(
            pl.col(SNP_col).is_in(clumped_variants[SNP_col]),
            ~pl.col(SNP_col).is_in(clumped_variants[f'{SNP_col}_lead']))
        join_columns = SNP_col, chrom_col, bp_col, ref_col, alt_col
        lead_sumstats = clumped_variants\
            .filter('is_lead')\
            .select(join_columns)\
            .join(sumstats, on=join_columns, how='left')
        plt.scatter(non_clump_sumstats['x'], non_clump_sumstats['y'],
                    c=non_clump_sumstats['color'], s=1, rasterized=True)
        plt.scatter(clump_sumstats['x'], clump_sumstats['y'], c='#DBA756', s=1)
        plt.scatter(lead_sumstats['x'], lead_sumstats['y'],
                    c='#DBA756', marker='D', edgecolors='k', s=10)
        if genome_build is not None:
            lead_sumstats = lead_sumstats\
                .pipe(lambda df: df.join(
                    df.filter(pl.col(p_col) < max_p_to_label)
                    .pipe(get_nearest_gene, genome_build=genome_build,
                          SNP_col=SNP_col, chrom_col=chrom_col, bp_col=bp_col,
                          gene_annotation_dir=gene_annotation_dir),
                    on=(SNP_col, chrom_col, bp_col), how='left'))\
                .drop_nulls('gene')
            for genes, x, y in zip(lead_sumstats['gene'], lead_sumstats['x'],
                                   lead_sumstats['y']):
                genes = ', '.join(genes)
                plt.text(*plt.gca().transData.inverted().transform(
                    plt.gca().transData.transform([x, y]) + 
                    np.array(text_padding_overrides.get(genes, text_padding))),
                    s=text_overrides.get(genes, genes),
                    fontsize=text_size_overrides.get(genes, text_size))
    # Put each chrom's label in the chrom's center; hide chr17/19/21 labels
    plt.xticks((cumulative_bp['start'] + cumulative_bp['end']) / 2,
               cumulative_bp[chrom_col]
               .replace({'17': '', '19': '', '21': ''}))
    # Set x and y labels and limits, resize
    plt.xlabel('Chromosome')
    plt.ylabel('-log$_{10}$(p)')
    padding = 10_000_000  # avoid clipping SNPs at the far left or right
    plt.xlim(sumstats['x'].min() - padding, sumstats['x'].max() + padding)
    plt.ylim(bottom=-np.log10(max_p))


def qqplot(ps, *, equal_aspect=True, **kwargs):
    """
    Generates a quantile-quantile (Q-Q) plot of p-values.
    Inspired by qqplot() from the qmplot package.
    
    Args:
        ps: a polars Series or 1D NumPy array of p-values to plot
        equal_aspect: if True, calls ax.set_aspect('equal') so the x and y axes
                      use the same number of inches per unit increase in x and
                      y
        **kwargs: passed to ax.scatter()
    """
    #
    import matplotlib.pyplot as plt
    import numpy as np
    if ps.is_empty():
        raise ValueError(f'ps is empty!')
    ax = kwargs.pop('ax') if 'ax' in kwargs else plt.gca()
    rasterized = kwargs.pop('rasterized') if 'rasterized' in kwargs else True
    ppoints = lambda n, a=0.5: (np.arange(n) + 1 - a) / (n + 1 - 2 * a)
    ax.scatter(-np.log10(ppoints(len(ps))),
               -np.log10(np.sort(ps).clip(5e-324)),
               edgecolors='none', rasterized=rasterized, **kwargs)
    xmax = ax.get_xlim()[1]
    ax.plot([0, xmax], [0, xmax], c='k', zorder=-1)
    if equal_aspect:
        ax.set_aspect('equal')
    ax.set_xlim(0)
    ax.set_ylim(0)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


###############################################################################
# [5] Statistics
###############################################################################


def standardize(array):
    """
    Standardize a polars Series, DataFrame or expression or a NumPy array to
    zero mean and unit variance (using 1 delta degree of freedom), columnwise.
    
    Args:
        array: the Series, DataFrame, expression or array to standardize

    Returns:
        The standardized Series, DataFrame, expression or array.
    """
    if isinstance(array, (pl.Series, pl.Expr)):
        return (array - array.mean()) / array.std()
    if isinstance(array, pl.DataFrame):
        return array.with_columns((pl.selectors.numeric() -
                                   pl.selectors.numeric().mean()) /
                                  pl.selectors.numeric().std())
    import numpy as np
    if not isinstance(array, np.ndarray):
        raise ValueError(f'array must be a polars Series, DataFrame or '
                         f'expression or a NumPy array!')
    return (array - np.mean(array, axis=0)) / np.std(array, ddof=1, axis=0)


def bonferroni(pvalues):
    """
    Performs Bonferroni correction on a polars Series, DataFrame, expression,
    or NumPy array of p-values.
    
    Args:
        pvalues: a polars Series, DataFrame, expression, or NumPy array of
                 p-values; may contain missing data

    Returns:
        Bonferroni-corrected p-values; same type as pvalues.
    """
    if isinstance(pvalues, (pl.Series, pl.Expr)):
        return (pvalues * (pvalues.len() - pvalues.null_count()))\
            .clip(upper_bound=1)
    import numpy as np
    if not isinstance(pvalues, np.ndarray):
        raise ValueError('pvalues must be a polars Series or expression or a '
                         'NumPy array!')
    return np.minimum(pvalues * (~np.isnan(pvalues).sum()), 1)


def fdr(pvalues):
    """
    Performs FDR correction on a polars Series or expression or 1D NumPy array
    of p-values.
    
    Args:
        pvalues: a polars Series or expression or 1D NumPy array of p-values;
                 may contain missing data

    Returns:
        FDR q-values; same type as pvalues.
    """
    if isinstance(pvalues, (pl.Series, pl.Expr)):
        # Couple of tricky things here:
        # 1) pvalues.arg_sort(descending=True, null_last=True) puts nulls at 
        #    the end. The pl.int_range(...) / num_non_null part gives the same 
        #    result as the np.linspace() call below for the first num_not_null 
        #    elements, but keeps going down for the remaining 
        #    pvalues.null_count() elements. These last pvalues.null_count() 
        #    elements are garbage data, but they're just there to make sure 
        #    pvalues.arg_sort(...) is the same size as pl.int_range(...). Once
        #    you divide pvalues.arg_sort(...) by pl.int_range(...), this 
        #    garbage data goes away because the corresponding elements of
        #    pvalues.arg_sort(descending=True) are null, so the quotient is
        #    null.
        # 2) gather(reverse_order.arg_sort()) is the inverse of gather(
        #    reverse_order): argsorting an argsort inverts the argsort!
        # 3) pvalues.null_count() is an unsigned int for expressions (but not
        #    Series), which leads to overflow if you take the negative of it. 
        #    To avoid this, cast it to pl.Int64.
        eager = isinstance(pvalues, pl.Series)
        num_null = pvalues.null_count() \
            if eager else pvalues.null_count().cast(pl.Int64)
        num_non_null = pvalues.len() - pvalues.null_count()
        reverse_order = pvalues.arg_sort(descending=True, nulls_last=True)
        return (pvalues.gather(reverse_order) /
                (pl.int_range(num_non_null, -num_null, -1, eager=eager) /
                 num_non_null))\
            .cum_min()\
            .gather(reverse_order.arg_sort())
    import numpy as np
    if not isinstance(pvalues, np.ndarray) or pvalues.ndim != 1:
        raise ValueError('pvalues must be a polars Series or expression or a '
                         '1D NumPy array!')
    qvalues = np.empty_like(pvalues)
    missing = np.isnan(pvalues)
    any_missing = missing.any()
    if any_missing:
        non_missing = ~missing
        pvalues = pvalues[non_missing]
        qvalues[missing] = np.nan
        reverse_order = np.argsort(-pvalues)
        # q[non_missing][reverse_order] doesn't work due to chained assignment
        qvalues[np.flatnonzero(non_missing)[reverse_order]] = \
            np.minimum.accumulate(pvalues[reverse_order] / np.linspace(
                1, 1 / len(pvalues), len(pvalues)))
    else:
        reverse_order = np.argsort(-pvalues)
        qvalues[reverse_order] = np.minimum.accumulate(
            pvalues[reverse_order] / np.linspace(
                1, 1 / len(pvalues), len(pvalues)))
    return qvalues


def fisher(table):
    """
    Computes the Fisher p-value from a 2 × 2 contingency table,
    using R's fisher.test function.

    Note that fisher.test differs from scipy.stats.fisher_exact.
    fisher.test gives the conditional maximum likelihood estimate of the OR,
    while fisher_exact gives the unconditional maximum likelihood estimate:
    docs.scipy.org/doc/scipy/reference/generated/scipy.stats.fisher_exact.html
    
    Example from Wikipedia: say you have two gene lists and want to test
    whether they are enriched for each other, i.e. whether being in one list
    is predictive of being in the other list.
    
    The four entries of the contingency table would be:
    1. The number of genes that are in both lists
    2. The number of genes that are in the first list and not the second
    3. The number of genes that are in the second list and not the first
    4. The number of genes that are not in either list
    
    Note that the four entries are non-overlapping and add up to the total
    number of genes.

    Args:
        table: a 2 × 2 Polars DataFrame, NumPy array, list of lists, tuple of
               tuples, etc. of p-values

    Returns:
        A 4-element namedtuple of floats with fields for the odds ratio (OR),
        its lower and upper 95% confidence intervals (lower_CI and upper_CI),
        and the Fisher p-value (p).
    """
    from collections import namedtuple
    import numpy as np
    table = np.asarray(table)
    assert table.shape == (2, 2)
    from ryp import to_py, to_r
    to_r(table, 'table', format='matrix')
    result = to_py('fisher.test(table)')
    OR = result['estimate']
    lower_CI, upper_CI = result['conf.int']
    p = result['p.value']
    FisherResult = namedtuple('FisherResult',
                              ('OR', 'lower_CI', 'upper_CI', 'p'))
    return FisherResult(OR=OR, lower_CI=lower_CI, upper_CI=upper_CI, p=p)


def auprc_score(labels, predictions):
    """
    Computes the area under the precision-recall curve.
    
    Args:
        labels: the ground-truth labels
        predictions: the predicted labels

    Returns:
        The area under the precision-recall curve
    """
    from sklearn.metrics import auc, precision_recall_curve
    precision, recall, thresholds = precision_recall_curve(labels, predictions)
    auprc = auc(recall, precision)
    return auprc


def ACAT(pvalues):
    """
    Calculates the ACAT combined p-value from a vector of p-values.
    
    Args:
        pvalues: a polars Series or NumPy array of p-values; may contain
                 missing data

    Returns:
        The ACAT combined p-value, as a single floating-point number.
    """
    import numpy as np
    from mpmath import mp
    if not (isinstance(pvalues, pl.Series) or (
            isinstance(pvalues, np.ndarray) and pvalues.ndim == 1)):
        raise ValueError('pvalues must be a polars Series or expression '
                         'or a 1D NumPy array!')
    mp.dps = np.ceil(-np.log10(pvalues.min())) + 100  # 100 = fudge factor
    x = np.array([mp.tan((0.5 - mp.mpf(pvalue)) * np.pi)
                  for pvalue in pvalues if pvalue is not None])
    return float(0.5 - mp.atan(x.mean()) / np.pi)


def inflation_factor(pvalues):
    """
    Calculates the genomic inflation factor from a vector of p-values.
    
    Args:
        pvalues: a polars Series or expression or NumPy array of p-values

    Returns:
        A single floating-point number with the genomic inflation factor.
    """
    from scipy.special import chdtri  # chdtri(df, p) == chi2.isf(p, df)
    if isinstance(pvalues, (pl.Series, pl.Expr)):
        return chdtri(1, pvalues.median()) / chdtri(1, 0.5)
    import numpy as np
    if not isinstance(pvalues, np.ndarray):
        raise ValueError('pvalues must be a polars Series or expression or a '
                         'NumPy array!')
    return chdtri(1, np.median(pvalues)) / chdtri(1, 0.5)


def inverse_normal_transform(values, *, c=3 / 8):
    """
    Calculates the rank-based inverse normal transform of a polars Series or
    expression or 1D NumPy array.
    
    Args:
        values: a polars Series or expression or 1D NumPy array.
        c: a parameter of the transformation: a fractional shift to the ranks
           that's applied before transforming. By default, we use the Blom
           transform (c = 3/8). For more details on c, see:
           ncbi.nlm.nih.gov/pmc/articles/PMC2921808/#S2title

    Returns:
        The rank-based inverse normal transform of values.
    """
    from scipy.special import ndtri  # ntdri(x) == norm.ppf(x)
    if isinstance(values, (pl.Series, pl.Expr)):
        rank = values.rank()
        transformed_rank = (rank - c) / \
                           (rank.len() - rank.null_count() - 2 * c + 1)
    else:
        import numpy as np
        from scipy.stats import rankdata
        if not isinstance(values, np.ndarray) or values.ndim != 1:
            raise ValueError('values must be a polars Series or expression or '
                             'a 1D NumPy array!')
        rank = rankdata(values)
        transformed_rank = (rank - c) / ((~np.isnan(rank)).sum() - 2 * c + 1)
    return ndtri(transformed_rank)


@polars_numpy_autoconvert()
def z_to_p(z_scores, *, high_precision=False):
    """
    Converts a polars Series or NumPy array of z-scores to p-values
    
    Args:
        z_scores: the polars Series or NumPy array of z-scores
        high_precision: if True, uses R's pnorm function for high-precision
                        output - important for very small p-values

    Returns:
        The corresponding p-values; same type as z_scores.
    """
    import numpy as np
    if high_precision:
        from ryp import to_py, to_r
        to_r(np.abs(z_scores), 'abs.z')
        pvalues = 2 * np.exp(np.float128(to_py(
            'pnorm(abs.z, lower.tail=False, log.p=True)')))
    else:
        from scipy.special import ndtr
        pvalues = 2 * ndtr(-np.abs(z_scores))  # ndtr(-x) == norm.sf(x)
    return pvalues


def t_to_p(t_scores, df):
    """
    Converts a polars Series or expression or NumPy array of t-scores to
    p-values
    
    Args:
        t_scores: the polars Series or expression or NumPy array of t-scores
        df: the number of degrees of freedom

    Returns:
        The corresponding p-values; same type as t_scores.
    """
    from scipy.special import stdtr
    if isinstance(t_scores, (pl.Series, pl.Expr)):
        return 2 * stdtr(df, -t_scores.abs())
    import numpy as np
    if not isinstance(t_scores, np.ndarray):
        raise ValueError('t_scores must be a polars Series or expression or a '
                         'NumPy array!')
    pvalues = 2 * stdtr(df, -np.abs(t_scores))  # stdtr(df, -x) == t.sf(x, df)
    return pvalues


def p_to_abs_z(pvalues):
    """
    Converts a polars Series or expression or NumPy array of p-values to
    abs(z-scores)
    
    Args:
        pvalues: the polars Series or expression or NumPy array of p-values

    Returns:
        The corresponding abs(z-scores); same type as pvalues.
    """
    from scipy.special import ndtri
    abs_z_scores = -ndtri(pvalues / 2)  # -ndtri(x) == norm.isf(x)
    return abs_z_scores


@polars_numpy_autoconvert(use_columns_from=['Y', 'X'])
def cov(X, Y=None, *, rowvar=False, ddof=1):
    """
    np.cov() but with proper support for Y, and with rowvar=False by default.

    Calculates the covariance between each pair of columns of X, or (if Y is
    not None), each column of X with each column of Y.

    Args:
        X: a polars DataFrame or 2D NumPy array for which to calculate the
           covariance, either with Y (if not None) or with itself
        Y: a polars DataFrame or 2D NumPy array for which to calculate the
           covariance with X
        rowvar: if True, calculates covariance across rows instead of down
                columns
        ddof: the number of delta degrees of freedom

    Returns:
        A covariance matrix; same type as X.
    """
    if Y is not None:
        assert type(X) is type(Y), \
            f'type(X) = {type(X).__name__}, type(Y) = {type(Y).__name__}'
        assert len(X.shape) == len(Y.shape) == 2, \
            f'X.shape = {X.shape:,}, Y.shape = {Y.shape:,}'
        assert X.shape[rowvar] == Y.shape[rowvar], \
            f'X.shape = {X.shape:,}, Y.shape = {Y.shape:,}'
    # a) Transpose if rowvar=False
    if not rowvar:
        X = X.T
        if Y is not None:
            Y = Y.T
    # b) De-mean
    X = X - X.mean(axis=1)[:, None]
    if Y is not None:
        Y = Y - Y.mean(axis=1)[:, None]
    # c) Calculate covariance
    if Y is None:
        c = (X @ X.T) / (X.shape[1] - ddof)
    else:
        c = (X @ Y.T) / (X.shape[1] - ddof)
    return c


@polars_numpy_autoconvert(use_columns_from=['Y', 'X'])
def cor(X, Y=None, *, rowvar=False, return_p=False, ddof=1):
    """
    np.corrcoef() but with proper support for Y, the ability to return
    p-values with return_p=True, and with rowvar=False by default.

    Correlates each pair of columns of X, or (if Y is not None), each column
    of X with each column of Y.

    Args:
        X: a polars DataFrame or 2D NumPy array to correlate, either with Y
           (if not None) or with itself
        Y: a polars DataFrame or 2D NumPy array to correlate with X
        rowvar: if True, correlates across rows instead of down columns
        return_p: if True, returns the p-value as well as the correlation
        ddof: the number of delta degrees of freedom

    Returns:
        The correlation matrix, or if return_p=True the correlation and p-value
        matrices as a two-element tuple. Returned matrices have the same type
        as X.
    """
    import numpy as np
    if Y is not None:
        assert type(X) is type(Y), \
            f'type(X) = {type(X).__name__}, type(Y) = {type(Y).__name__}'
        assert len(X.shape) == len(Y.shape) == 2, \
            f'X.shape = {X.shape:,}, Y.shape = {Y.shape:,}'
        assert X.shape[rowvar] == Y.shape[rowvar], \
            f'X.shape = {X.shape:,}, Y.shape = {Y.shape:,}'
    # a) Calculate covariance
    c = cov(X, Y, rowvar=rowvar, ddof=ddof)
    # b) Convert to correlation
    if Y is None:
        stddev = np.sqrt(np.diag(c))
        c /= stddev[:, None]
        c /= stddev[None, :]
    else:
        c /= X.std(axis=int(rowvar), ddof=ddof)[:, None]
        c /= Y.std(axis=int(rowvar), ddof=ddof)
    # c) Clip to [-1, 1]
    np.clip(c, -1, 1, out=c)
    # d) Calculate p-value if return_p=True.
    # Use np.errstate(divide='ignore') to ignore division-by-zero errors when
    # calculating p-values for the correlation of each column (or row, if
    # rowvar=True) of X with itself, in the case where Y is None.
    if return_p:
        from scipy.special import stdtr
        df = X.shape[rowvar] - 2
        with np.errstate(divide='ignore'):
            # stdtr(df, -x) == t.sf(x, df)
            p = 2 * stdtr(df, -np.abs(c * np.sqrt(df / (1 - c * c))))
        return c, p
    else:
        return c


###############################################################################
# [6] Matrix operations
###############################################################################

def upper_triangle(array, *, include_diagonal=True):
    """
    Returns the upper triangle of a square 2D NumPy array as a flat 1D array.
    If the NumPy array has shape (k, k), then the returned array will have
    length k(k+1)/2 if include_diagonal=True or k(k-1)/2 otherwise.
    
    Also accepts inputs of N > 2 dimensions, in which case the upper triangle
    is flattened along the last two axes, returning an (N-1)-dimensional array.

    For example:
    - upper_triangle(np.ones((100, 100))).shape == (100 * 99 // 2 + 100,)
                                                == (5050,)
    - upper_triangle(np.ones((70, 100, 100))).shape == (70, 5050)
    - upper_triangle(np.ones((30, 70, 100, 100))).shape == (30, 70, 5050)
    If include_diagonal=False, the last dimension would be 4950 not 5050
    
    Args:
        array: an N-dimensional NumPy array
        include_diagonal: if True, includes the diagonal along with the upper
                          triangular elements

    Returns:
        The upper triangle along the last two dimensions.
    """
    import numpy as np
    assert isinstance(array, np.ndarray), type(array)
    assert array.ndim >= 2, array.ndim
    assert array.shape[-2] == array.shape[-1], array.shape
    upper_triangle_mask = np.triu(np.ones(array.shape[-2:], dtype=bool),
                                  k=0 if include_diagonal else 1)
    return array[..., upper_triangle_mask]


def upper_triangle_to_square(array, *, include_diagonal=True, mirror=False):
    """
    Inverts the operation of upper_triangle(), converting a 1D NumPy array of
    upper triangular elements to a 2D array where the elements occupy the upper
    triangle (higher dimensions are not supported).
    
    Args:
        array: a 1D NumPy array of upper diagonal elements
        include_diagonal: if True, array is assumed to contain the diagonal
                          elements as well as the upper triangular elements;
                          if False, diagonal elements will be set to nan
        mirror: whether to mirror the result to the lower triangle; if False,
                lower triangular elements will be set to nan

    Returns:
        A 2D array where the elements of array occupy the upper triangle (and
        diagonal if include_diagonal=True, and lower triangle if mirror=True).
    """
    import numpy as np
    assert isinstance(array, np.ndarray), type(array)
    assert array.ndim == 1, array.ndim
    # Infer the size of the square matrix we're trying to reconstruct
    # If U is the size of the upper triangle array we got as input, and N was
    # the size of the original square matrix, then:
    # 1) include_diagonal=True: U = N(N - 1) / 2 + N = 0.5N^2 + 0.5N
    #    --> by the quadratic formula, N = (sqrt(8U + 1) - 1) / 2
    # 2) include_diagonal=False: U = N(N - 1) / 2 = 0.5N^2 - 0.5N
    #    --> by the quadratic formula, N = (sqrt(8U + 1) + 1) / 2
    N = (np.sqrt(8 * len(array) + 1) + (-1 if include_diagonal else 1)) / 2
    assert N == int(N), f'Array is the wrong length ({len(array):,}) to be ' \
                        f'the upper diagonal portion of a square matrix'
    N = int(N)
    square = np.full((N, N), np.nan)
    row_indices, col_indices = np.triu_indices_from(
        square, k=0 if include_diagonal else 1)
    square[row_indices, col_indices] = array
    if mirror:
        square[col_indices, row_indices] = array
    return square


def estimate_covariance(X):
    """
    Estimate the covariance matrix of X's observations via the Oracle
    Approximating Shrinkage (OAS) method. For details on the method, see:
    scikit-learn.org/stable/modules/generated/sklearn.covariance.OAS.html
    
    OAS was chosen because it is much faster than Scikit-learn's implementation
    of Ledoit-Wolf covariance estimation (sklearn.covariance.ledoit_wolf) or
    R's corpcor(), e.g. 100s vs 540s vs 1702s for np.random.RandomState(0)
    .random(size=(10000, 10000)). Which you choose has very little effect on
    downstream analyses, but any is more robust than using the raw covariance.
    
    Args:
        X: a polars DataFrame or 2D NumPy array (rows = observations, columns =
           features) to estimate the covariance of
    
    Returns:
        The estimated covariance matrix of X's observations, of shape
        (observations × observations).
    """
    from sklearn.covariance import oas
    covariance = oas(X.T)[0]
    return covariance


@polars_numpy_autoconvert()
def covariance_to_whitening_matrix(covariance, method='Cholesky'):
    """
    Calculates a whitening matrix from a covariance matrix.
    
    Based on rdrr.io/cran/whitening/src/R/getPhiPsiW.R.
    
    Args:
        covariance: a covariance matrix (2D NumPy array or polars DataFrame)
        method: one of 'Cholesky', 'ZCA', 'ZCA-cor', 'PCA', 'PCA-cor'

    Returns:
        The corresponding whitening matrix.
    """
    import numpy as np
    if not isinstance(covariance, np.ndarray) or covariance.ndim != 2:
        raise ValueError('covariance must be a 2D NumPy array or polars '
                         'DataFrame!')
    if method not in ('Cholesky', 'ZCA', 'ZCA-cor', 'PCA', 'PCA-cor'):
        raise ValueError(f"Unknown whitening method {method} (valid methods: "
                         f"'Cholesky', 'ZCA', 'ZCA-cor', 'PCA', 'PCA-cor')")
    if method == 'Cholesky':
        # dtrtri(chol, lower=True) can give an order of magnitude speed-up over
        # inv(chol) by exploiting the fact that the Cholesky matrix is lower
        # triangular.
        from scipy.linalg import cholesky
        from scipy.linalg.lapack import dtrtri
        whitening_matrix, info = dtrtri(cholesky(
            covariance, lower=True, check_finite=False),
            lower=True, overwrite_c=True)
        if info > 0:
            raise np.linalg.LinAlgError('Singular matrix')
        elif info < 0:
            raise ValueError(f'Invalid input to dtrtri (info = {info})')
    else:
        # Calculate eigenvalues and eigenvectors of the covariance (for ZCA and
        # PCA) or correlation (for ZCA-cor and PCA-cor). Use scipy.linalg.eigh
        # since np.linalg.eigh will rarely, and incorrectly, return all-zero
        # eigenvalues - possibly because the eigenvalue computation runs out of
        # memory?
        from scipy.linalg import eigh
        if method == 'ZCA' or method == 'PCA':
            eigenvalues, eigenvectors = eigh(covariance, check_finite=False)
        else:
            stddev = np.sqrt(np.diag(covariance))
            correlation = covariance / np.outer(stddev, stddev)
            eigenvalues, eigenvectors = eigh(correlation, check_finite=False)
        # For ZCA, the whitening matrix is U * S^-0.5 * V; for PCA, it's just
        # S^-0.5 * V. For PCA and PCA-cor, the eigenvectors' sign ambiguity
        # affects the final result, so multiply eigenvectors by -1 as needed to
        # ensure that the eigenvector matrix has a positive diagonal.
        if method == 'ZCA' or method == 'ZCA-cor':
            whitening_matrix = eigenvectors @ \
                               np.diag(1 / np.sqrt(eigenvalues)) @ \
                               eigenvectors.T
        else:
            makePosDiag = lambda matrix: matrix * np.sign(np.diag(matrix))
            whitening_matrix = np.diag(1 / np.sqrt(eigenvalues)) @ \
                makePosDiag(eigenvectors).T
        # For ZCA-cor and PCA-cor, right-multiply by a diagonal matrix of the
        # inverse standard deviations
        if method == 'ZCA-cor' or method == 'PCA-cor':
            whitening_matrix = whitening_matrix @ \
                               np.diag(1 / np.sqrt(np.diag(covariance)))
    return whitening_matrix


@polars_numpy_autoconvert()
def get_whitening_matrix(X, method='Cholesky'):
    """
    Given a polars DataFrame or 2D NumPy array (rows = observations, columns =
    features), estimates the covariance of X's observations, and uses this
    covariance to compute a whitening matrix according to the specified
    whitening method.
    
    Based on rdrr.io/cran/whitening/src/R/getPhiPsiW.R.
    
    Args:
        X: a polars DataFrame or 2D NumPy array
        method: one of 'Cholesky', 'ZCA', 'ZCA-cor', 'PCA', 'PCA-cor'

    Returns:
        The whitened version of X.
    """
    import numpy as np
    if not isinstance(X, np.ndarray) or X.ndim != 2:
        raise ValueError('X must be a polars DataFrame or 2D NumPy array!')
    if method not in ('Cholesky', 'ZCA', 'ZCA-cor', 'PCA', 'PCA-cor'):
        raise ValueError(f"Unknown whitening method {method} (valid methods: "
                         f"'Cholesky', 'ZCA', 'ZCA-cor', 'PCA', 'PCA-cor')")
    # Compute covariance matrix
    covariance = estimate_covariance(X)
    # Compute whitening matrix from the covariance matrix
    whitening_matrix = covariance_to_whitening_matrix(covariance, method)
    return whitening_matrix


@polars_numpy_autoconvert()
def whiten(X, method='Cholesky'):
    """
    Given a polars DataFrame or 2D NumPy array (rows = observations, columns =
    features), estimates the covariance of X's observations, uses this
    covariance to compute a whitening matrix according to the specified
    whitening method, and whitens the data.
    
    Based on rdrr.io/cran/whitening/src/R/getPhiPsiW.R.
    
    Args:
        X: a polars DataFrame or 2D NumPy array
        method: one of 'Cholesky', 'ZCA', 'ZCA-cor', 'PCA', 'PCA-cor'

    Returns:
        The whitened version of X.
    """
    whitened_X = get_whitening_matrix(X, method=method) @ X
    return whitened_X


def sparse_matrix_vector_op(sparse_matrix, op, vector, *, axis, inplace=False,
                            return_dtype=None, upcast_float32=True,
                            num_threads=1):
    """
    Apply an elementwise arithmetic operation to a sparse array/matrix and a
    vector, either rowwise or columnwise.
    
    Args:
        sparse_matrix: a CSR or CSC sparse array or matrix
        op: the operation to perform, as a string; must be one of:
            '+', '-', '*', '/', '//', '**', '%'
        vector: a 1D numeric vector
        axis: the axis to perform the operation across: columnwise (`axis=0`)
              or rowwise (`axis=1`)
        inplace: whether to perform the operation in place; will give an error
                 if `sparse_matrix` has integer dtype and the result of the
                 operation would be floating-point
        return_dtype: the data type of the output matrix; if None, will be
                      inferred automatically
        upcast_float32: whether to automatically upcast the output to float64
                        when it would otherwise be float32 (i.e. when
                        `sparse_matrix` and/or 'vector' are float32 and op is
                        not `'//'`) and inplace=False; has no effect otherwise.
                        Overridden by `dtype`, if specified.
        num_threads: the number of threads to use for the operation. Set
                     `num_threads=None` to use all available cores (as
                     determined by `os.cpu_count()`).
    
    Returns:
        A sparse matrix with the result of the operation, or None if
        `inplace=True`.
    """
    import numpy as np
    from scipy.sparse import csc_array, csr_array, csc_matrix, csr_matrix
    # Check inputs
    check_type(sparse_matrix, 'sparse_matrix',
               (csc_array, csr_array, csc_matrix, csr_matrix),
               'a sparse matrix')
    valid_ops = '+', '-', '*', '/', '//', '**', '%'
    if op not in valid_ops:
        check_type(op, 'op', str, 'a string representing a Python operator')
        error_message = \
            f'op must be one of these operators: {", ".join(valid_ops)}'
        raise ValueError(error_message)
    check_type(vector, 'vector', np.ndarray, 'a NumPy array')
    if vector.ndim != 1:
        error_message = \
            f'vector must be 1-dimensional, not {vector.ndim:,}-dimensional'
        raise ValueError(error_message)
    if axis != 0 and axis != 1:
        check_type(axis, 'axis', int, 'one of the integers 0 or 1')
        error_message = f'axis must be 0 or 1, not {axis}'
        raise ValueError(error_message)
    check_type(inplace, 'inplace', bool, 'boolean')
    try:
        np.dtype(return_dtype)
    except TypeError as e:
        error_message = 'return_dtype is not a valid NumPy data type'
        raise TypeError(error_message) from e
    check_type(upcast_float32, 'upcast_float32', bool, 'boolean')
    if num_threads is None:
        num_threads = os.cpu_count()
    else:
        check_type(num_threads, 'num_threads', int, 'a positive integer')
        check_bounds(num_threads, 'num_threads', 1)
    # Check dimensions
    if len(vector) != sparse_matrix.shape[axis]:
        error_message = (
            f'shapes do not match: len(vector) is {len(vector):,}, but '
            f'sparse_matrix.shape[{axis}] is {sparse_matrix.shape[axis]:,}')
        if len(vector) == sparse_matrix.shape[1 - axis]:
            error_message += (
                f'; however, sparse_matrix.shape[{1 - axis}] is '
                f'{sparse_matrix.shape[1 - axis]:,}, so perhaps you meant to '
                f'specify axis={1 - axis} instead of axis={axis}?')
        raise ValueError(error_message)
    # Get the result's data type, if `return_dtype` is not specified
    # (note: // always returns an integer dtype)
    if return_dtype is None:
        return_dtype = np.result_type(sparse_matrix.dtype, vector.dtype)
        if op == '//' and np.issubdtype(return_dtype, np.floating):
            return_dtype = np.int64
        elif return_dtype == np.float32 and not inplace and upcast_float32:
            return_dtype = np.float64
    # Get `sparse_matrix`'s data, indices and indptr
    data = sparse_matrix.data
    indices = sparse_matrix.indices
    indptr = sparse_matrix.indptr
    prange_import = \
        'from cython.parallel cimport prange' if num_threads > 1 else ''
    # If in-place...
    if inplace:
        # Raise an error if the result's dtype would be different from
        # `sparse_matrix`'s
        if inplace and return_dtype != data.dtype:
            error_message = (
                f'sparse_matrix has data type {data.dtype!r}, but '
                f'`sparse_matrix {op} vector` would have data type '
                f'{return_dtype!r}, so inplace=True cannot be specified')
            raise TypeError(error_message)
        # Run the operation with Cython
        if isinstance(sparse_matrix, (csc_array, csc_matrix)) == axis:
            # Rowwise for CSR, or columnwise for CSC
            cython_inline(f'''
                {prange_import}
                def run_op_inplace(
                        {cython_type(data.dtype)}[::1] data,
                        const {cython_type(indices.dtype)}[::1] indices,
                        const {cython_type(indptr.dtype)}[::1] indptr,
                        const {cython_type(vector.dtype)}[::1] vector,
                        const unsigned int num_threads):
                    cdef Py_ssize_t i, j
                    for j in {prange('indptr.shape[0] - 1', num_threads)}:
                        for i in range(indptr[j], indptr[j + 1]):
                            data[i] {op}= vector[j]
                ''')['run_op_inplace'](
                    data, indices, indptr, vector, num_threads)
        else:
            # Rowwise for CSC, or columnwise for CSR
            cython_inline(f'''
                {prange_import}
                def run_op_inplace(
                        {cython_type(data.dtype)}[::1] data,
                        const {cython_type(indices.dtype)}[::1] indices,
                        const {cython_type(indptr.dtype)}[::1] indptr,
                        const {cython_type(vector.dtype)}[::1] vector,
                        const unsigned int num_threads):
                    cdef Py_ssize_t i
                    for i in {prange('data.shape[0]', num_threads)}:
                        data[i] {op}= vector[indices[i]]
                ''')['run_op_inplace'](
                    data, indices, indptr, vector, num_threads)
    # If not in-place...
    else:
        # Allocate the output sparse matrix/array's data array
        result = np.empty_like(sparse_matrix.data, dtype=return_dtype)
        # Get the Cython type of the result
        result_cython_type = cython_type(return_dtype)
        # Run the operation with Cython
        if isinstance(sparse_matrix, (csc_array, csc_matrix)) == axis:
            # Rowwise for CSR, or columnwise for CSC
            cython_inline(rf'''
                {prange_import}
                def run_op(
                        {cython_type(data.dtype)}[::1] data,
                        const {cython_type(indices.dtype)}[::1] indices,
                        const {cython_type(indptr.dtype)}[::1] indptr,
                        const {cython_type(vector.dtype)}[::1] vector,
                        {result_cython_type}[::1] result,
                        const unsigned int num_threads):
                    cdef Py_ssize_t i, j
                    for j in {prange('indptr.shape[0] - 1', num_threads)}:
                        for i in range(indptr[j], indptr[j + 1]):
                            result[i] = \
                                <{result_cython_type}> data[i] {op} vector[j]
                ''')['run_op'](
                    data, indices, indptr, vector, result, num_threads)
        else:
            # Rowwise for CSC, or columnwise for CSR
            cython_inline(rf'''
                {prange_import}
                def run_op(
                        const {cython_type(data.dtype)}[::1] data,
                        const {cython_type(indices.dtype)}[::1] indices,
                        const {cython_type(indptr.dtype)}[::1] indptr,
                        const {cython_type(vector.dtype)}[::1] vector,
                        {result_cython_type}[::1] result,
                        const unsigned int num_threads):
                    cdef Py_ssize_t i
                    for i in {prange('data.shape[0]', num_threads)}:
                        result[i] = <{result_cython_type}> data[i] {op} \
                            vector[indices[i]]
                ''')['run_op'](
                    data, indices, indptr, vector, result, num_threads)
        # Return a new sparse matrix/array with the result
        return type(sparse_matrix)((result, indices, indptr),
                                   shape=sparse_matrix.shape)

    
def bincount(x, *, num_bins, num_threads=1):
    """
    A faster version of `numpy.bincount` that uses a fixed number of bins,
    lacks weights and supports multithreading.
    
    Args:
        x: a 1D NumPy array
        num_bins: the number of bins
        num_threads: the number of threads to use when counting

    Returns:
        A 1D NumPy array of length `num_bins`, containing the bin counts.
    """
    import numpy as np
    dtype = x.dtype
    cython_dtype = cython_type(dtype)
    counts = np.zeros(num_bins, dtype=np.int32)
    prange_import = \
        'from cython.parallel cimport prange' if num_threads > 1 else ''
    cython_inline(f'''
        {prange_import}
        def bincount(const {cython_dtype}[::1] arr,
                     int[::1] counts,
                     const unsigned num_threads):
            cdef long i
            for i in {prange('arr.shape[0]', num_threads=num_threads)}:
                counts[arr[i]] += 1
        ''')['bincount'](arr=x, counts=counts, num_threads=num_threads)
    return counts


def getnnz(sparse_array, axis=None, num_threads=1):
    """
    Count the number of stored values in a sparse array, including explicitly
    stored zeros. Matches the behavior of the now-deprecated getnnz() function
    for scipy sparse arrays. To count the number of non-zero elements (i.e.
    excluding explicitly stored zeros), use sparse_array.count_nonzero().
    
    Args:
        sparse_array: a csr_array or csc_array
        axis: whether to count the number of stored values in the whole array
              (axis=None), within each column (axis=0), or within each row
              (axis=1)
        num_threads: the number of threads to use when counting; only used for
                     csr when `axis` is 0 and for csc when `axis` is 1
    
    Returns:
        The number of stored values, either as a scalar (if axis=None) or a 1D
        array (if axis=0 or axis=1).
    """
    import numpy as np
    from scipy.sparse import csr_array, csc_array
    check_type(sparse_array, 'sparse_array', (csr_array, csc_array),
               'a csr_array or csc_array')
    if axis is None:
        return sparse_array.nnz
    elif axis != 0 and axis != 1:
        error_message = 'axis must be 0, 1, or None'
        raise ValueError(error_message)
    is_csr = isinstance(sparse_array, csr_array)
    if axis == is_csr:
        return np.diff(sparse_array.indptr)
    else:
        return bincount(sparse_array.indices,
                        num_bins=sparse_array.shape[is_csr],
                        num_threads=num_threads)

    
###############################################################################
# [7] Regression
###############################################################################


def linear_regression(X, y, *, add_intercept=True, report_intercept=False,
                      return_significance=True, return_residuals=False):
    """
    A fast linear regression implementation that calculate CIs and p-values.
    
    Based on the regression implementations of pingouin at
    github.com/raphaelvallat/pingouin/blob/master/pingouin/regression.py
    
    Args:
        X: the feature matrix (rows are observations, columns are features);
           may be a polars DataFrame or 2D NumPy array; must be numeric
        y: the outcome vector (one element per observation); may be a polars
           Series or 1D NumPy array; must be numeric
        add_intercept: if True, adds a column of ones to X before running the
                       regression; if False, assumes the user has added one
        report_intercept: if True, reports the intercept column added by
                          add_intercept in the returned results.
                          Can only be True if add_intercept=True.
        return_significance: if True, also returns the SEs, CIs and p-values,
                             not just the betas
        return_residuals: if True, also returns the residuals
    Returns:
        A polars DataFrame or NumPy array (whichever type X is) with one
        row per feature (column of X) and the following columns:
        - 'Feature': the names of the columns of X, only added when X is a
                     polars DataFrame
        - 'beta': the regression coefficients
        - 'SE': the standard errors of the regression coefficients
        - 'lower_CI' and 'upper_CI': the 95% confidence intervals of the
          regression coefficients
        - 'p': the p-values of the regression coefficients.
        If return_residuals=True, also returns a polars Series or 1D NumPy
        array (depending on the type of X) of length len(X) with the residuals
        for each of the observations.
    """
    # Tested against statsmodels:
    # >>> import numpy as np, statsmodels.api as sm; from timeit import timeit
    # >>> rng = np.random.default_rng(0)
    # >>> X = rng.standard_normal(size=(100000, 1000))
    # >>> y = rng.standard_normal(size=(100000))
    # >>> timeit(lambda: linear_regression(X, y), number=1)
    # 1.6827743009780534
    # >>> timeit(lambda: sm.OLS(y, sm.add_constant(X)).fit(), number=1)
    # 4.280662438017316
    # >>> from scipy.stats import pearsonr
    # >>> pearsonr(linear_regression(X, y)['p'], sm.OLS(y, X).fit().pvalues)
    # PearsonRResult(statistic=0.9998688713472692, pvalue=0.0)
    return GLM(X, y, add_intercept=add_intercept,
               report_intercept=report_intercept,
               return_significance=return_significance,
               return_residuals=return_residuals)


def logistic_regression(X, y, *, add_intercept=True, report_intercept=False,
                        return_significance=True, return_residuals=False,
                        max_iter=100, num_threads=1):
    """
    A fast logistic regression implementation that calculate CIs and p-values.
    
    Based on the regression implementations of pingouin at
    github.com/raphaelvallat/pingouin/blob/master/pingouin/regression.py
    
    Args:
        X: the feature matrix (observations are rows, features are columns);
           may be a polars DataFrame or 2D NumPy array; must be numeric
        y: the outcome vector (one element per observation); may be a polars
           Series or 1D NumPy array; must be boolean (dtype=bool)
        add_intercept: if True, adds a column of ones to X before running the
                       regression; if False, assumes the user has added one
        report_intercept: if True, reports the intercept column added by
                          add_intercept in the returned results.
                          Can only be True if add_intercept=True.
        return_significance: if True, also returns the SEs, CIs and p-values,
                             not just the betas
        return_residuals: if True, also returns the residuals
        max_iter: the maximum number of iterations for logistic regression
        num_threads: the number of threads used to compute losses and gradients
                     for logistic regression; set `num_threads=None` to use all
                     available cores (as determined by `os.cpu_count()`)
    Returns:
        A polars DataFrame or NumPy array (whichever type X is) with one
        row per feature (column of X) and the following columns:
        - 'Feature': the names of the columns of X, only added when X is a
                     polars DataFrame
        - 'OR': the odds ratios, i.e. exp(regression coefficients)
        - 'SE': the standard errors of the regression coefficients
        - 'lower_CI' and 'upper_CI': the 95% confidence intervals of the odds
                                     ratios
        - 'p': the p-values of the odds ratios.
        If return_residuals=True, also returns a polars Series or 1D NumPy
        array (depending on the type of X) of length len(X) with the residuals
        for each of the observations.
    """
    # Tested against statsmodels:
    # >>> import numpy as np, statsmodels.api as sm; from timeit import timeit
    # >>> rng = np.random.default_rng(0)
    # >>> X = rng.standard_normal(size=(100000, 1000))
    # >>> y = rng.binomial(1, 0.5, size=(100000)).astype(bool)
    # >>> timeit(lambda: logistic_regression(X, y), number=1)
    # 1.9138163239695132
    # >>> timeit(lambda: sm.Logit(y, sm.add_constant(X)).fit(method='lbfgs'),
    #            number=1)
    # 4.680954776005819
    # >>> timeit(lambda: sm.Logit(y, sm.add_constant(X)).fit(), number=1)
    # 6.372486169973854
    # >>> from scipy.stats import pearsonr
    # >>> pearsonr(logistic_regression(X, y)['p'],
    #              sm.Logit(y, sm.add_constant(X)).fit().pvalues[1:])
    # PearsonRResult(statistic=0.9999999997299595, pvalue=0.0)
    if num_threads is None:
        num_threads = os.cpu_count()
    return GLM(X, y, logistic=True, add_intercept=add_intercept,
               report_intercept=report_intercept,
               return_significance=return_significance,
               return_residuals=return_residuals, max_iter=max_iter,
               num_threads=num_threads)


def GLM(X, y, *, logistic=False, add_intercept=True, report_intercept=False,
        return_significance=True, return_residuals=False, max_iter=100,
        num_threads=1):
    """
    A fast linear and logistic regression implementation that calculate CIs and
    p-values.
    
    Based on the regression implementations of pingouin at
    github.com/raphaelvallat/pingouin/blob/master/pingouin/regression.py
    
    Args:
        X: the feature matrix (observations are rows, features are columns);
           may be a polars DataFrame or 2D NumPy array; must be numeric
        y: the outcome vector (one element per observation); may be a polars
           Series or 1D NumPy array; must be numeric for linear regression and
           boolean (dtype=bool) for logistic regression
        logistic: if True, do logistic regression; if False, do linear
        add_intercept: if True, adds a column of ones to X before running the
                       regression; if False, assumes the user has added one
        report_intercept: if True, reports the intercept column added by
                          add_intercept in the returned results.
                          Can only be True if add_intercept=True.
        return_significance: if True, also returns the SEs, CIs and p-values,
                             not just the betas
        return_residuals: if True, also returns the residuals
        max_iter: the maximum number of iterations for logistic regression
        num_threads: the number of threads used to compute losses and gradients
                     for logistic regression; set `num_threads=None` to use all
                     available cores (as determined by `os.cpu_count()`)
    Returns:
        A polars DataFrame with one row per feature (column of X) and the
        following columns:
        - 'Feature': the names of the columns of X, only added when X is a
                     polars DataFrame
        - 'beta' (linear regression only): the regression coefficients
        - 'OR' (logistic regression only): the odds ratios, i.e. exp(beta)
        - 'SE': the standard errors of the regression coefficients
        - 'lower_CI' and 'upper_CI': the 95% confidence intervals of the
          regression coefficients (for linear regression) or odds ratios (for
          logistic regression)
        - 'p': the p-values of the regression coefficients/odds ratios.
        If return_residuals=True, also returns a polars Series or 1D NumPy
        array (depending on the type of X) of length len(X) with the residuals
        for each of the observations.
    """
    import numpy as np
    from scipy.special import ndtr, ndtri, stdtr
    assert isinstance(X, (np.ndarray, pl.DataFrame)), type(X)
    assert isinstance(y, (np.ndarray, pl.Series, pl.DataFrame)), type(y)
    feature_names = None
    if not isinstance(X, np.ndarray):
        feature_names = X.columns
        X = X.to_numpy()
    if not isinstance(y, np.ndarray):
        y = y.to_numpy()
    assert np.issubdtype(X.dtype, np.number), X.dtypes
    if logistic:
        assert y.dtype == bool, y.dtype
    else:
        assert np.issubdtype(y.dtype, np.number), y.dtype
    assert X.ndim == 2, X.ndim
    assert y.ndim == 1, y.ndim
    assert len(X) == len(y), (X.shape, y.shape)
    assert not np.isnan(X).any()
    assert not np.isnan(y).any()
    assert np.linalg.matrix_rank(X) == min(X.shape), \
        f'np.linalg.matrix_rank(X) ({np.linalg.matrix_rank(X)}) < ' \
        f'min(X.shape) ({min(X.shape)})'
    if report_intercept and not add_intercept:
        raise ValueError("add_intercept is False, so report_intercept must "
                         "also be False since there's no intercept to report!")
    if add_intercept:
        X = np.column_stack((np.ones(len(X), dtype=X.dtype), X))
    if logistic:
        from scipy.optimize import minimize
        from sklearn._loss._loss import CyHalfBinomialLoss
        base_loss_and_gradient = CyHalfBinomialLoss().loss_gradient
        loss_out = np.empty(len(y))
        gradient_out = np.empty(len(y))
        y = y.astype(float)
        if num_threads is None:
            num_threads = os.cpu_count()
        
        def loss_and_gradient(beta, X, y):
            loss, grad_pointwise = base_loss_and_gradient(
                y_true=y, raw_prediction=X @ beta, sample_weight=None,
                loss_out=loss_out, gradient_out=gradient_out,
                n_threads=num_threads)
            # This was `return loss.sum(), X.T @ grad_pointwise` before scikit-
            # learn 1.4's github.com/scikit-learn/scikit-learn/pull/26721
            return loss.mean(), X.T @ (grad_pointwise / len(y))
        
        minimize_result = minimize(loss_and_gradient, x0=np.zeros(X.shape[1]),
                                   method='L-BFGS-B', jac=True, args=(X, y),
                                   options=dict(gtol=1e-4, maxiter=max_iter))
        if not minimize_result.success:
            raise RuntimeError(f'Logistic regression failed to converge in '
                               f'{max_iter} iterations; increase max_iter or '
                               f'scale your data to zero mean and unit '
                               f'variance')
        beta = minimize_result.x
        if return_significance:
            denom = np.tile(2 * (1 + np.cosh(X @ beta)), (X.shape[1], 1)).T
            SE = np.sqrt(np.diag(np.linalg.pinv((X / denom).T @ X)))
            # Use ndtr as a ufunc alternative to norm.sf:
            # ndtr(-x) == norm.sf(x)
            p = 2 * ndtr(-np.abs(beta / SE))
    else:
        beta, residuals, rank, singular_values = lstsq(X, y)
        if return_significance:
            df = len(X) - rank
            if df == 0:
                assert rank == len(X) == X.shape[1]
                raise ValueError('X is a full-rank square matrix, so the line '
                                 'perfectly fits the data and p-values are '
                                 'not defined; re-run with return_significance'
                                 '=False if you just want the betas')
            SE = np.sqrt(np.diag(np.linalg.pinv(X.T @ X) * residuals / df))
            # Use stdtr as a ufunc alternative to t.sf:
            # stdtr(df, -x) == t.sf(x, df)
            p = 2 * stdtr(df, -np.abs(beta / SE))
    if return_significance:
        # noinspection PyUnboundLocalVariable
        CI_width = -ndtri(0.025) * SE  # -ndtri(x) == norm.isf(x)
        lower_CI = beta - CI_width
        upper_CI = beta + CI_width
    if logistic:
        beta = np.exp(beta)
        if return_significance:
            # noinspection PyUnboundLocalVariable
            lower_CI = np.exp(lower_CI)
            # noinspection PyUnboundLocalVariable
            upper_CI = np.exp(upper_CI)
    if return_residuals:
        predictions = X @ beta
        if logistic:
            predictions = np.where(  # numerically stable sigmoid
                predictions >= 0, 1 / (1 + np.exp(-predictions)),
                np.exp(predictions) / (1 + np.exp(predictions)))
        residuals = y - predictions
        if feature_names is not None:
            residuals = pl.Series(residuals)
    if report_intercept:
        if feature_names is not None:
            feature_names = ['Intercept'] + feature_names
    elif add_intercept:
        if return_significance:
            # noinspection PyUnboundLocalVariable
            beta, SE, lower_CI, upper_CI, p = \
                beta[1:], SE[1:], lower_CI[1:], upper_CI[1:], p[1:]
        else:
            beta = beta[1:]
    # noinspection PyUnboundLocalVariable
    results = pl.DataFrame(
        ({'Feature': feature_names} if feature_names is not None else {}) |
        {'OR' if logistic else 'beta': beta} |
        ({'SE': SE, 'lower_CI': lower_CI, 'upper_CI': upper_CI, 'p': p}
         if return_significance else {}))
    # noinspection PyUnboundLocalVariable
    return (results, residuals) if return_residuals else results


@polars_numpy_autoconvert(use_columns_from='Y')
def linear_regressions(X, Y, *, add_intercept=True, report_intercept=False,
                       covariance=None, whitening_matrix=None,
                       return_significance=False, return_residuals=False):
    """
    Performs Y.shape[1] independent linear regressions; for each regression,
    one column of Y is regressed on all the features in X.
    
    Based on the regression implementations of pingouin at
    github.com/raphaelvallat/pingouin/blob/master/pingouin/regression.py
    
    If covariance or whitening_matrix is not None, use GLS instead of OLS.
    
    Args:
        X: the feature matrix (observations are rows, features are columns);
           may be a polars DataFrame or 2D NumPy array; must be numeric
        Y: the outcome matrix (observations are rows, labels are columns,
           with each column treated as an independent regression problem);
           may be a polars DataFrame or 2D NumPy array; must be numeric
        add_intercept: if True, adds a column of ones to X before running the
                       regression; if False, assumes the user has added one
        report_intercept: if True, reports the betas, p-values etc. for the
                          intercept column added by add_intercept as the first
                          row of each returned matrix. Can only be True if
                          add_intercept=True.
        covariance: a precomputed covariance matrix; when not None, enables GLS
        whitening_matrix: a precomputed Cholesky whitening matrix; when not
                          None, enables GLS
        return_significance: if True, also returns the SEs, CIs and p-values,
                             not just the betas
        return_residuals: if True, also returns the residuals
        
    Returns:
         A namedtuple containing 5 polars DataFrames or NumPy arrays (whichever
         type X is), each with one row per feature (column of X) and one column
         per regression (column of Y):
        - 'beta': the regression coefficients
        - 'SE': the standard errors of the regression coefficients
        - 'lower_CI' and 'upper_CI': the 95% confidence intervals of the
          regression coefficients
        - 'p': the p-values of the regression coefficients.
        If return_significance=False, returns just the betas, as a single
        polars DataFrame or NumPy array.
        If return_residuals=True, also return the residuals, as a 6th polars
        DataFrame or NumPy array.
    """
    import numpy as np
    assert isinstance(X, np.ndarray), type(X)
    assert isinstance(Y, np.ndarray), type(Y)
    assert np.issubdtype(X.dtype, np.number), X.dtype
    assert np.issubdtype(Y.dtype, np.number), Y.dtype
    assert X.ndim == Y.ndim == 2, (X.ndim, Y.ndim)
    assert len(X) == len(Y), (X.shape, Y.shape)
    assert not np.isnan(X).any()
    assert not np.isnan(Y).any()
    assert np.linalg.matrix_rank(X) == min(X.shape), \
        f'np.linalg.matrix_rank(X) ({np.linalg.matrix_rank(X)}) < ' \
        f'min(X.shape) ({min(X.shape)})'
    assert covariance is None or whitening_matrix is None
    if report_intercept and not add_intercept:
        raise ValueError("add_intercept is False, so report_intercept must "
                         "also be False since there's no intercept to report!")
    if add_intercept:
        X = np.column_stack((np.ones(len(X), X.dtype), X))
    if covariance is not None:
        if isinstance(covariance, pl.DataFrame):
            covariance = covariance.to_numpy()
        assert isinstance(covariance, np.ndarray), type(covariance)
        assert np.issubdtype(covariance.dtype, np.number), covariance.dtype
        assert covariance.shape == (len(X), len(X)), covariance.shape
        # noinspection PyUnresolvedReferences
        assert (covariance == covariance.T).all(), \
            'Covariance is not symmetric'
        assert not np.isnan(covariance).any()
        whitening_matrix = covariance_to_whitening_matrix(covariance)
    elif whitening_matrix is not None:
        if isinstance(whitening_matrix, pl.DataFrame):
            whitening_matrix = whitening_matrix.to_numpy()
        assert isinstance(whitening_matrix, np.ndarray), type(whitening_matrix)
        assert np.issubdtype(whitening_matrix.dtype, np.number), \
            whitening_matrix.dtype
        assert whitening_matrix.shape == (len(X), len(X)), \
            whitening_matrix.shape
        assert not np.isnan(whitening_matrix).any()
    if whitening_matrix is not None:
        X = whitening_matrix @ X
        Y = whitening_matrix @ Y
    beta, residuals, rank, singular_values = lstsq(X, Y)
    if not return_significance:
        return beta[1:] if add_intercept and not report_intercept else beta
    else:
        from collections import namedtuple
        from scipy.special import ndtri, stdtr
        df = len(X) - rank
        if df == 0:
            assert rank == len(X) == X.shape[1]
            raise ValueError('X is a full-rank square matrix, so the line '
                             'perfectly fits the data and p-values are not '
                             'defined; re-run with return_significance=False '
                             'if you just want the betas')
        SE = np.sqrt(np.outer(np.linalg.pinv(X.T @ X).diagonal(),
                              residuals / df))
        CI_width = -ndtri(0.025) * SE  # -ndtri(x) == norm.isf(x)
        lower_CI = beta - CI_width
        upper_CI = beta + CI_width
        p = 2 * stdtr(df, -np.abs(beta / SE))  # stdtr(df, -x) == t.sf(x, df)
        if return_residuals:
            residuals = Y - X @ beta
        if add_intercept and not report_intercept:
            beta, SE, lower_CI, upper_CI, p = \
                beta[1:], SE[1:], lower_CI[1:], upper_CI[1:], p[1:]
        if return_residuals:
            LinearRegressionsResult = namedtuple('LinearRegressionsResult', (
                'beta', 'SE', 'lower_CI', 'upper_CI', 'p', 'residuals'))
            return LinearRegressionsResult(beta=beta, SE=SE, lower_CI=lower_CI,
                                       upper_CI=upper_CI, p=p,
                                           residuals=residuals)
        else:
            LinearRegressionsResult = namedtuple('LinearRegressionsResult', (
                'beta', 'SE', 'lower_CI', 'upper_CI', 'p'))
            return LinearRegressionsResult(beta=beta, SE=SE, lower_CI=lower_CI,
                                       upper_CI=upper_CI, p=p)


def lstsq(X, Y):
    """
    A wrapper for np.linalg.lstsq that:
     - raises one of three different errors when rank < X.shape[1], depending
       on the value of len(X) relative to rank and X.shape[1]
     - sets rcond=None to silence warnings without changing lstsq's behavior
    
    Args:
        X: the X input to np.linalg.lstsq
           (rows are observations, columns are features)
        Y: the Y input to np.linalg.lstsq
           (rows are observations, columns are labels)

    Returns:
        The same return values as np.linalg.lstsq:
        numpy.org/doc/stable/reference/generated/numpy.linalg.lstsq.html
    """
    import numpy as np
    beta, residuals, rank, singular_values = np.linalg.lstsq(X, Y, rcond=None)
    if len(residuals) == 0:
        assert rank == len(X) == X.shape[1] or rank < X.shape[1]
        if rank == len(X) == X.shape[1]:
            # manually calculate residuals
            residuals = ((Y - X @ beta) ** 2).sum(axis=0)
        else:
            assert rank < X.shape[1]
            message = (
                f'X has {X.shape[1]:,} features but only {rank:,} are '
                f'linearly independent; this means one or more of the columns '
                f'in X are a linear combination of one or more of the other '
                f'columns. ')
            if rank == len(X) < X.shape[1]:
                message += (
                    f'In your case, this is entirely because there are only '
                    f'{len(X):,} observations, so there are fewer '
                    f'observations than features! Did you mix up which '
                    f'dimension contains your features and which contains '
                    f'your observations? You may need to transpose X.')
            elif rank < len(X) < X.shape[1]:
                message += (
                    f'In your case, this is partly because there are only '
                    f'{len(X):,} observations, so there are fewer '
                    f'observations than features! Did you mix up which '
                    f'dimension contains your features and which contains '
                    f'your observations? You may need to transpose X. '
                    f'However, there are even fewer independent features '
                    f'({rank:,}) than observations ({len(X):,}). Did you '
                    f'include the same column twice? Two intercepts? Did you '
                    f'forget to drop one of the levels of a one-hot encoded '
                    f'categorical variable, say by calling '
                    f'pl.DataFrame.to_dummies() without drop_first=True?')
            else:
                assert len(X) < rank < X.shape[1] or \
                       rank < len(X) == X.shape[1]
                message += (
                    f'In your case, this is not explained by having fewer '
                    f'observations ({len(X):,}) than features '
                    f'({X.shape[1]:,}). Are your data points constant (and '
                    f'therefore collinear with the intercept)? Did you '
                    f'include the same column twice? Two intercepts? Did you '
                    f'forget to drop one of the levels of a one-hot encoded '
                    f'categorical variable, say by calling '
                    f'pl.DataFrame.to_dummies() without drop_first=True?')
            raise ValueError(message)
    return beta, residuals, rank, singular_values


@polars_numpy_autoconvert(use_columns_from='data')
def regress_out(data, covariates, *, covariance=None, whitening_matrix=None):
    """
    Regress out covariates from each column of data.
    
    According to the Frisch-Waugh-Lovell (FWL) Theorem (hbs.edu/research-
    computing-services/Shared%20Documents/Training/fwltheorem.pdf), regressing
    out covariates from both X and Y and then running Y ~ X is equivalent to
    running Y ~ X + covariates.
    
    Do NOT include an intercept: this should be a null-op if the covariates are
    uncorrelated with the data; it wouldn't be if you included an intercept.
    (Also, to use the FWL trick, including an intercept here means you'd have
    to run Y ~ 0 + X with no intercept after, which feels strange.)
    
    If covariance or whitening_matrix are specified, regresses out using GLS
    instead, then back-transforms the residuals to the original space.
    
    Args:
        data: a numeric polars DataFrame or NumPy matrix
        covariates: a numeric polars DataFrame or NumPy matrix with the same
                    number of rows as data
        covariance: a precomputed covariance matrix; when not None, enables GLS
        whitening_matrix: a precomputed Cholesky whitening matrix; when not
                          None, enables GLS

    Returns:
        A matrix of the same size and type as data, containing the residuals
        of data after regressing out covariates.
    """
    if len(data) != len(covariates):
        raise ValueError('Data and covariates have different lengths!')
    if covariance is not None and whitening_matrix is not None:
        raise ValueError('Cannot specify both covariance and '
                         'whitening_matrix!')
    if covariance is not None:
        whitening_matrix = covariance_to_whitening_matrix(covariance)
    if whitening_matrix is not None:
        # noinspection PyUnresolvedReferences
        covariates = whitening_matrix @ covariates
        # noinspection PyUnresolvedReferences
        data = whitening_matrix @ data
    beta = linear_regressions(covariates, data, add_intercept=False,
                              return_significance=False)
    residuals = data - covariates @ beta
    if whitening_matrix is not None:
        # Back-transform GLS betas to the original space, by multiplying by the
        # inverse of the whitening matrix (which is just the Cholesky
        # decomposition of the covariance matrix)
        import numpy as np
        from scipy.linalg.lapack import dtrtri
        inverse_whitening_matrix, info = dtrtri(whitening_matrix, lower=True)
        if info > 0:
            raise np.linalg.LinAlgError('Singular matrix')
        elif info < 0:
            raise ValueError(f'Invalid input to dtrtri (info = {info})')
        # noinspection PyUnresolvedReferences
        residuals = inverse_whitening_matrix @ residuals
    return residuals


@polars_numpy_autoconvert(use_columns_from=['Y', 'X'])
def correlate(X, Y=None, *, rank=False, GLS=False, covariates=None,
              include_PCs=False):
    """
    Correlates each pair of columns of X, or (if Y is not None), each column
    of X with each column of Y.
    
    To test:
    
    import numpy as np
    import statsmodels.api as sm
    from utils import correlate, estimate_covariance
    rng = np.random.default_rng(0)
    X = rng.standard_normal(size=(1000, 100))
    y = rng.standard_normal(size=(1000, 50))
    covariates = rng.standard_normal(size=(1000, 20))
    a = sm.GLS(X[:, 0], sm.add_constant(X[:, 1]),
               sigma=estimate_covariance(X)).fit()
    b = correlate(X, GLS=True)
    assert np.allclose(a.pvalues[1], b.p[0, 1])
    
    Args:
        X: a polars DataFrame or 2D NumPy array to correlate, either with Y
           (if not None) or with itself
        Y: a polars DataFrame or 2D NumPy array to correlate with X
        rank: if True, rank-transforms X and Y
        GLS: if True, performs GLS instead of OLS, estimating the covariance
             matrix using the OAS method (sklearn.covariance.oas)
        covariates: if not None, must be a matrix of covariates; will be
                    regressed out of both X and Y (with GLS, if GLS=True)
                    before correlating
        include_PCs: if True, include the top N principal components (PCs) of
                     X as covariates (in addition to any covariates specified
                     via the covariates argument), where N is the optimal
                     number of PCs according to PCAforQTL's runElbow(). The Y
                     argument cannot currently be used with include_PCs. See
                     biorxiv.org/content/10.1101/2022.03.09.483661v1.full and
                     github.com/heatherjzhou/PCAForQTL.
    Returns:
        A namedtuple containing two polars DataFrames or NumPy arrays
        (depending on the type of X and Y) of shape (X.shape[1] ×
        Y.shape[1]), or (X.shape[1] × X.shape[1]) if Y is None:
        - correlation: correlations (or for GLS, pseudocorrelations; see below
                       for definition)
        - p: p-values for each correlation; when Y is None, p-values along the
             diagonal will be 1
    """
    # Definition of pseudocorrelation:
    #
    # For GLS, gls(y ~ x) and gls(x ~ y) have the same p-value but different
    # effect sizes beta_1 and beta_2. Use sqrt(beta_1 * beta_2) * sign(beta_1)
    # as GLS's correlation coefficient ("pseudocorrelation"?). sign(beta_2)
    # could be used instead of sign(beta_1) in that formula: they're the same.
    #
    # Why sqrt(beta_1 * beta_2) * sign(beta_1)? It's symmetric in beta_1 and
    # beta_2, and it's equal to the Pearson correlation for OLS. Here's why:
    #
    #     If beta_1 is the beta from lm(y ~ x), then:
    # (1) pearson_1 = beta_1 * sd(x) / sd(y)
    #
    #     and if beta_2 is the beta from lm(x ~ y), then:
    # (2) pearson_2 = beta_2 * sd(y) / sd(x)
    #
    #     Multiply (1) and (2):
    # (3) pearson_1 * pearson_2 = beta_1 * beta_2
    #
    #     But pearson_1 == pearson_2 since Pearson is symmetric:
    # (4) pearson^2 = beta_1 * beta_2
    #
    #     Take the square root
    # (5) |pearson| = sqrt(beta_1 * beta_2)
    #
    #     But beta_1 and beta_2 and pearson all have the same sign, so:
    # (6) pearson = sqrt(beta_1 * beta_2) * sign(beta_1)  # or beta_2
    import numpy as np
    from scipy.special import stdtr
    # Check inputs
    assert isinstance(X, np.ndarray), type(X)
    assert np.issubdtype(X.dtype, np.number), X.dtype
    assert X.ndim == 2, X.ndim
    assert not np.isnan(X).any()
    assert np.linalg.matrix_rank(X) == min(X.shape), \
        f'np.linalg.matrix_rank(X) ({np.linalg.matrix_rank(X)}) < ' \
        f'min(X.shape) ({min(X.shape)})'
    if Y is not None:
        assert isinstance(Y, np.ndarray), type(Y)
        assert np.issubdtype(Y.dtype, np.number), Y.dtype
        assert Y.ndim == 2, Y.ndim
        assert len(Y) == len(X), (len(Y), len(X))
        assert not np.isnan(Y).any()
    if covariates is not None:
        assert isinstance(covariates, np.ndarray), type(covariates)
        assert np.issubdtype(covariates.dtype, np.number), covariates.dtype
        assert covariates.ndim == 2, covariates.ndim
        assert len(covariates) == len(X), (len(covariates), len(X))
        assert not np.isnan(covariates).any()
        assert np.linalg.matrix_rank(covariates) == min(covariates.shape), \
            f'np.linalg.matrix_rank(covariates) ' \
            f'({np.linalg.matrix_rank(covariates)}) < min(covariates.shape) ' \
            f'({min(covariates.shape)})'
    # Define namedtuple to store this function's return values
    from collections import namedtuple
    CorrelationResult = namedtuple('CorrelationResult', ('correlation', 'p'))
    # If rank=True, rank the data
    if rank:
        from scipy.stats import rankdata
        X = rankdata(X, axis=0)
        if Y is not None:
            Y = rankdata(Y, axis=0)
    # If include_PCs=True, include the top N PCs as covariates
    if include_PCs:
        assert Y is None  # not implemented
        from ryp import r, to_py, to_r
        # r('source("https://raw.githubusercontent.com/heatherjzhou/PCAForQTL/"
        #           "master/R/22.01.04_main1.1_runElbow.R")')
        r('runElbow<-function(X=NULL,prcompResult=NULL){if(is.null(prcompResul'
          't)){if(is.null(X)){stop("Please input X or prcompResult.")}else{cat'
          '("Running PCA...\\n");prcompResult<-prcomp(X,center=TRUE,scale.=TRU'
          'E)}}else{if(class(prcompResult)!="prcomp"){stop("prcompResult must '
          'be a prcomp object returned by the function prcomp().")}};importanc'
          'eTable<-summary(prcompResult)$importance;x<-1:ncol(importanceTable)'
          ';y<-importanceTable[2,];x1<-x[1];y1<-y[1];x2<-x[length(x)];y2<-y[le'
          'ngth(y)];x0<-x;y0<-y;distancesDenominator<-sqrt((x2-x1)^2+(y2-y1)^2'
          ');distancesNumerator<-abs((x2-x1)*(y1-y0)-(x1-x0)*(y2-y1));distance'
          's<-distancesNumerator/distancesDenominator;numOfPCsChosen<-which.ma'
          'x(distances);names(numOfPCsChosen)<-NULL;return(numOfPCsChosen)}')
        to_r(X, 'X')
        r('prcompResult = prcomp(X, center=F)')
        numOfPCsChosen = to_py('runElbow(X, prcompResult)')
        print(f'PCA: {numOfPCsChosen} PCs chosen')
        PCs = to_py('prcompResult$x')[:, :numOfPCsChosen]
        if covariates is not None:
            covariates = np.concatenate((covariates, PCs), axis=1)
        else:
            covariates = PCs
    # If GLS, compute the whitening matrix
    whitening_matrix = get_whitening_matrix(X) if GLS else None
    # If there are any covariates, regress them out of both X and Y. If
    # GLS=True, specify whitening_matrix to use GLS to regress them out.
    if covariates is not None:
        X = regress_out(X, covariates, whitening_matrix=whitening_matrix)
        if Y is not None:
            Y = regress_out(Y, covariates, whitening_matrix=whitening_matrix)
    # If not GLS, just take the Pearson/Spearman correlation. This works even
    # if there are covariates, since we regressed them out!
    if not GLS:
        correlation, p = cor(X, Y, return_p=True)
        if Y is None:
            np.fill_diagonal(p, 1)
        return CorrelationResult(correlation=correlation, p=p)
    # Whiten everything
    intercept = whitening_matrix.sum(axis=1)
    X = whitening_matrix @ X
    if Y is not None:
        Y = whitening_matrix @ Y
    # Run linear regression on the whitened data; get pseudocorrelation.

    def all_pairs_linear_regression(X, Y, intercept):
        # Make a design matrix with the intercept, plus an empty slot in the
        # first position (0-based) to copy each row of data to.
        df = len(X) - 2
        design_matrix = np.stack((intercept, np.empty(len(intercept))), axis=1)
        beta, SE = (np.empty((X.shape[1], Y.shape[1])) for _ in range(2))
        for col_index in range(X.shape[1]):
            design_matrix[:, 1] = X[:, col_index]
            coef, residuals = np.linalg.lstsq(design_matrix, Y, rcond=None)[:2]
            beta[col_index] = coef[1]
            SE[col_index] = np.sqrt(np.linalg.pinv(
                design_matrix.T @ design_matrix)[1, 1] * residuals / df)
        return beta, SE
    
    # Run linear regression between all pairs of variables
    if Y is None:
        beta, SE = all_pairs_linear_regression(X, X, intercept)
        correlation = np.sqrt(beta * beta.T) * np.sign(beta)
        # Use np.errstate(divide='ignore') to ignore division-by-zero errors
        # when calculating p-values for the correlation of each column of X
        # with itself.
        with np.errstate(divide='ignore'):
            # stdtr(df, -x) == t.sf(x, df)
            p = 2 * stdtr(len(X) - 2, -np.abs(beta / SE))
    else:
        beta_X, SE_X = all_pairs_linear_regression(X, Y, intercept)
        beta_Y, SE_Y = all_pairs_linear_regression(Y, X, intercept)
        correlation = np.sqrt(beta_X * beta_Y.T) * np.sign(beta_X)
        p = 2 * stdtr(len(X) - 2, -np.abs(beta_X / SE_X))
    # Fill the diagonals of the p-values with 1s (they'd be 0s otherwise)
    if Y is None:
        np.fill_diagonal(p, 1)
    return CorrelationResult(correlation=correlation, p=p)


###############################################################################
# [8] Clustering
###############################################################################


def get_n_clusters(affinity_matrix, *, max_n_clusters=10):
    """
    A heuristic for estimating the number of clusters for spectral clustering
    based on eigengaps: gaps between eigenvalues.
    
    See docs.scipy.org/doc/scipy/reference/tutorial/arpack.html and
    github.com/mingmingyang/auto_spectral_clustering/blob/master/autosp.py.

    Args:
        affinity_matrix: the affinity matrix for spectral clustering, e.g.
                         the k-nearest neighbors graph
        max_n_clusters: the maximum number of clusters to consider; setting to
                        higher numbers is more runtime-intensive

    Returns:
        The heuristically estimated best and second-best number of clusters.
    """
    import numpy as np
    from scipy.sparse.csgraph import laplacian
    from scipy.sparse.linalg import eigsh
    lap = laplacian(affinity_matrix, normed=True)
    eigenvalues, eigenvectors = eigsh(-lap, k=max_n_clusters + 1, sigma=1)
    eigenvalues = -eigenvalues[::-1]
    eigengap = np.diff(eigenvalues)
    n_clusters_ranking = 2 + np.argsort(-eigengap[1:])
    best, second_best = n_clusters_ranking[:2]
    return best, second_best


@polars_numpy_autoconvert()
def spectral_cluster(features, *, n_neighbors=20, max_n_clusters=10,
                     n_clusters=None):
    """
    Performs spectral clustering on a matrix of features
    
    Args:
        features: a polars DataFrame or 2D NumPy array of features
        n_neighbors: the number of nearest neighbors used to construct the
                     k-nearest neighbors graph
        max_n_clusters: the maximum number of clusters to consider in
                        get_n_clusters(); setting to higher numbers is more
                        runtime-intensive
        n_clusters: if not None, use this number of clusters rather than
                    using get_n_clusters() to heuristically estimate the number
                    of clusters

    Returns:
        A 1D NumPy array or polars Series of cluster labels, depending on the
        input type.
    """
    from sklearn.cluster import spectral_clustering
    from sklearn.neighbors import kneighbors_graph
    features = standardize(features)
    affinity_matrix = kneighbors_graph(features, n_neighbors=n_neighbors,
                                       include_self=True)
    affinity_matrix = (affinity_matrix + affinity_matrix.T) / 2
    if n_clusters is None:
        best, second_best = get_n_clusters(affinity_matrix,
                                           max_n_clusters=max_n_clusters)
        n_clusters = best
    cluster_labels = pl.Series(spectral_clustering(
        affinity_matrix, n_clusters=n_clusters, random_state=0))
    return cluster_labels


###############################################################################
# [9] General human genetics
###############################################################################


def standardize_chromosomes(chromosome, *, return_numeric=False,
                            omit_chr_prefix=False):
    """
    Standardizes the chromosome name(s) in chromosome, which can be an integer,
    string, polars Series or polars expression.
    
    Args:
        chromosome: the chromosome name(s) to standardize
        return_numeric: whether to return chromosomes as numbers
        omit_chr_prefix: whether to leave out the 'chr' prefix. Mutually
                         exclusive with return_numeric.

    Returns:
        The corresponding standardized chromosome name(s):
        - 'chr1' to 'chr22': return as-is
        - 1 to 22 or '1' to '22': convert to 'chr1' to 'chr22'
        - 23, '23', 'chr23', 'X': convert to 'chrX'
        - 24, '24', 'chr24', 'Y': convert to 'chrY'
        - 'M' or 'MT': convert to 'chrM'
        - 25, '25', 'chr25': disallow; could refer to either chrM or chrXY (the
                             pseudoautosomal region of the X and Y chromosomes)
        Or, if return_numeric=True, disallow chrM and its aliases and convert
        everything to a number between 1 and 24, inclusive.
        Or, if omit_chr_prefix=True, leave out the 'chr' prefix.
    """
    if return_numeric and omit_chr_prefix:
        raise ValueError(f'Only one of return_numeric and omit_chr_prefix can '
                         f'be True!')
    check_type(chromosome, 'chromosome', (int, str, pl.Series, pl.Expr),
               'an int, str, or polars Series or expression')
    if isinstance(chromosome, pl.Expr):
        return chromosome.map_batches(lambda col: standardize_chromosomes(
            col, return_numeric=return_numeric,
            omit_chr_prefix=omit_chr_prefix))
    elif isinstance(chromosome, str) or is_integer(chromosome):
        if isinstance(chromosome, str):
            original_chromosome = chromosome
            chromosome = chromosome.removeprefix('chr').replace('X', '23')\
                .replace('Y', '24').replace('MT', 'M')
            if not ((chromosome.isdigit() and 1 <= int(chromosome) <= 24) or
                    (chromosome == 'M' and not return_numeric)):
                raise ValueError(f'Invalid chromosome '
                                 f'{original_chromosome!r}!')
        else:
            if chromosome not in range(1, 25):
                raise ValueError(f'chromosome == "{chromosome}" but should be '
                                 f'between 1 and 24 (chrY) inclusive!')
        if return_numeric:
            return int(chromosome)
        else:
            chromosome = str(chromosome).replace('23', 'X').replace('24', 'Y')
            return chromosome if omit_chr_prefix else f'chr{chromosome}'
    else:
        check_dtype(chromosome, 'chromosome',
                    (pl.String, pl.Categorical, pl.Enum, 'integer'))
        if chromosome.dtype in pl.INTEGER_DTYPES:
            mask = chromosome.is_in(range(1, 25))
            if not mask.all():
                error_message = (
                    'chroms were specified as integers, but some were not '
                    'between 1 and 24 (chrY) inclusive!')
                raise ValueError(error_message)
            chromosome = chromosome.cast(pl.String)
        else:
            # Remove this if statement once polars allows Categorical
            # .str.len_bytes(): github.com/pola-rs/polars/issues/9773
            if chromosome.dtype == pl.Categorical or \
                    chromosome.dtype == pl.Enum:
                chromosome = chromosome.cast(pl.String)
            chromosome = chromosome.str.replace('^chr', '')\
                .str.replace('X', '23').str.replace('Y', '24')\
                .str.replace('MT', 'M')
            valid = chromosome.cast(pl.Int8, strict=False).is_in(range(1, 25))
            if not return_numeric:
                valid |= chromosome == 'M'
            if not valid.all():
                raise ValueError('Some chromosomes are invalid!')
        if return_numeric:
            return chromosome.cast(pl.Int8)
        else:
            chromosome = chromosome.replace({'23': 'X', '24': 'Y'})
            return chromosome if omit_chr_prefix else 'chr' + chromosome


def load_alias_to_gene_map(gene_annotation_dir=f'{get_base_data_directory()}/'
                                               f'gene-annotations'):
    """
    Load a map from old gene names ("aliases") to their current gene names.
    
    Remove "ambiguous" aliases that map to multiple current gene names (like
    ACSM2, which maps to both ACSM2A and ACSM2B). This is done with
    pl.col.alias.is_unique() below.
    
    Also remove aliases that are gene names themselves, e.g. TTLL5 is an alias
    of TTLL10, but TTLL5 is also a gene name. This can happen when a gene
    "splits" into two genes. This is done with the two ~pl.col.alias.is_in
    conditions below.
    
    Args:
        gene_annotation_dir: the directory where the alias-to-gene map will be
                             cached. Must be run on the login node to generate
                             this cache, if it doesn't exist. Because
                             generating the cache can take a long time, you
                             will probably want to leave this argument at its
                             default value.

    Returns:
        A two-column DataFrame with old gene names in the "alias" column and
        their current gene names in the "gene" column.
    """
    gene_alias_file = f'{gene_annotation_dir}/gene_aliases.tsv'
    if not os.path.exists(gene_alias_file):
        raise_error_if_on_compute_node()
        os.makedirs(gene_annotation_dir, exist_ok=True)
        run(f'curl -fsSL https://ftp.ebi.ac.uk/pub/databases/genenames/new/'
            f'tsv/locus_groups/protein-coding_gene.txt | cut -f2,9,11 > '
            f'{gene_alias_file}')
    Ensembl_genes = get_Ensembl_map(gene_annotation_dir=gene_annotation_dir,
                                    no_unaliasing=True)\
        .filter(pl.col.most_recent_Ensembl_version ==
                pl.col.most_recent_Ensembl_version.first())\
        ['gene_name']
    alias_to_gene = pl.scan_csv(gene_alias_file, separator='\t')\
        .select(alias=pl.col.alias_symbol.str.split('|')
                .list.concat(pl.col.prev_symbol), gene='symbol')\
        .drop_nulls()\
        .explode('alias')\
        .filter(pl.col.alias.is_unique(), ~pl.col.alias.is_in(pl.col.gene),
                ~pl.col.alias.is_in(Ensembl_genes))\
        .collect()
    return alias_to_gene


def unalias(df, gene_col, *,
            gene_annotation_dir=f'{get_base_data_directory()}/'
                                f'gene-annotations'):
    """
    "Unaliases" the genes in df's gene_col by mapping old gene names to their
    current gene names, according to the map returned by
    load_alias_to_gene_map().
    
    Before matching a gene list from a third-party dataset
    Before matching gene lists from two third-party datasets to each other, run
    them both through this function. You don't have to run this on the result
    of Ensembl_to_gene(), though.
    
    Args:
        df: a polars DataFrame
        gene_col: the name of the column with the gene names to be unaliased
        gene_annotation_dir: the directory where the alias-to-gene map returned
                             by load_alias_to_gene_map() will be cached. Must
                             be run on the login node to generate this cache,
                             if it doesn't exist. Because generating the cache
                             can take a long time, you will probably want to
                             leave this argument at its default value.

    Returns:
        df with each matching gene name in gene_col mapped to its alias. gene
        names not matching any of the old gene names in
        load_alias_to_gene_map() are left as-is.
    """
    return map_df(df, gene_col, load_alias_to_gene_map(
        gene_annotation_dir=gene_annotation_dir),
                  key_col='alias', value_col='gene', retain_missing=True)


def get_Ensembl_map(ENSP=False,
                    gene_annotation_dir=f'{get_base_data_directory()}/'
                                        f'gene-annotations',
                    no_unaliasing=False):
    """
    Gets a map from ENSGs (or ENSPs, if ENSP=True) to their most recent gene
    symbols in the Ensembl database, according to the map returned by
    get_Ensembl_map(). The returned gene names are unaliased, so you don't have
    to run unalias() on them.
    
    Ensembl IDs retired on or before release-42 (in 2006) are not included
    since these releases lack GTF files on the Ensembl website, and these IDs
    would be unlikely to be encountered in modern genomics data anyway.
    
    Args:
        ENSP: whether to convert ENSPs (protein IDs) instead of ENSGs (gene
              IDs)
        gene_annotation_dir: the directory where the Ensembl to gene symbol map
                             will be cached. Must be run on the login node to
                             generate this cache, if it doesn't exist. Because
                             generating the cache can take a long time, you
                             will probably want to leave this argument at its
                             default value.
        no_unaliasing: if True, gene names will not be unalised
    Returns:
        A two-column polars DataFrame with an "Ensembl_ID" column containing
        each of the Ensembl IDs in the Ensembl database and a "gene_name"
        column containing their gene names.
    """
    Ensembl_ID_type = 'ENSP' if ENSP else 'ENSG'
    mapping_file = os.path.join(gene_annotation_dir,
                                f'{Ensembl_ID_type}_to_gene_name.tsv')
    if not os.path.exists(mapping_file):
        print(f'Generating "{mapping_file}"...')
        raise_error_if_on_compute_node()
        os.makedirs(gene_annotation_dir, exist_ok=True)
        latest_Ensembl_version = max(map(int, run(
            'curl -fsSL https://ftp.ensembl.org/pub/current_gtf/'
            'homo_sapiens/ | grep -oP "GRCh38\\.\\K[0-9]{3}" | uniq',
            stdout=subprocess.PIPE).stdout.split()))
        get_build = lambda i: "GRCh38" if i in range(76, 82) else "GRCh37" \
            if i in range(55, 76) else "NCBI36"
        GTF_basenames = [
            basename
            for i in range(latest_Ensembl_version, 81, -1) for basename in (
                [f'release-{i}/gtf/homo_sapiens/'
                 f'Homo_sapiens.GRCh38.{i}.chr_patch_hapl_scaff',
                 f'grch37/release-{i}/gtf/homo_sapiens/'
                 f'Homo_sapiens.GRCh37.{i}.chr_patch_hapl_scaff']
                if i in (87, 85, 82) else
                [f'release-{i}/gtf/homo_sapiens/'
                 f'Homo_sapiens.GRCh38.{i}.chr_patch_hapl_scaff'])
        ] + [
            f'release-{i}/gtf/homo_sapiens/Homo_sapiens.{get_build(i)}.{i}'
            for i in range(81, 47, -1)
        ] + [
            'release-47/gtf/Homo_sapiens.NCBI36.47',
            'release-46/homo_sapiens_46_36h/data/gtf/Homo_sapiens.NCBI36.46',
            'release-44/homo_sapiens_44_36f/data/gtf/Homo_sapiens.NCBI36.44',
            'release-43/homo_sapiens_43_36e/data/gtf/Homo_sapiens.NCBI36.43']
        sed_command = \
            's/.*gene_name "(\\S*)".*protein_id "(\\S*)".*/\\2\\t\\1/p' \
            if ENSP else \
            's/.*gene_id "(\\S*)".*gene_name "(\\S*)".*/\\1\\t\\2/p'
        cache_dir = os.path.join(gene_annotation_dir,
                                 f'{Ensembl_ID_type}_by_version')
        os.makedirs(cache_dir, exist_ok=True)
        version_mapping_files = [os.path.join(
            cache_dir, f'{basename.split("/")[1].replace("-", "_")}_grch37.tsv'
                       if basename.startswith('grch37/') else
                       f'{basename.split("/")[0].replace("-", "_")}.tsv')
            for basename in GTF_basenames]
        for basename, version_mapping_file in zip(GTF_basenames,
                                                  version_mapping_files):
            if not os.path.exists(version_mapping_file):
                print(f'Generating "{version_mapping_file}"...')
                run(f'curl -fsSL https://ftp.ensembl.org/pub/'
                    f'{basename}.gtf.gz | zcat | sed -nr \'{sed_command}\' | '
                    f'awk \'!seen[$1]++\' > {version_mapping_file}')
        run(f'awk \'!seen[$1]++ {{print $1, $2, gensub(/.*\\/([^/]+)\\..*/, '
            f'"\\\\1", "", FILENAME)}}\' OFS="\t" '
            f'{" ".join(version_mapping_files)} > {mapping_file}')
    return pl.read_csv(mapping_file, separator='\t', has_header=False,
                       new_columns=['Ensembl_ID', 'gene_name',
                                    'most_recent_Ensembl_version'],
                       comment_prefix='#')\
        .pipe(lambda df: df if no_unaliasing else unalias(
            df, 'gene_name', gene_annotation_dir=gene_annotation_dir))


def Ensembl_to_gene(df, Ensembl_ID_col, *, ENSP=False,
                    gene_annotation_dir=f'{get_base_data_directory()}/'
                                        f'gene-annotations'):
    """
    Converts a polars DataFrame with a column Ensembl_ID_col of Ensembl IDs
    (ENSGs or ENSPs) to their most recent gene symbols in the Ensembl database.
    Ensembl IDs may contain dots followed by version numbers; everything after
    the first dot is stripped before mapping to gene symbols. The returned gene
    names are unaliased, so you don't have to run unalias() on them.
    
    Args:
        df: a polars DataFrame
        Ensembl_ID_col: the name of the column with the Ensembl IDs
        ENSP: whether Ensembl_ID_col contains ENSPs (protein IDs) instead of
              ENSGs (gene IDs)
        gene_annotation_dir: the directory where the Ensembl to gene symbol map
                             returned by get_Ensembl_map() will be cached. Must
                             be run on the login node to generate this cache,
                             if it doesn't exist. Because generating the cache
                             can take a long time, you will probably want to
                             leave this argument at its default value.

    Returns:
        df with the Ensembl IDs in Ensembl_ID_col replaced with gene symbols;
        Ensembl IDs without gene symbols are set to null.
    """
    Ensembl_map = get_Ensembl_map(ENSP=ENSP,
                                  gene_annotation_dir=gene_annotation_dir)
    return df\
        .lazy()\
        .with_columns(pl.col(Ensembl_ID_col).str.split_exact('.', 1)
                      .struct.field('field_0').alias(Ensembl_ID_col))\
        .pipe(map_df, Ensembl_ID_col, Ensembl_map.lazy(), key_col='Ensembl_ID',
              value_col='gene_name')\
        .collect()


def load_GO_terms(GO_term_dir=f'{get_base_data_directory()}/GO-terms',
                  coding_genes_genome_build=None):
    """
    Loads a DataFrame of GO terms and the genes in each term.
    
    Args:
        GO_term_dir: directory to cache intermediates and the final result
        coding_genes_genome_build: must be None, 'hg19' or 'hg38'. If not None,
                                   unalias and subset to genes in that genome
                                   build, i.e. the genes returned by
                                   get_coding_genes(coding_genes_genome_build).

    Returns:
        A polars DataFrame with one row per GO term, with three columns:
        - name: the GO term name, e.g. biological process:taxis; all names
                start with biological_process, molecular_function, or
                cellular_component
        - ID: the GO term ID, e.g. GO:0042330
        - genes: a sorted list of genes in the GO term
        DataFrame rows will be sorted alphabetically by the name column.
    """
    GO_term_cache = f'{GO_term_dir}/GO_terms.parquet'
    if os.path.exists(GO_term_cache):
        # Load GO terms from cache, if already cached
        GO_terms = pl.read_parquet(GO_term_cache)
    else:
        import networkx as nx
        import obonet
        import re
        from io import StringIO
        # Load mapping of GO term (e.g. GO:0000001) to name (e.g. mitochondrion
        # inheritance) and namespace
        os.makedirs(GO_term_dir, exist_ok=True)
        obo_file = f'{GO_term_dir}/go-basic.obo'
        if not os.path.exists(obo_file):
            raise_error_if_on_compute_node()
            run(f'wget geneontology.org/ontology/go-basic.obo '
                f'-P {GO_term_dir}')
        condensed_GO_to_name_file = f'{GO_term_dir}/go-condensed.txt'
        if not os.path.exists(condensed_GO_to_name_file):
            run(f"grep '\\[Term\\]' -A3 {obo_file} | grep -v '\\[Term\\]' | "
                f"grep -v '^--' > {condensed_GO_to_name_file}")
        # Load mapping of GO IDs (e.g. GO:0000001) to names (e.g.
        # biological process:mitochondrion inheritance)
        GO_to_name = pl.scan_csv(condensed_GO_to_name_file, separator='\t',
                                 has_header=False)\
            .with_columns(pl.col.column_1.str.splitn(': ', 2)
                          .struct.rename_fields(['columns', 'values']))\
            .unnest('column_1')\
            .with_row_index()\
            .with_columns(pl.col.index // 3)\
            .collect()\
            .pivot(index='index', columns='columns', values='values')\
            .lazy()\
            .with_columns(name=pl.col.namespace.str.replace(
                '_', ' ', literal=True) + ':' + pl.col.name)\
            .drop('index', 'namespace')\
            .rename({'id': 'ID'})\
            .collect()
        # Load mapping of genes to GO IDs
        # Drop a tiny number of entries where the gene is missing
        # Map terms listed under alternative IDs to their canonical ID
        GO_file = f'{GO_term_dir}/goa_human.gaf.gz'
        if not os.path.exists(GO_file):
            raise_error_if_on_compute_node()
            run(f'wget geneontology.org/gene-associations/goa_human.gaf.gz '
                f'-P {GO_term_dir}')
        GO_terms = pl.read_csv(GO_file, separator='\t', has_header=False,
                               columns=[2, 4], new_columns=['gene', 'ID'],
                               comment_prefix='!').drop_nulls()
        obo_file_contents = open(obo_file).read()
        alt_ID_to_ID_map = {}
        for line in obo_file_contents.splitlines():
            if line.startswith('id: '):
                current_ID = line.removeprefix('id: ')
            elif line.startswith('alt_id: '):
                # noinspection PyUnboundLocalVariable
                alt_ID_to_ID_map[line.removeprefix('alt_id: ')] = current_ID
        GO_terms = GO_terms.with_columns(pl.col.ID.replace(alt_ID_to_ID_map))
        # Remove obsolete GO terms, if any
        obsolete_terms = set()
        for line in obo_file_contents.splitlines():
            if line.startswith('id: '):
                current_ID = line.removeprefix('id: ')
            elif line == 'is_obsolete: true':
                obsolete_terms.add(current_ID)
        GO_terms = GO_terms.filter(~pl.col.ID.is_in(obsolete_terms))
        # If a gene has a GO term, also make it part of all parents of the GO
        # term. part_of, regulates, positively_regulates and
        # negatively_regulates are all indicators of parent-child relationships
        # (just like is_a), but are not handled by obonet.read_obo(), so we
        # need to handle them ourselves.
        obo_file_contents = re.sub(
            'relationship: (part_of|regulates|positively_regulates|'
            'negatively_regulates)', 'is_a:', obo_file_contents)
        obo_graph = obonet.read_obo(StringIO(obo_file_contents))
        get_parents = cache(lambda GO_term: [GO_term] + list(
            nx.descendants(obo_graph, GO_term)))
        GO_terms = GO_terms.with_columns(pl.col.ID.map_elements(
            get_parents, return_dtype=pl.List(pl.String))).explode('ID')
        # Join with names
        GO_terms = GO_terms.join(GO_to_name, on='ID')
        # Drop duplicates
        GO_terms = GO_terms.unique()
        # Sort by GO term, then gene
        GO_terms = GO_terms.sort('name',  'gene')
        # Aggregate the genes for each GO term into a list
        GO_terms = GO_terms.group_by('name', 'ID').agg('gene')
        # Save
        GO_terms.write_parquet(GO_term_cache)
    # If coding_genes_genome_build is not None, unalias and subset to genes in
    # that genome build.
    if coding_genes_genome_build is not None:
        check_valid_genome_build(coding_genes_genome_build)
        # noinspection PyTypeChecker
        GO_terms = GO_terms\
            .pipe(unalias, gene_col='gene')\
            .filter(pl.col.gene.is_in(
                get_coding_genes(coding_genes_genome_build)['gene']))
    return GO_terms


def check_valid_genome_build(genome_build):
    """
    Checks whether a genome build is valid (for the functions here)
    
    Args:
        genome_build: the genome build; must be hg19 or hg38
    """
    assert genome_build == 'hg19' or genome_build == 'hg38', genome_build


def liftover(infile, outfile, *, input_genome_build='hg38',
             output_genome_build='hg19', log_file='/dev/null',
             bash_input=False,
             gene_annotation_dir=f'{get_base_data_directory()}/'
                                 f'gene-annotations'):
    """
    Converts an input BED file (with optional extra columns) in one genome
    build to an output BED file (with the same columns) in a different build.
    
    Args:
        infile: the input BED file; may have optional extra columns. If
                bash_input=True, must be a bash command that outputs to stdout.
        outfile: the output BED file; will have the same extra columns
        input_genome_build: the input file's genome build
        output_genome_build: the output file's genome build
        log_file: an optional log file location
        bash_input: if True, infile must be a bash command that outputs to
                    stdout, instead of a file.
        gene_annotation_dir: the directory where the liftOver chain file will
                             be cached. Must be run on the login node to
                             generate this cache, if it doesn't exist. Because
                             generating the cache can take a long time, you
                             will probably want to leave this argument at its
                             default value.
    """
    check_valid_genome_build(input_genome_build)
    check_valid_genome_build(output_genome_build)
    chain_file = os.path.join(
        gene_annotation_dir,
        f'{input_genome_build}To{output_genome_build.capitalize()}.'
        f'over.chain.gz')
    if not os.path.exists(chain_file):
        raise_error_if_on_compute_node()
        run(f'rsync -za rsync://hgdownload.cse.ucsc.edu/goldenPath/'
            f'{input_genome_build}/liftOver/{os.path.basename(chain_file)} '
            f'-O {chain_file}')
    infile = f"<({infile})" if bash_input else f"\"{infile}\""
    outfile = f'"{outfile}"'
    run(f'liftOver -bedPlus=3 {infile} {chain_file} {outfile} {log_file}')


def get_gencode_version():
    """
    Get the most recent available Gencode version. Must be run on a login node.
    
    Returns: the most recent available Gencode version.
    """
    import inspect
    calling_function = inspect.currentframe().f_back.f_code.co_name
    raise_error_if_on_compute_node(
        f'get_gencode_version() needs internet access; run on the login node! '
        f'Then, when running on the compute node, set gencode_version in '
        f'{calling_function} to the gencode version number it prints out.')
    gencode_version = int(run(
        'curl -fsSL https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human'
        '/latest_release/ | grep -oP "gencode\\.v\\K[0-9]{2,}" | uniq',
        stdout=subprocess.PIPE).stdout)
    print(f'Most recent Gencode version: {gencode_version}')
    return gencode_version


def get_coding_genes(genome_build='hg38', gencode_version=46,
                     gene_annotation_dir=f'{get_base_data_directory()}/'
                                         f'gene-annotations',
                     return_file=False):
    """
    Get the coordinates, gene symbols and Ensembl IDs of coding genes on
    autosomes, chrX and chrY.
    
    Args:
        genome_build: the genome build (hg19 or hg38) to get coordinates for
        gencode_version: the Gencode version to take coding genes from
        gene_annotation_dir: the directory where a BED file of the coding genes
                             will be cached. Must be run on the login node to
                             generate this cache, if it doesn't exist. Because
                             generating the cache can take a long time, you
                             will probably want to leave this argument at its
                             default value.
        return_file: If True, return a BED file path instead of a DataFrame.
                     Note: the start coordinates in the BED file are one less
                     than in the returned DataFrame, because BED files are
                     0-based while the returned DataFrame (and most sumstats)
                     are 1-based.
    Returns:
        A DataFrame with chrom, start, end, strand, gene, and Ensembl_IDs as
        columns, or a BED file of the same if return_file=True. The gene names
        are unaliased with unalias(), so you don't have to run unalias() again.
        Genes in the pseudoautosomal regions have chrX as their chromosome.
        Ensembl_IDs is a list column since gene symbols very occasionally have
        multiple Ensembl IDs.
    """
    check_valid_genome_build(genome_build)
    coding_genes_file = f'{gene_annotation_dir}/coding_genes_{genome_build}_' \
                        f'gencode_v{gencode_version}.bed'
    if os.path.exists(coding_genes_file):
        return coding_genes_file if return_file else pl.read_csv(
            coding_genes_file, separator='\t', has_header=False,
            new_columns=['chrom', 'start', 'end', 'strand', 'gene',
                         'Ensembl_IDs'])\
            .with_columns(pl.col.start + 1,  # convert BED to one-based
                          pl.col.Ensembl_IDs.str.split(','))
    os.makedirs(gene_annotation_dir, exist_ok=True)
    coding_genes_intermediate_file = \
        coding_genes_file.removesuffix('.bed') + '.intermediate.bed'
    if not os.path.exists(coding_genes_intermediate_file):
        raise_error_if_on_compute_node()
        gencode_URL = (
            f'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/'
            f'release_{gencode_version}'
            f'{"/GRCh37_mapping" if genome_build == "hg19" else ""}/gencode.v'
            f'{gencode_version}{"lift37" if genome_build == "hg19" else ""}.'
            f'annotation.gtf.gz')
        # Subtract 1 from the start coordinate because BED is 0-based, whereas
        # GTF is 1-based
        run(f'curl -fsSL {gencode_URL} | zcat | tr -d ";\\"" | awk \'$3 == '
            f'"gene" && $0 ~ /protein_coding|IG_.*_gene|TR_.*_gene/ {{for '
            f'(x = 1; x <= NF; x++) {{ if ($x == "gene_name") gene_name = '
            f'$(x + 1); else if ($x == "gene_id") gene_id = $(x + 1); }} '
            f'print $1, $4 - 1, $5, $7, gene_name, gene_id}}\' OFS="\t" | '
            f'sort -k1,1V -k2,3n > {coding_genes_intermediate_file}')
    coding_genes = pl.scan_csv(
        coding_genes_intermediate_file, separator='\t', has_header=False,
        new_columns=['chrom', 'start', 'end', 'strand', 'gene',
                     'Ensembl_IDs'])\
        .filter(pl.col.chrom != 'chrM')\
        .with_columns(Ensembl_IDs=pl.col.Ensembl_IDs.str.split_exact('.', 1)
                      .struct.field('field_0'))\
        .collect()
    # For hg19, 7 coding genes were subsequently merged into another gene
    # (PRAMEF21 -> PRAMEF20, NBPF16 -> NBPF15, FOXD4L2 -> FOXD4L4,
    # MRC1L1 -> MRC1, ANXA8L2 -> ANXA8L1, ASAH2C -> ASAH2B, CT45A4 -> CT45A3);
    # to avoid the same gene being present in two different locations, do NOT
    # unalias. Also, for both hg19 and hg38, ZNF475 is an alias of ZFP1, but
    # also a gene in its own right; we can similarly ignore the alias.
    assert len(coding_genes
               .with_columns(unaliased_gene='gene')
               .pipe(unalias, 'unaliased_gene',
                     gene_annotation_dir=gene_annotation_dir)
               .filter(pl.col.unaliased_gene != pl.col.gene)) == \
           (1 if genome_build == 'hg38' else 8)
    # No genes should have more than one chromosome, except those in the
    # pseudoautosomal regions
    chrX_PAR1_end = 2_781_479 if genome_build == 'hg38' else 2_699_520
    chrX_PAR2_start = 155_701_383 if genome_build == 'hg38' else 154_931_044
    chrY_PAR1_end = 2_781_479 if genome_build == 'hg38' else 2_649_520
    chrY_PAR2_start = 56_887_903 if genome_build == 'hg38' else 59_034_050
    non_PAR_coding_genes = coding_genes.filter(
        ~((pl.col.chrom == 'chrX') & (pl.col.end < chrX_PAR1_end)),
        ~((pl.col.chrom == 'chrX') & (pl.col.start >= chrX_PAR2_start)),
        ~((pl.col.chrom == 'chrY') & (pl.col.end < chrY_PAR1_end)),
        ~((pl.col.chrom == 'chrY') & (pl.col.end >= chrY_PAR2_start)))
    assert len(non_PAR_coding_genes.filter(
        pl.col.chrom.n_unique().over('gene') != 1)) == 0
    # No genes should have more than one strand
    assert len(coding_genes.filter(
        pl.col.strand.n_unique().over('gene') != 1)) == 0
    # Occasionally, non-pseudoautosomal genes may appear more than once, always
    # under different Ensembl IDs. Confirm that these genes are always
    # overlapping, then merge them, and make a list of their Ensembl IDs.
    duplicated_non_PAR_genes = \
        non_PAR_coding_genes.filter(pl.col.gene.is_duplicated())
    assert len(duplicated_non_PAR_genes.filter(
        pl.col.start.max().over('gene') >= pl.col.end.min().over('gene'))) == 0
    assert not duplicated_non_PAR_genes['Ensembl_IDs'].is_duplicated().any()
    coding_genes = coding_genes\
        .lazy()\
        .group_by('gene', maintain_order=True)\
        .agg(pl.first('chrom'), pl.min('start'), pl.max('end'),
             pl.first('strand'), 'Ensembl_IDs')\
        .with_columns(pl.col.Ensembl_IDs.list.unique().list.sort())\
        .select('chrom', 'start', 'end', 'strand', 'gene', 'Ensembl_IDs')\
        .collect()
    assert not coding_genes['gene'].is_duplicated().any()
    # Write to a file, and return
    coding_genes\
        .with_columns(pl.col.Ensembl_IDs.list.join(','))\
        .write_csv(coding_genes_file, separator='\t', include_header=False)
    run(f"rm '{coding_genes_intermediate_file}'")
    return coding_genes_file if return_file else \
        coding_genes.with_columns(pl.col.start + 1)


def get_coding_TSSs(genome_build='hg38', gencode_version=46,
                    gene_annotation_dir=f'{get_base_data_directory()}/'
                                        f'gene-annotations',
                    return_file=False):
    """
    Get the coordinates, gene symbols and Ensembl IDs of transcription start
    sites (TSSs) for coding genes on autosomes, chrX and chrY. Genes may have
    multiple TSSs, one per transcript.
    
    Args:
        genome_build: the genome build (hg19 or hg38) to get coordinates for
        gencode_version: the Gencode version to take coding genes from
        gene_annotation_dir: the directory where a BED file of the coding TSSs
                             will be cached. Must be run on the login node to
                             generate this cache, if it doesn't exist. Because
                             generating the cache can take a long time, you
                             will probably want to leave this argument at its
                             default value.
        return_file: If True, return a BED file path instead of a DataFrame.
                     The BED file will have start = bp - 1 and end = bp, since
                     BED files are 0-based while the returned DataFrame (and
                     most sumstats) are 1-based.
    Returns:
        A DataFrame with chrom, bp, strand, gene, and Ensembl_ID as columns, or
        a BED file of the same if return_file=True. The gene names are
        unaliased with unalias(), so you don't have to run unalias() again.
        Genes in the pseudoautosomal regions have chrX as their chromosome.
    """
    check_valid_genome_build(genome_build)
    coding_TSSs_file = f'{gene_annotation_dir}/coding_TSSs_{genome_build}_' \
                       f'gencode_v{gencode_version}.bed'
    if os.path.exists(coding_TSSs_file):
        return coding_TSSs_file if return_file else pl.read_csv(
            coding_TSSs_file, separator='\t', has_header=False,
            new_columns=['chrom', 'start', 'end', 'strand', 'gene',
                         'Ensembl_ID'])\
            .drop('start')\
            .rename({'end': 'bp'})
    os.makedirs(gene_annotation_dir, exist_ok=True)
    coding_TSSs_intermediate_file = \
        coding_TSSs_file.removesuffix('.bed') + '.intermediate.bed'
    if not os.path.exists(coding_TSSs_intermediate_file):
        raise_error_if_on_compute_node()
        gencode_URL = (
            f'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/'
            f'release_{gencode_version}'
            f'{"/GRCh37_mapping" if genome_build == "hg19" else ""}/gencode.v'
            f'{gencode_version}{"lift37" if genome_build == "hg19" else ""}.'
            f'annotation.gtf.gz')
        # Subtract 1 from the start coordinate because BED is 0-based, whereas
        # GTF is 1-based; also uniquify (sort -u) since multiple transcripts of
        # the same gene may share a TSS
        run(f'curl -fsSL {gencode_URL} | zcat | tr -d ";\\"" | awk \'$3 == '
            f'"transcript" && $0 ~ /protein_coding|IG_.*_gene|TR_.*_gene/ '
            f'{{for (x = 1; x <= NF; x++) {{ if ($x == "gene_name") '
            f'gene_name = $(x + 1); else if ($x == "gene_id") gene_id = '
            f'$(x + 1); }} print $1, $7 == "+" ? $4 - 1 : $5 - 1, '
            f'$7 == "+" ? $4 : $5, $7, gene_name, gene_id}}\' OFS="\t" | '
            f'sort -k1,1V -k2,3n -u > {coding_TSSs_intermediate_file}')
    coding_TSSs = pl.scan_csv(
        coding_TSSs_intermediate_file, separator='\t', has_header=False,
        new_columns=['chrom', 'start', 'end', 'strand', 'gene', 'Ensembl_ID'])\
        .filter(pl.col.chrom != 'chrM')\
        .with_columns(Ensembl_IDs=pl.col.Ensembl_IDs.str.split_exact('.', 1)
                      .struct.field('field_0'))\
        .collect()
    # For hg19, 7 coding genes were subsequently merged into another gene
    # (PRAMEF21 -> PRAMEF20, NBPF16 -> NBPF15, FOXD4L2 -> FOXD4L4,
    # MRC1L1 -> MRC1, ANXA8L2 -> ANXA8L1, ASAH2C -> ASAH2B, CT45A4 -> CT45A3);
    # to avoid the same gene being present in two different locations, do NOT
    # unalias. Also, for both hg19 and hg38, ZNF475 is an alias of ZFP1, but
    # also a gene in its own right; we can similarly ignore the alias.
    assert len(coding_TSSs
               .with_columns(unaliased_gene='gene')
               .pipe(unalias, 'unaliased_gene',
                     gene_annotation_dir=gene_annotation_dir)
               .filter(pl.col.unaliased_gene != pl.col.gene)
               .unique('gene')) == (1 if genome_build == 'hg38' else 8)
    # No TSSs should have more than one chromosome, except those in the
    # pseudoautosomal regions
    chrX_PAR1_end = 2_781_479 if genome_build == 'hg38' else 2_699_520
    chrX_PAR2_start = 155_701_383 if genome_build == 'hg38' else 154_931_044
    chrY_PAR1_end = 2_781_479 if genome_build == 'hg38' else 2_649_520
    chrY_PAR2_start = 56_887_903 if genome_build == 'hg38' else 59_034_050
    non_PAR_coding_TSSs = coding_TSSs.filter(
        ~((pl.col.chrom == 'chrX') & (pl.col.end < chrX_PAR1_end)),
        ~((pl.col.chrom == 'chrX') & (pl.col.start >= chrX_PAR2_start)),
        ~((pl.col.chrom == 'chrY') & (pl.col.end < chrY_PAR1_end)),
        ~((pl.col.chrom == 'chrY') & (pl.col.end >= chrY_PAR2_start)))
    assert len(non_PAR_coding_TSSs.filter(
        pl.col.chrom.n_unique().over('gene') != 1)) == 0
    # No TSSs should have more than one strand
    assert len(coding_TSSs.filter(
        pl.col.strand.n_unique().over('gene') != 1)) == 0
    # Make sure no genes share a TSS
    assert len(coding_TSSs.filter(
        pl.struct('chrom', 'start', 'end').is_duplicated())) == 0
    # Write to a file, and return
    coding_TSSs\
        .write_csv(coding_TSSs_file, separator='\t', include_header=False)
    run(f"rm '{coding_TSSs_intermediate_file}'")
    return coding_TSSs_file if return_file else \
        coding_TSSs.drop('start').rename({'end': 'bp'})


###############################################################################
# [10] GWAS
###############################################################################


def get_nearest_gene(sumstats, genome_build, *, gencode_version=46,
                     SNP_col='SNP', chrom_col='CHROM', bp_col='BP',
                     gene_col='gene', gene_distance_col='gene_distance',
                     TSS_col='TSS', TSS_distance_col='TSS_distance', TSS=None,
                     gene_annotation_dir=f'{get_base_data_directory()}/'
                                         f'gene-annotations',
                     include_details=False):
    """
    Gets the nearest coding gene(s) and/or TSS(s) for each variant in sumstats.
    
    Given a DataFrame of sumstats with columns SNP_col, chrom_col and bp_col,
    returns a DataFrame with those three columns plus (by default) four others,
    described below: gene_col, gene_distance_col, TSS_col, and
    TSS_distance_col. If include_details=True, include additional columns with
    more details (see below).
    
    Args:
        sumstats: the summary statistics, as a polars DataFrame
        genome_build: the genome build of sumstats
        gencode_version: the Gencode version to take coding genes from
        SNP_col: the name of the variant ID column in sumstats
        chrom_col: the name of the chromosome column in sumstats
        bp_col: the name of the base-pair column in sumstats
        gene_col: the name of the column in the output DataFrame containing the
                  gene symbol(s) of the nearest gene(s) to each variant, as a
                  list
        gene_distance_col: the name of the column in the output DataFrame
                           containing the distance from each variant to (each
                           of) the nearest gene(s) in gene_col, as a list. The
                           distance is 0 if the variant is inside the gene
                           (i.e. between the TSS and TES), and the distance to
                           whichever of the TSS and TES is closer, otherwise.
        TSS_col: the name of the column in the output DataFrame containing the
                 gene symbol(s) of the gene(s) with the nearest TSS(s) to each
                 variant, as a list. For genes with multiple transcripts, only
                 the transcript(s) with the nearest TSS are counted.
        TSS_distance_col: the name of the column in the output DataFrame
                          containing the distance from each variant to the
                          TSS(s) of (each of) the gene(s) in TSS_col, as a list
        TSS: if None (the default), include both gene_col/gene_distance_col and
             TSS_col/TSS_distance_col in the output. If False, include only
             gene_col/gene_distance_col. If True, include only
             TSS_col/TSS_distance_col.
        gene_annotation_dir: the directory where the coding gene locations
                             returned by get_coding_genes() will be cached.
                             Must be run on the login node to generate this
                             cache, if it doesn't exist. Because generating the
                             cache can take a long time, you will probably want
                             to leave this argument at its default value.
        include_details: whether to include additional columns in the output:
            For nearest genes (TSS=False or TSS=None):
            - previous_gene_boundary: the base-pair position of the previous
                                      gene boundary
            - previous_genes: a list of the previous gene(s)
            - next_gene_boundary: the base-pair position of the next gene
                                  boundary
            - next_genes: a list of the next gene(s)
            - genes_inside: a list of gene(s) containing the variant
            For nearest TSSs (TSS=True or TSS=None):
            - previous_TSS: the base-pair position of the previous TSS
            - previous_TSS_genes: a list of the previous TSS's gene(s)
            - next_TSS: the base-pair position of the next TSS
            - next_TSS_genes: a list of the next TSS's gene(s)
    
    Returns:
        A DataFrame with chrom_col, bp_col, gene_col/distance_col and/or
        TSS_col/TSS_distance_col columns, and optionally the other columns
        above if include_details=True.
    """
    # Sanity-check inputs
    check_valid_genome_build(genome_build)
    if sumstats.is_empty():
        raise ValueError(f'sumstats is empty!')
    if SNP_col not in sumstats:
        raise ValueError(f'{SNP_col!r} is not in sumstats; specify SNP_col')
    if chrom_col not in sumstats:
        raise ValueError(f'{chrom_col!r} is not in sumstats; specify '
                         f'chrom_col')
    if bp_col not in sumstats:
        raise ValueError(f'{bp_col!r} is not in sumstats; specify bp_col')
    if TSS is not True and TSS is not False and TSS is not None:
        raise ValueError(f'TSS must be True, False, or None')
    include_nearest_gene = TSS is not True
    include_nearest_TSS = TSS is not False
    if include_nearest_gene and gene_col in sumstats:
        raise ValueError(f'{gene_col!r} is already in sumstats; specify a '
                         f'different gene_col')
    if include_nearest_gene and gene_distance_col in sumstats:
        raise ValueError(f'{gene_distance_col!r} is already in sumstats; '
                         f'specify a different gene_distance_col')
    if include_nearest_TSS and TSS_col in sumstats:
        raise ValueError(f'{TSS_col!r} is already in sumstats; specify a '
                         f'different TSS_col')
    if include_nearest_TSS and TSS_distance_col in sumstats:
        raise ValueError(f'{TSS_distance_col!r} is already in sumstats; '
                         f'specify a different TSS_distance_col')
    if include_details and TSS is True:
        raise ValueError(f'include_details must be False when TSS=True')
    # Standardize sumstats' chromosomes and sort
    standardized_sorted_sumstats = sumstats\
        .lazy()\
        .select(SNP_col,
                pl.col(chrom_col).pipe(standardize_chromosomes),
                pl.col(chrom_col).alias('original_chrom'),
                bp_col)\
        .sort(chrom_col, bp_col)
    if include_nearest_gene:
        # Get coding genes
        coding_genes = get_coding_genes(
            genome_build=genome_build, gencode_version=gencode_version,
            gene_annotation_dir=gene_annotation_dir)
        # Get coding gene boundaries (starts and ends). Occasionally genes
        # share the same boundary, so make a list of the genes that share each
        # boundary. Almost all of these lists have just one gene, though.
        gene_boundaries = pl.concat([coding_genes.select('chrom', key, 'gene')
                                    .rename({'chrom': chrom_col, key: bp_col})
                                     for key in ('start', 'end')])\
            .group_by(chrom_col, bp_col).agg('gene')\
            .sort(chrom_col, bp_col)
        # If genes are nested (e.g. GPR52 is entirely contained within
        # RABGAP1L) or partially overlapping such that a gene boundary is
        # contained within another gene, include the containing gene in the
        # list as well. This is so that we can correctly infer which genes
        # contain each variant below.
        gene_boundaries = gene_boundaries\
            .lazy()\
            .drop('gene')\
            .with_row_index()\
            .join(gene_boundaries
                  .lazy()
                  .with_row_index()
                  .explode('gene')
                  .group_by('gene')
                  .agg(pl.col.index.min().alias('min_index'),
                       pl.col.index.max().alias('max_index'))
                  .with_columns(pl.int_ranges(
                      'min_index', pl.col.max_index + 1,
                      dtype=pl.UInt32).alias('index'))
                  .drop('min_index', 'max_index')
                  .explode('index')
                  .group_by('index')
                  .agg('gene'),
                  on='index')\
            .sort('index')\
            .drop('index')\
            .collect()
        # join_asof both forward, to the start of the next gene boundary, and
        # backward, to the end of the previous gene boundary.
        #
        # If any of the next and previous gene(s) are the same, then those are
        # the nearest genes, and since the variant is inside them, the distance
        # is 0.
        #
        # Otherwise, get the distances to both the next and previous genes:
        # when equal (very rare), take all the next and all the previous genes
        # as nearest genes; otherwise, just take the ones on the closer side.
        #
        # Deduplicate and sort each variant's list of nearest genes at the end.
        nearest_genes = standardized_sorted_sumstats\
            .join_asof(gene_boundaries.lazy()
                       .rename({bp_col: 'previous_gene_boundary',
                                'gene': 'previous_genes'}),
                       left_on=bp_col, right_on='previous_gene_boundary',
                       by=chrom_col, strategy='backward')\
            .join_asof(gene_boundaries.lazy()
                       .rename({bp_col: 'next_gene_boundary',
                                'gene': 'next_genes'}),
                       left_on=bp_col, right_on='next_gene_boundary',
                       by=chrom_col, strategy='forward')\
            .with_columns(pl.col('previous_genes', 'next_genes')
                          .fill_null([]))\
            .with_columns(distance_to_previous=pl.col(bp_col) -
                                               pl.col.previous_gene_boundary,
                          distance_to_next=pl.col.next_gene_boundary -
                                           pl.col(bp_col),
                          genes_inside=pl.col.previous_genes
                              .list.set_intersection(pl.col.next_genes))\
            .with_columns(pl.when(pl.col.genes_inside.list.len() > 0)
                          .then(pl.col.genes_inside)
                          .when(pl.col.distance_to_previous ==
                                pl.col.distance_to_next)
                          .then(pl.concat_list('previous_genes', 'next_genes'))
                          .when((pl.col.distance_to_previous <
                                 pl.col.distance_to_next) |
                                pl.col.distance_to_next.is_null())
                          .then(pl.col.previous_genes)
                          .otherwise(pl.col.next_genes)
                          .alias(gene_col),
                          pl.when(pl.col.genes_inside.list.len() > 0)
                          .then(0)
                          .when((pl.col.distance_to_previous <
                                 pl.col.distance_to_next) |
                                pl.col.distance_to_next.is_null())
                          .then(pl.col.distance_to_previous)
                          .otherwise(pl.col.distance_to_next)
                          .alias(gene_distance_col))\
            .with_columns(pl.col(gene_col).list.unique().list.sort())\
            .drop(chrom_col)\
            .rename({'original_chrom': chrom_col})\
            .collect()
        if include_details:
            nearest_genes = nearest_genes\
                .select(SNP_col, chrom_col, bp_col, gene_col,
                        gene_distance_col, 'previous_gene_boundary',
                        'previous_genes', 'next_gene_boundary', 'next_genes',
                        'genes_inside')
        else:
            nearest_genes = nearest_genes\
                .select(SNP_col, chrom_col, bp_col, gene_col,
                        gene_distance_col)
    if include_nearest_TSS:
        # Get coding TSSs; note that since no genes share a TSS, we do not need
        # to group_by(chrom_col, bp_col).agg('gene') as we did for
        # gene_boundaries. Instead, just box each gene as a one-element list.
        coding_TSSs = get_coding_TSSs(
            genome_build=genome_build, gencode_version=gencode_version,
            gene_annotation_dir=gene_annotation_dir)\
            .with_columns(pl.col.gene.reshape((-1, 1))
                          .cast(pl.List(pl.String)))
        # The logic for finding the nearest TSS is much simpler: just do two
        # asof joins to find the nearest TSS in each direction, then take the
        # closer of the two (or both, if equidistant)
        nearest_TSSs = standardized_sorted_sumstats\
            .join_asof(coding_TSSs.lazy()
                       .rename({'chrom': chrom_col, 'bp': 'previous_TSS',
                                'gene': 'previous_TSS_genes'}),
                       left_on=bp_col, right_on='previous_TSS', by=chrom_col,
                       strategy='backward')\
            .join_asof(coding_TSSs.lazy()
                       .rename({'chrom': chrom_col, 'bp': 'next_TSS',
                                'gene': 'next_TSS_genes'}),
                       left_on=bp_col, right_on='next_TSS', by=chrom_col,
                       strategy='forward')\
            .with_columns(pl.col('previous_TSS_genes', 'next_TSS_genes')
                          .fill_null([]))\
            .with_columns(distance_to_previous=pl.col(bp_col) -
                                               pl.col.previous_TSS,
                          distance_to_next=pl.col.next_TSS - pl.col(bp_col))\
            .with_columns(pl.when(pl.col.distance_to_previous ==
                                  pl.col.distance_to_next)
                          .then(pl.concat_list('previous_TSS_genes',
                                               'next_TSS_genes'))
                          .when((pl.col.distance_to_previous <
                                 pl.col.distance_to_next) |
                                pl.col.distance_to_next.is_null())
                          .then(pl.col.previous_TSS_genes)
                          .otherwise(pl.col.next_TSS_genes)
                          .alias(TSS_col),
                          pl.when((pl.col.distance_to_previous <
                                   pl.col.distance_to_next) |
                                  pl.col.distance_to_next.is_null())
                          .then(pl.col.distance_to_previous)
                          .otherwise(pl.col.distance_to_next)
                          .alias(TSS_distance_col))\
            .with_columns(pl.col(TSS_col).list.unique().list.sort())\
            .drop(chrom_col)\
            .rename({'original_chrom': chrom_col})\
            .collect()
        if include_details:
            nearest_TSSs = nearest_TSSs\
                .select(SNP_col, chrom_col, bp_col, TSS_col, TSS_distance_col,
                        'previous_TSS', 'previous_TSS_genes', 'next_TSS',
                        'next_TSS_genes')
        else:
            nearest_TSSs = nearest_TSSs\
                .select(SNP_col, chrom_col, bp_col, TSS_col, TSS_distance_col)
    if TSS is None:
        # noinspection PyUnboundLocalVariable
        return pl.concat([nearest_genes,
                          nearest_TSSs.drop(SNP_col, chrom_col, bp_col)],
                         how='horizontal')
    elif TSS is False:
        # noinspection PyUnboundLocalVariable
        return nearest_genes
    else:
        # noinspection PyUnboundLocalVariable
        return nearest_TSSs


def get_nearby_genes(sumstats, genome_build, cache_prefix, *,
                     gencode_version=46, max_distance=500_000,
                     SNP_col='SNP', chrom_col='CHROM', bp_col='BP',
                     gene_col='gene', gene_distance_col='gene_distance',
                     gene_annotation_dir=f'{get_base_data_directory()}/'
                                         f'gene-annotations'):
    """
    Gets all coding genes within max_distance for each variant in sumstats, via
    bedtools window.
    
    Given a DataFrame of sumstats with columns SNP_col, chrom_col and bp_col,
    returns a DataFrame with these three columns plus two others, described
    below: gene_col and gene_distance_col.
        
    Args:
        sumstats: the summary statistics, as a polars DataFrame. Sumstats must
                   be sorted (this is not checked).
        genome_build: the genome build of sumstats
        cache_prefix: the prefix of the location for the results of
                      `bedtools window` (e.g. `'file_prefix'` or
                      `'directory/file_prefix'`)
        gencode_version: the Gencode version to take coding genes from
        max_distance: the maximum distance from a variant to get nearby genes
        SNP_col: the name of the variant ID column in sumstats
        chrom_col: the name of the chromosome column in sumstats
        bp_col: the name of the base-pair column in sumstats
        gene_col: the name of the column in the output DataFrame containing the
                  gene symbol of each nearby gene, as a list
        gene_distance_col: the name of the column in the output DataFrame
                           containing the distance from the variant to each of
                           the nearby genes in gene_col, as a list. The
                           distance is 0 if the variant is inside the gene
                           (i.e. between the TSS and TES), and the distance to
                           whichever of the TSS and TES is closer, otherwise.
        gene_annotation_dir: the directory where the coding gene locations
                             returned by get_coding_genes() will be cached.
                             Must be run on the login node to generate this
                             cache, if it doesn't exist. Because generating the
                             cache can take a long time, you will probably want
                             to leave this argument at its default value.
    
    Returns:
        A DataFrame with SNP_col, chrom_col, bp_col, gene_col,
        gene_distance_col, and TSS_distance_col columns.
    """
    # Sanity-check inputs
    check_valid_genome_build(genome_build)
    if sumstats.is_empty():
        raise ValueError(f'sumstats is empty!')
    if SNP_col not in sumstats:
        raise ValueError(f'{SNP_col!r} not in sumstats; specify SNP_col')
    if chrom_col not in sumstats:
        raise ValueError(f'{chrom_col!r} not in sumstats; specify chrom_col')
    if bp_col not in sumstats:
        raise ValueError(f'{bp_col!r} not in sumstats; specify bp_col')
    # Get coding genes BED file
    coding_genes_bed_file = get_coding_genes(
        genome_build=genome_build, gencode_version=gencode_version,
        gene_annotation_dir=gene_annotation_dir, return_file=True)
    # Make sumstats BED file
    sumstats_bed_file = f'{cache_prefix}.sumstats.bed'
    sumstats\
        .with_columns(pl.col(chrom_col).pipe(standardize_chromosomes),
                      original_chrom=pl.col(chrom_col))\
        .select(chrom_col, pl.col(bp_col).sub(1).alias('start'),
                pl.col(bp_col).alias('end'), pl.col(SNP_col),
                'original_chrom')\
        .write_csv(sumstats_bed_file, separator='\t', include_header=False)
    # Run bedtools window
    nearby_genes_bed_file = f'{cache_prefix}.nearby_genes.bed'
    run(f'bedtools window -a {sumstats_bed_file} -b {coding_genes_bed_file} '
        f'-w {max_distance} > {nearby_genes_bed_file}')
    # Load output; remember to add 1 to the gene start to convert to 1-based
    nearby_genes = pl.scan_csv(
        nearby_genes_bed_file, separator='\t', has_header=False,
        new_columns=['__get_nearby_genes_variant_chrom',
                     '__get_nearby_genes_variant_start', bp_col, SNP_col,
                     chrom_col, '__get_nearby_genes_gene_chrom',
                     '__get_nearby_genes_gene_start',
                     '__get_nearby_genes_gene_end',
                     '__get_nearby_genes_gene_strand', gene_col,
                     '__get_nearby_genes_gene_Ensembl_ID'],
        schema_overrides={chrom_col: sumstats[chrom_col].dtype})\
        .select(SNP_col, chrom_col, bp_col, gene_col,
                pl.col.__get_nearby_genes_gene_start.add(1),
                '__get_nearby_genes_gene_end')\
        .with_columns(pl.when(pl.col(bp_col)
                              .is_between('__get_nearby_genes_gene_start',
                                          '__get_nearby_genes_gene_end'))
                      .then(0)
                      .otherwise(pl.min_horizontal(
                          pl.col('__get_nearby_genes_gene_start',
                                 '__get_nearby_genes_gene_end')
                          .sub(pl.col(bp_col)).abs()))
                      .alias(gene_distance_col))\
        .with_columns(pl.col(gene_distance_col).sort()
                      .over(SNP_col, chrom_col, bp_col))\
        .group_by(SNP_col, chrom_col, bp_col, maintain_order=True)\
        .agg(pl.col(gene_col, gene_distance_col))\
        .collect()
    # Align to sumstats, filling nulls with [] for variants with no genes
    # within max_distance.
    nearby_genes = sumstats\
        .select(SNP_col, chrom_col, bp_col)\
        .join(nearby_genes, on=(SNP_col, chrom_col, bp_col), how='left')\
        .with_columns(pl.col(gene_col, gene_distance_col).fill_null([]))
    # Clean up - but only at the end, so users can inspect what went wrong in
    # case of an error
    run(f"rm '{sumstats_bed_file}' '{nearby_genes_bed_file}'")
    return nearby_genes


def get_nearest_variant(genes, sumstats, cache_prefix, *, gene_col='gene',
                        gene_chrom_col='chrom', gene_start_col='start',
                        gene_end_col='end', SNP_col='SNP',
                        sumstats_chrom_col='CHROM', sumstats_bp_col='BP',
                        distance_col='distance'):
    """
    The reverse of `get_nearest_gene()`: for each gene in `genes`, gets the
    distance to the nearest variant in `sumstats` via `bedtools closest -d`.
    Based on Caleb Ji's `get_gene_distances` function.
    
    Given a DataFrame of sumstats with columns `SNP_col`, `chrom_col` and
    `bp_col`, and a DataFrame of genes with columns `gene_col`,
    `gene_chrom_col`, `gene_start_col`, and `gene_end_col`, returns a DataFrame
    with the four columns from the genes DataFrame plus two others: `SNP_col`,
    the nearest variant, and `distance_col`, the distance from the gene to this
    nearest variant.
    
    Args:
        genes: the list of genes, as a polars DataFrame. Must be sorted (this
               is not checked). To use all protein-coding genes, set
               `genes=get_coding_genes(genome_build)`.
        sumstats: the summary statistics, as a polars DataFrame. Must be sorted
                  (this is not checked), and `genes` and `sumstats` must use
                  the same genome build (this is also not checked).
        cache_prefix: the prefix of the location for the results of
                      `bedtools closest -d` (e.g. `'file_prefix'` or
                      `'directory/file_prefix'`)
        gene_col: the name of the gene symbol column in `genes`
        gene_chrom_col: the name of the chromosome column in `genes`
        gene_start_col: the name of the start column in `genes`
        gene_end_col: the name of the end column in `genes`
        SNP_col: the name of the variant ID column in `sumstats`.
        sumstats_chrom_col: the name of the chromosome column in `sumstats`
        sumstats_bp_col: the name of the base-pair column in `sumstats`
        distance_col: the name of the column in the output DataFrame containing
                      the distance from each gene to its nearest variant. The
                      distance is 0 if the variant is inside the gene (i.e.
                      between the TSS and TES), and the distance to whichever
                      of the TSS and TES is closer, otherwise.
    
    Returns:
        A DataFrame with `gene_col`, `gene_chrom_col`, `gene_start_col`,
        `gene_end_col`, `SNP_col`, and `distance_col` columns. `SNP_col` and
        `distance_col` will be `null` if the gene has no variants in `sumstats`
        on the same chromosome.
    """
    # Sanity-check inputs
    if genes.is_empty():
        error_message = 'genes is empty!'
        raise ValueError(error_message)
    if gene_col not in genes:
        error_message = f'{gene_col!r} not in genes; specify gene_col'
        raise ValueError(error_message)
    if gene_chrom_col not in genes:
        error_message = \
            f'{gene_chrom_col!r} not in genes; specify gene_chrom_col'
        raise ValueError(error_message)
    if gene_start_col not in genes:
        error_message = \
            f'{gene_start_col!r} not in genes; specify gene_start_col'
        raise ValueError(error_message)
    if gene_end_col not in genes:
        error_message = f'{gene_end_col!r} not in genes; specify gene_end_col'
        raise ValueError(error_message)
    if sumstats.is_empty():
        error_message = 'sumstats is empty!'
        raise ValueError(error_message)
    if SNP_col not in sumstats:
        error_message = f'{SNP_col!r} not in sumstats; specify SNP_col'
        raise ValueError(error_message)
    if sumstats_chrom_col not in sumstats:
        error_message = (
            f'{sumstats_chrom_col!r} not in sumstats; specify '
            f'sumstats_chrom_col')
        raise ValueError(error_message)
    if sumstats_bp_col not in sumstats:
        error_message = \
            f'{sumstats_bp_col!r} not in sumstats; specify sumstats_bp_col'
        raise ValueError(error_message)
    # Make gene BED file
    gene_bed_file = f'{cache_prefix}.genes.bed'
    genes\
        .select(pl.col(gene_chrom_col).pipe(standardize_chromosomes),
                gene_start_col, gene_end_col, gene_col)\
        .write_csv(gene_bed_file, separator='\t', include_header=False)
    # Make sumstats BED file
    sumstats_bed_file = f'{cache_prefix}.sumstats.bed'
    sumstats\
        .with_columns(pl.col(sumstats_chrom_col).pipe(standardize_chromosomes),
                      original_chrom=pl.col(sumstats_chrom_col))\
        .select(sumstats_chrom_col,
                pl.col(sumstats_bp_col).sub(1).alias('start'),
                pl.col(sumstats_bp_col).alias('end'), pl.col(SNP_col),
                'original_chrom')\
        .write_csv(sumstats_bed_file, separator='\t', include_header=False)
    # Run `bedtools closest -d`
    distance_bed_file = f'{cache_prefix}.distance.bed'
    run(f'bedtools closest -d -a {gene_bed_file} -b {sumstats_bed_file} > '
        f'{distance_bed_file}')
    nearest_variants = pl.scan_csv(
        distance_bed_file, separator='\t', has_header=False,
        new_columns=['__get_nearest_variant_gene_chrom', gene_start_col,
                     gene_end_col, '__get_nearest_variant_gene_strand',
                     gene_col, '__get_nearest_variant_gene_Ensembl_ID',
                     '__get_nearest_variant_variant_chrom',
                     '__get_nearest_variant_variant_start',
                     '__get_nearest_variant_variant_end', SNP_col,
                     gene_chrom_col, distance_col],
        schema_overrides={gene_chrom_col: genes[gene_chrom_col].dtype})\
        .filter(~pl.col.distance.eq(-1))\
        .select([gene_col, gene_chrom_col, gene_start_col, gene_end_col,
                 SNP_col, distance_col])\
        .collect()
    # Align with `genes`
    nearest_variants = genes\
        .select(gene_col, gene_chrom_col, gene_start_col, gene_end_col)\
        .join(nearest_variants,
              on=(gene_col, gene_chrom_col, gene_start_col, gene_end_col),
              how='left')
    # Clean up - but only at the end, so users can inspect what went wrong in
    # case of an error
    run(f"rm '{gene_bed_file}' '{sumstats_bed_file}' '{distance_bed_file}''")
    return nearest_variants


def get_bim_or_pvar_file_type(bim_or_pvar_file):
    """
    Gets  whether a file is .bim or .pvar
    
    Args:
        bim_or_pvar_file: the bim or pvar file

    Returns:
        True if ends in .pvar, False if ends in .bim, raises an error otherwise
    """
    is_pvar = bim_or_pvar_file.endswith('.pvar')
    if not is_pvar and not bim_or_pvar_file.endswith('.bim'):
        raise ValueError(f'bim_or_pvar_file "{bim_or_pvar_file}" must end '
                         f'with .bim or .pvar!')
    return is_pvar


def read_bim_or_pvar(bim_or_pvar_file, categorical=False):
    """
    Reads a plink variant file (.bim or .pvar) as a polars DataFrame.
    
    Args:
        bim_or_pvar_file: the .bim or .pvar file to read from; must end with
                          .bim or .pvar, and the extension determines which
                          file type it's parsed as
        categorical: whether to read in the REF and ALT columns as Categorical
                     instead of String

    Returns:
        A polars DataFrame with the contents of the .bim or .pvar file,
        with columns CHROM, SNP, CM (if present), BP, REF, ALT
    """
    if not os.path.exists(bim_or_pvar_file):
        raise FileNotFoundError(f'No such file or directory: '
                                f'{bim_or_pvar_file}')
    is_pvar = get_bim_or_pvar_file_type(bim_or_pvar_file)
    # If pvar, get the number of lines at the beginning starting with ##, and
    # the number of lines after that starting with (should be 0 or 1)
    if is_pvar:
        num_double_comment_lines = int(run(
            f"sed -n '/^##/!q;p' {bim_or_pvar_file} | wc -l",
            stdout=subprocess.PIPE).stdout)
        num_single_comment_lines = int(run(
            f"sed -n '/^##/!{{/^#/p;q}}' {bim_or_pvar_file} | wc -l",
            stdout=subprocess.PIPE).stdout)
        if num_single_comment_lines > 1:
            raise ValueError(f'pvar file "{bim_or_pvar_file}" has multiple '
                             f'header lines!')
        skip_rows = num_double_comment_lines
        has_header = bool(num_single_comment_lines)
    else:
        skip_rows = 0
        has_header = False
    dtypes = {'#CHROM' if has_header else 'column_1': pl.String}
    if categorical:
        if has_header:
            ref_col = 'REF'
            alt_col = 'ALT'
        else:
            first_line = pl.read_csv(bim_or_pvar_file, separator='\t',
                                     has_header=has_header, n_rows=0)
            if first_line.width not in (5, 6):
                raise ValueError(f'bim_or_pvar_file "{bim_or_pvar_file}" has '
                                 f'{first_line.width} columns, but should '
                                 f'have 5 or 6!')
            if first_line.width == 6:
                ref_col = 'column_6'
                alt_col = 'column_5'
            else:
                ref_col = 'column_5'
                alt_col = 'column_4'
        dtypes |= {ref_col: pl.Categorical('lexical'),
                   alt_col: pl.Categorical('lexical')}
    variants = pl.read_csv(bim_or_pvar_file, separator='\t',
                           has_header=has_header, skip_rows=skip_rows,
                           schema_overrides=dtypes)
    if has_header:
        variants = variants.rename({'#CHROM': 'CHROM', 'POS': 'BP',
                                    'ID': 'SNP'})
    else:
        if variants.width not in (5, 6):
            raise ValueError(f'bim_or_pvar_file "{bim_or_pvar_file}" has '
                             f'{variants.width} columns, but should have 5 or '
                             f'6!')
        if variants.width == 6:
            variants = variants.rename({'column_1': 'CHROM', 'column_2': 'SNP',
                                        'column_3': 'CM', 'column_4': 'BP',
                                        'column_5': 'ALT', 'column_6': 'REF'})
        else:
            variants = variants.rename({'column_1': 'CHROM', 'column_2': 'SNP',
                                        'column_3': 'BP', 'column_4': 'ALT',
                                        'column_5': 'REF'})
    return variants


def write_bim_or_pvar(variants, bim_or_pvar_file):
    """
    Writes a polars DataFrame to bim or pvar format. If the filename ends in
    .pvar, a header is added. If the centimorgan (CM) column isn't specified,
    it's set to 0 (if bim) or omitted (if pvar).
    
    Args:
        variants: a polars DataFrame of variants with columns CHROM, SNP, CM,
                  BP, ALT, REF (CM is optional)
        bim_or_pvar_file: the bim or pvar file to write to; must end with .bim
                          or .pvar, and the extension determines which file
                          type it's written as
    """
    if variants.is_empty():
        raise ValueError(f'variants is empty!')
    is_pvar = get_bim_or_pvar_file_type(bim_or_pvar_file)
    if is_pvar:
        if 'CM' in variants:
            columns = 'CHROM', 'BP', 'SNP', 'REF', 'ALT', 'CM'
        else:
            columns = 'CHROM', 'BP', 'SNP', 'REF', 'ALT'
    else:
        columns = 'CHROM', 'SNP', 'CM', 'BP', 'ALT', 'REF'
        if 'CM' not in variants:
            variants = variants.with_columns(CM=0)
    variants = variants.select(columns)\
        .rename({'CHROM': '#CHROM', 'BP': 'POS', 'SNP': 'ID'})
    variants.write_csv(bim_or_pvar_file, separator='\t',
                       include_header=is_pvar)


def make_bim_or_pvar_IDs_unique(bim_or_pvar_file, new_bim_or_pvar_file, 
                                separator='~'):
    """
    Make the variant IDs of a bim or pvar file unique when there are
    multiallelic variants, by replacing each SNP ID with separator.join(
    (SNP, CHROM, REF, ALT)). Write the result to a new bim or pvar file.

    On the off-chance a variant ID already contains the default separator, '~',
    you will get an error and will need to specify a different string via the 
    separator argument.
    
    Why include CHROM and not just ID/REF/ALT? Because the same variant may map 
    to multiple genomic locations (ncbi.nlm.nih.gov/snp/docs/rs_multi_mapping). 
    
    Why not include BP? Because it introduces a dependence on genome build. In
    any case, it's vanishingly unlikely that a multi-mapping variant would have
    the same base-pair position on two different chromosomes. 

    Args:
        bim_or_pvar_file: the bim or pvar file to read from; must end with .bim
                          or .pvar, and the extension determines which file
                          type it's parsed and written as
        new_bim_or_pvar_file: the bim or pvar file to write to after
                              uniquifying the variant IDs
        separator: the separator to join the SNP, CHROM, REF, and ALT columns
                   with when uniquifying; no variant IDs may contain the
                   separator
    """
    bim_or_pvar = read_bim_or_pvar(bim_or_pvar_file)
    if bim_or_pvar['SNP'].str.contains(separator).any():
        suffix = bim_or_pvar_file.split('.')[-1]
        raise ValueError(f'Some variant IDs in the {suffix} file '
                         f'{bim_or_pvar_file} contain the separator '
                         f'{separator!r}; specify a different separator via '
                         f'the separator argument!')
    bim_or_pvar\
        .with_columns(pl.concat_str('SNP', 'CHROM', 'REF', 'ALT',
                                    separator=separator))\
        .pipe(write_bim_or_pvar, new_bim_or_pvar_file)


def reverse_make_bim_or_pvar_IDs_unique(bim_or_pvar_file, new_bim_or_pvar_file,
                                        separator='~'):
    """
    Reverse the transformation of make_bim_or_pvar_IDs_unique() by removing the
    CHROM, REF and ALT added to the SNP ID by that function. Write the result
    to a new bim or pvar file.
    Args:
        bim_or_pvar_file: the bim or pvar file to read from; must end with .bim
                          or .pvar, and the extension determines which file
                          type it's parsed and written as
        new_bim_or_pvar_file: the bim or pvar file to write to after removing
                              the CHROM, REF and ALT from the variant IDs
        separator: the separator that the SNP, CHROM, REF, and ALT columns
                   were joined with when uniquifying
    """
    bim_or_pvar = read_bim_or_pvar(bim_or_pvar_file)
    if not bim_or_pvar['SNP'].str.contains(separator).all():
        suffix = bim_or_pvar_file.split('.')[-1]
        raise ValueError(f'not all variant IDs in the {suffix} file '
                         f'{bim_or_pvar_file} contain the separator '
                         f'{separator!r}; specify a different separator via '
                         f'the separator argument!')
    bim_or_pvar\
        .with_columns(SNP=pl.col.SNP.str.split_exact(separator, 1)
                      .struct.field('field_0'))\
        .pipe(write_bim_or_pvar, new_bim_or_pvar_file)


def merge_bfiles(bfiles, merged_bfile):
    """
    Merge multiple bed/bim/fam filesets (bfiles) into one big fileset
    (merged_bfile).
    
    As an optimization, just concatenates the raw bytes of the bed files rather
    than calling plink. As a result, all filesets in bfiles must have the same
    sample info (.fam file) and non-overlapping variants, though this is not
    checked for speed. This is almost always true when merging, e.g. in the
    very common case where you want to merge genetic data with one fileset per
    chromosome into a single fileset.

    Supports multiallelic variants, unlike plink 1.9's --merge-list and plink
    2.0's --pmerge-list. (More precisely, plink 2.0 doesn't support
    --pmerge-list on files with "split" multiallelic variants, but also doesn't
    yet, as of September 2023, support merging multiallelic variants with
    --make-pgen multiallelics=+. Even once this is supported, merge_bfiles()
    should still be much faster than --pmerge-list.)

    Args:
        bfiles: a list/tuple/etc. of prefixes of bed/bim/fam filesets to merge
        merged_bfile: the prefix of the merged bed/bim/fam fileset to be
                      created; {merged_bfile}.{bed,bim,fam} must not exist yet
    """
    if isinstance(bfiles, str):
        raise ValueError(f'bfiles must be a list/tuple/etc. of filesets, not '
                         f'a single string!')
    if len(bfiles) < 2:
        raise ValueError(f'bfiles must contain two or more filesets, not '
                         f'{len(bfiles):,}!')
    for suffix in 'bed', 'bim', 'fam':
        output_plink_file = f'{merged_bfile}.{suffix}'
        if os.path.exists(output_plink_file):
            raise FileExistsError(f'{output_plink_file} already exists!')
        for bfile in bfiles:
            input_plink_file = f'{bfile}.{suffix}'
            if not os.path.exists(input_plink_file):
                raise FileNotFoundError(f'No such file or directory: '
                                        f'{input_plink_file}')
    # Concatenate bed files: the first three bytes are the header, which is
    # always 01101100 00011011 00000001, and the rest is the data. Take the
    # header from the first file, and just the part after the header (tail
    # -qc+4) for the remaining files.
    run(f"(cat '{bfiles[0]}.bed'; tail -qc+4 " + ' '.join(
        f"'{bfile}.bed'" for bfile in bfiles[1:]) +
        f") > '{merged_bfile}.bed'")
    # Concatenate bim files
    run(f"cat " + ' '.join(f"'{bfile}.bim'" for bfile in bfiles) +
        f" > '{merged_bfile}.bim'")
    # Take the first fam file (since they're assumed to be all the same)
    run(f"cp '{bfiles[0]}.fam' '{merged_bfile}.fam'")


def merge_pfiles(pfiles, merged_pfile, temp_dir=os.environ.get('SCRATCH', '.'),
                 num_threads=1, memory=2000):
    """
    Merge multiple pgen/pvar/psam filesets (pfiles) into one big fileset
    (merged_pfile) using plink 2.0's --pmerge-list.
    
    --pmerge-list doesn't support multiallelic variants. (More precisely,
    plink2 doesn't support --pmerge-list on files with "split" multiallelic
    variants, but also doesn't yet, as of September 2023, support merging
    multiallelic variants with --make-pgen multiallelics=+.) We circumvent this
    by making the pvar IDs unique with make_bim_or_pvar_IDs_unique() before
    merging, and reset them after with reverse_make_bim_or_pvar_IDs_unique().

    Args:
        pfiles: a list/tuple/etc. of prefixes of pgen/pvar/psam filesets to
                merge
        merged_pfile: the prefix of the merged pgen/pvar/psam fileset to be
                      created; {merged_bfile}.{pgen,pvar,psam} must not exist
                      yet
        temp_dir: a directory to store the temporary uniquified pvars created
                  by make_bim_or_pvar_IDs_unique(), which will be deleted at
                  the end of the run if successful
        num_threads: the number of threads to use for merging; if None, use all
                     available cores
        memory: the number of megabytes of memory for plink to reserve during
                merging via the --memory flag
    """
    temp_pfiles = [os.path.join(temp_dir, f'{os.path.basename(pfile)}.tmp')
                   for pfile in pfiles]
    for pfile, temp_pfile in zip(pfiles, temp_pfiles):
        make_bim_or_pvar_IDs_unique(f'{pfile}.pvar', f'{temp_pfile}.pvar')
    pmerge_list = ';'.join(f'echo {pfile}.pgen {temp_pfile}.pvar {pfile}.psam'
                           for pfile, temp_pfile in zip(pfiles, temp_pfiles))
    run(f'plink2 '
        f'{f"--threads {num_threads} " if num_threads is not None else ""}'
        f'--memory {memory} '
        f'--pmerge-list <({pmerge_list}) '
        f'--make-pgen '
        f'--out {merged_pfile}')
    reverse_make_bim_or_pvar_IDs_unique(f'{merged_pfile}.pvar',
                                        f'{merged_pfile}.pvar')
    for temp_pfile in temp_pfiles:
        run(f'rm "{temp_pfile}.pvar"')
    run(f'rm "{merged_pfile}-merge.pgen" "{merged_pfile}-merge.pvar" '
        f'"{merged_pfile}-merge.psam"')


def flip_alleles(sumstats, *, flip_col='FLIP', ref_col='REF', alt_col='ALT',
                 beta_col='BETA', OR_col='OR', Z_col='Z', AAF_col='AAF',
                 AAF_cases_col='AAF_CASES', AAF_controls_col='AAF_controls'):
    """
    Flips alleles in sumstats where flip_col is True, and also flips effect
    sizes, odds ratios, Z-score, and allele frequencies for these variants.
    
    flip_col, ref_col, alt_col, and at least one of beta_col, OR_col and Z_col
    are required. AAF_col, AAF_cases_col, AAF_controls_col, and the remainder
    of beta_col, OR_col and Z_col are optional: if a column name doesn't appear
    in sumstats, it won't be flipped, but there won't be an error.
    
    Args:
        sumstats: a polars DataFrame of summary statistics
        flip_col: the name of a boolean column in sumstats saying which alleles
                  to flip
        ref_col: the name of the reference allele column in sumstats
        alt_col: the name of the alternate allele column in sumstats
        beta_col: the name of the effect size column in sumstats
        OR_col: the name of the odds ratio column in sumstats
        Z_col: the name of the Z-score column in sumstats
        AAF_col: the name of the alternate allele frequency column in sumstats
        AAF_cases_col: the name of the case alternate allele frequency column
                       in sumstats
        AAF_controls_col: the name of the control alternate allele frequency
                          column in sumstats
    Returns:
        sumstats with alleles flipped where flip_col is True.
    """
    if ref_col not in sumstats:
        raise ValueError(f'ref_col "{ref_col}" not in sumstats! Did you '
                         f'forget to specify a custom ref_col?')
    if alt_col not in sumstats:
        raise ValueError(f'alt_col "{alt_col}" not in sumstats! Did you '
                         f'forget to specify a custom alt_col?')
    if beta_col not in sumstats and OR_col not in sumstats and \
            Z_col not in sumstats:
        raise ValueError(f'none of beta_col "{beta_col}", OR_col "{OR_col}", '
                         f'and Z_col "{Z_col}" are in sumstats! Did you '
                         f'forget to specify a custom beta_col, OR_col or '
                         f'Z_col?')
    flip_transformations = {
        ref_col: pl.col(alt_col), alt_col: pl.col(ref_col),
        beta_col: -pl.col(beta_col), OR_col: 1 / pl.col(OR_col),
        Z_col: -pl.col(Z_col), AAF_col: 1 - pl.col(AAF_col),
        AAF_cases_col: 1 - pl.col(AAF_cases_col),
        AAF_controls_col: 1 - pl.col(AAF_controls_col)}
    return sumstats.with_columns(**{
        column_name: pl.when(pl.col(flip_col)).then(transformation)
                     .otherwise(pl.col(column_name))
        for column_name, transformation in flip_transformations.items()
        if column_name in sumstats})


def harmonize_sumstats_to_bim_or_pvar(
        sumstats, bim_or_pvar, *, SNP_col='SNP', ref_col='REF', alt_col='ALT',
        beta_col='BETA', OR_col='OR', Z_col='Z', AAF_col='AAF',
        AAF_cases_col='AAF_CASES', AAF_controls_col='AAF_controls',
        chrom_col='CHROM', convert_chromosomes=False):
    """
    Harmonizes sumstats to a bim or pvar file by flipping single nucleotide
    variants' alleles with flip_alleles() as necessary to match the bim/pvar.
    Indels that don't match, or single-nucleotide variants that don't match
    even with flipping, are removed. If convert_chromosomes=True, converts
    chromosomes to plink format ('1', '2', ..., '22', 'X', 'Y').
    
    flip_col, ref_col, alt_col, and at least one of beta_col, OR_col and Z_col
    are required. AAF_col, AAF_cases_col, AAF_controls_col, and the remainder
    of beta_col, OR_col and Z_col are optional: if a column name doesn't appear
    in sumstats, it won't be flipped, but there won't be an error.
    
    Args:
        sumstats: a polars DataFrame of sumstats
        bim_or_pvar: a polars DataFrame of the contents of a bim or pvar file
                     created with get_rs_numbers_bim_or_pvar() and loaded into
                     memory with read_bim_or_pvar()
        SNP_col: the name of the variant ID column in sumstats
        ref_col: the name of the reference allele column in sumstats
        alt_col: the name of the alternate allele column in sumstats
        beta_col: the name of the effect size column in sumstats
        OR_col: the name of the odds ratio column in sumstats
        Z_col: the name of the Z-score column in sumstats
        AAF_col: the name of the alternate allele frequency column in sumstats
        AAF_cases_col: the name of the case alternate allele frequency column
                       in sumstats
        AAF_controls_col: the name of the control alternate allele frequency
                          column in sumstats
        chrom_col: the name of the chromosome column in sumstats; only used if
                   convert_chromosomes=True
        convert_chromosomes: whether to convert the chromosomes in chrom_col to
                             plink format

    Returns:
        Sumstats with alleles flipped and chromosomes optionally converted to
        plink format.
    """
    # Use 0 as a dummy value just to check whether it's not null after joining.
    return sumstats\
        .lazy()\
        .join(bim_or_pvar.lazy().select('SNP', 'REF', 'ALT',
                                        matches_without_flips=0),
              left_on=(SNP_col, ref_col, alt_col),
              right_on=('SNP', 'REF', 'ALT'), how='left')\
        .with_columns(pl.col.matches_without_flips.is_not_null())\
        .join(bim_or_pvar.lazy().select('SNP', 'REF', 'ALT',
                                        matches_with_flips=0),
              left_on=(SNP_col, ref_col, alt_col),
              right_on=('SNP', 'ALT', 'REF'), how='left')\
        .with_columns(pl.col.matches_with_flips.is_not_null() &
                      pl.col(ref_col).str.len_bytes().eq(1) &
                      pl.col(alt_col).str.len_bytes().eq(1))\
        .pipe(flip_alleles, flip_col='matches_with_flips', ref_col=ref_col,
              alt_col=alt_col, beta_col=beta_col, OR_col=OR_col, Z_col=Z_col,
              AAF_col=AAF_col, AAF_cases_col=AAF_cases_col,
              AAF_controls_col=AAF_controls_col)\
        .filter(pl.col.matches_without_flips | pl.col.matches_with_flips)\
        .drop('matches_without_flips', 'matches_with_flips')\
        .pipe(lambda df: df.with_columns(pl.col(chrom_col).pipe(
            standardize_chromosomes, omit_chr_prefix=True))
            if convert_chromosomes else df)\
        .collect()


def ld_clump(sumstats, cache_prefix, *, pfile=None, bfile=None,
             clump_p1=5e-8, clump_p2=0.001, clump_r2=0.001, clump_kb=5000,
             SNP_col='SNP', chrom_col='CHROM', bp_col='BP', ref_col='REF',
             alt_col='ALT', p_col='P', separator='~', num_threads=1,
             memory=2000, verbose=True):
    """
    Performs linkage disequilibrium (LD) clumping on summary statistics using
    plink's --clump (cog-genomics.org/plink/2.0/postproc#clump) and the LD info
    from a plink pgen/pvar/psam ("pfile") or bed/bim/bam ("bfile") fileset.
    
    Specify exactly one of pfile or bfile. pfile/bfile can be a list/tuple,
    e.g. if there's 1 fileset per chromosome. 
    
    For example, try (on Niagara):
    import os
    import polars as pl
    from utils import ld_clump
    sumstats_file = '/scratch/w/wainberg/wainberg/sumstats/daner_MDDwoBP_' \
                    '20201001_2015iR15iex_HRC_MDDwoBP_iPSYCH2015i_' \
                    'UKBtransformed_Wray_FinnGen_MVPaf_2_HRC_MAF01.gz'
    sumstats = pl.read_csv(sumstats_file, separator=' ', null_values='-')
    clumped_variants = ld_clump(
        sumstats, cache_prefix=f'{os.environ["SCRATCH"]}/LD_pruning_test',
        pfile='scratch/wainberg/1000G/European_autosomal_chrX_hg38',
        chrom_col='CHR', ref_col='A2', alt_col='A1')
    
    For information on LD clumping, see:
    - cog-genomics.org/plink/1.9/postproc#clump
    - zzz.bwh.harvard.edu/plink/clump.shtml
    - explodecomputer.github.io/EEPE_2016/worksheets_win/bioinformatics.html
    
    sumstats will be harmonized to pfile/bfile via 
    harmonize_sumstats_to_bim_or_pvar() and then written to the temp file 
    f'{cache_prefix}.sumstats.tmp' prior to clumping. This means that it's okay
    if some ref and alt alleles are flipped between the sumstats file and the
    pfile/bfile, and it does not matter which column is ref_col and which is
    alt_col (though you may as well specify ref_col='REF', alt_col='ALT' or
    ref_col='A2', alt_col='A1' for clarity).
    
    ld_clump() ensures all variant IDs are unique by setting IDs for both 
    sumstats and pfile/bfile to ID + '~' + CHROM + '~' + REF + '~' + ALT. This
    requires making a temporary pvar file (or files, if pfile/bfile is a list).
    (On the off-chance a variant ID already contains '~', you will get an error
    and will need to specify a different string via the separator argument.)
    
    Why include CHROM and not just ID/REF/ALT? Because the same variant may map 
    to multiple genomic locations (ncbi.nlm.nih.gov/snp/docs/rs_multi_mapping). 
    
    Why not include BP? Because it would require the genome builds to match, 
    and it's vanishingly unlikely that a multi-mapping variant would have the 
    same base-pair position on two different chromosomes. 
    
    Args:
        sumstats: the summary statistics, as a polars DataFrame
        cache_prefix: the prefix of the location for the temporary sumstats
                      and pvar files and --clump's results
        pfile: the plink 2.x pgen/pvar/psam fileset LD info is taken from; can
               be a string (assumed to be a single bfile for all chromosomes)
               or a list/tuple of bfiles (assumed to be one per chromosome).
               Exactly one of pfile and bfile must be specified.
        bfile: the plink 1.x bed/bim/bam fileset LD info is taken from; can be
               a string (assumed to be a single bfile for all chromosomes) or a
               list/tuple of bfiles (assumed to be one per chromosome).
        clump_p1: the value of --clump-p1 passed to --clump; for definitions of
                  --clump args, see cog-genomics.org/plink/2.0/postproc#clump
        clump_p2: the value of --clump-p2 passed to --clump
        clump_r2: the value of --clump-r2 passed to --clump
        clump_kb: the value of --clump-kb passed to --clump
        SNP_col: the name of the variant ID column in sumstats
        chrom_col: the name of the chromosome column in sumstats
        bp_col: the name of the base-pair column in sumstats; does NOT have to
                be in the same genome build as pfile/bfile
        ref_col: the name of the reference allele column in sumstats
        alt_col: the name of the alternate column in sumstats
        p_col: the name of the p-value column in sumstats
        separator: the string to place between the variant ID, chromosome,
                   reference allele, and alternate allele when creating the
                   temporary sumstats and pvar files; no variant IDs may
                   contain the separator
        num_threads: the number of threads to use for LD clumping; if None, use
                     all available cores
        memory: the number of megabytes of memory for plink to reserve during
                LD-clumping via the --memory flag
        verbose: whether to print details of the LD clumping process

    Returns:
        A polars DataFrame with one row per variant in an LD clump, with the
        following columns:
        - SNP_col: the variant's ID
        - ref_col: the variant's reference allele
        - alt_col: the variant's alternate allele
        - chrom_col: the variant's chromosome, in the same format as sumstats
        - bp_col: the variant's base-pair position, in the same genome build as
                  sumstats
        - f'{SNP_col}_lead': the ID of the lead variant, i.e. the lowest
                             p-value variant in the clump
        - f'{ref_col}_lead': the lead variant's reference allele
        - f'{alt_col}_lead': the lead variant's alternate allele
        - f'{bp_col}_lead': the lead variant's base-pair position, in the same
                            genome build as sumstats
        - 'is_lead': whether the variant is the lead variant at its locus
                     (equivalent to testing whether SNP_col ==
                     f'{SNP_col}_lead' and similarly for ref_col, alt_col and
                     bp_col)
        - 'clump_start': the base-pair position of the start of the clump, in
                         the same genome build as sumstats
        - 'clump_end': the base-pair position of the end of the clump, in the
                       same genome build as sumstats
    """
    if sumstats.is_empty():
        raise ValueError(f'sumstats is empty!')
    if SNP_col not in sumstats:
        raise ValueError(f'SNP_col "{SNP_col}" not in sumstats! Did you '
                         f'forget to specify a custom SNP_col?')
    if chrom_col not in sumstats:
        raise ValueError(f'chrom_col "{chrom_col}" not in sumstats! Did you '
                         f'forget to specify a custom chrom_col?')
    if bp_col not in sumstats:
        raise ValueError(f'bp_col "{bp_col}" not in sumstats! Did you '
                         f'forget to specify a custom bp_col?')
    if ref_col not in sumstats:
        raise ValueError(f'ref_col "{ref_col}" not in sumstats! Did you '
                         f'forget to specify a custom ref_col?')
    if alt_col not in sumstats:
        raise ValueError(f'alt_col "{alt_col}" not in sumstats! Did you '
                         f'forget to specify a custom alt_col?')
    if p_col not in sumstats:
        raise ValueError(f'p_col "{p_col}" not in sumstats!')
    if sumstats[p_col].min() > clump_p1:
        raise ValueError('All p-values in sumstats are > clump_p1')
    if sumstats[SNP_col].str.contains(separator).any():
        raise ValueError(f'Some variant IDs in sumstats contain the '
                         f'separator "{separator}"; specify a different '
                         f'separator via the separator argument!')
    if pfile is None and bfile is None:
        raise ValueError('Must specify either pfile or bfile (but not both)')
    if pfile is not None and bfile is not None:
        raise ValueError('Do not specify both pfile and bfile')
    # Since only one of pfile and bfile is specified, refer to it with a single
    # variable, filesets. If there's only one fileset, box it in a tuple
    filesets = pfile if pfile is not None else bfile
    if isinstance(filesets, str):
        filesets = (filesets,)
    suffix = 'pvar' if pfile is not None else 'bim'
    # Load bim/pvar files from each fileset in pfile/bfile
    bim_or_pvars = {fileset: read_bim_or_pvar(f'{fileset}.{suffix}')
                    for fileset in filesets}
    for fileset, bim_or_pvar in bim_or_pvars.items():
        if bim_or_pvar['SNP'].str.contains(separator).any():
            raise ValueError(f'Some variant IDs in the {suffix} file '
                             f'{fileset}.{suffix} contain the separator '
                             f'{separator!r}; specify a different separator '
                             f'via the separator argument!')
    # Harmonize sumstats to bim/pvar
    old_length = len(sumstats)
    sumstats = sumstats\
        .with_columns(original_ref=ref_col, original_alt=alt_col)\
        .pipe(harmonize_sumstats_to_bim_or_pvar,
              pl.concat(bim_or_pvars.values()),
              SNP_col=SNP_col, ref_col=ref_col, alt_col=alt_col,
              chrom_col=chrom_col, convert_chromosomes=True)
    num_filtered = old_length - len(sumstats)
    if verbose:
        print(f'Removing {num_filtered:,} {plural("variant", num_filtered)} '
              f'in the sumstats '
              f'({100 * num_filtered / old_length:.2f}%) '
              f'that {"was" if num_filtered == 1 else "were"} not found in '
              f'the {suffix} file')
    # Create temporary sumstats file; map chromosomes to numbers to match plink
    # Escape commas in the SNP column (i.e. when SNPs have multiple rs numbers)
    # with two separators (~~ by default) because --clump delimits its own
    # output with commas
    temp_sumstats_file = f'{cache_prefix}.sumstats.tmp'
    sumstats\
        .with_columns(pl.col(SNP_col).str.replace(',', separator + separator))\
        .with_columns(pl.concat_str(SNP_col, chrom_col, ref_col, alt_col,
                                    separator=separator))\
        .write_csv(temp_sumstats_file, separator='\t')
    # Run --clump on each fileset. The default settings are:
    # --clump-p1 5e-8: start with genome-wide-significant SNPs as lead SNPs
    # --clump-p2 0.001: clumps will include all SNPs with p < 0.001...
    # --clump-r2 0.01: ...r2 > 0.01 with the lead SNP...
    # --clump-kb 5000: ...and within 5 MB of the lead SNP
    # Create a temporary pvar file for each fileset first, with ID set to
    # SNP + '_' + CHROM + '_' + REF + '_' + ALT.
    for file_index, (fileset, bim_or_pvar) in enumerate(
            bim_or_pvars.items(), start=1):
        temp_pvar_file = f'{cache_prefix}.pvar' if len(filesets) == 1 else \
            f'{cache_prefix}_{file_index}.pvar'
        bim_or_pvar\
            .with_columns(pl.concat_str('SNP', 'CHROM', 'REF', 'ALT',
                                        separator=separator))\
            .pipe(write_bim_or_pvar, temp_pvar_file)
        run(f'plink2 '
            f'{f"--threads {num_threads} " if num_threads is not None else ""}'
            f'--memory {memory} '
            f'--seed 0 ' +
            (f'--pfile {fileset.removesuffix(".pgen")} '
             if pfile is not None else
             f'--bfile {fileset.removesuffix(".bed")} ') +
            f'--pvar {temp_pvar_file} '
            f'--no-psam-pheno '  # optimization: avoids loading phenotypes
            f'--clump {temp_sumstats_file} '
            f'--clump-p1 {clump_p1} '
            f'--clump-p2 {clump_p2} '
            f'--clump-r2 {clump_r2} '
            f'--clump-kb {clump_kb} '
            f'--clump-snp-field {SNP_col} '
            f'--clump-field {p_col} '
            f'--out {cache_prefix}' +
            ('' if len(filesets) == 1 else f'_{file_index}'))
    # Aggregate clumping results, if more than one fileset
    clumping_results_file = f'{cache_prefix}.clumps'
    if len(filesets) > 1:
        run(f"(cat '{cache_prefix}_1.clumps'; tail -qn+2 " + ' '.join(
            f"'{cache_prefix}_{file_index}.clumps'"
            for file_index in range(2, len(filesets) + 1)) +
            f") > '{clumping_results_file}'")
    # Load clumping results; map each clumped variant to its lead variant
    # (also remember to add a mapping from each lead variant to itself, and
    # unescape the double-separator back to comma at the right moment)
    if not os.path.exists(clumping_results_file):
        raise RuntimeError(f'Clumping results file {clumping_results_file} is '
                           f'empty - this is a bug in ld_clump()!')
    if os.path.exists(f'{clumping_results_file}.missing_id'):
        raise RuntimeError(f'Some variants in sumstats were missing from '
                           f'{"pfile" if pfile is not None else "bfile"} - '
                           f'this is a bug in ld_clump()')
    clumped_variants = pl.read_csv(clumping_results_file, separator='\t',
                                   columns=['ID', 'SP2'])\
        .with_columns(pl.when(pl.col.SP2 != '.').then(pl.col.SP2)
                      .str.split(','))\
        .explode('SP2')\
        .select(pl.all().str.replace(separator + separator, ','))\
        .pipe(lambda df: pl.concat([
            df.filter(pl.col.SP2 != 'NONE')
            .with_columns(SP2=pl.col.SP2.str.split_exact('(', 1)
                          .struct.field('field_0')),
            # add a mapping from each lead variant to itself
            pl.DataFrame({'ID': df['ID'].unique()}).with_columns(SP2='ID')]))\
        .rename({'ID': 'lead_variant', 'SP2': SNP_col})\
        .select(SNP_col, 'lead_variant')\
        .with_columns(pl.col(SNP_col).str.split_exact('~', 3).struct
                      .rename_fields([SNP_col, chrom_col, ref_col, alt_col]),
                      pl.col('lead_variant').str.split_exact('~', 3).struct
                      .rename_fields([f'{SNP_col}_lead', f'{chrom_col}_lead',
                                      f'{ref_col}_lead', f'{alt_col}_lead']))\
        .unnest(SNP_col, 'lead_variant')\
        .pipe(lambda df: df if sumstats[chrom_col].dtype == pl.String else df
              .with_columns(pl.col(chrom_col, f'{chrom_col}_lead').cast(int)))\
        .drop(f'{chrom_col}_lead')
    assert not clumped_variants.select(SNP_col, chrom_col, ref_col, alt_col)\
        .is_duplicated().any()
    # Add the base-pair extent of each clump
    clumped_variants = clumped_variants\
        .join(sumstats.select(SNP_col, chrom_col, ref_col, alt_col, bp_col,
                              'original_ref', 'original_alt'),
              on=(SNP_col, chrom_col, ref_col, alt_col), how='left')\
        .join(sumstats.select(SNP_col, chrom_col, ref_col, alt_col, bp_col,
                              'original_ref', 'original_alt')
              .rename({bp_col: f'{bp_col}_lead',
                       'original_ref': 'original_ref_lead',
                       'original_alt': 'original_alt_lead'}),
              left_on=(f'{SNP_col}_lead', chrom_col, f'{ref_col}_lead',
                       f'{alt_col}_lead'),
              right_on=(SNP_col, chrom_col, ref_col, alt_col), how='left')\
        .with_columns(clump_start=pl.min(bp_col).over(f'{SNP_col}_lead'),
                      clump_end=pl.max(bp_col).over(f'{SNP_col}_lead'),
                      is_lead=pl.col(SNP_col).eq(pl.col(f'{SNP_col}_lead')) &
                              pl.col(ref_col).eq(pl.col(f'{ref_col}_lead')) &
                              pl.col(alt_col).eq(pl.col(f'{alt_col}_lead')) &
                              pl.col(bp_col).eq(pl.col(f'{bp_col}_lead')))\
        .with_columns(pl.col('original_ref').alias(ref_col),
                      pl.col('original_alt').alias(alt_col),
                      pl.col('original_ref_lead').alias(f'{ref_col}_lead'),
                      pl.col('original_alt_lead').alias(f'{alt_col}_lead'))\
        .select(SNP_col, ref_col, alt_col, chrom_col, bp_col,
                f'{SNP_col}_lead', f'{ref_col}_lead', f'{alt_col}_lead',
                f'{bp_col}_lead', 'is_lead', 'clump_start', 'clump_end')
    assert clumped_variants.null_count().sum_horizontal().item() == 0
    # Sort by chromosome and base-pair position
    clumped_variants = clumped_variants.pipe(
        sort_sumstats, chrom_col=chrom_col, bp_col=bp_col)
    # Clean up - but only at the end, so users can inspect what went wrong in
    # case of an error
    run(f"rm '{temp_sumstats_file}'")
    if len(filesets) == 1:
        run(f"rm '{cache_prefix}.pvar' '{cache_prefix}.clumps' "
            f"'{cache_prefix}.log'")
    else:
        run(f"rm '{cache_prefix}.clumps' " +
            ' '.join(f"'{cache_prefix}_{file_index}.{suffix}'"
                     for file_index in range(1, len(filesets) + 1)
                     for suffix in ('pvar', 'clumps', 'log')))
    return clumped_variants


def sort_sumstats(sumstats, *, chrom_col='CHROM', bp_col='BP'):
    """
    Sorts summary statistics by chromosome and base-pair position.
    
    Args:
        sumstats: a polars DataFrame of summary statistics
        chrom_col: the name of the chromosome column in sumstats
        bp_col: the name of the base-pair position column in sumstats

    Returns:
        sumstats, sorted by chromosome and base-pair position.
    """
    assert 'numeric_chrom' not in sumstats
    return sumstats\
        .with_columns(numeric_chrom=pl.col(chrom_col).pipe(
            standardize_chromosomes, return_numeric=True))\
        .sort('numeric_chrom', bp_col)\
        .drop('numeric_chrom')


# noinspection PyShadowingBuiltins
def munge_sumstats(raw_sumstats_files, munged_sumstats_file, REF, ALT, P, *,
                   SNP=None, CHROM=None, BP=None, BETA=None, OR=None, SE=None,
                   N=None, N_CASES=None, N_CONTROLS=None, NEFF=None, AAF=None,
                   AAF_CASES=None, AAF_CONTROLS=None, INFO=None,
                   quantitative=False, dbSNP=None, preamble=None,
                   separator='\t', verbose=True, **read_csv_kwargs):
    """
    "Munges" the sumstats file(s) raw_sumstats_files into a consistent format,
    saving to munged_sumstats_file (if not None) and returning the sumstats.
    
    If dbSNP is not None, infers rs numbers. You should always do this if CHROM
    and BP are available, even if the sumstats already came with rs numbers!
    
    Make sure to set schema_overrides={'CHROM': str} explicitly when the
    sumstats include sex chromosomes and there is no chr prefix (e.g. when it's
    '1' and 'X' rather than 'chr1' and 'chrX').
    
    Performs the following steps, in order, on each file in raw_sumstats_files:
    1. Reads in the file with pl.read_csv(). Use separator to specify the
       separator (tab by default, not comma!) and **kwargs to specify any other
       arguments you want. Use separator='whitespace' to use any number of
       consecutive whitespace characters as the separator.
    2. Creates columns for each of the arguments from SNP to INFO. Specify a
       single column name or polars expression of other columns, e.g.:
       - AAF='allele_frequency'
       - AAF=(pl.col.FCAS * pl.col.NCAS + pl.col.FCON * pl.col.NCON) /
             (pl.col.NCAS + pl.col.NCON)
       - AAF=pl.when(pl.col.A1 == pl.col.MINOR_ALLELE).then(pl.col.MAF)
             .otherwise(1 - pl.col.MAF)
       Certain columns are required: REF, ALT, P, CHROM + BP (if dbSNP is not
       None), SNP (if dbSNP is None), and either BETA or OR (but not both). N
       is required for quantitative traits (quantitative=True) and at least one
       of NEFF, AAF, and N_CASES + N_CONTROLS is required for case-control
       traits. Optional columns that are None won't be included in the output
       unless they can be inferred from columns that are specified.
    3. QC: removes variants with missing data in any column, non-ACGT alleles,
       P outside (0, 1], OR <= 0, SE <= 0, N/N_CASES/N_CONTROLS/NEFF <= 0, AAF/
       AAF_CASES/AAF_CONTROLS outside (0, 1), INFO outside (0, 1].
       If verbose=True, prints stats on how many were removed.
    4. Standardizes chromosomes to chr1, chr2, ... chr22, chrX, chrY
    5. Converts variants to their minimal representations with
       get_minimal_representations().
    6. If dbSNP is not None, infers rs numbers based on CHROM and BP. If
       verbose=True, prints the number and % of variants that had rs numbers in
       dbSNP and the number and % that had to be flipped in order to match; if
       SNP is not None as well, prints the number and % of rs numbers (variant
       IDs starting with rs) in the SNP column that matched the ones from
       dbSNP. These messages are useful for catching errors: for instance, you
       may have specified a different genome build of dbSNP than the chrom/bp
       positions in raw_sumstats_files, in which case most variants won't have
       rs numbers, or you may have erroneously specified the alt columns as ref
       and vice versa.
    7. If dbSNP is not None, flips alleles for single-nucleotide variants if
       necessary to match dbSNP, and sets BETA = -BETA, OR = 1 / OR, AAF = 1 -
       AAF, AAF_CASES = 1 - AAF_CASES, and AAF_CONTROLS = 1 - AAF_CONTROLS for
       these variants. Only REF/ALT flips are considered, NOT strand flips,
       because modern GWAS datasets don't have strand flips - e.g. there isn't
       a single non-ambiguous variant (i.e. not A/T, T/A, C/G, or G/C) in the
       Als et al. 2023 depression GWAS or the Watson et al. 2019 anorexia GWAS
       that matches dbSNP with a strand flip.
    8. Removes variants without rs numbers in dbSNP and variants that aren't
       unique (based on SNP + REF + ALT, + CHROM/BP if they're not None).
    9. Sorts by CHROM and then BP, if those are not None.
    10. Saves to munged_sumstats_file, unless munged_sumstats_file is None. If
        munged_sumstats_file ends in .gz, sumstats will be block-gzipped with
        bgzip for convenience.
    11. Returns the sumstats.

    For example, let's munge the depression sumstats on Niagara at
    /scratch/w/wainberg/wainberg/sumstats/daner_MDDwoBP_20201001_2015iR15iex_
    HRC_MDDwoBP_iPSYCH2015i_UKBtransformed_Wray_FinnGen_MVPaf_2_HRC_MAF01.gz.
    
    The header line is:
    CHR SNP BP A1 A2 FRQ_A_294322 FRQ_U_741438 INFO OR SE P ngt Direction \
        HetISqt HetDf HetPVa Nca Nco Neff_half
    
    So, filling in all the matching columns from left to right, we have:
    munge_sumstats(
        raw_sumstats_files='sumstats/daner_MDDwoBP_20201001_2015iR15iex_HRC_'
                           'MDDwoBP_iPSYCH2015i_UKBtransformed_Wray_FinnGen_'
                           'MVPaf_2_HRC_MAF01.gz',
        munged_sumstats_file='MD.gz',
        CHROM='CHR', SNP='SNP', BP='BP', ALT='A1', REF='A2', INFO='INFO',
        OR='OR', SE='SE', P='P', N_CASES='Nca', N_CONTROLS='Nco', ...
    
    To complete the picture, we must also specify AAF, which is a function of
    the case and control Ns and frequencies. We can also specify NEFF, which
    is two times the Neff_half column. We must also specify that the input
    sumstats are space-delimited:
        ..., AAF=(pl.col.FRQ_A_294322 * pl.col.Nca + pl.col.FRQ_U_741438 *
                  pl.col.Nco) / (pl.col.Nca + pl.col.Nco)',
        NEFF=2 * pl.col.Neff_half, separator=' ')

    Args:
        raw_sumstats_files: the filename of the input sumstats, or a
                            list/tuple/etc. of filenames to be concatenated
        munged_sumstats_file: if not None, the output sumstats filename to save
                              to
        REF: the reference allele column/expression; mandatory
        ALT: the alternate allele column/expression; mandatory
        P: the p-value column/expression; mandatory
        SNP: The variant ID column in raw_sumstats_files, or a polars
             expression to generate a variant ID column; mandatory when dbSNP
             is None.
        CHROM: the chromosome column/expression; mandatory when dbSNP is not
               None or BP is not None
        BP: the base-pair position column/expression; mandatory when dbSNP is
            not None or CHROM is not None
        BETA: The effect size column/expression. Specify exactly one of BETA
              and OR; if OR is specified, BETA will be calculated as log(OR).
        OR: The odds ratio column/expression. Specify exactly one of BETA/OR.
            Not explicitly disallowed for quantitative traits!
        SE: The standard error column/expression. If None, will be back-
            calculated from BETA and P: SE = |BETA| / p_to_abs_z(P).
            p_to_abs_z() is an awk implementation of utils.py's p_to_abs_z()
            function. log(OR) is substituted for BETA if BETA is None.
        N: The sample size column/expression. Mandatory for quantitative traits
           (quantitative=True). For case-control traits (quantitative=False),
            will be calculated as N_CASES + N_CONTROLS if both N_CASES and
            N_CONTROLS are not None. Can be a number, in which case N is
            assumed to be that number for all variants.
        N_CASES: The number of cases column/expression. Must be None for
                 quantitative traits (quantitative=True). Can be a number, in
                 which case N_CASES is assumed to be the same for all variants.
        N_CONTROLS: The number of controls column/expression. Must be None for
                    quantitative traits (quantitative=True). Can be a number,
                    in which case N_CONTROLS is assumed to be the same for all
                    variants.
        NEFF: The effective sample size column/expression. For quantitative
              traits (quantitative=True), will be set to N if missing. For
              case-control traits (quantitative=False), will be set to
              ((4 / (2 * AAF * (1 - AAF) * INFO)) - BETA^2) / SE^2 if AAF is
              not None (substituting log(OR) for BETA if BETA is None, and
              skipping the "* INFO" if INFO is None), and
              4 / (1 / N_CASES + 1 / N_CONTROLS) if both N_CASES and N_CONTROLS
              are not None; see comment below for justification. As a result,
              for case-control traits, either NEFF or AAF or both N_CASES and
              N_CONTROLS are mandatory. Can be a number, in which case NEFF is
              assumed to be the same for all variants.
        AAF: the alternate allele frequency column/expression
        AAF_CASES: the case alternate allele frequency column/expression
        AAF_CONTROLS: the control alternate allele frequency column/expression
        INFO: the imputation INFO score column/expression
        quantitative: is the trait quantitative (True) or case-control (False)?
        dbSNP: a DataFrame returned by load_dbSNP(); must match the genome
               build of raw_sumstats_files' BP column! If not specified, do not
               infer rs numbers.
        preamble: an optional function to apply immediately after loading each
                  file in raw_sumstats_files (so use the original column names
                  in the expression, not the new column names). For instance,
                  to filter to variants with minor allele frequencies between
                  0.01 and 0.99, use preamble=lambda df: df.filter(
                  pl.col.all_maf.between(0.01, 0.99)).
        separator: the separator to use when parsing each raw sumstats file
                   with pl.read_csv(), or 'whitespace' to use any number of
                   consecutive whitespace characters as the separator, like
                   delim_whitespace=True in pandas.read_table().
        verbose: whether to print details of the munging process
        **read_csv_kwargs: keyword arguments to pl.read_csv(). By default,
                           null_values='NA' is set.
    
    Returns:
        A DataFrame of the munged sumstats.
    """
    # Check inputs
    if REF is None:
        raise ValueError('REF is always mandatory')
    if ALT is None:
        raise ValueError('ALT is always mandatory')
    if P is None:
        raise ValueError('P is always mandatory')
    if dbSNP is None and SNP is None:
        raise ValueError('SNP is mandatory when dbSNP is None')
    if CHROM is None and BP is not None:
        raise ValueError('CHROM is mandatory when BP is not None')
    if CHROM is not None and BP is None:
        raise ValueError('BP is mandatory when CHROM is not None')
    if dbSNP is not None and BP is None:
        raise ValueError('CHROM and BP are mandatory when dbSNP is not None')
    if BETA is None and OR is None:
        raise ValueError('Neither BETA nor OR specified; specify exactly one')
    if BETA is not None and OR is not None:
        raise ValueError('Both BETA and OR specified; specify exactly one')
    if quantitative and N is None:
        raise ValueError('N is mandatory when quantitative=True')
    if quantitative and N_CASES is not None:
        raise ValueError('N_CASES cannot be specified when quantitative=True')
    if quantitative and N_CONTROLS is not None:
        raise ValueError('N_CONTROLS cannot be specified when '
                         'quantitative=True')
    if not quantitative and NEFF is None and AAF is None and \
            (N_CASES is None or N_CONTROLS is None):
        raise ValueError('quantitative=True and NEFF is None, so either AAF '
                         'or both of N_CASES and N_CONTROLS need to be '
                         'specified to infer it')
    # Wrap string arguments in pl.col()
    if isinstance(SNP, str):
        SNP = pl.col(SNP)
    if isinstance(REF, str):
        REF = pl.col(REF)
    if isinstance(ALT, str):
        ALT = pl.col(ALT)
    if isinstance(P, str):
        P = pl.col(P)
    if isinstance(CHROM, str):
        CHROM = pl.col(CHROM)
    if isinstance(BP, str):
        BP = pl.col(BP)
    if isinstance(BETA, str):
        BETA = pl.col(BETA)
    if isinstance(OR, str):
        OR = pl.col(OR)
    if isinstance(SE, str):
        SE = pl.col(SE)
    if isinstance(N, str):
        N = pl.col(N)
    if isinstance(N_CASES, str):
        N_CASES = pl.col(N_CASES)
    if isinstance(N_CONTROLS, str):
        N_CONTROLS = pl.col(N_CONTROLS)
    if isinstance(NEFF, str):
        NEFF = pl.col(NEFF)
    if isinstance(AAF, str):
        AAF = pl.col(AAF)
    if isinstance(AAF_CASES, str):
        AAF_CASES = pl.col(AAF_CASES)
    if isinstance(AAF_CONTROLS, str):
        AAF_CONTROLS = pl.col(AAF_CONTROLS)
    if isinstance(INFO, str):
        INFO = pl.col(INFO)
    # If SE is missing, back-calculate it from BETA and P: |Z| = |BETA| / SE,
    # so SE = |BETA| / |Z|. We can get |Z| from P via p_to_abs_z().
    if SE is None:
        # noinspection PyUnresolvedReferences
        SE = (BETA if BETA is not None else OR.log()).abs() / p_to_abs_z(P)
    # Try to infer N and NEFF if not specified
    #
    # For quantitative traits, NEFF = N. For binary traits, NEFF can be given
    # in two mathematically equivalent ways: 4/(1/N_CASES + 1/N_CONTROLS)
    # (cell.com/ajhg/pdfExtended/S0002-9297(21)00145-2) and 4v(1-v)N where v =
    # N_CASES/N (medrxiv.org/content/10.1101/2021.09.22.21263909v1.full-text).
    #
    # But this doesn't account for varying case-control ratios across the
    # cohorts in a meta-analysis, which biases estimates that depend on N
    # (medrxiv.org/content/10.1101/2021.09.22.21263909v1.full).
    #
    # Instead, biorxiv.org/content/10.1101/2021.03.29.437510v4.full suggests
    # a per-variant Neff = (4 / (2 * MAF * (1 - MAF) * INFO) - BETA^2) / SE^2,
    # bounded to between 0.5 and 1.1 times "the total (effective) sample size",
    # which seems to be max(4/(1/N_CASES + 1/N_CONTROLS)) across SNPs:
    # github.com/privefl/paper-misspec/blob/main/code/investigate-misspec-N.R
    # To get a global Neff, they take the 80th %ile of the per-variant Neffs.
    # This formula is based on Equation 1 of the "New formula used in LDpred2"
    # section of sciencedirect.com/science/article/abs/pii/S0002929721004201.
    #
    # github.com/GenomicSEM/GenomicSEM/wiki/2.1-Calculating-Sum-of-Effective-
    # Sample-Size-and-Preparing-GWAS-Summary-Statistics drops the BETA^2 and
    # INFO, i.e. Neff = 4 / (2 * MAF * (1 - MAF)) / SE^2, bounded to between
    # 0.5 and 1.1 times 4/(1/N_CASES + 1/N_CONTROLS), where N_CASES/N_CONTROLS
    # are again global rather than per-SNP. This formula is based on an old
    # version of the above formula that was used in the original LDpred2 paper.
    #
    # Here, we use Neff = (4 / (2 * AAF * (1 - AAF) * INFO) - BETA^2) / SE^2
    # (dropping the INFO part if not available), without bounding (because it
    # seems like a hack). If MAF not available, but N_CASES and N_CONTROLS are
    # available, fall back to using 4/(1/N_CASES + 1/N_CONTROLS). Note that
    # even if AAF = 1 - MAF instead of MAF, AAF * (1 - AAF) == MAF * (1 - MAF).
    # "((4 / (2 * ..." has been simplified to "((2 / (..." below.
    if N is None and N_CASES is not None and N_CONTROLS is not None:
        N = N_CASES + N_CONTROLS
    if NEFF is None:
        if quantitative and N is not None:
            NEFF = N
        elif AAF is not None:
            if INFO is not None:
                NEFF = (2 / (AAF * (1 - AAF) * INFO) - (
                    BETA if BETA is not None else OR.log()) ** 2) / SE ** 2
            else:
                NEFF = (2 / (AAF * (1 - AAF)) - (
                    BETA if BETA is not None else OR.log()) ** 2) / SE ** 2
        elif N_CASES is not None and N_CONTROLS is not None:
            NEFF = 4 / (1 / N_CASES + 1 / N_CONTROLS)
        else:
            raise ValueError('NEFF is unexpectedly None; internal error')
    # For case-control traits, output ORs even if input has betas, & vice versa
    if not quantitative and OR is None:
        OR = BETA.exp()
    if quantitative and BETA is None:
        BETA = OR.log()
    # Ensure REF and ALT are capitalized
    REF = REF.str.to_uppercase()
    ALT = ALT.str.to_uppercase()
    # Enumerate columns to include, and the formula for each column
    column_formulas = {
        column_name: formula for column_name, formula in {
            'SNP': SNP, 'CHROM': CHROM, 'BP': BP, 'REF': REF, 'ALT': ALT,
            'AAF': AAF, 'AAF_CASES': AAF_CASES, 'AAF_CONTROLS': AAF_CONTROLS,
            'INFO': INFO, 'BETA' if quantitative else 'OR':
                BETA if quantitative else OR, 'SE': SE, 'P': P, 'N': N,
            'N_CASES': N_CASES, 'N_CONTROLS': N_CONTROLS, 'NEFF': NEFF}.items()
        if formula is not None}
    # Load files in raw_sumstats_files; apply the preamble function (if not
    # None) and select columns
    raw_sumstats_files = (raw_sumstats_files,) \
        if isinstance(raw_sumstats_files, str) else tuple(raw_sumstats_files)
    columns = tuple(set.union(*(set(formula.meta.root_names())
                                if isinstance(formula, pl.Expr) else formula
                                for formula in column_formulas.values()
                                if not isinstance(formula, (int, float)))))
    default_read_csv_kwargs = dict(null_values='NA')
    read_csv_kwargs = default_read_csv_kwargs | read_csv_kwargs \
        if read_csv_kwargs is not None else default_read_csv_kwargs
    sumstats = pl.concat((
        read_csv_delim_whitespace(raw_sumstats_file, columns=columns,
                                  **read_csv_kwargs)
        if separator == 'whitespace' else
        # Use scan_csv() instead of read_csv() when possible, i.e. when the raw
        # sumstats file isn't gzipped and delim_whitespace=False.
        pl.read_csv(raw_sumstats_file, columns=columns, separator=separator,
                    **read_csv_kwargs)
        if raw_sumstats_file.endswith('.gz') else
        pl.scan_csv(raw_sumstats_file, separator=separator, **read_csv_kwargs)
        .select(columns))
        .lazy()
        .pipe(preamble if preamble is not None else lambda df: df)
        .select(**column_formulas)
        for raw_sumstats_file in raw_sumstats_files)\
        .collect()
    # QC: remove variants with missing data in any column, non-ACGT alleles, P
    # outside (0, 1], OR <= 0, SE <= 0, N/N_CASES/N_CONTROLS/NEFF <= 0, AAF/
    # AAF_CASES/AAF_CONTROLS outside (0, 1), INFO outside (0, 1].
    # If verbose=True, print stats on how many were removed.
    if verbose:
        num_initial_variants = len(sumstats)
    # noinspection PyUnresolvedReferences
    filters = {filter_name: filter for filter_name, filter in ({
        'non-ACGT alleles': sumstats['REF'].str.contains('[^ACGT]') |
                            sumstats['ALT'].str.contains('[^ACGT]'),
        'P outside (0, 1]': ~sumstats['P'].is_between(0, 1, closed='right'),
        'OR <= 0': sumstats['OR'] <= 0 if 'OR' in sumstats else None,
        'SE <= 0': sumstats['SE'] <= 0 if 'SE' in sumstats else None,
        'N <= 0': sumstats['N'] <= 0 if 'N' in sumstats else None,
        'N_CASES <= 0': sumstats['N_CASES'] <= 0
                         if 'N_CASES' in sumstats else None,
        'N_CONTROLS <= 0': sumstats['N_CONTROLS'] <= 0
                           if 'N_CONTROLS' in sumstats else None,
        'NEFF <= 0': sumstats['NEFF'] <= 0 if 'NEFF' in sumstats else None,
        'AAF outside (0, 1)': ~sumstats['AAF'].is_between(0, 1, closed='none')
                              if 'AAF' in sumstats else None,
        'AAF_CASES outside (0, 1)': ~sumstats['AAF_CASES'].is_between(
            0, 1, closed='none') if 'AAF_CASES' in sumstats else None,
        'AAF_CONTROLS outside (0, 1)': ~sumstats['AAF_CONTROLS'].is_between(
            0, 1, closed='none') if 'AAF_CONTROLS' in sumstats else None,
        'INFO <= 0': sumstats['INFO'] <= 0 if 'INFO' in sumstats else None} | {
        f'null {column_name}': sumstats[column_name].is_null()
        for column_name in column_formulas}).items()
        if filter is not None and filter.any()}
    if len(filters) > 0:
        union_of_filters = reduce(lambda a, b: a | b, filters.values())
        total_num_filtered = union_of_filters.sum()
        if total_num_filtered == len(sumstats):
            raise ValueError('All variants would be filtered out!')
        if total_num_filtered > 0:
            if verbose:
                last_filter_name = tuple(filters)[-1] \
                    if len(filters) > 1 else None
                print(f'Removing ' + ', '.join(
                    f'{"and " if filter_name == last_filter_name else ""}'
                    f'{num_filtered:,} {plural("variant", num_filtered)} with '
                    f'{filter_name}'
                    for filter_name, filter, num_filtered in
                    ((filter_name, filter, filter.sum())
                     for filter_name, filter in filters.items())) +
                      f', for a total of {total_num_filtered:,} '
                      f'{plural("variant", total_num_filtered)} '
                      f'({100 * total_num_filtered / len(sumstats):.2f}%)')
            sumstats = sumstats.filter(~union_of_filters)
    elif verbose:
        print('All variants pass initial QC filters!')
    # Standardize chromosomes
    if CHROM is not None:
        sumstats = sumstats\
            .with_columns(pl.col.CHROM.pipe(standardize_chromosomes))
    # Convert variants to their minimal representations
    sumstats = sumstats.pipe(get_minimal_representations,
                             bp_col='BP' if BP is not None else None)
    # If dbSNP is not None, infer rs numbers and flip variants to match dbSNP;
    # if verbose=True, print stats
    if dbSNP is not None:
        if verbose and SNP is not None:
            sumstats = sumstats.rename({'SNP': 'original_SNP'})
        sumstats = sumstats\
            .pipe(get_rs_numbers, dbSNP=dbSNP, verbose=verbose)\
            .select('SNP', pl.exclude('SNP'))\
            .pipe(flip_alleles)
        if verbose:
            num_nonmissing_rs = sumstats['SNP'].is_not_null().sum()
            percent_nonmissing_rs = 100 * num_nonmissing_rs / len(sumstats)
            num_flipped = sumstats['FLIP'].sum()
            percent_flipped = 100 * num_flipped / num_nonmissing_rs
            print(f'{num_nonmissing_rs:,} variants '
                  f'({percent_nonmissing_rs:.2f}%) had rs numbers in dbSNP, '
                  f'of which {num_flipped:,} ({percent_flipped:.2f}%) had to '
                  f'be flipped in order to match')
            if SNP is not None:
                matching_rs = sumstats\
                    .select('original_SNP', 'SNP')\
                    .filter(pl.col.original_SNP.str.starts_with('rs'),
                            pl.col.SNP.is_not_null())\
                    .with_columns(pl.col.SNP.str.split(','))\
                    .with_columns(match=pl.col.SNP.list.contains(
                        pl.col.original_SNP))
                num_matching_rs = matching_rs['match'].sum()
                percent_matching_rs = 100 * num_matching_rs / len(matching_rs)
                likely_merges = matching_rs\
                    .filter(~pl.col.match, ~pl.col.SNP.list.len().eq(1))\
                    .select(pl.col.SNP.list.get(0).str.slice(2).cast(int) <
                            pl.col.original_SNP.str.slice(2).cast(int))\
                    .to_series()
                num_likely_merges = likely_merges.sum()
                percent_likely_merges = 100 * num_likely_merges / \
                                        len(likely_merges)
                print(f'{len(matching_rs):,} variant IDs in the SNP column/'
                      f'expression you specified start with "rs" and have at '
                      f'least one rs number in dbSNP; {num_matching_rs:,} '
                      f'({percent_matching_rs:.2f}%) of these match dbSNP. Of '
                      f'those that didn\'t, {len(likely_merges):,} have '
                      f'exactly one matching rs number, and for '
                      f'{num_likely_merges:,} ({percent_likely_merges:.2f}%) '
                      f'of these, the matching rs number is numerically '
                      f'smaller than the original, which suggests the two rs '
                      f'numbers might have been merged.')
                sumstats = sumstats.drop('original_SNP')
        sumstats = sumstats.drop('FLIP')
    # Remove variants without rs numbers in dbSNP
    variants_with_rs_numbers = sumstats['SNP'].is_not_null()
    num_without_rs_numbers = len(sumstats) - variants_with_rs_numbers.sum()
    if num_without_rs_numbers > 0:
        if verbose:
            print(f'Removing {num_without_rs_numbers:,} '
                  f'{plural("variant", num_without_rs_numbers)} '
                  f'({100 * num_without_rs_numbers / len(sumstats):.2f}%) '
                  f'without rs numbers in dbSNP')
        sumstats = sumstats.filter(variants_with_rs_numbers)
    # Remove non-unique variants
    variant_columns = ['SNP', 'REF', 'ALT'] + \
                      (['CHROM', 'BP'] if BP is not None else [])
    unique_mask = sumstats.select(variant_columns).is_unique()
    num_non_unique = len(sumstats) - unique_mask.sum()
    if num_non_unique > 0:
        if verbose:
            print(f'Removing {num_non_unique:,} non-unique '
                  f'{plural("variant", num_non_unique)} '
                  f'({100 * num_non_unique / len(sumstats):.2f}%)')
        sumstats = sumstats.filter(unique_mask)
    # Print how many variants were retained
    if verbose:
        # noinspection PyUnboundLocalVariable
        print(f'{len(sumstats):,} of {num_initial_variants:,} variants '
              f'({100 * len(sumstats) / num_initial_variants:.2f}%) were '
              f'retained')
    # Sort
    if BP is not None:
        sumstats = sumstats.pipe(sort_sumstats)
    # Save, if munged_sumstats_file is not None; block-gzip if
    # munged_sumstats_file ends in .gz
    if munged_sumstats_file is not None:
        sumstats\
            .with_columns(pl.selectors.float()
                          .map_elements('{:.12g}'.format,
                                        return_dtype=pl.String))\
            .write_csv(munged_sumstats_file.removesuffix('.gz'),
                       separator='\t')
        if munged_sumstats_file.endswith('.gz'):
            run(f'bgzip -f {munged_sumstats_file.removesuffix(".gz")}')
    # Return the sumstats, even if saving
    return sumstats
    

def munge_regenie_sumstats(raw_sumstats_files, munged_sumstats_file, *,
                           quantitative=False, dbSNP=None):
    """
    A wrapper for munge_sumstats() when all sumstats are in regenie format.
    
    Args:
        raw_sumstats_files: a sumstats file (or list thereof) in regenie format
        munged_sumstats_file: a file path where the munged sumstats file will
                              be output
        quantitative: are the sumstats files for quantitative traits?
        dbSNP: a DataFrame returned by load_dbSNP(); must match the genome
               build of raw_sumstats_files' BP column! If not specified, do not
               infer rs numbers.
    """
    munge_sumstats(raw_sumstats_files, munged_sumstats_file,
                   CHROM=pl.col.CHROM.cast(pl.String).replace({'23': 'X'}),
                   BP='GENPOS', SNP='ID', REF='ALLELE0', ALT='ALLELE1',
                   AAF='A1FREQ',
                   AAF_CASES=None if quantitative else 'A1FREQ_CASES',
                   AAF_CONTROLS=None if quantitative else 'A1FREQ_CONTROLS',
                   INFO='INFO', N='N',
                   N_CASES=None if quantitative else 'N_CASES',
                   N_CONTROLS=None if quantitative else 'N_CONTROLS',
                   BETA='BETA', SE='SE', P=10 ** -pl.col.LOG10P,
                   quantitative=quantitative,
                   separator=' ', dbSNP=dbSNP)


###############################################################################
# [11] rs number mapping
###############################################################################


def get_minimal_representations_numpy(ref_vector, alt_vector, bp_vector):
    """
    Converts variants - represented as 1D NumPy arrays ref_vector, alt_vector,
    and optionally bp_vector - to their minimal representations.

    Modifies bp_vector, ref_vector and alt_vector in-place.
    
    cureffi.org/2014/04/24/converting-genetic-variants-to-their-minimal-
    representation explains what a minimal representation is.
    
    github.com/ericminikel/minimal_representation/blob/master/
    minimal_representation.py contains the code this function is based on.
    
    This operation only affects indels: by definition, single-nucleotide
    variants are already in their minimal representation.
    
    Args:
        ref_vector: the reference alleles of the variants
        alt_vector: the alternate alleles of the variants
        bp_vector: the base pairs of the variants; may be None
    """
    import numpy as np
    assert isinstance(ref_vector, np.ndarray) and ref_vector.ndim == 1
    assert isinstance(alt_vector, np.ndarray) and alt_vector.ndim == 1
    assert len(alt_vector) == len(ref_vector)
    if bp_vector is not None:
        assert isinstance(bp_vector, np.ndarray) and bp_vector.ndim == 1
        assert len(bp_vector) == len(ref_vector)
        assert bp_vector.dtype in ('int32', int), bp_vector.dtype
        int_type = 'int' if bp_vector.dtype == 'int32' else 'long'
    # noinspection PyUnboundLocalVariable
    cython_inline(rf'''
    def get_minimal_representations_cython(str[:] ref_vector, str[:] alt_vector
            {f', {int_type}[:] bp_vector' if bp_vector is not None else ''}):
        cdef size_t i, min_len, ref_len, alt_len, ref_start, ref_end, \
            alt_start, alt_end
        {f'cdef {int_type} bp' if bp_vector is not None else ''}
        cdef str ref, alt
        for i in range(ref_vector.shape[0]):
            ref = ref_vector[i]
            alt = alt_vector[i]
            ref_len = len(ref)
            alt_len = len(alt)
            min_len = ref_len if ref_len < alt_len else alt_len
            if min_len == 1:
                continue
            {f'bp = bp_vector[i]' if bp_vector is not None else ''}
            ref_start = 0
            alt_start = 0
            ref_end = ref_len - 1
            alt_end = alt_len - 1
            while ref[ref_end] == alt[alt_end]:
                ref_end -= 1
                alt_end -= 1
                min_len -= 1
                if min_len == 1: break
            else:
                while ref[ref_start] == alt[alt_start]:
                    ref_start += 1
                    alt_start += 1
                    {'bp += 1' if bp_vector is not None else ''}
                    min_len -= 1
                    if min_len == 1: break
            ref_vector[i] = ref[ref_start:ref_end + 1]
            alt_vector[i] = alt[alt_start:alt_end + 1]
            {'bp_vector[i] = bp' if bp_vector is not None else ''}
    ''')['get_minimal_representations_cython'](
        **{'ref_vector': ref_vector, 'alt_vector': alt_vector} |
          (({'bp_vector': bp_vector}) if bp_vector is not None else {}))


def get_minimal_representations(df, *, ref_col='REF', alt_col='ALT',
                                bp_col='BP'):
    """
    A wrapper for get_minimal_representations_numpy() for polars DataFrames.
    
    Args:
        df: a polars DataFrame
        ref_col: the name of the reference allele column in df
        alt_col: the name of the alternate allele column in df
        bp_col: the name of the base-pair column in df; may be None

    Returns:
        A same-sized df with each variant converted to its minimal
        representation.
    """
    import pyarrow as pa
    if not isinstance(df, pl.DataFrame):
        raise ValueError('df must be a polars DataFrame!')
    if df.is_empty():
        raise ValueError(f'df is empty!')
    if ref_col not in df:
        raise ValueError(f'{ref_col!r} not in df; specify ref_col')
    if alt_col not in df:
        raise ValueError(f'{alt_col!r} not in df; specify alt_col')
    if bp_col is not None and bp_col not in df:
        raise ValueError(f'{bp_col!r} not in df; specify bp_col')
    ref_vector = df[ref_col].to_numpy()
    alt_vector = df[alt_col].to_numpy()
    bp_vector = df[bp_col].to_numpy(writable=True) \
        if bp_col is not None else None
    get_minimal_representations_numpy(ref_vector, alt_vector, bp_vector)
    df = df.with_columns(pl.from_arrow(pa.array(
                            ref_vector, type=pa.large_utf8())).alias(ref_col),
                         pl.from_arrow(pa.array(
                            alt_vector, type=pa.large_utf8())).alias(alt_col),
                         **{bp_col: pl.from_numpy(bp_vector)[:, 0]}
                         if bp_col is not None else {})
    return df


def get_minimal_representations_awk(include_bp=True):
    """
    An awk implementation of get_minimal_representations(). Assumes the
    variant's reference allele, alternate allele and base-pair position are
    stored in the variables ref, alt and bp.
    
    Args:
        include_bp: if False, do not include the correction for base-pair
                    position (include_bp=False is useful when there's no
                    base-pair column)
    
    Returns:
        A code string to be integrated into a larger awk command.
    """
    import re
    return re.sub(r'\s+', ' ', f'''
        if (length(ref) > 1 || length(alt) > 1) {{
            while (length(alt) > 1 && length(ref) > 1 &&
                   substr(alt, length(alt)) == substr(ref, length(ref))) {{
                alt = substr(alt, 1, length(alt) - 1);
                ref = substr(ref, 1, length(ref) - 1);
            }}
            while (length(alt) > 1 && length(ref) > 1 &&
                   substr(alt, 1, 1) == substr(ref, 1, 1)) {{
                alt = substr(alt, 2);
                ref = substr(ref, 2);
                {"bp++;" if include_bp else ""}
            }}
        }}
        ''').strip().rstrip()


def load_dbSNP(genome_build, *,
               dbSNP_dir=f'{get_base_data_directory()}/dbSNP'):
    """
    Load all autosomal + chrX/Y/M variants in dbSNP, caching intermediate and
    final results in the cache directory dbSNP_dir. Multiallelic variants are
    split, then variants are converted to their minimal representations.
    Variants that map to multiple genomic locations
    (ncbi.nlm.nih.gov/snp/docs/rs_multi_mapping) are removed.
    
    Must be run on a node with a large amount of memory! 120 GiB should be
    sufficient. (You will need more - closer to 188 GiB, the size of a
    Niagara compute node - if the dbSNP cache has not been generated. Also,
    you will need to first run it on a login node, to download the files, then
    when it runs out of memory, switch to a compute node. It is hacky.)
    
    Args:
        genome_build: the genome build (hg38 or hg19)
        dbSNP_dir: the cache directory to store intermediate and final results.
                   Because generating the cache can take a long time, you will
                   probably want to leave this argument at its default value.

    Returns:
        A polars DataFrame with columns CHROM, BP, REF, ALT and RSID.
    """
    check_valid_genome_build(genome_build)
    dbSNP_file = os.path.join(dbSNP_dir, f'{genome_build}.tsv')
    read_csv_kwargs = dict(separator='\t', schema_overrides={
        'CHROM': pl.Categorical, 'BP': pl.Int32,
        'REF': pl.Categorical('lexical'), 'ALT': pl.Categorical('lexical')})
    if os.path.exists(dbSNP_file):
        return pl.read_csv(dbSNP_file, **read_csv_kwargs)
    else:
        full_dbSNP_file = os.path.join(
            dbSNP_dir, f'GCF_000001405.'
                       f'{40 if genome_build == "hg38" else 25}.gz')
        if not os.path.exists(full_dbSNP_file):
            raise_error_if_on_compute_node()
            run(f'mkdir -p {dbSNP_dir} && '
                f'wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/'
                f'{os.path.basename(full_dbSNP_file)} -O {full_dbSNP_file}')
        dbSNP_file_with_multimapping = os.path.join(
            dbSNP_dir, f'{genome_build}_with_multimapping.tsv')
        if not os.path.exists(dbSNP_file_with_multimapping):
            dbSNP_chromosome_IDs = {
                'NC_000001.11': 'chr1', 'NC_000002.12': 'chr2',
                'NC_000003.12': 'chr3', 'NC_000004.12': 'chr4',
                'NC_000005.10': 'chr5', 'NC_000006.12': 'chr6',
                'NC_000007.14': 'chr7', 'NC_000008.11': 'chr8',
                'NC_000009.12': 'chr9', 'NC_000010.11': 'chr10',
                'NC_000011.10': 'chr11', 'NC_000012.12': 'chr12',
                'NC_000013.11': 'chr13', 'NC_000014.9': 'chr14',
                'NC_000015.10': 'chr15', 'NC_000016.10': 'chr16',
                'NC_000017.11': 'chr17', 'NC_000018.10': 'chr18',
                'NC_000019.10': 'chr19', 'NC_000020.11': 'chr20',
                'NC_000021.9': 'chr21', 'NC_000022.11': 'chr22',
                'NC_000023.11': 'chrX', 'NC_000024.10': 'chrY',
                'NC_012920.1': 'chrM'
            } if genome_build == 'hg38' else {
                'NC_000001.10': 'chr1', 'NC_000002.11': 'chr2',
                'NC_000003.11': 'chr3', 'NC_000004.11': 'chr4',
                'NC_000005.9': 'chr5', 'NC_000006.11': 'chr6',
                'NC_000007.13': 'chr7', 'NC_000008.10': 'chr8',
                'NC_000009.11': 'chr9', 'NC_000010.10': 'chr10',
                'NC_000011.9': 'chr11', 'NC_000012.11': 'chr12',
                'NC_000013.10': 'chr13', 'NC_000014.8': 'chr14',
                'NC_000015.9': 'chr15', 'NC_000016.9': 'chr16',
                'NC_000017.10': 'chr17', 'NC_000018.9': 'chr18',
                'NC_000019.9': 'chr19', 'NC_000020.10': 'chr20',
                'NC_000021.8': 'chr21', 'NC_000022.10': 'chr22',
                'NC_000023.10': 'chrX', 'NC_000024.9': 'chrY',
                'NC_012920.1': 'chrM'}
            dbSNP_awk_array = ''.join(
                f'dbSNP_map["{ID}"] = "{chrom}"; '
                for ID, chrom in dbSNP_chromosome_IDs.items())
            # Memory optimization: seen is deleted after every chromosome
            # (assumes chromosomes are contiguous in dbSNP, which they are)
            run(f'awk \'BEGIN {{while ((getline line) && (line ~ /^##/)); '
                f'print "CHROM", "BP", "REF", "ALT", "RSID"; {dbSNP_awk_array}'
                f'}} $1 in dbSNP_map {{split($5, alts, ","); chrom = '
                f'dbSNP_map[$1]; if (chrom != prev_chrom) delete seen; '
                f'prev_chrom = chrom; for (i in alts) {{ bp = $2; ref = $4; '
                f'alt = alts[i]; {get_minimal_representations_awk()}; '
                f'if (!seen[bp":"ref":"alt":"$3]++) print chrom, bp, ref, '
                f'alt, $3}}}}\' OFS="\t" <(zcat {full_dbSNP_file}) > '
                f'{dbSNP_file_with_multimapping}')
        # Remove multi-mapping variants; this is equivalent (aside from
        # sorting) to dbSNP = dbSNP.unique(['RSID', 'REF', 'ALT'], keep='none')
        dbSNP = pl.read_csv(dbSNP_file_with_multimapping, **read_csv_kwargs)
        dbSNP = dbSNP\
            .cast({'CHROM': pl.Enum([f'chr{i}' for i in range(1, 23)] +
                                    ['chrX', 'chrY', 'chrM'])})\
            .sort('RSID', 'REF', 'ALT')\
            .with_columns(rle_id=pl.struct('RSID', 'REF', 'ALT').rle_id())\
            .filter(pl.col.rle_id.ne(pl.col.rle_id.shift(fill_value=-1)) &
                    pl.col.rle_id.ne(pl.col.rle_id.shift(-1, fill_value=-1)))\
            .drop('rle_id')\
            .sort('CHROM', 'BP', 'REF', 'ALT', 'RSID')
        dbSNP.write_csv(dbSNP_file, separator='\t')
        return dbSNP


def check_valid_dbSNP(dbSNP):
    """
    Checks if a dbSNP instance is valid.
    
    Args:
        dbSNP: a dbSNP instance
    """
    if not isinstance(dbSNP, pl.DataFrame):
        raise TypeError(f'dbSNP must be a DataFrame returned by '
                        f'load_dbSNP(), but has type {type(dbSNP).__name__}')
    if dbSNP.columns != ['CHROM', 'BP', 'REF', 'ALT', 'RSID']:
        raise ValueError(f"dbSNP must be a DataFrame returned by load_dbSNP() "
                         f"with columns ['CHROM', 'BP', 'REF', 'ALT', 'RSID']")
    if dbSNP.height < 1_000_000_000:
        raise ValueError(f'dbSNP must be a DataFrame returned by '
                         f'load_dbSNP() and should have at least a billion '
                         f'rows, but yours has only {dbSNP.height} rows')


def get_rs_numbers(df, dbSNP, *, chrom_col='CHROM', bp_col='BP', ref_col='REF',
                   alt_col='ALT', rs_col='SNP', flip_col='FLIP',
                   fall_back_to_old_IDs=False, verbose=True):
    """
    Given a DataFrame of variants with chrom_col, bp_col, ref_col, and alt_col
    columns, adds a column rs_col to the DataFrame with the rs numbers (joined
    with commas, in the rare case a variant has multiple). Also adds a column
    flip_col, saying which variants needed to have their ref and alt alleles
    flipped to match dbSNP; flipping is only attempted for single-nucleotide
    variants.
    
    If rs_col is already a column of df, variants present in dbSNP will have
    their IDs in rs_col overwritten with their rs numbers. rs numbers for
    variants not present in dbSNP (including multi-mapping variants; see
    load_dbSNP()) will be set to null, unless fall_back_to_old_IDs=True, in
    which case their original IDs will be retained.
    
    If rs_col is not already a column of df, variants missing from dbSNP will
    always have their rs numbers set to null, as there are no old IDs to fall
    back to.
    
    Note: this function does not include defunct rs numbers that have been
    merged into other rs numbers, since they don't appear in the dbSNP files
    used in load_dbSNP().

    Args:
        df: a DataFrame with chrom_col, bp_col, ref_col, and alt_col columns
        dbSNP: a DataFrame returned by load_dbSNP(); must match the genome
               build of df's bp_col!
        chrom_col: the name of the chromosome column in df
        bp_col: the name of the base-pair column in df
        ref_col: the name of the reference allele column in df
        alt_col: the name of the alternate allele column in df
        rs_col: the name of the rs number column to be added to df
        flip_col: the name of the flip column to be added to df. True where
                  alleles had to be flipped to match dbSNP, False where they
                  matched without flipping, null if the variant didn't match
                  dbSNP either with or without flipping
        fall_back_to_old_IDs: if True, variants missing from dbSNP will retain
                              their original variant IDs, instead of having
                              them set to null. Requires chrom_col and bp_col
                              to already be present in df.
        verbose: whether to print what's happening at each step

    Returns: df with two additional columns: rs_col, containing the rs numbers,
             and flip_col, containing which variants were flipped.
    """
    if df.is_empty():
        raise ValueError(f'df is empty!')
    if chrom_col not in df:
        raise ValueError(f'"{chrom_col}" not in df; specify chrom_col')
    if bp_col not in df:
        raise ValueError(f'"{bp_col}" not in df; specify bp_col')
    if ref_col not in df:
        raise ValueError(f'"{ref_col}" not in df; specify ref_col')
    if alt_col not in df:
        raise ValueError(f'"{alt_col}" not in df; specify alt_col')
    if flip_col in df:
        raise ValueError(f'"{flip_col}" already in df; rename it or specify '
                         f'a different column name for flip_col')
    if fall_back_to_old_IDs and rs_col not in df:
        raise ValueError(f'You specified fall_back_to_old_IDs=True, but '
                         f'rs_col "{rs_col}" is not in df; specify it')
    check_valid_dbSNP(dbSNP)
    # Construct the rs number column piecewise by chromosome for efficiency
    df = df.with_row_index()
    rs_numbers = None
    for df_chrom_ID in df[chrom_col].unique(maintain_order=True):
        try:
            chrom = standardize_chromosomes(df_chrom_ID)
        except ValueError:
            raise ValueError(f'df contains non-standard chromosome '
                             f'"{df_chrom_ID}"!')
        if verbose:
            print(f'Getting rs numbers for {chrom}...')
        # Subset to chromosome; convert df to minimal representations
        dbSNP_chrom = dbSNP.filter(pl.col.CHROM == chrom)\
            .drop('CHROM')\
            .rename({'BP': bp_col, 'REF': ref_col, 'ALT': alt_col,
                     'RSID': rs_col})
        df_chrom = df.filter(pl.col(chrom_col) == df_chrom_ID)\
            .select('index', bp_col, ref_col, alt_col)\
            .pipe(get_minimal_representations, bp_col=bp_col, ref_col=ref_col,
                  alt_col=alt_col)\
            .with_columns(pl.col(bp_col).cast(pl.Int32),
                          pl.col(ref_col, alt_col).cast(pl.Categorical))
        # Allow matches without ref/alt flips...
        matches_without_flips = df_chrom\
            .join(dbSNP_chrom, on=[bp_col, ref_col, alt_col], how='left',
                  coalesce=True)\
            .drop_nulls(rs_col)
        # ...or with ref/alt flips, but only for SNVs (since for indels, e.g.
        # "21:15847757:A:AG" and "21:15847757:AG:A" are different variants: the
        # first is an insertion and the second is a deletion)
        # Remove .cast(pl.String) once polars allows Categorical
        # .str.len_bytes(): github.com/pola-rs/polars/issues/9773
        matches_with_flips = df_chrom\
            .filter(pl.col(ref_col).cast(pl.String).str.len_bytes() == 1,
                    pl.col(alt_col).cast(pl.String).str.len_bytes() == 1)\
            .join(dbSNP_chrom, left_on=[bp_col, ref_col, alt_col],
                  right_on=[bp_col, alt_col, ref_col], how='left')\
            .drop_nulls(rs_col)
        # Ensure no variants match both with and without flips (theoretically
        # possible since dbSNP has lots of edge cases, but we don't support it)
        assert not matches_with_flips['index'].is_in(
            matches_without_flips['index']).any()
        # Merge matches with and without flips; in the rare case that a variant
        # has multiple rs numbers, report all of them as a comma-separated list
        chrom_rs_numbers = pl.concat([
            matches_without_flips.with_columns(pl.lit(False).alias(flip_col)),
            matches_with_flips.with_columns(pl.lit(True).alias(flip_col))])\
            .group_by('index')\
            .agg(pl.col(rs_col).sort_by(pl.col(rs_col).str.slice(2).cast(int))
                 .str.join(','), pl.first(flip_col))
        rs_numbers = chrom_rs_numbers if rs_numbers is None else \
            rs_numbers.extend(chrom_rs_numbers)
        del dbSNP_chrom, df_chrom, matches_without_flips, matches_with_flips, \
            chrom_rs_numbers
    if rs_col in df:
        df = df\
            .lazy()\
            .drop(rs_col)\
            .join(rs_numbers.lazy(), on='index', how='left')\
            .drop('index')\
            .collect()
    else:
        df = df\
            .join(rs_numbers, on='index', how='left')\
            .drop('index')
    # Print how many variants had matching rs numbers
    if verbose:
        num_mapped = len(df) - df[rs_col].null_count()
        print(f'{num_mapped:,} of {len(df):,} variants '
              f'({100 * num_mapped / len(df):.2f}%) had matching rs numbers')
    # Fall back to old IDs, if specified
    if rs_col in df and fall_back_to_old_IDs:
        df = df.with_columns(pl.col(rs_col).fill_null(df[rs_col]))
    return df


def get_positions(df, dbSNP, *, rs_col='SNP', ref_col='REF', alt_col='ALT',
                  chrom_col='CHROM', bp_col='BP', flip_col='FLIP',
                  fall_back_to_old_positions=False):
    """
    The reverse of get_rs_numbers(): given a DataFrame of variants with rs_col,
    ref_col, and alt_col columns, adds columns chrom_col and bp_col to the
    DataFrame giving the chromosome and base-pair positions of each variant.
    Also adds a column flip_col, saying which variants needed to have their ref
    and alt alleles flipped to match dbSNP; flipping is only attempted for
    single-nucleotide variants.
    
    Use this function to remap sumstats to a different genome build!
    
    If chrom_col and bp_col are already columns of df, variants present in
    dbSNP will have their chromosomes and base-pairs in chrom_col and bp_col
    overwritten with those from dbSNP. Variants not present in dbSNP (including
    multi-mapping variants; see load_dbSNP()) will have their chromosomes and
    base-pair positions set to null, unless fall_back_to_old_positions=True, in
    which case their original chromsomes/base-pair positions wil be retained.
    
    If chrom_col and bp_col are not already columns of df, variants missing
    from dbSNP will always have their chromosomes and base-pair positions set
    to null, as there are no old chromosomes and base-pair positions to fall
    back to.
    
    Note: this function does not include defunct rs numbers that have been
    merged into other rs numbers, since they don't appear in the dbSNP files
    used in load_dbSNP(). Unlike for get_rs_numbers(), which also has this
    behavior, here it is a limitation!
    
    Args:
        df: a DataFrame with rs_col, ref_col, and alt_col columns
        dbSNP: a DataFrame returned by load_dbSNP() for the genome build you
               want to get chrom/bp positions for
        rs_col: the name of the rs number column in df
        ref_col: the name of the reference allele column in df
        alt_col: the name of the alternate allele column in df
        chrom_col: the name of the chromosome column to be added to df
        bp_col: the name of the base-pair column to be added to df
        flip_col: the name of the flip column to be added to df. True where
                  alleles had to be flipped to match dbSNP, False where they
                  matched without flipping, null if the variant didn't match
                  dbSNP either with or without flipping
        fall_back_to_old_positions: if True, variants missing from dbSNP will
                                    retain their original chromosomes and
                                    base-pair positions, instead of having them
                                    set to null. Requires chrom_col and bp_col
                                    to already be present in df.

    Returns: df with three additional columns: chrom_col and bp_col, containing
             the chromosomes and base-pair positions, and flip_col, containing
             which variants were flipped.
    """
    if df.is_empty():
        raise ValueError(f'df is empty!')
    if rs_col not in df:
        raise ValueError(f'"{rs_col}" not in df; specify rs_col')
    if ref_col not in df:
        raise ValueError(f'"{ref_col}" not in df; specify ref_col')
    if alt_col not in df:
        raise ValueError(f'"{alt_col}" not in df; specify alt_col')
    if chrom_col in df and bp_col not in df:
        raise ValueError(f'chrom_col "{chrom_col}" is present in df but '
                         f'bp_col "{bp_col}" is not; either both must be '
                         f'present, or neither')
    if chrom_col not in df and bp_col in df:
        raise ValueError(f'bp_col "{bp_col}" is present in df but chrom_col '
                         f'"{chrom_col}" is not; either both must be present, '
                         f'or neither')
    if flip_col in df:
        raise ValueError(f'"{flip_col}" already in df; rename it or specify '
                         f'a different column name for flip_col')
    check_valid_dbSNP(dbSNP)
    # If fall_back_to_old_positions=True, save chromosomes and base-pair
    # positions for later; otherwise, drop them so that they don't mess up the
    # join with dbSNP
    if fall_back_to_old_positions:
        if chrom_col not in df:
            raise ValueError(f'You specified fall_back_to_old_positions=True, '
                             f'but chrom_col "{chrom_col}" and bp_col '
                             f'"{bp_col}" are not in df; specify them')
        df = df.rename({chrom_col: f'_GET_VARIANT_POSITIONS_{chrom_col}',
                        bp_col: f'_GET_VARIANT_POSITIONS_{bp_col}'})
    else:
        df = df.drop(chrom_col, bp_col)
    # Convert df's ref and alt to their minimal representations
    df = df\
        .pipe(get_minimal_representations, ref_col=ref_col, alt_col=alt_col,
              bp_col=None)\
        .with_columns(pl.col(ref_col, alt_col).cast(pl.Categorical))
    # Unlike for get_rs_numbers(), there isn't much of an efficiency gain in
    # constructing the chromosome and base-pair columns piecewise by chromosome
    # because the variant's chromosome isn't known a priori, so it needs to be
    # matched against the entirety of dbSNP. However, as an optimization,
    # subset dbSNP to just the rsIDs in df.
    dbSNP = dbSNP\
        .rename({'CHROM': chrom_col, 'BP': bp_col, 'REF': ref_col,
                 'ALT': alt_col, 'RSID': rs_col})\
        .filter(pl.col(rs_col).is_in(df[rs_col]))
    # Allow matches without ref/alt flips, or with ref/alt flips for SNVs only
    # (since for indels, e.g. "21:15847757:A:AG" and "21:15847757:AG:A" are
    # different variants: the first is an insertion, the second is a deletion)
    # Remove .cast(pl.String) once polars allows Categorical .str.len_bytes():
    # github.com/pola-rs/polars/issues/9773
    df = df\
        .lazy()\
        .join(dbSNP.lazy(), on=[rs_col, ref_col, alt_col], how='left')\
        .join(dbSNP.lazy(), left_on=[rs_col, ref_col, alt_col],
              right_on=[rs_col, alt_col, ref_col], how='left',
              suffix='_flipped')\
        .with_columns(pl.when(pl.col(ref_col).cast(pl.String)
                              .str.len_bytes() == 1,
                              pl.col(alt_col).cast(pl.String)
                              .str.len_bytes() == 1)
                      .then(pl.col(f'{chrom_col}_flipped',
                                   f'{bp_col}_flipped')))\
        .with_columns(pl.when(pl.col(chrom_col).is_not_null())
                      .then(False)
                      .when(pl.col(f'{chrom_col}_flipped').is_not_null())
                      .then(True)
                      .alias(flip_col))\
        .collect()
    # Ensure no variants match both with and without flips (theoretically
    # possible since dbSNP has lots of edge cases, but we don't support it)
    assert len(df.filter(pl.col(chrom_col).is_not_null(),
                         pl.col(f'{chrom_col}_flipped').is_not_null())) == 0
    # Merge matches with and without flips
    df = df\
        .with_columns(
            pl.col(chrom_col).fill_null(pl.col(f'{chrom_col}_flipped')),
            pl.col(bp_col).fill_null(pl.col(f'{bp_col}_flipped')))\
        .drop(f'{chrom_col}_flipped', f'{bp_col}_flipped')
    # If fall_back_to_old_positions=True, retain the original chromosomes
    # and base-pair positions for variants not found in dbSNP
    if fall_back_to_old_positions:
        df = df\
            .with_columns(pl.col(chrom_col).fill_null(
                              pl.col(f'_GET_VARIANT_POSITIONS_{chrom_col}')),
                          pl.col(bp_col).fill_null(
                              pl.col(f'_GET_VARIANT_POSITIONS_{bp_col}')))\
            .drop(f'_GET_VARIANT_POSITIONS_{chrom_col}',
                  f'_GET_VARIANT_POSITIONS_{bp_col}')
    return df


def get_rs_numbers_bim_or_pvar(bim_or_pvar_file, dbSNP, *, verbose=True):
    """
    Given a plink bim/pvar file, replace the IDs in the second column with
    rs numbers, saving the original file to {bim_or_pvar_file}.no_rs_numbers.
    Variants not in dbSNP (or multi-mapping variants; see load_dbSNP()) will
    retain their original variant IDs.
    
    The file must end with .bim or .pvar; file format will be inferred based on
    the extension. If .pvar, strips the header lines starting with ##.
    
    Args:
        bim_or_pvar_file: the bim or pvar file to read from and write to; must
                          end with .bim or .pvar, and the extension determines
                          which file type it's parsed and written as
        dbSNP: a DataFrame returned by load_dbSNP(); must match the genome
               build of bim_or_pvar_file!
        verbose: whether to print what's happening at each step
    """
    check_valid_dbSNP(dbSNP)
    is_pvar = get_bim_or_pvar_file_type(bim_or_pvar_file)
    file_type = 'pvar' if is_pvar else 'bim'
    old_bim_or_pvar_file = f'{bim_or_pvar_file}.no_rs_numbers'
    if verbose and os.path.exists(old_bim_or_pvar_file):
        print(f'old {file_type} file "{old_bim_or_pvar_file}" already exists; '
              f'returning without updating rs numbers, since rs numbers seem '
              f'to have already been added')
        return
    if verbose:
        print(f'Loading {file_type} file "{bim_or_pvar_file}"...')
    variants = read_bim_or_pvar(bim_or_pvar_file)
    if verbose:
        print(f'Getting rs numbers...')
    variants = get_rs_numbers(variants, dbSNP=dbSNP, verbose=verbose,
                              fall_back_to_old_IDs=True)
    if verbose:
        print(f'Moving "{bim_or_pvar_file}" to "{old_bim_or_pvar_file}"...')
    run(f'mv "{bim_or_pvar_file}" "{old_bim_or_pvar_file}"')
    if verbose:
        print(f'Saving to "{bim_or_pvar_file}"...')
    write_bim_or_pvar(variants, bim_or_pvar_file)


def get_positions_bim_or_pvar(bim_or_pvar_file, new_bim_or_pvar_file, dbSNP, *,
                              separator='~', verbose=True):
    """
    Given a plink bim/pvar file, replace its base-pair positions with those
    from dbSNP, remove variants not found in dbSNP (or multi-mapping variants;
    see load_dbSNP()), and save the result to a new bim/pvar file,
    new_bim_or_pvar_file. Then, create bed/fam or pgen/psam files with the same
    prefix as new_bim_or_pvar_file, with the same variants removed.
    
    Use this function to remap plink files to a different genome build!
    
    The file must end with .bim or .pvar; file format will be inferred based on
    the extension. If .pvar, strips the header lines starting with ##.
    
    Args:
        bim_or_pvar_file: the bim or pvar file to read from; must end with .bim
                          or .pvar, and the extension determines which file
                          type it's parsed as
        new_bim_or_pvar_file: the bim or pvar file to write to; must end with
                              .bim or .pvar, and the extension determines which
                              file type it's written as. The prefix of this
                              file will be used for the created bed/fam and
                              pgen/psam files.
        dbSNP: a DataFrame returned by load_dbSNP() for the genome build you
               want to get chrom/bp positions for
        separator: the string to place between the variant ID, reference
                   allele, and alternate allele when creating temporary
                   bim/pvar files (an internal step in this function); no
                   variant IDs may contain the separator
        verbose: whether to print what's happening at each step
    """
    check_valid_dbSNP(dbSNP)
    is_pvar = get_bim_or_pvar_file_type(bim_or_pvar_file)
    file_type = 'pvar' if is_pvar else 'bim'
    if os.path.exists(new_bim_or_pvar_file):
        raise ValueError(f'new {file_type} file "{new_bim_or_pvar_file}" '
                         f'already exists! Delete it before running '
                         f'get_positions_bim_or_pvar()')
    if verbose:
        print(f'Loading {file_type} file "{bim_or_pvar_file}"...')
    variants = read_bim_or_pvar(bim_or_pvar_file)
    if variants['SNP'].str.contains(separator).any():
        raise ValueError(f'Some variant IDs in the {file_type} file '
                         f'"{bim_or_pvar_file}" contain the separator '
                         f'{separator!r}; specify a different separator '
                         f'via the separator argument!')
    if verbose:
        print('Getting new positions...')
    # Modify this once polars allows Categorical .str.slice():
    # github.com/pola-rs/polars/issues/9773
    new_variants = get_positions(variants, dbSNP=dbSNP)\
        .drop_nulls(['CHROM', 'BP'])\
        .with_columns(pl.col.CHROM.cast(str).str.slice(3))\
        .drop('FLIP')
    # Make variant IDs unique, like make_bim_or_pvar_IDs_unique() but without
    # including chromosome, because variants can (very rarely) switch
    # chromosomes between builds. Also make a file of uniquified variant IDs
    # that were found in dbSNP.
    bim_or_pvar_file_tmp = \
        bim_or_pvar_file.replace(f'.{file_type}', f'.tmp.{file_type}')
    found_variants_file = \
        bim_or_pvar_file.replace(f'.{file_type}', '.found')
    if verbose:
        print(f'Saving a {file_type} with the original positions (with '
              f'uniquified variant IDs) to "{bim_or_pvar_file_tmp}"...')
    variants\
        .with_columns(pl.concat_str('SNP', 'REF', 'ALT', separator=separator))\
        .pipe(write_bim_or_pvar, bim_or_pvar_file_tmp)
    if verbose:
        print(f'Saving a list of variants found in dbSNP (with uniquified '
              f'variant IDs) to "{found_variants_file}"...')
    new_variants\
        .select(pl.concat_str('SNP', 'REF', 'ALT', separator=separator))\
        .write_csv(found_variants_file)
    prefix = bim_or_pvar_file.removesuffix(f'.{file_type}')
    new_prefix = new_bim_or_pvar_file.removesuffix(f'.{file_type}')
    if verbose:
        print(f'Subsetting to variants in dbSNP, saving to '
              f'"{new_prefix}.tmp.{{pgen,pvar,psam}}"...')
    run(f'plink2 '
        f'{"" if verbose else "--silent"} '
        f'--{"pfile" if is_pvar else "bfile"} "{prefix}" '
        f'--pvar "{bim_or_pvar_file_tmp}" '
        f'--extract "{found_variants_file}" '
        f'--make-pgen '
        f'--out "{new_prefix}.tmp"')
    if verbose:
        print(f'Updating "{new_prefix}.tmp.pvar" to the new positions...')
    temp_pvar_file = f'{new_prefix}.tmp.pvar'
    new_variant_order = read_bim_or_pvar(temp_pvar_file, categorical=True)\
        .with_columns(SNP=pl.col.SNP.str.split_exact(separator, 1)
                      .struct.field('field_0'))
    assert len(new_variants) == len(new_variant_order)
    new_variants = new_variant_order\
        .select('SNP', 'REF', 'ALT')\
        .join(new_variants, on=('SNP', 'REF', 'ALT'), how='left')
    write_bim_or_pvar(new_variants, temp_pvar_file)
    if verbose:
        print(f'Sorting by the new positions, saving to "{new_prefix}.'
              f'{{{"pgen,pvar,psam" if is_pvar else "bed,bim,fam"}}}"...')
    run(f'plink2 '
        f'{"" if verbose else "--silent"} '
        f'--pfile "{new_prefix}.tmp" '
        f'--sort-vars '
        f'--make-pgen '
        f'--out "{new_prefix}{"" if is_pvar else ".tmp2"}"')
    if not is_pvar:
        # Work around "Error: Fixed-width .bed/.pgen output doesn't support
        # sorting yet. Generate a regular sorted .pgen first, and then reformat
        # it.": groups.google.com/g/plink2-users/c/Vt-g0dblpw4
        run(f'plink2 '
            f'{"" if verbose else "--silent"} '
            f'--pfile "{new_prefix}.tmp2" '
            f'--make-bed '
            f'--out "{new_prefix}"')
    if verbose:
        print('Cleaning up...')
    run(f'rm "{bim_or_pvar_file_tmp}" "{found_variants_file}" '
        f'"{new_prefix}.log" "{new_prefix}.tmp.pgen" "{new_prefix}.tmp.pvar" '
        f'"{new_prefix}.tmp.psam" "{new_prefix}.tmp.log"')
    if not is_pvar:
        run(f'rm "{new_prefix}.tmp2.pgen" "{new_prefix}.tmp2.pvar" '
            f'"{new_prefix}.tmp2.psam" "{new_prefix}.tmp2.log"')
