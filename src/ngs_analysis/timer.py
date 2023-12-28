"""A timer that can be used in a "with" statement. The results
can be printed out or saved to a dictionary (repeat uses in the
same closure will be saved too).
"""
import time
import inspect
import uuid

_timings = {}
invocation_ids = {}

class Timer:
    def __init__(self, name=None, verbose=False):
        """A useful timer. See `example` below.

        :param name: used in timings dictionary, if None then not added to
            dictionary
        :param verbose: if True, print the name and time elapsed when finished; 
            if a string, print the string when starting, and time elapsed when 
            finished
        """
        self.name = name
        self.verbose = verbose

    def __enter__(self):
        self.start = time.time()
        self.closure = identify_this_function(1)
        if isinstance(self.verbose, str):
            print(self.verbose, end=' ', flush=True)
        return self

    def __exit__(self, type, value, traceback):
        end = time.time()
        elapsed_time = end - self.start
        if self.verbose is True:
            if self.name is None:
                print('Done in:  {elapsed_time:.3g} seconds')
                return
            else:
                print(f'{self.name}: {elapsed_time:.3g} seconds')
        elif isinstance(self.verbose, str):
            print(f'done in {elapsed_time:.3g} seconds')
        
        this = self.closure
        
        if this not in _timings:
            _timings[this] = {}

        if self.name in _timings[this]:
            if isinstance(_timings[this][self.name], list):
                _timings[this][self.name].append(elapsed_time)
            else:
                # convert to list
                _timings[this][self.name] = [_timings[this][self.name], elapsed_time]
        else:
            _timings[this][self.name] = elapsed_time


def get_all_timings():
    """Gets timings measured in the parent frame.
    """
    return _timings.get(identify_this_function(1), {})


def example():

    def f_0():
        with Timer('task1'):
            time.sleep(0.5)
        with Timer('task1'):
            time.sleep(0.2)
        print(get_all_timings())

    def f_1():
        with Timer('task2', verbose=True):
            time.sleep(1)
        print(get_all_timings())

    def f_2():
        for _ in range(5):
            with Timer('task3'):
                time.sleep(0.1)
        print(get_all_timings())

    def f_3():
        with Timer('task4', verbose='Running task4...'):
            time.sleep(3)
        print(get_all_timings())
        
    f_0()
    f_1()
    f_2()
    f_3()
    f_0()


def identify_this_function(go_up=0):
    """Get a unique ID specific to this invocation of the function.

    :param go_up: how many extra frames to go up in the stack
    """
    frame = inspect.currentframe()
    caller_frame = frame.f_back
    for _ in range(go_up):
        caller_frame = caller_frame.f_back

    frame_id = id(caller_frame)

    if frame_id not in invocation_ids:
        invocation_ids[frame_id] = str(uuid.uuid4())

    return invocation_ids[frame_id]
