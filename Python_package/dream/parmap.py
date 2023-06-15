from concurrent.futures import ProcessPoolExecutor

def parmap(func, iterable, max_workers=8):
    """
    Apply a given function to each element of an iterable in parallel using a
    pool of processes.

    Parameters
    ----------
    func : callable
        The function to apply to each element of the iterable.
    iterable : iterable
        The iterable containing the data to be processed.
    max_workers : int, optional (default=8)
        The maximum number of processes to use.

    Returns
    -------
    list
        A list containing the result of applying the function to each element
        of the iterable in parallel.
    """
    results = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        for i in executor.map(func, iterable):
            results.append(i)
    return results
