import doctest

def WindowSize(n):
    '''Returns the appropriate window size for sickle
    The size is 0.1 times the length of the read.
    If this length is less than 1, then the window is set
    to be equal to the length of the read.
    
    >>> WindowSize(1)
    1
    >>> WindowSize(10)
    1
    >>> WindowSize(100)
    10
    '''
    
    # check n is valid integer > 0
    if not isinstance(n, int) or n <= 0:
        raise Exception("The length of the read must be > 0")

    # calculate and return window size
    if n <= 10:
        return 1
    else:
        return int(n*0.1)
        
if __name__ == "__main__":
    doctest.testmod()

