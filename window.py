def slide_window()
    '''
    This function procceses the quality scores from the sequence
    Input: 
        scores:list with quality scores 
        window_size: size of window to calculate avg
        threshold: threshold to trimm sequence
    Output:
        5'end cut index
        3'end cut index
    '''
    
    
    
def calc_avg(scores)
    '''
    Calculate avg given list of scares
    Input:
        scores: list containing scores
    Output:
        avg: avg of scores
        
    >>> calc_avg([1])
    1.0
    >>> calc_avg ([1,1,1,1])
    1.0
    '''
    return sum(scores)/float(len(scores))
