def slide_window(scores,window_size,threshold)
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
    size = len (scores)
    #Averages of all sliding windows
    averages = []
    for x in range(size-window_size)
        averages.append(calc_avg(scores[x:x+window_size]))
    
    first = ''
    second = ''
    x = 0
    #If average started bad , slide until we get a good one
    if averages[0] < threshold:
        while averages[x] < threshold and x < len(averages):
            first  = x
            x+=1
        #Found good average or end of seq, keep sliding until gets bad again
        while averages[x] > threshold and x < len(averages):
            second = x
            x=+1
        #If the reads were never good return 0 0
        if second == '':
            return 0,0
        else:
            return first, second
                
    else:
        #If we started with a good average slide until we find a bad one
        first = x
        while averages[x] > threshold and x < len(averages):
            second = x
            x=+1
        
        return first, second
        
    
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
