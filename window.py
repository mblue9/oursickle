import doctest

def slide_window(scores,window_size,threshold):
    '''
    This function procceses the quality scores from the sequence
    Input: 
        scores:list with quality scores 
        window_size: size of window to calculate avg
        threshold: threshold to trimm sequence
    Output:
        5'end cut index
        3'end cut index
   
    #Case 1: read starts good ends bad
    >>> slide_window([30,30,30,30,30,5,5,5,5],4,15)
    (0, 6)
    
    #Case2: read starts bad ends good
    >>> slide_window([5,5,5,5,5,30,30,30,30,30],4,15)
    (3, 9)
     
     '''
    
    
    
  
    size = len(scores)
    #Averages of all sliding windows
    averages = []
    for x in range(size-window_size+1):
        averages.append(calc_avg(scores[x:x+window_size]))
    
    first = ''
    second = ''
    x = 0
    #If average started bad, slide until we get a good one
    if averages[0] < threshold:
        while x < len(averages) and averages[x] < threshold:
            x+=1
            first  = x           
        #Found good average or end of seq, keep sliding until gets bad again
        while x < len(averages) and averages[x] > threshold:
            second = x
            x+=1
            
        #If the read scores are all below threshold return 0 0
        if second == '':
            return -1,-1
        else:
            return first, second + window_size - 1
                
    else:
        #If we started with a good average, slide along read until we find a bad one
        first = x
        while x < len(averages) and averages[x] > threshold:
            second = x
            x+=1
            

        return first, second + window_size -1
        
    
def calc_avg(scores):
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
    
# only run this block if we're
# executing this script directly;
# don't run if importing!
if __name__ == "__main__":
    doctest.testmod()
    

