# https://groklearning.com/learn/unimelb-comp10001-2021-s2/proj1/0/

# PROJECT 1


def changes(ts):
    ''' 
    Takes in ts, a list of numbers represeting a time series
    Caclulates the changes between sequential number values in ts 
    Returns as a generator to reduce memory use
    '''
    
    # change =  final - initial
    return (ts[i + 1] - ts[i] for i in range(len(ts) - 1))


def non_increasing_tst(ts, delta):
    ''' 
    Takes in 
        - ts, a list of numbers represeting a time series
        - delta, the magnitude of change between values condidered noise
    calcaultes if there is not an increase greater than delta in ts
    Returns a Boolean
    '''
    
    return not any(c > delta for c in changes(ts))


####################

def changes(ts):
    ''' 
    Takes in ts, a list of numbers represeting a time series
    Caclulates the changes between sequential number values in ts 
    Returns as a generator to reduce memory use
    '''
    
    # change =  final - initial
    return (ts[i + 1] - ts[i] for i in range(len(ts) - 1))


def change_to_char(change, delta):
    ''' 
    Takes in 
        - change, the change between numbers
        - delta, the magnitude of change between values condidered noise
    calculates a character corresponding to the change being up, steady or down 
    Returns a single string character
    '''
    
    if abs(change) <= delta:
        return 's'
    if change > 0:
        return 'u'
    return 'd'
        

def contains_tst(ts, pattern, delta):
    '''
    Takes in
        - ts, a list of numbers representing a time series
        - pattern, a string representing a possible pattern in the time series
    calculates if the pattern is in the time series
    Returns a boolean
    '''
    
    # create a ts pattern string and check it contains 'pattern' arguement
    ts_pattern = ''.join(change_to_char(c, delta) for c in changes(ts))
    return pattern in ts_pattern



####################


from math import sqrt

def ts_distance(ts1, ts2):
    '''
    Takes in two time series
    calculates the distance between two time series as defined below
    https://groklearning.com/learn/unimelb-comp10001-2021-s2/proj1/7/
    Returns a float
    '''
    
    # name variables to improve readability
    n = len(ts1)
    
    # calculate and store variables which are used more than once
    ts_diffs = [abs(ts1[i] - ts2[i]) for i in range(n)]
    mean = sum(ts_diffs) / n
    
    return sqrt((1 / n) * sum((dif - mean) ** 2 for dif in ts_diffs))


def closest_tst(ts, ts_list):
    ''' 
    Takes in
        - ts, a list of numbers representing a time series
        - ts_list, a list of time series
    calculates the index of the time series in ts_list closest to ts 
    Returns an integer
    '''
    
    # calculate and store distances
    # return index of closest, which will have the min distance
    
    dists = [ts_distance(ts, cur_ts) for cur_ts in ts_list]
    return dists.index(min(dists))


####################



from math import sqrt

def mean_distance(i1, seqs):
    '''
    Takes in 
        - i1, an index 
        - seqs, a list of sequences (lists containing numbers)
    calculates the average vector distance of the sequence at i1 to others
    of the same length
    Returns a float
    '''
    
    # calculate the sum of distances
    dist_sum = 0
    s1 = seqs[i1]
    for i2 in range(len(seqs)):
        if i1 != i2:
            s2 = seqs[i2]
            dist_sum += sqrt(sum((s1[i] - s2[i]) ** 2 for i in range(len(s1))))

    # return the mean
    return dist_sum / (len(seqs) - 1) 


def anomaly_tst(ts, w, threshold):
    '''
    Takes in 
        - ts, a list of floats representing a time series
        - w, an integer 
        - threshold, a float
    calcaultes the index of the largest anomaly if any are in ts
    where an anomaly is when
    the mean distance of one subsequence to others
    is greater than the threshold value
    otherwise -1 is returned
    Returns an integer
    '''

    # calculate and store mean_distances
    subseqs = [ts[i:i + w] for i in range(len(ts) - w + 1)]
    mean_distances = [mean_distance(i1, subseqs) for i1 in range(len(subseqs))]
    
    # return mean_distances index if an anomaly occurs or -1 for no index
    max_mean = max(mean_distances) 
    if max_mean > threshold:
        return mean_distances.index(max_mean)
    return -1





















# PROJECT 2






def manhatan(loctn_a, loctn_b):
    ''' return manhatan distance between location a and b
    
    Args: 
        loctn_a / loctn_b (tuple): locations with integer axis values
        
    Returns: integer
    '''
    return sum(abs(a_val - b_val) for a_val, b_val in zip(loctn_a, loctn_b))


def execute1_lds(plan, start):
    ''' return distance to complete deliveries, 
        without considering capacity and using the route in the order of plan
    
    Args:
        plan (list): of deliveries (dlv)
            dlv[0] (integer): number of boxes being delivered
            dlv[1] (tuple): box pickup location 
            dlv[2] (tuple): box delivery location 
        start (tuple): truck start location
    
    Returns: integer
    '''
    
    dist = 0
    for dlv in plan:
        
        # add distance between deliveries
        dist += manhatan(start, dlv[1])  
        
        # add distance for deliveries
        dist += manhatan(dlv[1], dlv[2]) 
        
        # update start location
        start = dlv[2]  
    
    return dist




############



from math import ceil

def manhatan(loctn_a, loctn_b):
    ''' return manhatan distance between location a and b
    
    Args: 
        loctn_a / loctn_b (tuple): locations with integer axis values
        
    Returns: integer
    '''
    return sum(abs(a_val - b_val) for a_val, b_val in zip(loctn_a, loctn_b))


def execute2_lds(plan, start, capacity):
    ''' return distance to complete deliveries using route in the order of plan
    
    Args:
        plan (list): of deliveries (dlv)
            dlv[0] (integer): number of boxes being delivered
            dlv[1] (tuple): box pickup location 
            dlv[2] (tuple): box delivery location 
        start (tuple): truck start location
        capacity (integer): number of boxes truck can hold
    
    Returns: integer
    '''

    dist = 0
    for dlv in plan:
        
        # add distance between deliveries
        dist += manhatan(start, dlv[1])  
        
        # add distance for deliveries
        dist += manhatan(dlv[1], dlv[2]) * (ceil(dlv[0] / capacity) * 2 - 1) 
        
        # update start location
        start = dlv[2]  
    
    return dist





##########




from math import ceil

def manhatan(loctn_a, loctn_b):
    ''' return manhatan distance between location a and b
    
    Args: 
        loctn_a / loctn_b (tuple): locations with integer axis values
        
    Returns: integer
    '''
    return sum(abs(a_val - b_val) for a_val, b_val in zip(loctn_a, loctn_b))


def optimise1_lds(plan, start, capacity):
    ''' return distance to complete deliveries using greedy route 
    
    Args:
        plan (list): of deliveries (dlv)
            dlv[0] (integer): number of boxes being delivered
            dlv[1] (tuple): box pickup location 
            dlv[2] (tuple): box delivery location 
        start (tuple): truck start location
        capacity (integer): number of boxes truck can hold
    
    Returns: integer
    '''
    
    dist = 0
    
    # add distance for deliveries
    for dlv in plan:
        dist += manhatan(dlv[1], dlv[2]) * (ceil(dlv[0] / capacity) * 2 - 1)
    
    # add distance between deliveries
    while plan:
        
        # find min index
        dists = [manhatan(start, dlv[1]) for dlv in plan]
        min_indx = dists.index(min(dists))
   
        # update values
        dist += dists[min_indx]
        start = plan[min_indx][2]
        del plan[min_indx]

    return dist







###########






from math import ceil

def manhatan(loctn_a: tuple, loctn_b: tuple) -> int:
    ''' return manhatan distance between location a and b
    
    Args: 
        loctn_a / loctn_b (tuple): locations with integer axis values
        
    Returns: integer
    '''
    return sum(abs(a_val - b_val) for a_val, b_val in zip(loctn_a, loctn_b))


def min_dis_between_delveries(plan, start):
    ''' return minimum distance between deliveries
    
    Args: 
        plan (list): of deliveries (dlv)
            dlv[0] (integer): number of boxes being delivered
            dlv[1] (tuple): box pickup location 
            dlv[2] (tuple): box delivery location 
        start (tuple): truck start location
        
    Returns: integer
    '''
    if not plan:  # base case
        return 0

    # recursion to find minimum distance 
    return min(
        manhatan(start, dlv[1])
        + min_dis_between_delveries(plan[:i] + plan[i + 1:], dlv[2]) 
        for i, dlv in enumerate(plan)
    )


def optimise2_lds(plan, start, capacity):
    ''' return distance to complete deliveries using shortest route 
    
    Args:
        plan (list): of deliveries (dlv)
            dlv[0] (integer): number of boxes being delivered
            dlv[1] (tuple): box pickup location 
            dlv[2] (tuple): box delivery location 
        start (tuple): truck start location
        capacity (integer): number of boxes truck can hold
    
    Returns: integer
    '''
    
    # add distance between deliveries
    dist = min_dis_between_delveries(plan, start)
     
    # add distance for deliveries
    for dlv in plan:
        dist += manhatan(dlv[1], dlv[2]) * (ceil(dlv[0] / capacity) * 2 - 1)  

    return dist



