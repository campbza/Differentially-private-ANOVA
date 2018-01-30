import numpy as np

def datagen(means_list, v, n):
    '''
    Generate normally distributed data with means from means_list
    and standard deviation v.  For each group, sample n points.

    We truncate the datapoints to be within [0,1].
    '''
    #ARCOMMENT: means_list can now be any size.
    data = []
    for mean in means_list:
        data.append(np.array(np.random.normal(mean,v,n)).tolist())

    k = len(data)
    i = 0
    while i < k:
        l = len(data[i])
        j = 0
        while j < l:
            if data[i][j] > 1.0:
                data[i][j] = 1.0
            elif data[i][j] < 0.0:
                data[i][j] = 0.0
            j += 1
        i += 1
    return data

def sample(data, samplesize):
#input: data (format: list of three lists), and a sampling size
#output: datalist populated with 3 groups, each of size "samplesize", sampled from each list, no replacement
    sample1 = np.array(np.random.choice(data[0], samplesize, replace = False)).tolist()
    sample2 = np.array(np.random.choice(data[1], samplesize, replace = False)).tolist()
    sample3 = np.array(np.random.choice(data[2], samplesize, replace = False)).tolist()
    sample_data = [sample1, sample2, sample3]
    return sample_data

