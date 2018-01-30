import csv

def normalize(xs):
#input: list of lists
#output: normalized list of lists
    n = len(xs)
    maximum = max([max(xs[i]) for i in range(n)])
    ys = xs
    for i in range(n):
        l = len(xs[i])
        for j in range(l):
            ys[i][j] = xs[i][j] / float(maximum)
    return ys

def f(x, c1, c2):
#input: csv file x, column c1 (categorical variable - string), column c2 (response variable - real number)
#output: data list of lists
    with open(x, mode = 'r') as infile:
	    data = [(row[c1],row[c2]) for row in csv.DictReader(infile)]
    data_dict = {}
    for k,v in data:
	    if k in data_dict:
		    data_dict[k].append(v)
	    else:
		    data_dict[k] = [v]
    for i in data_dict:
        data_dict[i] = filter(None, data_dict[i])
        for j in range(len(data_dict[i])):
            data_dict[i][j] = float(data_dict[i][j])
    data_list = normalize([data_dict[i] for i in data_dict])
    return data_list

