#####################
#Ultimate ANOVA file#
#####################

import numpy as np
from csv_to_datalist import *
from datagen import *
import sys

def many_means(data):
#input: list of lists
#output: mean of each list in the data
    l = len(data)
    xs = []
    i = 0
    while i < l:
        xs.append(np.mean(data[i]))
        i += 1
    return xs

def overall_mean(data):
#input: data list
#output: mean over all records in the data
    l = len(data)
    flat_list = [item for sublist in data for item in sublist] #list-of-lists to list
    overall_mean = np.mean(flat_list)
    return overall_mean
    

def SSA(data, epsilon, total_size):
#input: data list, epsilon value (optional), total_size
#output: epsilon-differentially private SSA
    l = len(data)
    sizes = [len(data[i]) for i in range(l)]
    means = many_means(data)
    mean = overall_mean(data)
    ssa = 0
    i = 0
    while i < l:
        ssa += sizes[i] * ((means[i] - mean) ** 2)
        i += 1
    if epsilon == None:
        return ssa
    else:
        ssa += np.random.laplace(0.0, (9.0+5.0/total_size)/epsilon) #9+5./total_size is our provable upper-bound on the sensitivity of SSA
        return ssa

def SSE(data, epsilon):
#input: data list, epsilon value
#output: epsilon-differentially private SSE
    l = len(data)
    means = many_means(data)
    sse = 0
    i = 0
    while i < l:
        m = len(data[i])
        x = 0
        j = 0
        while j < m:
            x += (data[i][j] - means[i]) ** 2
            j += 1
        sse += x
        i += 1
    if epsilon == None:
        return sse
    else:
        sse += np.random.laplace(0.0, 7.0/epsilon) #7 is our provable upper-bound on the sensitivity of SSE 
        return sse

def fstar(n, dfa, dfe, mse, epsilon, total_size):
    '''
    We generate a sampling distribution for our test statistic. Since we are adding noise to 
    SSE and SSA, we cannot consult a true F-distribution. Instead, we must generate a sampling 
    distribution that takes into account the noise we are adding and consult this to get our 
    p-value.
    '''
#input: sample size, degrees of freedom, mse, and epsilon, total_size
#output: n random variables drawn from chisquare with noise added
    if epsilon != None:
        numerator = (mse*np.random.chisquare(dfa, n) + np.random.laplace(0.0, (9.0 + 5.0/total_size)/(epsilon/2.0), n))/(dfa)
        denominator = (mse*np.random.chisquare(dfe, n) + np.random.laplace(0.0, 7.0/(epsilon/2.0), n))/(dfe)
    else:
        numerator = (mse*np.random.chisquare(dfa, n))/(dfa)
        denominator = (mse*np.random.chisquare(dfe, n))/(dfe)
    return numerator/denominator

def anova(data, epsilon, fout, variance):
#input: normalized data (list of lists), epsilon, file to write to, variance (None to use estimate or a number to use if real.)
#output: returns a boolean value: True if p-value <0.05, False otherwise.
# This output is used to count how many significant p-values we have in anova_test() function.
# writes output to file fout

    number_of_groups = len(data)
    total_size = sum([len(data[i]) for i in range(len(data))])
    dfa = number_of_groups - 1
    dfe = total_size - number_of_groups

    ## Calculate SSE,SSA,MSE,MSA, and F
    if epsilon == None:
        sse = SSE(data, None)
        ssa = SSA(data, None,total_size)
    else:
        sse = SSE(data, epsilon / 2)
        ssa = SSA(data, epsilon / 2, total_size)
    mse = sse / dfe
    msa = ssa / dfa
    f = msa / mse

    ## Generate null distribution.  Checks the variance input variable to see if we
    ## should use the real variance or MSE as the estimate of variance.
    ## we use 100K samples to generate this distribution
    if variance != None: ## use the value of variance passed as input
        fstarsim = fstar(100000, dfa, dfe, variance, epsilon, total_size)
    else: # use the estimate (MSE) as variance.
        fstarsim = fstar(100000, dfa, dfe, mse, epsilon, total_size)

    ## Determine the proportion of samples that are greater than the f-ratio.
    pval = np.mean(fstarsim > f)

    ## write a line to the outfile.
    fout.write(str(total_size) + ',' + str(sse) + ',' + str(ssa) + ',' + str(epsilon) + ',' + str(pval) + ',' + str(f) + '\n')
    
    if pval < 0.05: ## hard-coded p-value
        return True
    return False

def anova_test(num_runs, epsilon_vals, filename, means_list, stddev, group_counts,realvar):
    fout = open(filename, 'w') #open filename
    i = 0
    
    while i < len(epsilon_vals):
        print('EPSILON VAL',epsilon_vals[i])
        j = 0
        allsig = False
        while j < len(group_counts) and not allsig:
            print('  GROUP COUNT',group_counts[j])
            k = 0
            numsig = 0
            while k < num_runs:
                if k % 20 == 0:
                    print('  run %d of %d' % (k,num_runs))
                data = datagen(means_list,stddev,group_counts[j])
                ## IF/ELSE handles whether we pass in the real variance
                if realvar:
                    issig = anova(data, epsilon_vals[i], fout, stddev**2) 
                else:
                    issig = anova(data, epsilon_vals[i], fout, None) 
                if issig:
                    numsig+=1
                k += 1
            if numsig == num_runs: 
                print('  all values are significant; not running any larger groups.')
                allsig = True
            else:
                print('  %d of %d runs were significant' % (numsig,num_runs))
            j += 1
        i += 1
    fout.close() #close file
    print('Wrote data to ' + filename)
    return


if __name__ == '__main__':
    options = ['estimate','realvar','noisy']
    if len(sys.argv) != 2 or sys.argv[1].lower() not in options:
        print('Usage: python3 anova.py <experiment>\n\t\t<experiment> is one of "estimate","realvar", or "noisy"')
        sys.exit()

    num_runs = 1000
    epsilon_vals = [None,1,.5,.1,.01]
    group_counts = [30,300,3000]
    

    experiment = sys.argv[1].lower()
    if experiment == 'estimate':
        m1 = 0.35
        m2 = 0.5
        m3 = 0.65
        stddev = 0.15
        filename = 'test_estimate_%druns_%.2fm1_%.2fm2_%.2fm3_%.2fvar_%depsilons_%dcounts.csv' % \
            (num_runs,m1,m2,m3,stddev,len(epsilon_vals),len(group_counts))
        anova_test(num_runs,epsilon_vals,filename,[m1,m2,m3],stddev,group_counts,realvar=False)
    elif experiment == 'realvar':
        m1 = 0.35
        m2 = 0.5
        m3 = 0.65
        stddev = 0.15
        filename = 'test_realvar_%druns_%.2fm1_%.2fm2_%.2fm3_%.2fvar_%depsilons_%dcounts.csv' % \
          (num_runs,m1,m2,m3,stddev,len(epsilon_vals),len(group_counts))
        anova_test(num_runs,epsilon_vals,filename,[m1,m2,m3],stddev,group_counts,realvar=True)
    elif experiment == 'noisy':
        means_list = [0.5,0.5,0.5,0.6,0.4,0.43]
        stddev = 0.2
        filename = 'test_noisy_%druns_%dmeans_%.2fvar_%depsilons_%dcounts.csv' % \
            (num_runs,len(means_list),stddev,len(epsilon_vals),len(group_counts))
        anova_test(num_runs,epsilon_vals,filename,means_list,stddev,group_counts,realvar=False)

