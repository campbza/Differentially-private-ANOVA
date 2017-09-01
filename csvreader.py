import csv
import matplotlib.pyplot as plt
import sys
import os

def plot_lines(sizes, eps_vals, significance, outfile, threshold, title):
    x_vals = sizes
    y_vals = significance
    labels = [eps_vals[i] for i in range(len(eps_vals))]
    for i in range(len(labels)):
        labels[i] = 'epsilon = ' + str(labels[i])
    colors = ['r','g','b','y','m','c','k']
    shapes = ['o','s','d','^','x','<','v']
    fig=plt.figure(figsize=(8,5))
    ax=plt.subplot(1,1,1)
    for i in range(len(y_vals)):
        style_code = '-'+colors[i]+shapes[i]
        ax.plot(x_vals,y_vals[i],style_code,label=labels[i],lw=2,ms=8)
    ax.set_xscale('log')
    ax.set_xlim([0,max(x_vals) + 1])
    plt.title(title)
    plt.xlabel('Database size (log scale)')
    plt.ylabel('Percent significance at ' + str(threshold))
    plt.legend(loc='lower right')

    plt.tight_layout()
    plt.savefig(outfile)
    print('Wrote to ' + outfile)
    return

def plot_lines_allow_truncated(x_vals,y_vals,labels, outfile, threshold, title):
    all_xvals = []
    for x in x_vals:
        all_xvals += x
    all_xvals = list(set(all_xvals))
    all_xvals.sort()

    new_yvals = []
    for i in range(len(y_vals)):
        y = y_vals[i]
        x = x_vals[i]
        val_dict = {x[j]:y[j] for j in range(len(y))}
        for x in all_xvals:
            if x not in val_dict:
                val_dict[x] = 1.0
        x_vals[i] = all_xvals
        y_vals[i] = [val_dict[s] for s in all_xvals]

    for i in range(len(labels)):
        if labels[i] == 'None':
            labels[i] = 'epsilon = $\infty$'
        else:
            labels[i] = 'epsilon = ' + str(labels[i])
    colors = ['r','g','b','y','m','c','k']
    shapes = ['o','s','d','^','x','<','v']
    fig=plt.figure(figsize=(8,5))
    ax=plt.subplot(1,1,1)
    for i in range(len(y_vals)):
        style_code = '-'+colors[i]+shapes[i]
        ax.plot(x_vals[i],y_vals[i],style_code,label=labels[i],lw=2,ms=8)
        print(labels[i])
        print('   n    proportion')
        for j in range(len(x_vals[i])):
            print('  ',x_vals[i][j],y_vals[i][j])
    ax.set_xscale('log')
    ax.set_xlim([0,max([max(x) for x in x_vals]) + 1])
    ax.set_ylim([0,1.05])
    plt.title(title)
    plt.xlabel('Database size (log scale)')
    plt.ylabel('Percent significance at ' + str(threshold))
    plt.legend(loc='best')

    plt.tight_layout()
    plt.savefig(outfile)
    print('Wrote to ' + outfile)
    if 'pdf' in outfile:
        os.system('pdfcrop %s %s' % (outfile,outfile))
    return

def pvals_significance(infile, outfile, graphtitle, threshold):
#input: csv file, outfile name (.png), graphtitle, and threshold for significant pvals
#output: graph of infile data
    file_contents = open(infile, 'r')
    reader = csv.reader(file_contents)
    pval_dict = {}
    sizes = []
    epsilons = []
    for line in reader:
        size = int(line[0])
        if size not in sizes:
            sizes.append(size) #gather list of sizes (these are total size of the database)
        epsilon = line[3]
        if epsilon not in epsilons:
            epsilons.append(epsilon) #gather list of epsilon values
        key = (size,epsilon)
        if key in pval_dict: #populate dictionary keyed by size,epsilon tuples
            pval_dict[key].append(float(line[4])) #line[4] in csv is the pval
        else:
            pval_dict[key] = [float(line[4])]
    for item in pval_dict: #replace the lists with a percent of how many in that list are significant
        l = len(pval_dict[item])
        sig_pvals = [x for x in pval_dict[item] if x < threshold]
        m = len(sig_pvals)
        pval_dict[item] = float(m) / float(l)
    sizes.sort()
    epsilon_dict = {}
    labels = []
    significance = [] #list of list of percent significance, each list corresponding to epsilon
    for e in epsilons:
        epsilon_dict[e] = []
    for s in sizes:
        for e in epsilons:
            epsilon_dict[e].append(pval_dict[(s,e)])
            if epsilon_dict[e] not in significance:
                significance.append(epsilon_dict[e])
    plot_lines(sizes, epsilons, significance, outfile, threshold, graphtitle)
    return 

def pvals_significance_allow_truncated_lines(infile, outfile, graphtitle, threshold):
#input: csv file, outfile name (.png), graphtitle, and threshold for significant pvals
#output: graph of infile data
    file_contents = open(infile, 'r')
    reader = csv.reader(file_contents)
    epsilon_dict = {} # dictionary of (epsilon_val: size_dictionary) pairs, where each value is (size,pval_list).
    for line in reader:
        if len(line)<5: # line is truncated; skip
            continue
        size = int(line[0])
        epsilon = line[3]
        if epsilon != 'None':
            epsilon = float(epsilon)
        pval = float(line[4])
        if epsilon not in epsilon_dict:
            epsilon_dict[epsilon] = {}
        if size not in epsilon_dict[epsilon]:
            epsilon_dict[epsilon][size] = []
        epsilon_dict[epsilon][size].append(pval)
    for epsilon in epsilon_dict:
        for size in epsilon_dict[epsilon]:
            num_runs = len(epsilon_dict[epsilon][size])
            epsilon_dict[epsilon][size] = sum([1 for x in epsilon_dict[epsilon][size] if x < threshold])/float(num_runs)

    # sort epsilons
    epsilon_keys = epsilon_dict.keys()
    if 'None' in epsilon_keys:
            extra = ['None']
    else:
        extra = []
    to_sort = [s for s in epsilon_keys if s != 'None']
    to_sort.sort(reverse=True)
    epsilon_keys = extra + [s for s in to_sort]

    x_vals = []
    y_vals = []
    labels = []
    for epsilon in epsilon_keys:
        sizes = [key for key in epsilon_dict[epsilon].keys()]
        sizes.sort()

        values = [epsilon_dict[epsilon][size] for size in sizes]
        x_vals.append(sizes)
        y_vals.append(values)
        labels.append(epsilon)

    plot_lines_allow_truncated(x_vals,y_vals,labels, outfile, threshold, graphtitle)
    return 

import sys
if __name__ == '__main__':
    if len(sys.argv) != 4:
        print('ERROR! <input.csv> <output.png> <graph_title> expected.')
        sys.exit()
    inputfile = sys.argv[1]
    outputfile = sys.argv[2]
    graphtitle = sys.argv[3]
    threshold = 0.05
    print(graphtitle)
    pvals_significance_allow_truncated_lines(inputfile,outputfile,graphtitle,threshold)

