import matplotlib.pyplot as plt
import numpy as np
import csv, sys

def parse(e_file):
    '''
    Returns x<float> and y<float> coordinates [xs], [ys] from a given excel file as produced by ROCCurve.py - save_csv(name,points)
    '''
    xs =[]
    ys =[]
    if not 'base' in str(e_file):
        with open(str(e_file), 'r',newline='') as csvfile:
            spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
            i = 0
            for row in spamreader:
                i+=1
                if i %10 ==0:
                    xs.append(float(row[0]))
                    ys.append(float(row[1]))
    else:
        with open(str(e_file), 'r',newline='') as csvfile:
            spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
            i = 0
            for row in spamreader:
                i+=1
                if i %100000 == 0:
                    xs.append(float(row[0]))
                    ys.append(float(row[1]))                        
    return xs,ys

def parse_all(prefix, peaks=True):
    '''
    Parses all files with a given prefix <str> in one batch
    peaks(bool) - True if parsing peak output files
    return: all [([<float>,<float>,..],[<float>,<float>,..]),([<float>,<float>,..],[<float>,<float>,..])...]
    '''
    all = []
    i = 0

    while True:
        if peaks:
            try:
                file = prefix+'_'+str(i)+'Peaks.csv'       
                all.append(parse(file))
                i+=1
            except:
                break
        else:
            try:
                file = prefix+'_'+str(i)+'Bases.csv'       
                all.append(parse(file))
                i+=1
            except:
                break

    return all

def get_averages_and_errors(parsed):
    '''
    Averaging function for combining replicates of the same run
    parsed <list> - output of parse_all(prefix)
    returns:
    Xs <[float,float,..]> horizontal steps to graph, either no. peaks or no. bases
    average_Ys <[float,float,..]> corresponding mean values 
    high_bounds <[float,float,..]> corresponding maximum y value for each x value
    low_bounds <[float,float,..]> corresponding maximum x value for each y value
    '''
    Xs = []
    average_Ys = []
    high_bounds = []
    low_bounds = []

    replicates = len(parsed)

    for i in range(len(parsed[0][0])):
        Xs.append(parsed[0][0][i])
        sumy = 0
        low = 1
        high = 0
        for y in range(replicates):
            pointy = parsed[y][1][i]
            sumy+=pointy
            if pointy < low:
                low = pointy
            if pointy > high:
                high = pointy

        average_Ys.append(sumy/replicates)
        low_bounds.append(low)
        high_bounds.append(high)
    
    return Xs, average_Ys, high_bounds, low_bounds


def make_arrays(xs,average,low_bounds,high_bounds):
    '''
    Transforms output of get_averages_and_errors(parsed) from lists to numpy arrays for use in ribbon(args,names)
    '''
    x = np.array(xs, dtype=float)
    y = np.array(average, dtype=float)
    low =  np.array(low_bounds, dtype=float)
    high = np.array(high_bounds, dtype=float)

    return x,y,low,high

def ribbon(args,names):
    '''
    Graphs ROC curves (ribbon plots) 
    args [<str>,<str>,..] the prefixes of each experimental group
    names [<str>,<str>,..] the full names of each experimental group to display
    '''
    averages = []

    # Peaks
    
    for arg in args:
        to_plot = parse_all(str(arg),peaks=True)
        averages.append(get_averages_and_errors(to_plot))
    
    print('Graphing at peak level')
    fig, ax = plt.subplots()
    ax.set_ylabel('True Positive Rate (Sensitivity)')
    ax.set_xlabel('False Positive Rate (1 - Specificity)')
    for group in averages:
        x,y,low,high = make_arrays(group[0], group[1],group[2],group[3])
        ax.plot(x, y, '-')
        ax.fill_between(x, high, low, alpha=0.2)
    
    title = input('Please set a title for the graph - at the Peak Level')
    ax.set_title(title)
    fig.canvas.set_window_title(title)
    plt.legend(names)
    plt.show()
    
    # Bases

    averages = []

    print('Graphing at base level')
    for arg in args:
        to_plot = parse_all(str(arg),peaks=False)
        averages.append(get_averages_and_errors(to_plot))


    fig, ax = plt.subplots()
    ax.set_ylabel('True Positive Rate (Sensitivity)')
    ax.set_xlabel('False Positive Rate (1 - Specificity)')
    for group in averages:
        x,y,low,high = make_arrays(group[0], group[1],group[2],group[3])
        ax.plot(x, y, '-')
        ax.fill_between(x, high, low, alpha=0.2)
    
    title = input('Please set a title for the graph at the Base Level')
    ax.set_title(title)
    fig.canvas.set_window_title(title)
    plt.legend(names)

    plt.show()
  
