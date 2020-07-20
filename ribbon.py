import matplotlib.pyplot as plt
import numpy as np
import csv, sys

def parse(e_file):
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
                if i %10000 == 0:
                    xs.append(float(row[0]))
                    ys.append(float(row[1]))                        
    return xs,ys

def parse_all(prefix):
    all = []
    i = 0

    while True:
        try:
            file = prefix+str(i)+'.csv'       
            all.append(parse(file))
            print('woop', len(all[i][0]))
            i+=1
        except:
            break

    return all

def get_averages_and_errors(parsed):
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
    
    x = np.array(xs, dtype=float)
    y = np.array(average, dtype=float)
    low =  np.array(low_bounds, dtype=float)
    high = np.array(high_bounds, dtype=float)

    return x,y,low,high

def ribbon(args):
    averages = []

    for arg in args:
        to_plot = parse_all(str(arg))
        averages.append(get_averages_and_errors(to_plot))

    fig, ax = plt.subplots()
    for group in averages:
        x,y,low,high = make_arrays(group[0], group[1],group[2],group[3])
        ax.plot(x, y, '-')
        ax.fill_between(x, high, low, alpha=0.2)

    plt.show()

def produce_ribbons(arguments):
    prefixes = []
    for arg in arguments:
        prefixes.append(str(arg))
    
    ribbon(prefixes)
    return
    
