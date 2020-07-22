'''
ROC curve plotter 

(c) Richard Pienaar - 2020
University of Queensland  

'''
import ribbon
import os
import bed, sys, ival
import csv

def get_prefixes():
    '''
    Returns list [str,str..] of prefix names - experimental group

    Works by going into the Test_data subfolder
    '''
    file_names = os.listdir(".\Test_data")
    prefixes = []
    for filename in file_names:
        if filename.endswith('.bed'):
            prefixes.append(filename[:-4])
        else:
            prefixes.append(filename)
    return prefixes

def total_bases(bfile):
    '''
    returns count(int) of number of bases in a bedfile (bed.Bedfile)
    '''
    count = 0
    try:
        for entry in bfile:
            interval = entry.getInterval()
            count+= len(interval)
    except:
        return 0
    return count

def count_entries(P,entry,T=True,metric='pValue'):
    '''
    Returns list of all overlapping regions of bedEntry (entry) in bedfile (P) as:
     (Interval(interest of P and entry),metric,T)
     where: 
     metric = p.value or score
     and T is whether its a True(T=True) or False(T=False) Positive entry   
    '''
    OlP = P.getOverlap(entry)
    out = []
    if metric == 'pValue':
        if OlP != []:
            for region in OlP:
                out.append((ival.isect(region.getInterval(),entry.getInterval()), entry.pValue, T))
    else:
        if OlP != []:
            for region in OlP:
                out.append((ival.isect(region.getInterval(),entry.getInterval()), entry.score, T))
    return out

def get_points_to_plot(TP,TN,FP,FN,test,all,name2,bybase=True,precise=True):
    '''
    calculates points for ROC curve plotting and returns them as xs, ys [float,float,...]
    TP (bfile) --> True Positives
    TN (bfile) --> True Negatives
    FP (bfile) --> False Positives
    FN (bfile) --> False Negatives
    all (bfile) --> All test data
    name2 (str) --> name of the experimental group
    bybase (bool) --> calc by base or by peak
    precise (bool) --> Reports every single base when True, otherwise collapses:
                - e.g if the region 80-86 is all Tp, precise True returns x.y values corresponding to 
                - 80,81,82,83,84,85,86 whereas i precise is False only the vx,y values of 80 and 86 are returned (the inbetween can be inferred) 
    '''
    empty = total_bases(FN) # check for fully encompassing test data
    if bybase:
        YDENOM = total_bases(TP)+total_bases(FN) 
        XDENOM = total_bases(FP)+total_bases(TN)
    else:
        YDENOM=0
        XDENOM=0

        for entry in all:
            if TP.getOverlap(entry) != []:
                YDENOM+=1
            elif empty != 0:
                if FN.getOverlap(entry) != []:
                    YDENOM+=1
                else:
                    XDENOM+=1
            else:
                XDENOM+=1
    data = [] # each entry - presort
    # for graphing
    ysum = 0 
    xsum = 0
    ys = []
    xs = []

    if bybase:
        for entry in test: # calculating values to graph
            TPs = count_entries(TP,entry,metric='pValue')
            FPs = count_entries(FP,entry,False,metric='pValue')
            for ol in TPs:
                data.append(ol)
            for ol in FPs:
                data.append(ol)
        
        # sort by T/F 
        data.sort(key = lambda data: data[2])
        
        # Correction for overlaps overlapping due to the way ival.py works
        i =0
        while True:
            
            try:
                inter = ival.isect(data[i][0],data[i+1][0])
                if inter != None :
                    data.append((ival.isect(data[i][0],data[i+1][0]),data[i][1],data[i][2]))
                    data.remove(data[i])
                    data.remove(data[i+1])
                else:
                    i+=1
            except:
                break
    else:
        for entry in test: # calculating values to graph
            data.append((entry, entry.pValue))
            
    if 'idr' in name2.lower():
        data.sort(key = lambda data: data[1], reverse=True) 
    else:
        data.sort(key = lambda data: data[1])

    if bybase:
        if not precise:

            for zone in data:
                if zone[2]: # True Positive
                    ysum+= len(zone[0])
                else: # False Positive
                    xsum+= len(zone[0])
                ys.append(ysum / YDENOM)
                xs.append(xsum / XDENOM)    
        
        else: # EVERY DATA POINT
            for zone in data: 
                if zone[2]: # True Positive
                    oldsum = ysum
                    while ysum < (oldsum+ len(zone[0])):
                        ysum+= 1
                        xs.append(xsum / XDENOM)
                        ys.append(ysum / YDENOM)
                else: # False Positive
                    oldsum = xsum
                    while xsum < (oldsum + len(zone[0])):
                        xsum+=1
                        xs.append(xsum / XDENOM)
                        ys.append(ysum / YDENOM)
    else:
        for zone in data:
            if TP.getOverlap(zone[0]) != []:
                ysum+=1
            else:
                xsum+=1
            ys.append(ysum / YDENOM)
            xs.append(xsum / XDENOM) 
    
    return xs,ys

def save_csv(name,points):
    '''
    Saves ROC points (from get_points_to_plot()) as csv for parsing later
    FIXME: Maybe make not to the general directory
    '''
    with open(name+'.csv', 'w',newline='') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=',',dialect='excel')
        for e in range(len(points[0])):
            spamwriter.writerow([points[0][e],points[1][e]])
        csvfile.close() 

def assert_data(True_Pos,True_Neg,All,False_pos,False_Neg):
    '''
    True if bedtools operations have been performed successfully
    True_Pos,True_Neg,All,False_pos,False_Neg - <bed.Bedfile>
    '''
    # checks files are correct (no base is labelled twice)
    if total_bases(True_Pos)+total_bases(True_Neg)+total_bases(False_pos)+total_bases(False_Neg)-total_bases(All) == 0:
        return True
    else:
        # debug
        print(total_bases(True_Pos), total_bases(True_Neg), total_bases(False_pos), total_bases(False_Neg), total_bases(All))
        print(total_bases(True_Pos)+total_bases(True_Neg)+total_bases(False_pos)+total_bases(False_Neg)-total_bases(All))
        return False

def largest_test_set(prefixes):
    '''
    Returns <bed.BedFile> the largest test data set, must be all encompasaing i.e every point is labelled
    '''
    winner = 0
    winnerfile = None
    for expgroup in prefixes: 
        try:
            file = bed.BedFile(".\Test_data\\"+str(expgroup)+".bed",format='Peaks')
        except:
            file = bed.BedFile(".\Test_data\\"+str(expgroup)+".bed",format='idr')
        
        winningcount = total_bases(file)
        if winningcount > winner:
            winnerfile = file
            winner = winningcount

    return winnerfile

def name_groups(prefixes):
    '''
    Request user input on the exact name of each experimental group, for graphing
    returns [str,str,str..] - the names 
    '''
    names =[]
    for pref in prefixes:
        names.append(input('Please specify this groups full name '+str(pref)+" : "))
    return names

def main():
    '''
    Generates ROC curves from bedfiles (.bed)
    Assuming correct folder management, should return comparitive ribbon plots of all the .bed files of each experimental group in the required subdirectories
    ''' 
    # file management stuff
    prefixes = get_prefixes()
    names = name_groups(prefixes) # full names for graphing 
    runs = os.listdir(".\True_Positives")
    Directorynames= [".\True_Positives",".\False_Positives",".\True_Negatives",".\False_Negatives"]
    all = [] # used to infer the appropriate grouping for each run 

    # check for existence of csv files already
    try:
        ribbon.ribbon(prefixes,names)
    
    except: 
        print('loading bed files - General')
        for prefix in prefixes:
            try:
                all.append(bed.BedFile('.\Test_Data\\'+str(prefix)+'.bed',format='Peaks'))
            except:
                all.append(bed.BedFile('.\Test_Data\\'+str(prefix)+'.bed',format='idr'))

        # loading file steps
        for expgroup in prefixes:
            TPs,FPs,TNs,FNs = [],[],[],[]
            print('loading bed files - '+str(expgroup))
            for subdir in Directorynames:
                for bfile in runs: 
                    if str(expgroup) == str(bfile)[0:len(str(expgroup))] and str(bfile)[len(str(expgroup))] == '_':
                        print(str(subdir)+"\\"+str(bfile))
                        if 'idr' not in str(expgroup).lower():
                            if subdir == Directorynames[0]:
                                TPs.append(bed.BedFile(str(subdir)+"\\"+str(bfile),format='Peaks'))
                            elif subdir == Directorynames[1]:
                                FPs.append(bed.BedFile(str(subdir)+"\\"+str(bfile),format='Peaks'))
                            elif subdir == Directorynames[2]:
                                TNs.append(bed.BedFile(str(subdir)+"\\"+str(bfile),format='Peaks'))
                            else:
                                FNs.append(bed.BedFile(str(subdir)+"\\"+str(bfile),format='Peaks'))
                        else:
                            if subdir == Directorynames[0]:
                                TPs.append(bed.BedFile(str(subdir)+"\\"+str(bfile),format='idr'))
                            elif subdir == Directorynames[1]:
                                FPs.append(bed.BedFile(str(subdir)+"\\"+str(bfile),format='idr'))
                            elif subdir == Directorynames[2]:
                                try:
                                    TNs.append(bed.BedFile(str(subdir)+"\\"+str(bfile),format='idr')) # check for basic data - resulting from bedtools
                                except:
                                    TNs.append(bed.BedFile(str(subdir)+"\\"+str(bfile),format='simple'))
                            else:
                                try:
                                    FNs.append(bed.BedFile(str(subdir)+"\\"+str(bfile),format='idr'))              
                                except:
                                    FNs.append(bed.BedFile(str(subdir)+"\\"+str(bfile),format='simple'))              

            try:
                test = bed.BedFile(".\Test_data\\"+str(expgroup)+".bed",format='Peaks')
            except:
                test = bed.BedFile(".\Test_data\\"+str(expgroup)+".bed",format='idr')
            
            # Calculation steps
            for i in range(len(TPs)):
                for possibility in all:
                    if assert_data(TPs[i],TNs[i],possibility,FPs[i],FNs[i]):
                        print('calculating points - Peak level', str(expgroup), str(i))
                        points = get_points_to_plot(TPs[i],TNs[i],FPs[i],FNs[i],test,possibility,expgroup,bybase=False) # Peak level
                        save_csv(str(expgroup)+"_"+str(i)+'Peaks',points) # write output to csv
                        
                        print('calculating points - Base level', str(expgroup), str(i))
                        points = get_points_to_plot(TPs[i],TNs[i],FPs[i],FNs[i],test,possibility,expgroup,bybase=True,precise=True) # Base level
                        save_csv(str(expgroup)+"_"+str(i)+'Bases',points) # write output to csv

        ribbon.ribbon(prefixes,names) # Graph curves with newly generated data


if __name__ == "__main__": 
    main()
