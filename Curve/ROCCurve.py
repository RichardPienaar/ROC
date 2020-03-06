'''
ROC curve plotter 

To run:
--> ROCCurve.py param1 param2
params (min1):
 - all (ChIP-R)
 - idr (idr)

Generates:
- Labelled pvalue plot (Yaxis 0-->1)
- Unlabelled pvalue plot (Yaxis0-->maxpvalue)
- ROC curve By Peak
- ROC curve By Base

'''
import bed, sys, ival
import matplotlib.pyplot as plt
import math

def total_bases(bfile):
    '''
    returns count(int) of number of bases in a bedfile
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

def get_points_to_plot(TP,TN,FP,FN,test,all,name2,bybase=True):
    '''
    calculates points for plotting and returns them as xs, ys [float,float,...]
    TP (bfile) --> True Positives
    TN (bfile) --> True Negatives
    FP (bfile) --> False Positives
    FN (bfile) --> False Negatives
    all (bfile) --> All test data
    name2 (str) --> idr / all
    bybase (bool) --> calc by base or by peak
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
    else:
        for entry in test: # calculating values to graph
            data.append((entry, entry.pValue))
            
    if name2 == 'idr':
        data.sort(key = lambda data: data[1], reverse=True) 
    else:
        data.sort(key = lambda data: data[1])

    if bybase:
        for zone in data:
            if zone[2]: # True Positive
                ysum+= len(zone[0])
            else: # False Positive
                xsum+= len(zone[0])
            ys.append(ysum / YDENOM)
            xs.append(xsum / XDENOM)    
    else:
        for zone in data:
            if TP.getOverlap(zone[0]) != []:
                ysum+=1
            else:
                xsum+=1
            ys.append(ysum / YDENOM)
            xs.append(xsum / XDENOM) 
    
    return xs,ys

def draw_ROC_curve(xs,ys,name2,bybase=True):
    '''
    draws ROC curve

    params:
    - xs, ys [float,float,..] -> From get_points_to_plot(TP,TN,FP,FN,test,all,name2,bybase)
    - name2 (str) -> idr / all
    - bybase (bool) -> points generated at peak or base level

    '''
       
    # ROC graph
    axes = plt.gca() 
    axes.set_ylabel('True Positive Rate (Sensitivity)')
    axes.set_xlabel('False Positive Rate (1 - Specificity)')
    axes.set_xlim(0,1)
    axes.set_ylim(0,1)

    fig = plt.gcf()

    if name2.upper() == 'IDR':
        if bybase:
            axes.set_title('ROC Curve - IDR - Base Pair Level')
            
        else:
            axes.set_title('ROC Curve - IDR - Peak Level')
        fig.canvas.set_window_title('ROC Curve - IDR')

    else:
        if bybase:
            axes.set_title('ROC Curve - ChIP-R - Base Pair Level')
        else:
            axes.set_title('ROC Curve - ChIP-R - Peak Level')
        fig.canvas.set_window_title('ROC Curve - ChIP-R')
    
    colours = ('blue','olive','red','orange')
    for i in range(len(xs)):
        plt.plot(xs[i],ys[i],color=colours[i])
        
    plt.legend(('71','72','73','74'), loc='upper left')
    plt.show()   

def graph_Pvalues(test,name2,labelled=True):
    '''
    plots pvalues
    - Labelled (bool) --> labelled axes or not
    - name2 (str) --> idr / all

    NOTE: retransforms idr pvalues before plotting
    '''
    data = [] # to append each entry - presort
    # for graphing
    xsum = 0
    ys = []
    xs = []

    for entry in test: 
        if name2 == 'idr':
            data.append((entry, (10**-entry.pValue))) # transform if idr 
        else:
            data.append((entry, entry.pValue))

    data.sort(key = lambda data: data[1]) # ordering by pvalue
   
    for zone in data:
        xsum+=1
        xs.append(xsum)   
        ys.append(zone[1]) 

    # Graphing 
    axes = plt.gca() 
    axes.set_ylabel('pValue')
    axes.set_xlabel('count')
    
    if not labelled:
        plt.xticks([])
        plt.yticks([])
    
    else:
        axes.set_ylim(0,1)

    fig = plt.gcf()

    if name2.upper() == 'IDR':
        axes.set_title('IDR pValues')
        fig.canvas.set_window_title('IDR pValues') 

    else:
        axes.set_title('ChIP-R pvalues')
        fig.canvas.set_window_title('ChIP-R pvalues') 
       
    plt.plot(xs,ys,'bo')
    plt.show()   

def assert_data(True_Pos,True_Neg,All,False_pos,False_Neg):
    '''
    True if bedtools operations have been performed successfully
    '''
    # checks files are correct (no base is labelled twice)
    if total_bases(True_Pos)+total_bases(True_Neg)+total_bases(False_pos)+total_bases(False_Neg)-total_bases(All) == 0:
        return True
    else:
        # debug
        print(total_bases(True_Pos), total_bases(True_Neg), total_bases(False_pos), total_bases(False_Neg), total_bases(All))
        print(total_bases(True_Pos)+total_bases(True_Neg)+total_bases(False_pos)+total_bases(False_Neg)-total_bases(All))
        return False

def main():
    # accepts 'IDR' and 'all' (or both) as arguments. File names need to follow the formatting. 
    # 71,2,3,4 refer to SR code endings

    print('loading bed files - general')

    for arg in sys.argv[1: ]:

        if str(arg) == 'all': #ChIP-R
            OT = bed.BedFile('m0_all.bed',  format='Peaks')
            all = bed.BedFile('m0_all.bed', format='Peaks')

        elif str(arg) == 'idr':
            OT = bed.BedFile('combined_low.bed',  format='idr') # retest_low.bed is IDR with score instead of pValue
            all = bed.BedFile('m0_all.bed', format='Peaks')

        elif str(arg) == 'MED':
            OT = bed.BedFile('med_m1_all.bed',  format='Peaks')
            all = bed.BedFile('med_m0_all.bed')

        graph_Pvalues(OT,str(arg)) # pvalues, first labelled
        graph_Pvalues(OT,str(arg), labelled=False) # next unlabelled
        xs = []
        ys = [] # by peak
        xs_bb = [] 
        ys_bb = [] # by base

        for i in range(1,5):
            print('loading bed files - specific')

            True_Pos = bed.BedFile('True_Positives_'+str(arg)+'_7'+str(i)+'.bed') 
            True_Neg = bed.BedFile('True_Negatives_'+str(arg)+'_7'+str(i)+'.bed')
            False_Pos = bed.BedFile('False_Positives_'+str(arg)+'_7'+str(i)+'.bed')
            False_Neg = bed.BedFile('False_Negatives_'+str(arg)+'_7'+str(i)+'.bed')

            print('checking data validity',i,'/4')
            if assert_data(True_Pos,True_Neg,all,False_Pos,False_Neg):
                print('calculating points - By Peak')
                nxs,nys = get_points_to_plot(True_Pos,True_Neg,False_Pos,False_Neg,OT,all,str(arg),bybase=False)
                print('calculating points - By Base')
                nxs1, nys1 = get_points_to_plot(True_Pos,True_Neg,False_Pos,False_Neg,OT,all,str(arg),bybase=True)
                xs.append(nxs)
                ys.append(nys)
                xs_bb.append(nxs1)
                ys_bb.append(nys1)
            else:
                raise ValueError

        # draw plots
        print('drawing ROC Curve - By Peak')
        draw_ROC_curve(xs,ys,str(arg),bybase=False)
        print('drawing ROC Curve - By Base')
        draw_ROC_curve(xs_bb,ys_bb,str(arg),bybase=True)

if __name__ == "__main__": 
    main()
