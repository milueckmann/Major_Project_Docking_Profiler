#!/usr/bin/python

#Author: Michael Lueckmann, PhD student, University of Copenhagen
#Purpose: Generate overview plots from a VLS-docking run with ICM, using a SDF-file containing ICM-VLS scores as input


import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import glob
import math
import re
import sys
from matplotlib.backends.backend_pdf import PdfPages
from operator import itemgetter

####SET PARAMETERS TO PLOT AT THE END OF THIS SCRIPT####

####SET INPUT FILENAME####
sdfFile = "Model_input_unsorted.sdf"

####SET THRESHOLD IN PERCENT####
threshold = 10

####SET FILENAME FOR PDF OUTPUT#### 
pp = PdfPages("Plot_Docking_Profile.pdf")


##########################################
####NOTHING SHOULD BE CHANGED FROM HERE###
##########################################

pattern_molecule = re.compile("\$\$\$\$") # regex search pattern for end of molecule in an SDF file

####CHECK IF INPUT IS IN SDF FILE FORMAT####
#check if "$$$$" pattern is in input file
def checkinput():
    with open(sdfFile) as f:
        found = False
        for line in f:
            if re.match(pattern_molecule,line):
                found = True

    if found == False:
        print "Error: input file must be in SDF format" 
        sys.exit()

checkinput()

####GENERATE VALUE LIST####
#parses all associated information from input SDF file into a value list ("scorelist.tsv")
paramlist = []
flength = 0 #counter for total number of data points

with open(sdfFile) as f, open("scorelist.tsv", "w") as b:
    for line in f:
        if ">" in line: #search pattern for associated parameters
            line = re.sub('<', '', line) # clean line
            headerline = line.strip("<>\n").split(" ") #clean line
            b.write("%s " % (headerline[1]))
            paramlist.append(headerline[1])
        if ("$$$$") in line: #search pattern for end of molecule
            b.write("\n")
            break
    for line in f:
        if ">" in line:
            nextline = next(f)
            valueline = nextline.strip(" \n\t")
            b.write("%s " % (valueline))
        if ("$$$$") in line:
            flength+=1 #counter for total nuber of data points, used for bars number in histogram plotter
            b.write("\n")

scorelist = "scorelist.tsv"

####GET "SCORE" COLUMN NUMBER####
#find column number of "Score" values in SDF file
def getscorecolumnnumber():
    for i,element in enumerate(paramlist):
        if element == "Score":
            return i

column_list = list(range(1,len(paramlist))) #generate list of numbers depending on number of associated parameters
int_list = map(int, column_list) #convert to integers, later used in the plotter functions

####SORT BY SCORES####
#sort all molecules by their docking score ("Score")
with open("scorelist.tsv", "r") as a, open("scorelist_sorted.tsv", "w") as b:
    header = " ".join(paramlist)
    b.write("%s\n" % (header))
    scorecolumnno = getscorecolumnnumber()
    next(a) #skip header line
    for line in a:
        lines = [line.split() for line in a]
    for element in lines: #convert to list of lists of mixed string (NAME) and float (all other values)
        for i,j in enumerate(element):
            try:
                element[i] = float(j)
            except ValueError:
                pass
    lines.sort(key=itemgetter(scorecolumnno)) #sort by docking scores
    str_lines = [[str(j) for j in element] for element in lines] #converts back to a list of lists of strings
    for element in str_lines:
        clean_line = " ".join(element)
        b.write("%s\n" % (clean_line))


####TAKE TOP X % BY SCORES####
#extracts top X (see threshold) percent by score out of all molecules 
with open("scorelist_sorted.tsv", "r") as a, open("scorelist_sorted_topX.tsv", "w") as b:
    header = " ".join(paramlist)
    b.write("%s\n" % (header))
    top_lines = (flength + threshold // 2) // threshold #python integer division with rounding
    next(a) #skip header line
    i=0 #counter for top X % lines
    for line in a:
        if i < top_lines:
            b.write("%s" % (line))
        i+=1

####CALCULATE AVERAGES####
#calculate averages for all associated parameters of top X % and of all molecules
with open("scorelist.tsv", "r") as a, open("scorelist_sorted_topX.tsv", "r") as b, open("averages.tsv", "w") as c:
    for line in a:
        av_lines = [line.split() for line in a]
        for element in av_lines: #delete "NAME" entries
            del element[0]
    int_lines = [[int(float(j)) for j in element] for element in av_lines] #convert to int/float
    whole_list = np.array(int_lines)
    whole_average = np.mean(whole_list, axis=0) #calculate averages of whole scorelist
    whole_average_line = " ".join(str(e) for e in whole_average) #convert to string
    av_header = " ".join(paramlist[1:])
    c.write("%s\n%s\n%s\n\n" % ("WHOLE LIST AVERAGES", av_header, whole_average_line))
    for line in b:
        top_av_lines = [line.split() for line in b]
        for element in top_av_lines: #delete names
            del element[0]
    top_int_lines = [[int(float(j)) for j in element] for element in top_av_lines] #convert to int/float
    top_list = np.array(top_int_lines)
    top_average = np.mean(top_list, axis=0) #calculate averages of top % scorelist
    top_average_line = " ".join(str(e) for e in top_average) #convert to string
    c.write("%s%s%s\n%s\n%s\n\n" % ("TOP ", threshold, " % LIST AVERAGES", av_header, top_average_line))

#########################
####PLOTTER FUNCTIONS####
#########################

####SCATTERPLOTTER####
#makes a simple scatter plot of any two parameters
def scatterplotter(xval,yval):
    load=np.loadtxt(scorelist, skiprows=1, usecols=(int_list))
    for i,element in enumerate(paramlist, start=-1):
        if element == "NAME":
           pass 
        else:
            data=load[:,int(i)]
            if xval == element:
                x=data
                xlabel=element
            if yval == element:
                y=data
                ylabel=element
                r_matr = np.ma.corrcoef(x,y) #calculates correlation cofficient
                r_list = r_matr.tolist() #converts matrix to list
    scat_fig = plt.figure()
    plt.plot(x,y,"bo")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.annotate("r =", xy=(0.02, 0.95), xycoords='axes fraction') #annotates correlation coefficient in top left corner
    plt.annotate(r_list[0][1], xy=(0.065, 0.95), xycoords='axes fraction')
    return scat_fig

####HISTOGRAM PLOTTER####
#makes a simple histogram of any parameter
def histoplotter(xval):
    load=np.loadtxt(scorelist, skiprows=1, usecols=(int_list))
    for i,element in enumerate(paramlist, start=-1):
        if element == "NAME":
           pass
        else:
            data=load[:,int(i)]
            if xval == element:
                x=data
                xlabel=element
    histo_fig = plt.figure()
    b = math.sqrt(flength) #bins number = square root of total number of data points
    n, bins, patches = plt.hist(x, b, normed=1, facecolor='green', alpha=0.75)
    plt.xlabel(xlabel)
    plt.ylabel("ratio")
    return histo_fig

####COMPARISON BAR PLOTTER####
#generates a comparative plot of average values of the whole list versus the list with the top X %
def average_compare():
    comp_fig = plt.figure() 
    whole_average_list = np.array(whole_average).tolist() #convert array to list
    top_average_list = np.array(top_average).tolist()  
    av_xlabel = paramlist[1:] #skipping the first "NAME" column
    N = len(whole_average_list)
    width = 0.3
    index = np.arange(N)
    rects1 = plt.bar(index+width, whole_average_list, width, color = "blue", label = "whole list averages", alpha = 0.4)
    rects2 = plt.bar(index+width*2, top_average_list, width, color = "red", label = "top ranked averages", alpha = 0.4)
    plt.xticks(index+width*1.5, (av_xlabel), fontsize=7)
    plt.xlabel("Parameter")
    plt.ylabel("Value")
    plt.legend()
    return comp_fig


####SET PARAMETERS TO PLOT BELOW####
plot0 = average_compare()
plot1 = scatterplotter("Score","Natom")
plot2 = scatterplotter("Score","Nflex")
plot3 = scatterplotter("Score","Hbond")
plot4 = scatterplotter("Score","Hphob")
plot5 = scatterplotter("Score","VwInt")
plot6 = scatterplotter("Score","Eintl")
plot7 = scatterplotter("Score","Dsolv")
plot8 = scatterplotter("Score","SolEl")
plot9 = scatterplotter("Score","mfScore")
plot10 = scatterplotter("Score","RecConf")
plot11 = scatterplotter("Score","molLogP")
plot12 = histoplotter("molLogP")

pp.savefig(plot0)
pp.savefig(plot1)
pp.savefig(plot2)
pp.savefig(plot3)
pp.savefig(plot4)
pp.savefig(plot5)
pp.savefig(plot6)
pp.savefig(plot7)
pp.savefig(plot8)
pp.savefig(plot9)
pp.savefig(plot10)
pp.savefig(plot11)
pp.savefig(plot12)

pp.close()
