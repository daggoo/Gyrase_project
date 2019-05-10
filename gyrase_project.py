# coding=utf-8

from fitter import Fitter
from os import listdir
from os.path import isfile, join
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio import SeqUtils as su
import scipy
from scipy import signal
from scipy.stats import poisson, gamma
import math
from math import exp
from math import factorial
from Parsers_lib import *
from Bio.SeqUtils import GC as GC_count
from scipy.optimize import curve_fit
from Bio import motifs
from Bio.Alphabet  import IUPAC 
from Bio.Seq import Seq
from random import shuffle
from scipy.stats.stats import pearsonr   
import random
from IPython.parallel import Client
import time
from matplotlib_venn import venn2,venn3, venn3_circles

#path = '/data/Gyrase/Pipe/Fragments_ends/'.decode("utf-8")

#path = '/data/Gyrase/Pipe/Fragments_ends/Data_second_mapping/case_6_AC_SW_1/'
path = '/data/Gyrase/Pipe/Fragments_ends/Data_second_mapping/Final_data/25Sept/'
#path = 'F:\Dropbox\_Работа\CRISPR\Girase\Fragments_ends\\fitting\Fragments_ends\\'.decode("utf-8")
#outfile_list_of_peaks = open("F:\Dropbox\_Работа\CRISPR\Girase\Fragments_ends\\fitting\Fragments_ends\\List_of_peaks.txt".decode("utf-8"), 'w+')

genome = open(u"/data/Gyrase/IGV_tracks/E_coli_w3110_G_Mu.fasta", 'r')
plasmidSC101 = open(u"/data/Gyrase/Pipe/pSC101.fasta", 'r')
plasmidBR322 = open(u"/data/Gyrase/Pipe/pBR322.fasta", 'r')
plasmidMu = open(u"/data/Gyrase/Pipe/Mu_bacteriophage.fasta", 'r')

plasmidUC19 = open(u"/data/Gyrase/Pipe/pUC19.gb", 'r')
genome_U00096_2 = open(u"/data/Gyrase/Pipe/Fragments_ends/Data_second_mapping/Transcription/U00096_2.fasta", 'r')
#genome = open(u"F:\Dropbox\_Работа\CRISPR\Girase\E_coli_w3110_G_Mu.fasta", 'r')

for record in SeqIO.parse(genome, "fasta"):
    genomefa = record.seq
    
for record in SeqIO.parse(genome_U00096_2, "fasta"):
    genomefa_2 = record.seq
genome2=SeqIO.read('/data/Gyrase/Pipe/Fragments_ends/Data_second_mapping/Transcription/U00096_2.gb','genbank')
genome_w3110=SeqIO.read('/data/Gyrase/Pipe/Fragments_ends/Data_second_mapping/Transcription/ncbi_w3110.gb','genbank')
print len(genomefa)
print len(genomefa_2)

    
for record in SeqIO.parse(plasmidSC101, "fasta"):
    pSC101 = record.seq
    
for record in SeqIO.parse(plasmidBR322, "fasta"):
    pBR322 = record.seq	
    
for record in SeqIO.parse(plasmidUC19, "genbank"):
    pUC19=record.seq
    
for record in SeqIO.parse(plasmidMu, "fasta"):
    pMu=record.seq

 
seq = pSC101[(4730):(4750)]
seq = pBR322[(980):(1000)]
#print seq
# onlyfiles = [f for f in listdir((path)) if isfile(join((path), f))]
files = listdir((path).encode(sys.getfilesystemencoding()))

wigfiles = []
for file in files:
    if file.endswith(".wig"):
        wigfiles.append(file)

xcoord=np.arange(0,4647999)

d1_b = 274750
d1_e = 371950
d2_b = 793850
d2_e = 807450
d3_b = 1199000
d3_e = 1214100

#win_width = 110
win_width_l = 66
win_width_r = 66
win_width=110

def get_value(i, ends):

    if i < 0:
        j = len(ends) + i
    elif i >= len(ends):
        j = i - len(ends)
    elif i >= d1_b and i <= d1_e:
        j = d1_e + i - d1_b + 1
    elif i >= d2_b and i <= d2_e:
        j = d2_e + i - d2_b + 1
    elif i >= d3_b and i <= d3_e:
        j = d3_e + i - d3_b + 1
    else:
        j=i

    return ends[j]

def return_averaged_line(ends):
    y_coords=[]
    mean = 0.0
    window=100000
    window_float=window*1.0

    mean=0.0
    for i in range(-1 * window, window):
        mean = mean + get_value(i,ends)
    mean=mean/(2*window_float)
    y_coords.append(mean)
    last=len(ends)

    for i in range(1, len(ends)):
	if (i >= d1_b and i <= d1_e) or (i >= d2_b and i <= d2_e) or (i >= d3_b and i <= d3_e) :
	    mean = mean + (get_value(i + window, ends) - get_value(i - window, ends)) / (2 * window_float)
	    y_coords.append(mean)
	else:
	    mean=mean + (get_value(i+window,ends) - get_value(i-window,ends))/(2*window_float)
	    y_coords.append(mean)	

    return y_coords

def peaks_distribution_over_genome():
	#peaks_coords=open(path+ "List_of_peaks.txt", 'r')
        peaks_coords=open(path+ "List_of_peaks.txt", 'r')
        fig=plt.figure(figsize=(20,3))  
        f,axarr = plt.subplots(9,1)  
	
	i=0;
	for line in peaks_coords:
            peaks=[]
	    line=line.rstrip()
	    
	    line=line.split('\t')
	    sample= str(line[0])
	    print sample
            for peak in line[1:len(line)]:
		peak=peak.split('-')
                peaks.append(float(peak[0]))
	    print len(peaks)           
 
            bins=np.linspace(0, 4647454, 500)
            axarr[i].set_xlim(0, 4647454)

            ticks1=[0, 500000, 1000000, 1500000, 2000000, 2500000, 3000000, 3500000, 4000000, 4500000]
            xticknames1=['', '500kb', '1000kb', '1500kb', '2000kb', '2500kb', '3000kb', '3500kb', '4000kb', '4500kb']
            ticks2=[223771, 2729394, 3428115, 3469603, 3600731, 3694454, 4212776]
            xticknames2=['rRNA H', 'rRNA G ', 'rRNA D', 'rRNA B', 'rRNA A', 'rRNA C', 'rRNA E']

	    if i==8 :
		axarr[i].set_xticks(ticks1, minor=False)
		axarr[i].set_xticklabels(xticknames1, rotation=0, fontsize=5)
		#plt.setp(axarr[i].set_xticklabels(xticknames1), rotation=0, fontsize=10)
		axarr[i].set_xticks(ticks2, minor=True)
		axarr[i].set_xticklabels(xticknames2, minor=True)
		plt.setp(axarr[i].set_xticklabels(xticknames2, minor=True), rotation=90, fontsize=5)
		axarr[i].tick_params(direction='out', pad=10, which='minor')
		axarr[i].xaxis.grid(True, which='minor', linewidth=1, linestyle='--', alpha=0.4)
		axarr[i].tick_params(axis='both', which='major', labelsize=5)
		axarr[i].tick_params(axis='both', which='major', labelsize=5)
	    else:
		axarr[i].axes.get_xaxis().set_visible(False)
		axarr[i].set_xticks(ticks2, minor=True)
		axarr[i].xaxis.grid(True, which='minor', linewidth=1, linestyle='--', alpha=0.4)	
		axarr[i].tick_params(axis='both', which='major', labelsize=5)
		
	    
            #axarr[i].set_title(sample, fontsize=9)
	    axarr[i].text(.5,.8,sample,
		    horizontalalignment='center',
		transform=axarr[i].transAxes,fontsize=5)	    

            

            xori_ter=0
            yter=1400000
            yori=3711828
            axarr[i].hist(peaks, bins, alpha=0.3)
            axarr[i].plot(yter, xori_ter, 'ro', markersize=5)
            axarr[i].plot(yori, xori_ter, 'go', markersize=5)
	    i=i+1 
       # plt.tight_layout()
	f.savefig(path+"Peak_distribution.png", dpi=600)        
	peaks_coords.close()

def peak_analysis_check_occucrence_in_samples():
	peaks_coords=open(path+ "List_of_peaks.txt", 'r')
        #peaks_coords=open(path+ "List_of_peaks.txt", 'r')
	outfile= open(path+ "List_of_peaks_score.txt".decode("utf-8"), 'w+')	
	list_of_peaks={}
	ordered_list_of_samples=[]

	for line in peaks_coords:
		peaks=[]
		line=line.rstrip()
		line=line.split('\t')
		sample= str(line[0])
		#print sample
		for peak in line[1:len(line)]:
		    peak=peak.split('-')
		    peaks.append(int(peak[0]))
		    list_of_peaks[sample]=peaks
		ordered_list_of_samples.append(sample)
	all_peaks=[]
	score={}
	total_score=0
	

	for peak_file in list_of_peaks:
	    print peak_file
	    for peak in list_of_peaks[peak_file]:
		if peak not in all_peaks:
			all_peaks.append(peak)

	for peak in all_peaks:
	    tmp_score={}
            for peak_file in list_of_peaks:
		
		if peak in list_of_peaks[peak_file]:
			tmp_score[peak_file]=1			
			score[peak]=tmp_score
		else:
			tmp_score[peak_file]=0			
			score[peak]=tmp_score			
	
	outfile.write("peak")		
	for i in range(0, len(ordered_list_of_samples)):
		outfile.write("\t" +ordered_list_of_samples[i])
	outfile.write("\n")
	#print len(ordered_list_of_samples)
	for peak in score:
		#print score[peak]
		total_score=0
		outfile.write(str(peak))
		for i in range(0, len(ordered_list_of_samples)):
			total_score=total_score+score[peak][ordered_list_of_samples[i]]
			outfile.write("\t" +str(score[peak][ordered_list_of_samples[i]]))		
		outfile.write("\t" + str(total_score)+ "\n")
	outfile.close()
	
def peaks_analysis2_create_files_common_peaks():
    peaks_coords=open(path+ "List_of_peaks.txt", 'r')
    '''half ready'''
    all_unique_peaks=[]
    all_unique_peaks_with_height={}
    for line in peaks_coords:
	peaks = {}
	line=line.rstrip()
	line = line.split('\t')
	sample = str(line[0])
	print sample
	for peak in line[1:len(line)]:
	    peak = peak.split('-')	    
	    peaks[int(peak[0])]=float((peak[1]))
		
	if "Microcin_IP_Mu_1" in sample:
	    peaks_micro_mu_1 = peaks
	if "Microcin_IP_Mu_3" in sample:
	    peaks_micro_mu_3 = peaks
	if "Microcin_IP_50mkm_4" in sample:
	    peaks_micro_4 = peaks	
	if "Cfx_IP_Mu_10_mkM_1" in sample:
	    peaks_cfx_mu_1 = peaks		
	if "Cfx_IP_Mu_10_mkM_3" in sample:
	    peaks_cfx_mu_3 = peaks	
	if "Cfx_IP_10mkm_3" in sample:
	    peaks_cfx_3 = peaks		    
	if "Oxo_IP_Mu_120_mkM_1" in sample:
	    peaks_oxo_mu_1 = peaks		
	if "Oxo_IP_120_mkM_1" in sample:
	    peaks_oxo_1 = peaks	
	if "Oxo_IP_120_mkM_2" in sample:
	    peaks_oxo_2 = peaks		
	    
	#print len(peaks)
    
    common_micro={}
    common_cfx={}
    common_oxo={}
    common_unique={}
    
    for p in peaks_micro_mu_1:
	if p in peaks_micro_mu_3:
	    common_micro[p]=np.mean([peaks_micro_mu_1[p],peaks_micro_mu_3[p]])
	    
    for p in peaks_cfx_mu_1:
	if p in peaks_cfx_mu_3:
	    common_cfx[p]=np.mean([peaks_cfx_mu_1[p],peaks_cfx_mu_3[p]])

    common_unique= dict(common_micro)    
    
    for p in common_cfx:
	if p not in common_unique:
	    common_unique[p]=common_cfx[p]
	else:
	    common_unique[p]=np.mean([common_micro[p],common_cfx[p]])
    print 'len common micro ' + str(len(common_micro))
    print 'len common cfx ' + str(len(common_cfx))
    print 'len common unique ' + str(len(common_unique))

    outfile= open(path+"List_of_peaks_MICROCIN_CFX_OXO.txt".decode("utf-8"), 'w+')
    
    outfile.write("MICROCIN")
    for peak in sorted(common_micro.iterkeys()):
	outfile.write('\t' + str(peak) + '-' + str(common_micro[peak]) )
    outfile.write("\nCFX")
    for peak in sorted(common_cfx.iterkeys()):
	outfile.write('\t' + str(peak) + '-' + str(common_cfx[peak] )    	)
    outfile.close()
    

def curve_fitting(wigfile):
    print wigfile
    pars_com = Parsers()
    wigin = open(path+wigfile, 'r')
    ends=pars_com.Wig_pars_return_normalized(wigin)
    fig=plt.figure(figsize=(16, 8), dpi=100)

    y_coords=[]
    mean = 0.0
    window=100000
    window_float=window*1.0

    mean=0.0
    for i in range(-1 * window, window):
        mean = mean + get_value(i,ends)
    mean=mean/(2*window_float)
    y_coords.append(mean)
    last=len(ends)

    for i in range(1, len(ends)):
	if (i >= d1_b and i <= d1_e) or (i >= d2_b and i <= d2_e) or (i >= d3_b and i <= d3_e) :
	    mean = mean + (get_value(i + window, ends) - get_value(i - window, ends)) / (2 * window_float)
	    y_coords.append(mean)
	else:
	    mean=mean + (get_value(i+window,ends) - get_value(i-window,ends))/(2*window_float)
	    y_coords.append(mean)

    return y_coords, ends

def plot_averaged_graph_one_sample (ends, line_after_division,line_un, wigfile): 

    mask_array=[]
    for k in range(0, len(ends),1):
        if (k >= d1_b and k <= d1_e):
            mask_array.append(True)
        elif (k >= d2_b and k <= d2_e):
            mask_array.append(True)
        elif (k >= d3_b and k <= d3_e):
            mask_array.append(True)
        else:
            mask_array.append(False)
    mc = np.ma.masked_array(return_averaged_line(ends), mask=mask_array)
    mc_after_difision=np.ma.masked_array(return_averaged_line(line_after_division), mask=mask_array)
    mc_un =np.ma.masked_array(return_averaged_line(line_un), mask=mask_array)

    plt.figure(figsize=(16, 8), dpi=100)
    plt.plot(xcoord,mc, '.', label='initial averaged 200kb', color='blue')
    plt.plot(xcoord,mc_after_difision, '.', label='after division', color='orange')	
    plt.plot(xcoord,mc_un, '-', label='UN line', color='green')
 
    plt.title(wigfile[:-4] + "_200kb")
    plt.legend()
    plt.savefig(path + wigfile[:-4] + '_200kb.png', dpi=100, figsize=(16, 8))
    plt.close()

def plot_averaged_graph_two_samples (ends_exp, line_after_division_exp, ends_contr, line_after_division_contr, line_after_division_exp_contr, name ):	
    mask_array=[]
    for k in range(0, len(ends_exp),1):
        if (k >= d1_b and k <= d1_e):
            mask_array.append(True)
        elif (k >= d2_b and k <= d2_e):
            mask_array.append(True)
        elif (k >= d3_b and k <= d3_e):
            mask_array.append(True)
        else:
            mask_array.append(False)

    mc = np.ma.masked_array(return_averaged_line(ends_exp), mask=mask_array)
    mc_after_division_exp=np.ma.masked_array(return_averaged_line(line_after_division_exp), mask=mask_array)
    mc_after_division_contr=np.ma.masked_array(return_averaged_line(line_after_division_contr), mask=mask_array)
    mc_after_division_exp_contr=np.ma.masked_array(return_averaged_line(line_after_division_exp_contr), mask=mask_array)
    mc_after_division_exp_contr_not_av=np.ma.masked_array(line_after_division_exp_contr, mask=mask_array)

    plt.figure(figsize=(16, 8), dpi=100)
    plt.plot(xcoord,mc, '.', label='initial averaged 200kb', color='blue')
    plt.plot(xcoord,mc_after_division_exp, '.', label='after division exp', color='orange')	
    plt.plot(xcoord,mc_after_division_contr, '.', label='after division contr' , color='green')	
    plt.plot(xcoord,mc_after_division_exp_contr, '.', label='after division exp by contr', color='red')	
    #plt.plot(xcoord,mc_after_division_exp_contr_not_av, '.', label='after division exp by contr not av', color='grey')	
    plt.title(name)
    plt.legend()
    plt.savefig(path + name , dpi=100, figsize=(16, 8))
    plt.close()
    #genome.close()

def return_AU_stat(x):
    AU_test=[5,7,9,11,12,14,16,17,19,20,22,23,24,26,27,28,30,31,32,34,35]
    
    AU_test20=20*1.75
    AU_test25=25*1.64
    AU_test30=30*1.60
    AU_test40=40*1.50
    AU_test50=50*1.44
    AU_test75=75*1.36
    AU_test100=100*1.30

    if len(AU_test)>x:
 	quantile=AU_test[int(x)]
    elif x>=20 and x<25:
	quantile=AU_test20
    elif x>=25 and x<30:
	quantile=AU_test25
    elif x>=30 and x<40:
	quantile=AU_test30
    elif x>=40 and x<50:
	quantile=AU_test40
    elif x>=50 and x<75:
	quantile=AU_test50
    elif x>=75 and x<100:
	quantile=AU_test75
    else:
	quantile=AU_test100
    return quantile

def create_gc_plot_fitting_AU(fname, peaks_coord, ends_exp,ends_cont):
    ends=ends_exp
    print fname
    # Peaks calling.
    # threshold - threshold height of reads ends peak upon specific position
    # win_width - width of area of interest under peak centre
    # Returns array "peaks", that contains coordinates of the first position of a gap

    win_width = 110
    peaks = []
    coord_list = ""
    aver = 0

   
    for i in range(len(ends) - 1 - 5):
        if ends[i] > return_AU_stat(ends_cont[i]):
            if ends[i + 5] >return_AU_stat(ends_cont[i+5]):
                #coord_list = coord_list + "\t" + str(i + 1) + "-" + str((ends[i] + ends[i + 5]) / 2)
		coord_list = coord_list + "\t" + str(i + 1) + "-" + str(min(ends[i] ,ends[i + 5]))
                peak = []
                peak.append(i + 1)
                peak.append(i + 1 - win_width + 1)
                peak.append(i + 1 + win_width + 1)
                peaks.append(peak)

    print 'len(peaks)=' + str(len(peaks))

    peaks_coord[fname] = coord_list[1:]

    if (len(peaks) > 25000):
        return

    
    # Open fasta file with genome
    # genome = open("D:/Bio_Lin_Share/Gyrase/Replic_1_Mu_SGS/Cfx_seq_data_G/E_coli_w3110_G_Mu.fasta", 'r')
    #genome = open(u"/data/Gyrase/IGV_tracks/E_coli_w3110_G_Mu.fasta", 'r')
    #for record in SeqIO.parse(genome, "fasta"):
    #    genomefa = record.seq

    # Returns "seqs" array contains sequences under peaks with coordinates from "peaks"
    k = 0
    seqs = []
    for i in range(len(peaks)):
        seq = genomefa[int(peaks[i][1]):int(peaks[i][2])]
        seqs.append(seq)
        k = k + 1
    print 'len(seqs)=' + str(len(seqs))

    if len(seqs) == 0:
        return

    create_motif_plot(fname,seqs)
    
    instances=[]
    for seq in seqs:
	instances.append(Seq(str(seq)))
    #return instances,seqs
    return instances,seqs
    

def create_motif_plot(fname,seqs):
    ############################

    ############################

    # PWM construction
    # Scans sequences stack by columns, counts the number of particular letters
    # Returns pwm - "matrix"
    matrix = []
    template = seqs[0]
    #if len(seqs[0]) < 220:
    #    template = seqs[10]

    # seqs[0]
    for i in range(len(template)):
        column = [0, 0, 0, 0]
        for j in range(len(seqs)):
            #if len(seqs[j]) < 220:
            #    continue
            if seqs[j][i] == str('A'):
                column[0] = column[0] + 1
            elif seqs[j][i] == str('T'):
                column[1] = column[1] + 1
            elif seqs[j][i] == str('G'):
                column[2] = column[2] + 1
            elif seqs[j][i] == str('C'):
                column[3] = column[3] + 1

        # ATGC_pwm.write(str(float(column[0])/(column[0]+column[1]+column[2]+column[3]))+'\t'+str(float(column[1])/(column[0]+column[1]+column[2]+column[3]))+'\t'+str(float(column[2])/(column[0]+column[1]+column[2]+column[3]))+'\t'+str(float(column[3])/(column[0]+column[1]+column[2]+column[3]))+'\n')
        matrix.append(column)

    print 'len(matrix)=' + str(len(matrix))

    # Counts different combinations of nucleotides, only GC is used subsequently by default
    # Returns reductive pwms
    GC_percent = []
    GA_percent = []
    GT_percent = []
    AT_percent = []
    CT_percent = []
    A_percent = []
    T_percent = []
    G_percent = []
    C_percent = []

    for i in range(len(matrix)):
        GC = float((int(matrix[i][2]) + int(matrix[i][3]))) / (
            int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
        GT = float((int(matrix[i][1]) + int(matrix[i][2]))) / (
            int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
        AT = float((int(matrix[i][0]) + int(matrix[i][1]))) / (
            int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
        GA = float((int(matrix[i][0]) + int(matrix[i][2]))) / (
            int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
        CT = float((int(matrix[i][1]) + int(matrix[i][3]))) / (
            int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
        A = float((int(matrix[i][0]))) / (int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
        T = float((int(matrix[i][1]))) / (int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
        G = float((int(matrix[i][2]))) / (int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
        C = float((int(matrix[i][3]))) / (int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
        GC_percent.append(GC)
        GT_percent.append(GT)
        AT_percent.append(AT)
        CT_percent.append(CT)
        GA_percent.append(GA)
        A_percent.append(A)
        T_percent.append(T)
        G_percent.append(G)
        C_percent.append(C)

    '''   
    print GC_percent
    print AT_percent
    print A_percent
    print T_percent
    print G_percent
    print C_percent
    '''

    ############################

    # GC statistics module
    # Counts average GC% over the whole genome
    GC_genome = GC_count(genomefa) / 100
    print 'GC_genome=' + str(GC_genome)

    # Counts GC% pValue in the particular pwm column
    # Returns pValue array and auxiliary Zero array for plotting
    alignment_thick = len(seqs)
    pValue = []
    Zero = []
    for i in range(len(GC_percent)):
        pValue.append(scipy.stats.binom_test(float(GC_percent[i]) * alignment_thick, n=alignment_thick, p=GC_genome))
        Zero.append(1)

    ############################

    # Plotting
    x_axis = []
    for i in range(len(GC_percent)):
        x_axis.append(-109 + i)
	#x_axis.append(-65 + i)
    print 'len(x_axis)=' + str(len(x_axis))

    ax_range = [-110, +110, 0.2, 1]

    # fig =
    plt.figure(figsize=(16, 8), dpi=100)
    # GC% pwm plotting
    #plt.suptitle(wigfile[:-4].split("\\")[len(wigfile[:-4].split("\\")) - 1], fontsize=20)
    plt.suptitle(fname[:-4], fontsize=20)
    plot1 = plt.subplot()
    plot1.plot(x_axis, GC_percent, color='green', linewidth=1)
    plot1.plot(x_axis, GC_percent, 'go', markersize=3)
    plot1.set_xticks(np.arange(-120, 112, 10))
    plot1.axis([-110, 110, 0.3, 1])
    plot1.set_xlim(-110, 110)
    plot1.annotate('GC%', xytext=(85, 0.68), xy=(40, 0.85), color='green', weight="bold", size=15)
    plot1.annotate('p-Value', xytext=(85, 0.63), xy=(-105, 0.64), color='cyan', weight="bold", size=15)
    plot1.set_xlabel('Position, nt', size=17)
    plot1.set_ylabel('GC%', size=17)
    # pValue plotting
    plot2 = plot1.twinx()
    plot2.plot(x_axis, pValue, 'k', linewidth=0.5, alpha=0.6)
    plot2.fill_between(x_axis, pValue, Zero, color='cyan', alpha=0.2)
    plot2.set_yticks(np.arange(0, 1.01, 0.01), minor=False)
    plot2.set_yscale('log')
    plot2.set_yticks([0.005], minor=True)
    plot2.yaxis.grid(True, which='minor', linewidth=2, linestyle='--', alpha=0.3)
    plot2.annotate('Confidence level = 0.005', xytext=(66, 0.0032), xy=(40, 0.8), color='black', size=15)
    plot2.set_ylim(0.0000001, 1.0)
    plot2.set_xlim(-110, 110)
    plot2.set_ylabel('p-Value, logarithmic scale', size=17)

    # plt.show()
    # plt.savefig("D:/Bio_Lin_Share/Gyrase/Replic_1_Mu_SGS/Cfx_seq_data_G/Over_90_peaks_motif.png", dpi=320, figsize=(80, 40))

    plt.savefig(path+fname, dpi=100, figsize=(16, 8))
    plt.close()
    genome.close()
    
    
def create_motif_plot_from_peaks(fname,peaks):
    ############################
    seqs = []
    for i in range(len(peaks)):
	seq = genomefa[int(peaks[i]- win_width + 1):int(peaks[i]+ win_width + 1)]
	seqs.append(seq)
    print 'len(seqs)=' + str(len(seqs))

    if len(seqs) == 0:
	return
    ############################

    ############################

    # PWM construction
    # Scans sequences stack by columns, counts the number of particular letters
    # Returns pwm - "matrix"
    matrix = []
    template = seqs[0]
    #if len(seqs[0]) < 220:
    #    template = seqs[10]

    # seqs[0]
    for i in range(len(template)):
        column = [0, 0, 0, 0]
        for j in range(len(seqs)):
            #if len(seqs[j]) < 220:
            #    continue
            if seqs[j][i] == str('A'):
                column[0] = column[0] + 1
            elif seqs[j][i] == str('T'):
                column[1] = column[1] + 1
            elif seqs[j][i] == str('G'):
                column[2] = column[2] + 1
            elif seqs[j][i] == str('C'):
                column[3] = column[3] + 1

        # ATGC_pwm.write(str(float(column[0])/(column[0]+column[1]+column[2]+column[3]))+'\t'+str(float(column[1])/(column[0]+column[1]+column[2]+column[3]))+'\t'+str(float(column[2])/(column[0]+column[1]+column[2]+column[3]))+'\t'+str(float(column[3])/(column[0]+column[1]+column[2]+column[3]))+'\n')
        matrix.append(column)

    print 'len(matrix)=' + str(len(matrix))

    # Counts different combinations of nucleotides, only GC is used subsequently by default
    # Returns reductive pwms
    GC_percent = []
    GA_percent = []
    GT_percent = []
    AT_percent = []
    CT_percent = []
    A_percent = []
    T_percent = []
    G_percent = []
    C_percent = []

    for i in range(len(matrix)):
        GC = float((int(matrix[i][2]) + int(matrix[i][3]))) / (
            int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
        GT = float((int(matrix[i][1]) + int(matrix[i][2]))) / (
            int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
        AT = float((int(matrix[i][0]) + int(matrix[i][1]))) / (
            int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
        GA = float((int(matrix[i][0]) + int(matrix[i][2]))) / (
            int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
        CT = float((int(matrix[i][1]) + int(matrix[i][3]))) / (
            int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
        A = float((int(matrix[i][0]))) / (int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
        T = float((int(matrix[i][1]))) / (int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
        G = float((int(matrix[i][2]))) / (int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
        C = float((int(matrix[i][3]))) / (int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
        GC_percent.append(GC)
        GT_percent.append(GT)
        AT_percent.append(AT)
        CT_percent.append(CT)
        GA_percent.append(GA)
        A_percent.append(A)
        T_percent.append(T)
        G_percent.append(G)
        C_percent.append(C)

    '''   
    print GC_percent
    print AT_percent
    print A_percent
    print T_percent
    print G_percent
    print C_percent
    '''

    ############################

    # GC statistics module
    # Counts average GC% over the whole genome
    GC_genome = GC_count(genomefa) / 100
    print 'GC_genome=' + str(GC_genome)

    # Counts GC% pValue in the particular pwm column
    # Returns pValue array and auxiliary Zero array for plotting
    alignment_thick = len(seqs)
    pValue = []
    Zero = []
    for i in range(len(GC_percent)):
        pValue.append(scipy.stats.binom_test(float(GC_percent[i]) * alignment_thick, n=alignment_thick, p=GC_genome))
        Zero.append(1)

    ############################

    # Plotting
    x_axis = []
    for i in range(len(GC_percent)):
        x_axis.append(-109 + i)
	#x_axis.append(-65 + i)
    print 'len(x_axis)=' + str(len(x_axis))

    ax_range = [-110, +110, 0.2, 1]

    # fig =
    plt.figure(figsize=(16, 8), dpi=100)
    # GC% pwm plotting
    #plt.suptitle(wigfile[:-4].split("\\")[len(wigfile[:-4].split("\\")) - 1], fontsize=20)
    plt.suptitle(fname[:-4], fontsize=20)
    plot1 = plt.subplot()
    plot1.plot(x_axis, GC_percent, color='green', linewidth=1)
    plot1.plot(x_axis, GC_percent, 'go', markersize=3)
    plot1.set_xticks(np.arange(-120, 112, 10))
    plot1.axis([-110, 110, 0.3, 1])
    plot1.set_xlim(-110, 110)
    plot1.annotate('GC%', xytext=(85, 0.68), xy=(40, 0.85), color='green', weight="bold", size=15)
    plot1.annotate('p-Value', xytext=(85, 0.63), xy=(-105, 0.64), color='cyan', weight="bold", size=15)
    plot1.set_xlabel('Position, nt', size=17)
    plot1.set_ylabel('GC%', size=17)
    # pValue plotting
    plot2 = plot1.twinx()
    plot2.plot(x_axis, pValue, 'k', linewidth=0.5, alpha=0.6)
    plot2.fill_between(x_axis, pValue, Zero, color='cyan', alpha=0.2)
    plot2.set_yticks(np.arange(0, 1.01, 0.01), minor=False)
    plot2.set_yscale('log')
    plot2.set_yticks([0.005], minor=True)
    plot2.yaxis.grid(True, which='minor', linewidth=2, linestyle='--', alpha=0.3)
    plot2.annotate('Confidence level = 0.005', xytext=(66, 0.0032), xy=(40, 0.8), color='black', size=15)
    plot2.set_ylim(0.0000001, 1.0)
    plot2.set_xlim(-110, 110)
    plot2.set_ylabel('p-Value, logarithmic scale', size=17)

    # plt.show()
    # plt.savefig("D:/Bio_Lin_Share/Gyrase/Replic_1_Mu_SGS/Cfx_seq_data_G/Over_90_peaks_motif.png", dpi=320, figsize=(80, 40))

    plt.savefig(path+fname, dpi=100, figsize=(16, 8))
    plt.close()
    genome.close()
    return GC_percent,len(peaks),GC_genome,seqs

def september_create_motif_pictures_and_files():
    peaks_coords=open(path+ "List_of_peaks_MICROCIN_CFX_OXO_2.txt", 'r')
    
    for line in peaks_coords:
	peaks=[]
	line=line.split('\t')
	sample=str(line[0])
	for peak in line[1:len(line)]:
	    peak = peak.split('-')
	    peaks.append(int(peak[0]))
	data=create_motif_plot_from_peaks(sample,peaks)
	
	outfile= open(path+"" + sample +"_data.txt".decode("utf-8"), 'w+')    
	outfile.write(sample)
	for GC in data[0]:
	    outfile.write('\t' + str(GC) )
	outfile.write("\n")
	outfile.write(str(data[1]) + "\n")
	outfile.write(str(data[2]) + "\n")    
	outfile.close()
	
	outfile= open(path+"" + sample +"_seqs.txt".decode("utf-8"), 'w+')    
	outfile.write(">_" + sample+ "\n")
	i=0
	for seq in data[3]:	    
	    i+=1
	    outfile.write(">_"+str(i) + "\n"+ str(seq) + "\n" )	  
	outfile.close()	
    
    peaks_coords.close()
	
	    
def create_motif_plot_two_samples(fname,seqs1,seqs2,seqs3):
    ############################
    ############################

    # PWM construction
    # Scans sequences stack by columns, counts the number of particular letters
    # Returns pwm - "matrix"
    both_seqs=[seqs1,seqs2,seqs3]
    clr=["green","blue","red"]
    label=["Micro","Cfx", "Oxo"]
    cl=0	
    fig=plt.figure(figsize=(16, 8), dpi=100)
    for seqs in both_seqs: 	
	    matrix = []
	    template = seqs[0]
	   # if len(seqs[0]) < 220:
	   #	template = seqs[10]

	    # seqs[0]
	    for i in range(len(template)):
		column = [0, 0, 0, 0]
		for j in range(len(seqs)):
		    #if len(seqs[j]) < 220:
		    #    continue
		    if seqs[j][i] == str('A'):
		        column[0] = column[0] + 1
		    elif seqs[j][i] == str('T'):
		        column[1] = column[1] + 1
		    elif seqs[j][i] == str('G'):
		        column[2] = column[2] + 1
		    elif seqs[j][i] == str('C'):
		        column[3] = column[3] + 1

		# ATGC_pwm.write(str(float(column[0])/(column[0]+column[1]+column[2]+column[3]))+'\t'+str(float(column[1])/(column[0]+column[1]+column[2]+column[3]))+'\t'+str(float(column[2])/(column[0]+column[1]+column[2]+column[3]))+'\t'+str(float(column[3])/(column[0]+column[1]+column[2]+column[3]))+'\n')
		matrix.append(column)

	    print 'len(matrix)=' + str(len(matrix))

	    # Counts different combinations of nucleotides, only GC is used subsequently by default
	    # Returns reductive pwms
	    GC_percent = []
	    GA_percent = []
	    GT_percent = []
	    AT_percent = []
	    CT_percent = []
	    A_percent = []
	    T_percent = []
	    G_percent = []
	    C_percent = []

	    for i in range(len(matrix)):
		GC = float((int(matrix[i][2]) + int(matrix[i][3]))) / (
		    int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
		GT = float((int(matrix[i][1]) + int(matrix[i][2]))) / (
		    int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
		AT = float((int(matrix[i][0]) + int(matrix[i][1]))) / (
		    int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
		GA = float((int(matrix[i][0]) + int(matrix[i][2]))) / (
		    int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
		CT = float((int(matrix[i][1]) + int(matrix[i][3]))) / (
		    int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
		A = float((int(matrix[i][0]))) / (int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
		T = float((int(matrix[i][1]))) / (int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
		G = float((int(matrix[i][2]))) / (int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
		C = float((int(matrix[i][3]))) / (int(matrix[i][0]) + int(matrix[i][1]) + int(matrix[i][2]) + int(matrix[i][3]))
		GC_percent.append(GC)
		GT_percent.append(GT)
		AT_percent.append(AT)
		CT_percent.append(CT)
		GA_percent.append(GA)
		A_percent.append(A)
		T_percent.append(T)
		G_percent.append(G)
		C_percent.append(C)

	    ############################

	    # GC statistics module
	    # Counts average GC% over the whole genome
	    GC_genome = GC_count(genomefa) / 100
	    print 'GC_genome=' + str(GC_genome)

	    # Counts GC% pValue in the particular pwm column
	    # Returns pValue array and auxiliary Zero array for plotting
	    alignment_thick = len(seqs)
	    pValue = []
	    Zero = []
	    for i in range(len(GC_percent)):
		pValue.append(scipy.stats.binom_test(float(GC_percent[i]) * alignment_thick, n=alignment_thick, p=GC_genome))
		Zero.append(1)

	    ############################

	    # Plotting
	    x_axis = []
	    for i in range(len(GC_percent)):
		#x_axis.append(-65 + i)
		x_axis.append(-109 + i)
	    print 'len(x_axis)=' + str(len(x_axis))

	    ax_range = [-110, +110, 0.2, 1]

	    # fig =

	    # GC% pwm plotting
	    #plt.suptitle(wigfile[:-4].split("\\")[len(wigfile[:-4].split("\\")) - 1], fontsize=20)
	    plt.suptitle(fname[:-4], fontsize=20)
	    plot1 = plt.subplot()
	    plot1.plot(x_axis, GC_percent,  label=label[cl], color=clr[cl], linewidth=1)
	    #plot1.plot(x_axis, GC_percent, 'go', markersize=3)
	    ax = fig.add_subplot(111)
    	    #for xy in zip(x_axis, GC_percent):                                       # <--
    		#ax.annotate('(%s,  %.1f)'  % xy, xy=xy, textcoords='data', fontsize=5)
		#ax.annotate('%s'  % xy[0], xy=xy, textcoords='data', fontsize=5)

	    plot1.set_xticks(np.arange(-120, 112, 10))
	    plot1.axis([-110, 110, 0.3, 1])
	    plot1.set_xlim(-110, 110)
	    plot1.annotate('GC%', xytext=(85, 0.68), xy=(40, 0.85), color='green', weight="bold", size=15)

	    #plot1.annotate('p-Value', xytext=(85, 0.63), xy=(-105, 0.64), color='cyan', weight="bold", size=15)
	    plot1.set_xlabel('Position, nt', size=17)
	    plot1.set_ylabel('GC%', size=17)
	    # pValue plotting
	    cl=cl+1
	    '''
	    plot2 = plot1.twinx()
	    plot2.plot(x_axis, pValue, 'k', linewidth=0.5, alpha=0.6)
	    plot2.fill_between(x_axis, pValue, Zero, color='cyan', alpha=0.2)
	    plot2.set_yticks(np.arange(0, 1.01, 0.01), minor=False)
	    plot2.set_yscale('log')
	    plot2.set_yticks([0.005], minor=True)
	    plot2.yaxis.grid(True, which='minor', linewidth=2, linestyle='--', alpha=0.3)
	    plot2.annotate('Confidence level = 0.005', xytext=(66, 0.0032), xy=(40, 0.8), color='black', size=15)
	    plot2.set_ylim(0.0000001, 1.0)
	    plot2.set_xlim(-110, 110)
	    plot2.set_ylabel('p-Value, logarithmic scale', size=17)
	    
	    '''
    plt.legend()
    plt.savefig(path+fname, dpi=100, figsize=(16, 8))
    genome.close()
    plt.close()

def fitting_divide_case_7_Mitya(ex_file,cont_file, out_file):
    
    control= curve_fitting(cont_file)
    experim=curve_fitting(ex_file)

    

    exper_ends = [1.0 * x +1.0  for x in experim[1]]
    control_ends= [1.0 * x+1.0  for x in control[0]]
    #un_exper_ends = [1.0 * x+1.0  for x in un_experim[0]]
    #un_control_ends= [1.0 * x +1.0 for x in un_control[0]]
	
 
    ends_divide=[]
    count_zero=0
    count_zero_exp=0
    count_zero_control=0

    for i in range (0, len(exper_ends)):
	if exper_ends[i]!=0 and control_ends[i]!=0:
		ends_divide.append(exper_ends[i]/ control_ends[i])
	else:
		ends_divide.append(0)
    

   # ends_divide_experim=np.divide(exper_ends,un_exper_ends)
   # ends_divide_control=np.divide(control_ends,un_control_ends)

    plot_averaged_graph_one_sample(exper_ends, ends_divide, control_ends, ex_file)

    return ends_divide

def case_7_Mitya():
    
	outfile_list_of_peaks= open(path + "List_of_peaks.txt".decode("utf-8"), 'w+')
    
	peaks_coord = {}
	control=[]
	seqs=[]
	all_microcin_seqs=[]
	all_cfx_seqs=[]
	all_oxo_seqs=[]
    
    
	un_cont_mu_file="Un_IN_Mu_1_ends.wig"
	un_ex_mu_file="Un_IP_Mu_1_ends.wig"
	out_file="Un_IP_Mu_1_vs_Un_IN_Mu_1"
	un_mu_ends_fitting=fitting_divide_case_7_Mitya(un_ex_mu_file,un_cont_mu_file,out_file)

	un_cont_mu_new_file="Un_IN_Mu_1_new_ends.wig"
	un_ex_mu_new_file="Un_IP_Mu_1_new_ends.wig"
	out_file="Un_IP_Mu_1_new_vs_Un_IN_Mu_1_new"
	un_mu_new_ends_fitting=fitting_divide_case_7_Mitya(un_ex_mu_new_file,un_cont_mu_new_file,out_file)	
	
	out_file="Un_IP_Mu_1_new_vs_Un_IP_Mu_1"	
	seqs=create_gc_plot_fitting_AU(out_file+"_case_7_Mitya.png", peaks_coord,un_mu_new_ends_fitting,un_mu_ends_fitting)
	
	out_file="Un_IP_Mu_1_vs_Un_IP_Mu_1_new"
	seqs=create_gc_plot_fitting_AU(out_file+"_case_7_Mitya.png", peaks_coord,un_mu_ends_fitting,un_mu_new_ends_fitting)
		
	un_cont_file="Un_IN_1_ends.wig"
	un_ex_file="Un_IN_2_ends.wig"
	out_file="Un_IN_2_vs_Un_IN_1"
	un_ends_fitting=fitting_divide_case_7_Mitya(un_ex_file,un_cont_file,out_file)
    
    
	''' OXO '''
    
	'''
	cont_file="Cfx_IN_Mu_10_mkM_3_ends.wig"
	ex_file="Oxo_IP_Mu_120_mkM_1_ends.wig"
	out_file="Oxo_IP_Mu_120_mkM_1_vs_Un_IP_Mu_1_new"
	ends_fitting=fitting_divide_case_7_Mitya(ex_file,cont_file,out_file)    
	seqs=create_gc_plot_fitting_AU(out_file+"_case_7_Mitya.png", peaks_coord,ends_fitting,un_mu_new_ends_fitting)
	all_oxo_seqs=all_oxo_seqs + seqs[1]	
    
    
    
	cont_file="Cfx_IN_10_mkM_3_ends.wig"
	ex_file="Oxo_IP_120_mkM_1_ends.wig"
	out_file="Oxo_IP_120_mkM_1_vs_Un_1"
	ends_fitting=fitting_divide_case_7_Mitya(ex_file,cont_file,out_file)    
	seqs=create_gc_plot_fitting_AU(out_file+"_case_7_Mitya.png", peaks_coord,ends_fitting,un_ends_fitting)
	all_oxo_seqs=all_oxo_seqs + seqs[1]
    
    
	cont_file="Cfx_IN_10_mkM_3_ends.wig"
	ex_file="Oxo_IP_120_mkM_2_ends.wig"
	out_file="Oxo_IP_120_mkM_2_vs_Un_1"
	ends_fitting=fitting_divide_case_7_Mitya(ex_file,cont_file,out_file)    
	seqs=create_gc_plot_fitting_AU(out_file+"_case_7_Mitya.png", peaks_coord,ends_fitting,un_ends_fitting)	
	all_oxo_seqs=all_oxo_seqs + seqs[1]		
    
	print 'len all oxo seq ' + str(len(all_oxo_seqs))
	create_motif_plot("all_oxo.png",all_oxo_seqs)
	'''
	
	''' MICRO '''
	'''
	cont_file="Microcin_IN_Mu_1_ends.wig"
	ex_file="Microcin_IP_Mu_1_ends.wig"
	out_file="Microcin_IP_Mu_1_vs_Un_Mu_1"
	ends_fitting=fitting_divide_case_7_Mitya(ex_file,cont_file,out_file)	
	seqs=create_gc_plot_fitting_AU(out_file+"_case_7_Mitya.png", peaks_coord,ends_fitting,un_mu_ends_fitting)
	all_microcin_seqs=all_microcin_seqs + seqs[1]
	
	cont_file="Cfx_IN_Mu_10_mkM_3_ends.wig"
	ex_file="Microcin_IP_Mu_3_ends.wig"
	out_file="Microcin_IP_Mu_3_vs_Un_Mu_1"
	ends_fitting=fitting_divide_case_7_Mitya(ex_file,cont_file,out_file)    
	seqs=create_gc_plot_fitting_AU(out_file+"_case_7_Mitya.png", peaks_coord,ends_fitting,un_mu_ends_fitting)
	all_microcin_seqs=all_microcin_seqs + seqs[1]

	cont_file="Cfx_IN_10_mkM_3_ends.wig"
	ex_file="Microcin_IP_50_mkM_4_ends.wig"
	out_file="Microcin_IP_50_mkM_4_vs_Un_1"
	ends_fitting=fitting_divide_case_7_Mitya(ex_file,cont_file,out_file)    
	seqs=create_gc_plot_fitting_AU(out_file+"_case_7_Mitya.png", peaks_coord,ends_fitting,un_ends_fitting)
	all_microcin_seqs=all_microcin_seqs + seqs[1]

	print 'len all microcin seq ' + str(len(all_microcin_seqs))
	create_motif_plot("all_microcin.png",all_microcin_seqs)
	
	'''
	''' CFX '''
	'''
	cont_file="Cfx_IN_Mu_1_ends.wig"
	ex_file="Cfx_IP_Mu_10_mkM_1_ends.wig"
	out_file="Cfx_IP_Mu_10_mkM_1_vs_Un_Mu_1"
	ends_fitting=fitting_divide_case_7_Mitya(ex_file,cont_file,out_file)
	seqs=create_gc_plot_fitting_AU(out_file+"_case_7_Mitya.png", peaks_coord,ends_fitting,un_mu_ends_fitting)
	all_cfx_seqs=all_cfx_seqs + seqs[1]

	cont_file="Cfx_IN_Mu_10_mkM_3_ends.wig"
	ex_file="Cfx_IP_Mu_10_mkM_3_ends.wig"
	out_file="Cfx_IP_Mu_10_mkM_3_vs_Un_Mu_1"
	ends_fitting=fitting_divide_case_7_Mitya(ex_file,cont_file,out_file)    
	seqs=create_gc_plot_fitting_AU(out_file+"_case_7_Mitya.png", peaks_coord,ends_fitting,un_mu_ends_fitting)
	all_cfx_seqs=all_cfx_seqs + seqs[1]
	

	cont_file="Cfx_IN_10_mkM_3_ends.wig"
	ex_file="Cfx_IP_10_mkM_3_ends.wig"
	out_file="Cfx_IP_10_mkM_3_vs_Un_1"
	ends_fitting=fitting_divide_case_7_Mitya(ex_file,cont_file,out_file)    
	seqs=create_gc_plot_fitting_AU(out_file+"_case_7_Mitya.png", peaks_coord,ends_fitting,un_ends_fitting)	
	all_cfx_seqs=all_cfx_seqs + seqs[1]		
	'''
	
	'''september'''
	
	cont_file="Cfx_IN_Mu_10_mkM_3_ends.wig"
	ex_file="Cfx_IP_Mu_900_nM_3_ends.wig"
	out_file="Cfx_IP_Mu_900_nM_3_vs_Un_Mu_1_new"
	ends_fitting=fitting_divide_case_7_Mitya(ex_file,cont_file,out_file)    
	seqs=create_gc_plot_fitting_AU(out_file+"_case_7_Mitya.png", peaks_coord,ends_fitting,un_mu_new_ends_fitting)	
	all_cfx_seqs=all_cfx_seqs + seqs[1]
	
	cont_file="Cfx_IN_10_mkM_3_ends.wig"
	ex_file="Cfx_IP_900_nM_3_ends.wig"
	out_file="Cfx_IP_900_nM_3_vs_Un_1"
	ends_fitting=fitting_divide_case_7_Mitya(ex_file,cont_file,out_file)    
	seqs=create_gc_plot_fitting_AU(out_file+"_case_7_Mitya.png", peaks_coord,ends_fitting,un_ends_fitting)	
	all_cfx_seqs=all_cfx_seqs + seqs[1]	
	
	print 'len all cfx seq ' + str(len(all_cfx_seqs))
	create_motif_plot("all_cfx.png",all_cfx_seqs)
	

	
	#create_motif_plot_two_samples("mircocin_cfx.png",all_microcin_seqs,all_cfx_seqs)
	
	
	
	
	
	
	
	'''
	cont_file="Un_IN_Mu_1_ends.wig"
	ex_file="Un_IP_Mu_1_ends.wig"
	out_file="Un_IP_Mu_1_vs_Un_IN_Mu_1"

	ends_fitting=fitting_divide(ex_file,cont_file,out_file)
	#create_gc_plot_fitting(out_file+"av.png", peaks_coord,ends_fitting[2],ends_fitting[3])
	create_gc_plot_fitting(out_file+"av.png", peaks_coord,ends_fitting[0],ends_fitting[3])
	create_gc_plot_fitting(out_file+"fe.png", peaks_coord,ends_fitting[1],ends_fitting[3])
	'''

	'''
	cont_file="Cfx_IN_10_mkM_3_ends.wig"
	ex_file="Cfx_IP_10_mkM_3_ends.wig"
	out_file="Cfx_IP_10_mkM_3_vs_Cfx_IN_10_mkM_3"

	ends_fitting=fitting_divide(ex_file,cont_file,un_ex_file,un_cont_file,out_file)
	#create_gc_plot_fitting(out_file+"av.png", peaks_coord,ends_fitting[2],ends_fitting[3])
	#create_gc_plot_fitting(out_file+"av.png", peaks_coord,ends_fitting[0])
	create_gc_plot_fitting(out_file+"fe_UN.png", peaks_coord,ends_fitting[4],ends_fitting[5])
	'''

	for peak_file in peaks_coord:
	    outfile_list_of_peaks.write(peak_file + "\t" +peaks_coord[peak_file]+"\n")
	    
	outfile_list_of_peaks.close()

def fitting_divide_case_6_AC_SW(ex_file,cont_file,un_ex_file, un_cont_file, out_file):
    
    control= curve_fitting(cont_file)
    experim=curve_fitting(ex_file)
    un_control= curve_fitting(un_cont_file)
    un_experim= curve_fitting(un_ex_file)

    

    exper_ends = [1.0 * x +1.0  for x in experim[1]]
    control_ends= [1.0 * x+1.0  for x in control[1]]
    un_exper_ends = [1.0 * x+1.0  for x in un_experim[0]]
    un_control_ends= [1.0 * x +1.0 for x in un_control[0]]
	
 
    ends_divide_experim=[]
    ends_divide_control=[]
    count_zero=0
    count_zero_exp=0
    count_zero_control=0

    for i in range (0, len(exper_ends)):
	if exper_ends[i]!=0 and un_exper_ends[i]!=0:
		ends_divide_experim.append(exper_ends[i]/ un_exper_ends[i])
	else:
		ends_divide_experim.append(0)

	if control_ends[i]!=0 and un_control_ends[i]!=0:
		ends_divide_control.append(control_ends[i]/ un_control_ends[i])
	else:
		ends_divide_control.append(0)
    

    

    ends_divide_experim=np.divide(exper_ends,un_exper_ends)
    ends_divide_control=np.divide(control_ends,un_control_ends)
    #ends_divide=np.divide(ends_divide_experim,ends_divide_control)

    plot_averaged_graph_one_sample(experim[1], ends_divide_experim, un_experim[1], ex_file)
    plot_averaged_graph_one_sample(control[1], ends_divide_control,un_control[1], cont_file)    
    #plot_averaged_graph_two_samples(experim[1], ends_divide_experim, control_ends, ends_divide_control, ends_divide, out_file)

    return ends_divide_experim, ends_divide_control

def case_6_AC_SW():
    
	peaks_coord = {}
	control=[]
	seqs=[]
	all_microcin_seqs=[]
	all_cfx_seqs=[]
	
	un_cont_mu_file="Un_IN_Mu_1_ends.wig"
	un_ex_mu_file="Un_IP_Mu_1_ends.wig"	
	un_cont_mu_new_file="Un_IN_Mu_1_new_ends.wig"
	un_ex_mu_new_file="Un_IP_Mu_1_new_ends.wig"	    
	un_cont_file="Un_IN_1_ends.wig"
	un_ex_file="Un_IN_2_ends.wig"

	cont_file="Microcin_IN_Mu_1_ends.wig"
	ex_file="Microcin_IP_Mu_1_ends.wig"
	out_file="Microcin_IP_Mu_1_vs_Microcin_IN_Mu_1"	
	ends_fitting=fitting_divide_case_6_AC_SW(ex_file,cont_file,un_ex_mu_file,un_cont_mu_file,out_file)
	seqs=create_gc_plot_fitting_AU(out_file+"_case_6_AC_SW.png", peaks_coord,ends_fitting[0],ends_fitting[1])
	all_microcin_seqs=all_microcin_seqs + seqs[1]
 
	
	cont_file="Cfx_IN_Mu_10_mkM_3_ends.wig"
	ex_file="Microcin_IP_Mu_3_ends.wig"
	out_file="Microcin_IP_Mu_3_vs_Cfx_IN_Mu_10_mkM_3"
	ends_fitting=fitting_divide_case_6_AC_SW(ex_file,cont_file,un_ex_mu_file,un_cont_mu_file,out_file)
	seqs=create_gc_plot_fitting_AU(out_file+"_case_6_AC_SW.png", peaks_coord,ends_fitting[0],ends_fitting[1])
	all_microcin_seqs=all_microcin_seqs + seqs[1]
	
	'''
	cont_file="Microcin_IN_Mu_1_ends.wig"
	ex_file="Microcin_IP_Mu_3_ends.wig"
	out_file="Microcin_IP_Mu_3_vs_Microcin_IN_Mu_1_ends"
	ends_fitting=fitting_divide_case_6_AC_SW(ex_file,cont_file,un_ex_mu_file,un_cont_mu_file,out_file)
	seqs=create_gc_plot_fitting_AU(out_file+"_case_6_AC_SW.png", peaks_coord,ends_fitting[0],ends_fitting[1])	
	all_microcin_seqs=all_microcin_seqs + seqs[1]
	
	
	cont_file="Cfx_IN_10_mkM_3_ends.wig"
	ex_file="Microcin_IP_3_ends.wig"
	out_file="Microcin_IP_3_vs_Cfx_IN_10_mkM_3"
	ends_fitting=fitting_divide_case_6_AC_SW(ex_file,cont_file,un_ex_file,un_cont_file,out_file)
	seqs=create_gc_plot_fitting_AU(out_file+"_case_6_AC_SW", peaks_coord,ends_fitting[0],ends_fitting[1])
	all_microcin_seqs=all_microcin_seqs + seqs[1]
	'''
	
	cont_file="Cfx_IN_10_mkM_3_ends.wig"
	ex_file="Microcin_IP_50_mkM_4_ends.wig"
	out_file="Microcin_IP_50_mkM_4_vs_Cfx_IN_10_mkM_3"
	ends_fitting=fitting_divide_case_6_AC_SW(ex_file,cont_file,un_ex_file,un_cont_file,out_file)
	seqs=create_gc_plot_fitting_AU(out_file+"_case_6_AC_SW.png", peaks_coord,ends_fitting[0],ends_fitting[1])
	all_microcin_seqs=all_microcin_seqs + seqs[1]
	

	print 'len all microcin seq ' + str(len(all_microcin_seqs))
	create_motif_plot("all_microcin.png",all_microcin_seqs)

	cont_file="Cfx_IN_Mu_1_ends.wig"
	ex_file="Cfx_IP_Mu_10_mkM_1_ends.wig"
	out_file="Cfx_IP_Mu_10_mkM_1_vs_Cfx_IN_Mu_1"
	ends_fitting=fitting_divide_case_6_AC_SW(ex_file,cont_file,un_ex_mu_file,un_cont_mu_file,out_file)
	seqs=create_gc_plot_fitting_AU(out_file+"_case_6_AC_SW.png", peaks_coord,ends_fitting[0],ends_fitting[1])
	all_cfx_seqs=all_cfx_seqs + seqs[1]

	cont_file="Cfx_IN_Mu_10_mkM_3_ends.wig"
	ex_file="Cfx_IP_Mu_10_mkM_3_ends.wig"
	out_file="Cfx_IP_Mu_10_mkM_3_vs_Cfx_IN_Mu_10_mkM_3"
	ends_fitting=fitting_divide_case_6_AC_SW(ex_file,cont_file,un_ex_mu_file,un_cont_mu_file,out_file)
	seqs=create_gc_plot_fitting_AU(out_file+"_case_6_AC_SW.png", peaks_coord,ends_fitting[0],ends_fitting[1])
	all_cfx_seqs=all_cfx_seqs + seqs[1]
        
	cont_file="Cfx_IN_10_mkM_3_ends.wig"
	ex_file="Cfx_IP_10_mkM_3_ends.wig"
	out_file="Cfx_IP_10_mkM_3_vs_Cfx_IN_10_mkM_3"
	ends_fitting=fitting_divide_case_6_AC_SW(ex_file,cont_file,un_ex_file,un_cont_file,out_file)
	seqs=create_gc_plot_fitting_AU(out_file+"_case_6_AC_SW.png", peaks_coord,ends_fitting[0],ends_fitting[1])
	all_cfx_seqs=all_cfx_seqs + seqs[1]

	print 'len all cfx seq ' + str(len(all_cfx_seqs))
	create_motif_plot("all_cfx.png",all_cfx_seqs)
	create_motif_plot_two_samples("mircocin_cfx.png",all_microcin_seqs,all_cfx_seqs)
	

		

	

	out_file="Un_IP_Mu_1_new_vs_Un_IN_Mu_1_new"
	ends_fitting=fitting_divide_case_6_AC_SW(un_ex_mu_new_file,un_cont_mu_new_file,un_ex_mu_file,un_cont_mu_file,out_file)
	seqs=create_gc_plot_fitting_AU(out_file+"_case_6_AC_SW", peaks_coord,ends_fitting[0],ends_fitting[1])	
	
	
	out_file="Un_IP_Mu_1_vs_Un_IN_Mu_1"
	ends_fitting=fitting_divide_case_6_AC_SW(un_ex_mu_file,un_cont_mu_file,un_ex_mu_new_file,un_cont_mu_new_file,out_file)
	seqs=create_gc_plot_fitting_AU(out_file+"_case_6_AC_SW", peaks_coord,ends_fitting[0],ends_fitting[1])		
	
	'''
	cont_file="Un_IN_Mu_1_ends.wig"
	ex_file="Un_IP_Mu_1_ends.wig"
	out_file="Un_IP_Mu_1_vs_Un_IN_Mu_1"

	ends_fitting=fitting_divide(ex_file,cont_file,out_file)
	#create_gc_plot_fitting(out_file+"av.png", peaks_coord,ends_fitting[2],ends_fitting[3])
	create_gc_plot_fitting(out_file+"av.png", peaks_coord,ends_fitting[0],ends_fitting[3])
	create_gc_plot_fitting(out_file+"fe.png", peaks_coord,ends_fitting[1],ends_fitting[3])
	'''

	'''
	cont_file="Cfx_IN_10_mkM_3_ends.wig"
	ex_file="Cfx_IP_10_mkM_3_ends.wig"
	out_file="Cfx_IP_10_mkM_3_vs_Cfx_IN_10_mkM_3"

	ends_fitting=fitting_divide(ex_file,cont_file,un_ex_file,un_cont_file,out_file)
	#create_gc_plot_fitting(out_file+"av.png", peaks_coord,ends_fitting[2],ends_fitting[3])
	#create_gc_plot_fitting(out_file+"av.png", peaks_coord,ends_fitting[0])
	create_gc_plot_fitting(out_file+"fe_UN.png", peaks_coord,ends_fitting[4],ends_fitting[5])
	'''

	for peak_file in peaks_coord:
	    outfile_list_of_peaks.write(peak_file + "\t" +peaks_coord[peak_file]+"\n")
	outfile_list_of_peaks.close()

def atgc_content():
	a_count=0
	t_count=0
	g_count=0
	c_count=0
	for nt in genomefa:
		#print nt
		if  nt == 'A':
			a_count=a_count+1
		if  nt == 'T':
			t_count=t_count+1
		if  nt == 'G':
			g_count=g_count+1
		if  nt =='C':
			c_count=c_count+1
	print 'len genomefa ' + str (len(genomefa))
	print 'A content ' + str(float(a_count)/len(genomefa))
	print 'T content ' + str(float(t_count)/len(genomefa))
	print 'G content ' + str(float(g_count)/len(genomefa))
	print 'C content ' + str(float(c_count)/len(genomefa))
	print 'GC content ' + str(float(c_count+g_count)/len(genomefa))

def venn_diagram():
    peaks_coords=open(path+ "List_of_peaks.txt", 'r')
    peaks_micro_mu_1={}
    peaks_micro_mu_3={}
    peaks_micro_4={}
    peaks_cfx_mu_1={}
    peaks_cfx_mu_3={}
    peaks_cfx_3={}
    peaks_oxo_mu_1={}
    peaks_oxo_1={}
    peaks_oxo_2={}
    
    for line in peaks_coords:
	peaks = {}
	line = line.split('\t')
	sample = str(line[0])
	#print sample
	for peak in line[1:len(line)]:
	    peak = peak.split('-')
	    peaks[int(peak[0])]=float(peak[1])

	    
	if "Microcin_IP_Mu_1" in sample:
	    peaks_micro_mu_1 = peaks
	if "Microcin_IP_Mu_3" in sample:
	    peaks_micro_mu_3 = peaks
	if "Microcin_IP_50_mkM_4" in sample:
	    peaks_micro_4 = peaks	    
	if "Cfx_IP_Mu_10_mkM_1" in sample:
	    peaks_cfx_mu_1 = peaks		
	if "Cfx_IP_Mu_10_mkM_3" in sample:
	    peaks_cfx_mu_3 = peaks
	if "Cfx_IP_10_mkM_3" in sample:
	    peaks_cfx_3 = peaks	
	if "Oxo_IP_Mu_120_mkM_1" in sample:
	    peaks_oxo_mu_1 = peaks
	if "Oxo_IP_120_mkM_1" in sample:
	    peaks_oxo_1 = peaks
	if "Oxo_IP_120_mkM_2" in sample:
	    peaks_oxo_2 = peaks		    
    
	
    ''' 
    print len (peaks_micro_mu_1)
    print len (peaks_micro_mu_3 )
    print len (peaks_micro_4 	)    
    print len (peaks_cfx_mu_1 	)	
    print len (peaks_cfx_mu_3 )
    print len (peaks_cfx_3 )
    print len (peaks_oxo_mu_1 )
    print len (peaks_oxo_1 )
    print len (peaks_oxo_2  ) 
    '''
    all_micro={}
    all_cfx={}
    all_oxo={}
    all_micro_2={}
    all_cfx_2={}
    all_oxo_2={}    
    '''venn1'''
    
    
    for peak in peaks_micro_mu_1.iterkeys():
	if peak not in all_micro:
	    all_micro[peak]=[]
	    all_micro[peak].append(peaks_micro_mu_1[peak])
	    
    for peak in  peaks_micro_mu_3.iterkeys() :
	if peak not in all_micro:
	    all_micro[peak]=[]
	    all_micro[peak].append(peaks_micro_mu_3[peak])
	else:
	    all_micro[peak].append(peaks_micro_mu_3[peak])	    
		    
    for peak in  peaks_micro_4.iterkeys():
	if peak not in all_micro:
	    all_micro[peak]=[]
	    all_micro[peak].append(peaks_micro_4[peak])   
	else:
	    all_micro[peak].append(peaks_micro_4[peak])   
    #print len (all_micro)    
    

    
    
    for peak in peaks_cfx_mu_1.iterkeys():
	if peak not in all_cfx:
	    all_cfx[peak]=[]
	    all_cfx[peak].append(peaks_cfx_mu_1[peak])
	    
    for peak in  peaks_cfx_mu_3.iterkeys() :
	if peak not in all_cfx:
	    all_cfx[peak]=[]
	    all_cfx[peak].append(peaks_cfx_mu_3[peak])
	else:
	    all_cfx[peak].append(peaks_cfx_mu_3[peak])
		    
    for peak in  peaks_cfx_3.iterkeys():
	if peak not in all_cfx:
	    all_cfx[peak]=[]
	    all_cfx[peak].append(peaks_cfx_3[peak])
	else:
	    all_cfx[peak].append(peaks_cfx_3[peak])	    
  
	    
    #print len (all_cfx)      
	    
    for peak in peaks_oxo_mu_1.iterkeys():
	if peak not in all_oxo:
	    all_oxo[peak]=[]
	    all_oxo[peak].append(peaks_oxo_mu_1[peak])
	else:
	    all_oxo[peak].append(peaks_oxo_mu_1[peak])	
	    
		    
    for peak in  peaks_oxo_1.iterkeys():
	if peak not in all_oxo:
	    all_oxo[peak]=[]
	    all_oxo[peak].append(peaks_oxo_1[peak])
	else:
	    all_oxo[peak].append(peaks_oxo_1[peak])	
	    
    for peak in  peaks_oxo_2.iterkeys():
	if peak not in all_oxo:
	    all_oxo[peak]=[]
	    all_oxo[peak].append(peaks_oxo_2[peak])
	else:
	    all_oxo[peak].append(peaks_oxo_2[peak])			    
    #print len (all_oxo)
	   
	   
	   
    
    seqs_m = []
    seqs_c = []
    seqs_o = []
    for pk in all_micro.iterkeys():
	seq = genomefa[int(pk- win_width + 1):int(pk+ win_width + 1)]
	seqs_m.append(seq)
    print 'len(seqs)=' + str(len(seqs_m))
    
    for pk in all_cfx.iterkeys():
	seq = genomefa[int(pk- win_width + 1):int(pk+ win_width + 1)]
	seqs_c.append(seq)
    print 'len(seqs)=' + str(len(seqs_c))
    
    for pk in all_oxo.iterkeys():
	seq = genomefa[int(pk- win_width + 1):int(pk+ win_width + 1)]
	seqs_o.append(seq)
    print 'len(seqs)=' + str(len(seqs_o))    
	   
    create_motif_plot_two_samples("all", seqs_m,seqs_c,seqs_o)
    
    plt.figure(figsize=(16, 8), dpi=100)    
    venn3([set(all_micro.iterkeys()), set(all_cfx.iterkeys()), set(all_oxo.iterkeys())], ('Micro = 2733 ', 'Cfx = 5502', 'Oxo = 7698'))
    #plt.show()
    plt.title("Venn1")
    plt.savefig(path +"venn1.png", dpi=100, figsize=(16, 8))
    plt.close()
    
    
    '''venn2'''
    for peak in all_micro:
	i=0
	height=[]
	if peak in peaks_micro_mu_1:
	    i+=1
	    height.append(peaks_micro_mu_1[peak])
	if peak in peaks_micro_mu_3:
	    i+=1
	    height.append(peaks_micro_mu_3[peak])	    
	if peak in peaks_micro_4:
	    i+=1
	    height.append(peaks_micro_4[peak])	    
	if i>=2:
	    all_micro_2[peak]=height 
	
	#if (peak in peaks_micro_mu_1 and peak in peaks_micro_mu_3) or (peak in peaks_micro_4 and peak in peaks_micro_mu_1) or (peak in peaks_micro_mu_3 and peak in peaks_micro_4):
	    #all_micro_2.append(peak)
    #print len (all_micro_2)     
    for peak in all_cfx:
	i=0
	height=[]
	if peak in peaks_cfx_mu_1:
	    i+=1
	    height.append(peaks_cfx_mu_1[peak])
	if peak in peaks_cfx_mu_3:
	    i+=1
	    height.append(peaks_cfx_mu_3[peak])	    
	if peak in peaks_cfx_3:
	    i+=1
	    height.append(peaks_cfx_3[peak])	    
	if i>=2:
	    all_cfx_2[peak]=height
	    
	#if peak in peaks_cfx_mu_1 and peak in peaks_cfx_mu_3   or peak in    peaks_cfx_mu_1 and peak in peaks_cfx_3   or peak in  peaks_cfx_3 and peak in peaks_cfx_mu_3:
	    #all_cfx_2.append(peak)	    
    #print len (all_cfx_2)
    
    for peak in all_oxo:
	i=0
	height=[]
	if peak in peaks_oxo_mu_1 or peak==657709:
	    i+=1
	    height.append(peaks_oxo_mu_1[peak])
	if peak in peaks_oxo_1:
	    i+=1
	    height.append(peaks_oxo_1[peak])	    
	if peak in  peaks_oxo_2:
	    i+=1
	    height.append( peaks_oxo_2[peak])	    
	if i>=2 or peak==657709:
	    all_oxo_2[peak]=height	
	
	
    outfile= open(path+"List_of_peaks_MICROCIN_CFX_OXO_2.txt".decode("utf-8"), 'w+')    
    outfile.write("MICROCIN_2")
    for peak in sorted(all_micro_2.iterkeys()):
	outfile.write('\t' + str(peak) + '-' + str(np.mean(all_micro_2[peak]) ))
    
    outfile.write("\nCFX_2")
    for peak in sorted(all_cfx_2.iterkeys()):
	outfile.write('\t' + str(peak) + '-' + str(np.mean(all_cfx_2[peak]) ))
    
    outfile.write("\nOXO_2")
    for peak in sorted(all_oxo_2.iterkeys()):
	outfile.write('\t' + str(peak) + '-' + str(np.mean(all_oxo_2[peak]) ))  
    outfile.close()  
    
    outfile= open(path+"List_of_peaks_MICROCIN_CFX_OXO_1.txt".decode("utf-8"), 'w+')    
    outfile.write("MICROCIN")
    for peak in sorted(all_micro.iterkeys()):
	outfile.write('\t' + str(peak) + '-' + str(np.mean(all_micro[peak]) ))
    
    outfile.write("\nCFX")
    for peak in sorted(all_cfx.iterkeys()):
	outfile.write('\t' + str(peak) + '-' + str(np.mean(all_cfx[peak]) ))
    
    outfile.write("\nOXO")
    for peak in sorted(all_oxo.iterkeys()):
	outfile.write('\t' + str(peak) + '-' + str(np.mean(all_oxo[peak]) ))  
    outfile.close()        
	
	#if peak in peaks_oxo_mu_1 and peak in peaks_oxo_1    or  peak in peaks_oxo_mu_1  and peak in peaks_oxo_1   or peak in   peaks_oxo_1 and  peak in  peaks_oxo_2:
	    #all_oxo_2.append(peak)		    
    #print len (all_oxo_2)   
	   
    plt.figure(figsize=(16, 8), dpi=100)    
    #venn3([set(all_micro_2.iterkeys()), set(all_cfx_2.iterkeys()), set(all_oxo_2.iterkeys())], ('Micro2 = 525', 'Cfx2 = 575', 'Oxo2= 4489'))
    venn2([set(all_micro_2.iterkeys()), set(all_cfx_2.iterkeys())], ('Micro2 = 525', 'Cfx2 = 575'))
    #plt.show()
    plt.title("Venn2")
    plt.savefig(path +"venn2_micro_cfx.png", dpi=100, figsize=(16, 8))
    plt.close() 
    
    plt.figure(figsize=(16, 8), dpi=100)    
    venn3([set(all_micro_2.iterkeys()), set(all_cfx_2.iterkeys()), set(all_oxo_2.iterkeys())], ('Micro2 = 525', 'Cfx2 = 575', 'Oxo2= 4490'))
    #plt.show()
    plt.title("Venn2")
    plt.savefig(path +"venn2.png", dpi=100, figsize=(16, 8))
    plt.close()     
    
    all_peaks1={}
    all_peaks2={}
    all_peaks3={}
    
    for peak in all_micro.iterkeys():
	if peak not in all_peaks1:
	    all_peaks1[peak]=[]
	    all_peaks1[peak].extend(all_micro[peak])
	    
    for peak in all_cfx.iterkeys():
	if peak not in all_peaks1:
	    all_peaks1[peak]=[]
	    all_peaks1[peak].extend(all_cfx[peak])
	else:
	    all_peaks1[peak].extend(all_cfx[peak])
	    
    for peak in all_oxo.iterkeys():
	if peak not in all_peaks1:
	    all_peaks1[peak]=[]
	    all_peaks1[peak].extend(all_oxo[peak])
	else:
	    all_peaks1[peak].extend(all_oxo[peak])	    
    
       
    
    for peak in all_micro_2:
	if peak not in all_peaks2:
	    all_peaks2[peak]=[]
	    all_peaks2[peak].extend(all_micro_2[peak])	    
	    
    for peak in all_cfx_2:
	if peak not in all_peaks2:
	    all_peaks2[peak]=[]
	    all_peaks2[peak].extend(all_cfx_2[peak])	
	else:
	    all_peaks2[peak].extend(all_cfx_2[peak])
	    
    for peak in all_oxo_2:
	if peak not in all_peaks2:
	    all_peaks2[peak]=[]
	    all_peaks2[peak].extend(all_oxo_2[peak])	
	else:
	    all_peaks2[peak].extend(all_oxo_2[peak])
    
       
    for peak in all_micro_2:
	if peak not in all_peaks3:
	    all_peaks3[peak]=[]
	    all_peaks3[peak].extend(all_micro_2[peak])	    
	    
    for peak in all_cfx_2:
	if peak not in all_peaks3:
	    all_peaks3[peak]=[]
	    all_peaks3[peak].extend(all_cfx_2[peak])	
	else:
	    all_peaks3[peak].extend(all_cfx_2[peak])
	    
  
    
    outfile_list_of_peaks_merged= open(path + "List_of_peaks_merged_3.txt".decode("utf-8"), 'w+')
    outfile_list_of_peaks_merged.write("all_peaks_3")
    for peak in all_peaks3:
	outfile_list_of_peaks_merged.write("\t" + str(peak) + "-" +str(np.mean(all_peaks3[peak])))
    outfile_list_of_peaks_merged.close()
    
    outfile_list_of_peaks_merged= open(path + "List_of_peaks_merged_2.txt".decode("utf-8"), 'w+')
    outfile_list_of_peaks_merged.write("all_peaks_2")
    for peak in all_peaks2:
	outfile_list_of_peaks_merged.write("\t" + str(peak) + "-" +str(np.mean(all_peaks2[peak])))
    outfile_list_of_peaks_merged.close()
    
    
    outfile_list_of_peaks_merged= open(path + "List_of_peaks_merged_1.txt".decode("utf-8"), 'w+')
    outfile_list_of_peaks_merged.write("all_peaks_1")
    for peak in all_peaks1:
	outfile_list_of_peaks_merged.write("\t" + str(peak) + "-" +str(np.mean(all_peaks1[peak])))
    outfile_list_of_peaks_merged.close()
	    
    return all_peaks1,all_peaks2,all_peaks3

def motif_analysis():
    peaks_coords=open(path+ "List_of_peaks_merged_3.txt", 'r')
    all_microcin_seqs = []
    all_cfx_seqs = []
    all_seqs = []
    all_seqs_new = []
    peaks_micro_1=[]
    peaks_micro_3=[]
    peaks_cfx_1=[]
    peaks_cfx_3=[]
    background = {'A': 0.245774783354, 'C': 0.2537191331, 'G': 0.254184334046, 'T': 0.246130246797}

    all_unique_peaks=[]
    all_unique_peaks_with_height={}
    for line in peaks_coords:
        peaks = []
        line = line.split('\t')
        sample = str(line[0])
        print sample
        
        for peak in line[1:len(line)]:
            peak = peak.split('-')
            peaks.append(int(peak[0]))
            if int(peak[0]) not in all_unique_peaks:
                all_unique_peaks.append(int(peak[0]))
		all_unique_peaks_with_height[int(peak[0])]=[]
		all_unique_peaks_with_height[int(peak[0])].append(float(peak[1]))
	    else:
		all_unique_peaks_with_height[int(peak[0])].append(float(peak[1]))
	    
	if "Microcin_IP_Mu_1" in sample:
	    peaks_micro_1 = peaks
	if "Microcin_IP_Mu_3" in sample:
	    peaks_micro_3 = peaks	
	if "Cfx_IP_Mu_10_mkM_1" in sample:
	    peaks_cfx_1 = peaks		
	if "Cfx_IP_Mu_10_mkM_3" in sample:
	    peaks_cfx_3 = peaks		
        #print len(peaks)

        # Returns "seqs" array contains sequences under peaks with coordinates from "peaks"
        k = 0
        seqs = []
        for i in range(len(peaks)):
            seq = genomefa[(int(peaks[i]) - win_width_l + 1):(int(peaks[i]) + win_width_r + 1)]
            seqs.append(seq)
        print 'len(seqs)=' + str(len(seqs))

        if len(seqs) == 0:
            return
        #create_motif_plot(sample, seqs)

        if "Microcin" in sample:
            all_microcin_seqs = all_microcin_seqs + seqs
        elif "Cfx" in sample:
            all_cfx_seqs = all_cfx_seqs + seqs

    common_micro=[]
    common_cfx=[]
    for p in peaks_micro_1:
	if p in peaks_micro_3:
	    common_micro.append(p)
    for p in peaks_cfx_1:
	if p in peaks_cfx_3:
	    common_cfx.append(p)
    common_unique=[]
    for p in common_micro+common_cfx:
	if p not in common_unique:
	    common_unique.append(p)
	    
    print 'len common micro ' + str(len(common_micro))
    print 'len common cfx ' + str(len(common_cfx))
    
    
    
    common_unique=peaks
    seqs=[]   
    for i in range(len(common_unique)):
	seq = genomefa[(int(common_unique[i]) - win_width_l + 1):(int(common_unique[i]) + win_width_r + 1)]
	seqs.append(seq)
    print 'len common unique ' + str(len(seqs))	       
	    
 
    
    #all_seqs = all_microcin_seqs + all_cfx_seqs
    #create_motif_plot("all_micro_cfx_seq.png", all_seqs)
    
    
    all_seqs=seqs
    create_motif_plot("all_intersection_seq.png", all_seqs)
    
    
    #hard_copy = ['A'] * 1626 + ['C'] * 1678 + ['G'] * 1681 + ['T'] * 1627
    hard_copy = ['A'] * (int(background['A']* len (all_seqs)) +3) + ['C'] * (int(background['C']* len (all_seqs))) +  ['G'] * (int(background['G']* len (all_seqs))) + ['T'] * (int(background['T']* len (all_seqs)))
    print len(hard_copy)
    print len (all_seqs)
    for k in range(65, 69):
        soft_copy=[]
        soft_copy.extend(hard_copy)
        shuffle(soft_copy)
        l = 0
        for j in range(0, len(all_seqs)):
            ls=list(all_seqs[j])
            ls[k]= soft_copy[l]
            l = l + 1
            all_seqs[j]="".join(ls)

    create_motif_plot("all_seq_edt.png", all_seqs)
    #create_motif_plot_two_samples("mircocin_cfx.png", all_microcin_seqs, all_cfx_seqs)
    
    instances = []
    for seq in all_seqs:
        instances.append(Seq(seq))
	

    m = motifs.create(instances)
    pwm = m.counts.normalize()
    pssm = pwm.log_odds(background)
    peaks_output=[[],[],[],[]]
    
    #print "Peak Av_height  Position : score "
    
    #for peak in all_unique_peaks:
    for peak in common_unique:
        test_seq = Seq(str(genomefa[peak - win_width_l + 1 :peak + win_width_r + 1 ]), IUPAC.unambiguous_dna)  #
	search_result=pssm.search(test_seq,threshold=-100)
	max_score=-100
	max_pos=0
	b=0
	for position, score in search_result:
	    if position==0:
	    #if score>max_score:
		max_score=score
		max_pos=position
		b=1
	if b==1:	
	    peaks_output[0].append(peak)
	    peaks_output[1].append(np.mean(all_unique_peaks_with_height[peak]))
	    peaks_output[2].append(max_pos)
	    peaks_output[3].append(float(max_score))	    
	    #print("%d  %5.3f  %d  %5.3f" % (peak, np.mean(all_unique_peaks_with_height[peak]), max_pos, max_score))    
    print 'cor coefficient ' + str(pearsonr(peaks_output[1],peaks_output[3]))
          
    
    test_seq = Seq(str(genomefa), IUPAC.unambiguous_dna)
    whole_genome_scores=pssm.calculate(test_seq)
    
    print 'mean height peaks ' + str(np.mean(peaks_output[1]))
    print 'mean motif score peaks ' + str(np.mean(peaks_output[3]))
    print 'mean motif score genome ' + str(np.mean(whole_genome_scores))
    
    plt.figure(figsize=(16, 8), dpi=100)
    f,axarr = plt.subplots(3,1)  
    ax=axarr[0]
    ends=[]
    ends.extend(peaks_output[1])
    #ax.hist(ends,label="height")
    ax.hist(ends, color='green', alpha = 0.5)  
    ax.set_xlabel('Height')
    ax.set_ylabel('Number of peaks') 
    ax.tick_params(axis='both', which='major', labelsize=9)
    #ax.legend()
        
    ax=axarr[1]
    ends=[]
    ends.extend(peaks_output[3])
    ax.hist(ends,label="peaks score")
    ax.hist(ends, color='blue', alpha = 0.1)    
    ax.set_xlim((-30,30))
    ax.set_xlabel('Score')
    ax.set_ylabel('Number of peaks')   
    ax.tick_params(axis='both', which='major', labelsize=9)
    #ax.legend()
    
    ax=axarr[2]
    #ax.hist(whole_genome_scores,label="wholegenome score")
    ax.hist(whole_genome_scores, color='red', alpha = 0.5)
    ax.set_xlim((-30,30))
    ax.set_yticks(range(0,2000000,500000))
    ax.set_xlabel('Score')
    ax.set_ylabel('Number of sites')    
    ax.tick_params(axis='both', which='major', labelsize=9)
    #ax.legend()
    
    plt.tight_layout()
    plt.savefig(path+'height_motif_score_distributuion_all_peaks_3', dpi=100, figsize=(16, 8))
    
    #if 1==1:
	#return
    
    outfile= open(path +"List_of_sites_th_0_all_peaks_3_calculate.txt".decode("utf-8"), 'w+')
    for i in range(0,len(whole_genome_scores)):
	outfile.write(str(whole_genome_scores[i]) + "\n")
    outfile.close()
    
    
        
    test_seq = Seq(str(genomefa), IUPAC.unambiguous_dna).reverse_complement()
    whole_genome_scores=pssm.calculate(test_seq)   
    outfile= open(path +"List_of_sites_th_0_all_peaks_3_calculate_rc.txt".decode("utf-8"), 'w+')
    for i in range(0,len(whole_genome_scores)):
	outfile.write(str(whole_genome_scores[i]) + "\n")
    outfile.close()    
    
    '''    
    test_seq = Seq(str(genomefa), IUPAC.unambiguous_dna)
    print 'len genomefa ' + str(len(str(genomefa)))
    print 'len genomefa ' + str(len(genomefa))    
    
    pwm = m.counts.normalize()
    pssm = pwm.log_odds(background)
    
    distribution=pssm.distribution(background=background, precision = 10**4)
    
    threshold = distribution.threshold_fpr(0.01)
    print 'threshold fpr = ' + str(threshold)
    
    threshold = distribution.threshold_fnr(0.01)
    print 'threshold fnr = ' + str(threshold)
    
    threshold=distribution.threshold_balanced(1000)
    print 'threshold balanced = ' + str(threshold)
    
    
    size=0
    outfile= open(path +"List_of_sites_th_0_all_peaks_2.txt".decode("utf-8"), 'w+')
    
    for i in range (0, len (genomefa)- (win_width_l + win_width_r +1)):
	if i%200000==0:
	    print i
	test_seq = Seq(str(genomefa[i  :i+ win_width_l  + win_width_r + 1 ]), IUPAC.unambiguous_dna)  #
	search_result=pssm.search(test_seq,threshold=0)
    
	#search_result=pssm.search(test_seq, threshold=0)
	for position, score in search_result:
	    size=size+1
	    if position>0:
		outfile.write("Position %d: score = %5.3f" % (position+i, score)+ "\n")
	    #print("Position %d: score = %5.3f" % (position, score))
	#print size 
    outfile.close()
    '''
    
def september_calculate_correlation():
    peaks_coords=open(path+ "List_of_peaks_MICROCIN_CFX_OXO_2.txt", 'r')
    micro={}
    cfx={}
    oxo={}
    micro_cfx={}
    micro_cfx_oxo={}
    
    
    

    all_unique_peaks_with_height={}
    for line in peaks_coords:
        peaks = {}
        line = line.split('\t')
        sample = str(line[0])
        print sample
        
        for peak in line[1:len(line)]:
            peak = peak.split('-')
            peaks[int(peak[0])]= float(peak[1])
	    
	    if 'MICROCIN' in sample or 'CFX' in sample:
		if int(peak[0]) not in micro_cfx:
		    micro_cfx[int(peak[0])]= float(peak[1])
		if int(peak[0]) in micro_cfx:
		    tmp=micro_cfx[int(peak[0])]
		    micro_cfx[int(peak[0])]=(1.0 * tmp +float(peak[1]))/2
	    if int(peak[0]) not in micro_cfx_oxo:
		micro_cfx_oxo[int(peak[0])]=[]
		micro_cfx_oxo[int(peak[0])].append(float(peak[1]))		
	    if int(peak[0]) in micro_cfx_oxo:
		micro_cfx_oxo[int(peak[0])].append(float(peak[1]))	    
				
   
	if "MICROCIN" in sample:
	    micro = peaks		
	if "CFX" in sample:
	    cfx = peaks		
	if "OXO" in sample:
	    oxo = peaks		

    for peak in micro_cfx_oxo:
	tmp=micro_cfx_oxo[peak]
	micro_cfx_oxo[peak]=np.mean(tmp)
	
    september_help_correlation("Micro", micro, 1 )
    september_help_correlation("Cfx", cfx, 2 )
    september_help_correlation("Oxo", oxo, 2 )
    september_help_correlation("Micro_Cfx", micro_cfx, 3 )
    september_help_correlation("Micro_Cfx_Oxo", micro_cfx_oxo, 2 )
    peaks_coords.close()
     
def september_help_correlation(name, peaks, delta ):
    print "\n" + name +"\n"
    background = {'A': 0.245774783354, 'C': 0.2537191331, 'G': 0.254184334046, 'T': 0.246130246797}
    common_unique=list(peaks.keys())
    seqs=[]   
    
    for i in range(len(common_unique)):
	seq = genomefa[(int(common_unique[i]) - win_width_l + 1):(int(common_unique[i]) + win_width_r + 1)]
	seqs.append(seq)
    print 'len common unique ' + str(len(seqs))	      
    all_seqs=seqs
    hard_copy = ['A'] * (int(background['A']* len (all_seqs)) +delta) + ['C'] * (int(background['C']* len (all_seqs))) +  ['G'] * (int(background['G']* len (all_seqs))) + ['T'] * (int(background['T']* len (all_seqs)))
    print len(hard_copy)
    print len (all_seqs)
    for k in range(65, 69):
        soft_copy=[]
        soft_copy.extend(hard_copy)
        shuffle(soft_copy)
        l = 0
        for j in range(0, len(all_seqs)):
            ls=list(all_seqs[j])
            ls[k]= soft_copy[l]
            l = l + 1
            all_seqs[j]="".join(ls)
 
    instances = []
    for seq in all_seqs:
        instances.append(Seq(seq))
	

    m = motifs.create(instances)
    pwm = m.counts.normalize()
    pssm = pwm.log_odds(background)
    peaks_output=[[],[],[],[]]
 
    for peak in common_unique:
        test_seq = Seq(str(genomefa[peak - win_width_l + 1 :peak + win_width_r + 1 ]), IUPAC.unambiguous_dna)  #
	search_result=pssm.search(test_seq,threshold=-100)
	max_score=-100
	max_pos=0
	b=0
	for position, score in search_result:
	    if position==0:
	    #if score>max_score:
		max_score=score
		max_pos=position
		b=1
	if b==1:	
	    peaks_output[0].append(peak)
	    peaks_output[1].append(peaks[peak])
	    peaks_output[2].append(max_pos)
	    peaks_output[3].append(float(max_score))	    
	    #print("%d  %5.3f  %d  %5.3f" % (peak, np.mean(all_unique_peaks_with_height[peak]), max_pos, max_score))    
    print 'cor coefficient ' + str(pearsonr(peaks_output[1],peaks_output[3]))
    
    
    operons=[[3597167, 3602272, '-'], [3466047, 3471144, '-'], [3690984, 3695995, '-'], [3424645, 3429656, '-'], [4212776, 4217871, '+'], [2725847, 2730930, '-'], [223771, 228875, '+']]
    
    peaks_operon_in_operon=[]
    peaks_operon_in_peak=[]
    peaks_operon_in_height=[]
    peaks_operon_in_score=[]
    peaks_operon_up_operon=[]
    peaks_operon_up_peak=[]
    peaks_operon_up_height=[]
    peaks_operon_up_score=[]
    peaks_operon_down_operon=[]
    peaks_operon_down_peak=[]
    peaks_operon_down_height=[]
    peaks_operon_down_score=[]
    width=5000
    
    for i in range(0, len(peaks_output[0])):
	peak=peaks_output[0][i]
	height=peaks_output[1][i]
	score=peaks_output[3][i]
	
	for operon in operons:
	    if peak >=operon[0] and peak <=operon[1]:
		peaks_operon_in_operon.append(operon)
		peaks_operon_in_peak.append(peak)
		peaks_operon_in_height.append(height)
		peaks_operon_in_score.append(score)		
		
	    if peak>=operon[0]-width:
		if peak<operon[0]:
		    if operon[2]=='-':
			peaks_operon_down_operon.append(operon)
			peaks_operon_down_peak.append(peak)
			peaks_operon_down_height.append(height)
			peaks_operon_down_score.append(score)
		    elif operon[2]=='+':
			peaks_operon_up_operon.append(operon)
			peaks_operon_up_peak.append(peak)
			peaks_operon_up_height.append(height)
			peaks_operon_up_score.append(score)
			
	    if peak>operon[1]:
		if peak<=operon[1]+width:
		    if operon[2]=='-':
			peaks_operon_up_operon.append(operon)
			peaks_operon_up_peak.append(peak)
			peaks_operon_up_height.append(height)
			peaks_operon_up_score.append(score)
		    elif operon[2]=='+':
			peaks_operon_down_operon.append(operon)
			peaks_operon_down_peak.append(peak)
			peaks_operon_down_height.append(height)
			peaks_operon_down_score.append(score)
			    
			    
    print str (len (peaks_operon_in_height) + len (peaks_operon_up_height) + len (peaks_operon_down_height))
    
    plt.figure(figsize=(16, 8), dpi=100)
   
    plot1 = plt.subplot()    
    grey=plot1.scatter(peaks_output[3], peaks_output[1],color='grey',alpha =0.5)
    plot1.set_ylim(0, np.max(peaks_output[1]) + 10)    
        

   
    

    azure= plot1.scatter(peaks_operon_down_score, peaks_operon_down_height,color='maroon')    
    blue=plot1.scatter(peaks_operon_in_score, peaks_operon_in_height,color='blue') 
    red= plot1.scatter(peaks_operon_up_score, peaks_operon_up_height,color='red') 
    plot1.legend(( blue, red, azure),
           ('Inside', 'Upstream', 'Downstream'),
           scatterpoints=1,
           loc='upper right',
           ncol=1,
           fontsize=10)
    plot1.set_xlabel('Score')
    plot1.set_ylabel('Height') 
    plot1.tick_params(axis='both', which='major', labelsize=9)
    if name=='Micro':
	title="Microcin"
    if name=='Cfx':
	title="Ciprofloxacin"
    if name=='Oxo':
	title="Oxolinic acid"
    if name=='Micro_Cfx':
	title="Microcin and Ciprofloxacin"
    if name=='Micro_Cfx_Oxo':
	title="Microcin and Ciprofloxacin and Oxolinic acid"
    plt.title(title)
    plt.savefig(path+name +"_scatter.png", dpi=100, figsize=(16, 8))
    
    plt.close()
    
    print  name
    operon= {}
    for i in range (0,len(operons)):
	operon[i]=[[],[],[]]
	for j in range (0, len (peaks_operon_in_operon)):
	    if peaks_operon_in_operon[j]==operons[i]:
		operon[i][1].append(j)
	for j in range (0, len (peaks_operon_up_operon)):	
	    if peaks_operon_up_operon[j]==operons[i]:
		operon[i][0].append(j)
	for j in range (0, len (peaks_operon_down_operon)):
	    if peaks_operon_down_operon[j]==operons[i]:
		operon[i][2].append(j)
		
    for i in range (0,len(operons)):
	print str(operons[i]) + "\t", 
	for index in operon[i][0]:
	    print str(peaks_operon_up_peak[index] ) + "-" + str (peaks_operon_up_height[index] ) + "-" + str (peaks_operon_up_score[index] ) + ";",
	print  "\t",
	for index in operon[i][1]:
	    print str(peaks_operon_in_peak[index] ) + "-" + str (peaks_operon_in_height[index] ) + "-" + str (peaks_operon_in_score[index] ) + ";",
	print  "\t",
	for index in operon[i][2]:
	    print str(peaks_operon_down_peak[index] ) + "-" + str (peaks_operon_down_height[index] ) + "-" + str (peaks_operon_down_score[index] ) + ";",
	print  "\n" 	
	
    
    
    
def motif_analysis_random():
    peaks_coords=open(path+ "List_of_peaks_merged_2.txt", 'r')
    all_microcin_seqs = []
    all_cfx_seqs = []
    all_seqs = []
    all_seqs_new = []
    peaks_micro_1=[]
    peaks_micro_3=[]
    peaks_cfx_1=[]
    peaks_cfx_3=[]
    background = {'A': 0.245774783354, 'C': 0.2537191331, 'G': 0.254184334046, 'T': 0.246130246797}

    all_unique_peaks=[]
    all_unique_peaks_with_height={}
    for line in peaks_coords:
        peaks = []
        line = line.split('\t')
        sample = str(line[0])
        print sample
        
        for peak in line[1:len(line)]:
            peak = peak.split('-')
            peaks.append(int(peak[0]))
            if int(peak[0]) not in all_unique_peaks:
                all_unique_peaks.append(int(peak[0]))
		all_unique_peaks_with_height[int(peak[0])]=[]
		all_unique_peaks_with_height[int(peak[0])].append(float(peak[1]))
	    else:
		all_unique_peaks_with_height[int(peak[0])].append(float(peak[1]))
	

    
    
    common_unique=peaks
    seqs=[]   
    for i in range(len(common_unique)):
	seq = genomefa[(int(common_unique[i]) - win_width_l + 1):(int(common_unique[i]) + win_width_r + 1)]
	seqs.append(seq)
    print 'len common unique ' + str(len(seqs))	       
	    
 
    
    all_seqs=seqs
    #create_motif_plot("all_intersection_seq.png", all_seqs)
  
    hard_copy = ['A'] * (int(background['A']* len (all_seqs)) +2) + ['C'] * (int(background['C']* len (all_seqs))) +  ['G'] * (int(background['G']* len (all_seqs))) + ['T'] * (int(background['T']* len (all_seqs)))
    print len(hard_copy)
    print len (all_seqs)
    for k in range(65, 69):
        soft_copy=[]
        soft_copy.extend(hard_copy)
        shuffle(soft_copy)
        l = 0
        for j in range(0, len(all_seqs)):
            ls=list(all_seqs[j])
            ls[k]= soft_copy[l]
            l = l + 1
            all_seqs[j]="".join(ls)

    #create_motif_plot("all_seq_edt.png", all_seqs)
    #create_motif_plot_two_samples("mircocin_cfx.png", all_microcin_seqs, all_cfx_seqs)
    
    instances = []
    for seq in all_seqs:
        instances.append(Seq(seq))
	

    m = motifs.create(instances)
    pwm = m.counts.normalize()
    print pwm
    perfect_seq=""
    for i in range(0,len(pwm[0])):	
	if pwm[0][i] >=pwm[1][i] and pwm[0][i] >=pwm[2][i]  and pwm[0][i] >=pwm[3][i]:
	    perfect_seq=perfect_seq+"a"
	    #print "A",
	elif pwm[1][i] >=pwm[0][i] and pwm[1][i] >=pwm[2][i]  and pwm[1][i] >=pwm[3][i]:
	    #print "C",
	    perfect_seq=perfect_seq+"c"
	elif pwm[2][i] >=pwm[0][i] and pwm[2][i] >=pwm[1][i]  and pwm[2][i] >=pwm[3][i]:
	    #print "G",
	    perfect_seq=perfect_seq+"g"
	elif pwm[3][i] >=pwm[0][i] and pwm[3][i] >=pwm[1][i]  and pwm[3][i] >=pwm[2][i]:
	    #print "T",	
	    perfect_seq=perfect_seq+"t"
    
    
    pssm = pwm.log_odds(background)
    peaks_output=[[],[],[],[]]   
    test_seq = Seq(str(perfect_seq), IUPAC.unambiguous_dna)
    whole_genome_scores=pssm.calculate(test_seq)   
    #outfile= open(path +"List_of_sites_th_0_all_peaks_2_calculate_rc.txt".decode("utf-8"), 'w+')
    
    print perfect_seq
    seq=Seq(perfect_seq)
    print su.GC(seq)         
    print whole_genome_scores  
      
    
    
    ls=list(perfect_seq)
    random.shuffle(ls)
    random_perfect_seq="".join(ls)
    seq=Seq(random_perfect_seq)
    print random_perfect_seq
    print su.GC(seq)    
    
    test_seq = Seq(str(random_perfect_seq), IUPAC.unambiguous_dna)
    whole_genome_scores=pssm.calculate(test_seq)   
    #outfile= open(path +"List_of_sites_th_0_all_peaks_2_calculate_rc.txt".decode("utf-8"), 'w+')
    print whole_genome_scores     
    
    ls=list(perfect_seq)
    random.shuffle(ls)
    random_perfect_seq="".join(ls)
    seq=Seq(random_perfect_seq)
    print random_perfect_seq
    print su.GC(seq)
    
    test_seq = Seq(str(random_perfect_seq), IUPAC.unambiguous_dna)
    whole_genome_scores=pssm.calculate(test_seq)   
    #outfile= open(path +"List_of_sites_th_0_all_peaks_2_calculate_rc.txt".decode("utf-8"), 'w+')
    print whole_genome_scores 
    
    
    ls=list(perfect_seq)
    random.shuffle(ls)
    random_perfect_seq="".join(ls)
    seq=Seq(random_perfect_seq)
    print random_perfect_seq
    print su.GC(seq)

    test_seq = Seq(str(random_perfect_seq), IUPAC.unambiguous_dna)
    whole_genome_scores=pssm.calculate(test_seq)   
    #outfile= open(path +"List_of_sites_th_0_all_peaks_2_calculate_rc.txt".decode("utf-8"), 'w+')
    print whole_genome_scores     
    
def motif_analysis_plasmid():
    peaks_coords=open(path+ "List_of_peaks_merged_3.txt", 'r')
    background = {'A': 0.245774783354, 'C': 0.2537191331, 'G': 0.254184334046, 'T': 0.246130246797}
    all_unique_peaks=[]
    all_unique_peaks_with_height={}
    for line in peaks_coords:
        peaks = []
        line = line.split('\t')
        sample = str(line[0])
        print sample
        
        for peak in line[1:len(line)]:
            peak = peak.split('-')
            peaks.append(int(peak[0]))
            
        # Returns "seqs" array contains sequences under peaks with coordinates from "peaks"
        k = 0
        seqs = []
        for i in range(len(peaks)):
            seq = genomefa[(int(peaks[i]) - win_width_l + 1):(int(peaks[i]) + win_width_r + 1)]
            seqs.append(seq)
        print 'len(seqs)=' + str(len(seqs))

        if len(seqs) == 0:
            return

    
    common_unique=peaks
    seqs=[]   
    for i in range(len(common_unique)):
	seq = genomefa[(int(common_unique[i]) - win_width_l + 1):(int(common_unique[i]) + win_width_r + 1)]
	seqs.append(seq)
    print 'len common unique ' + str(len(seqs))	       
	    

    all_seqs=seqs
   
    
    #hard_copy = ['A'] * 1626 + ['C'] * 1678 + ['G'] * 1681 + ['T'] * 1627
    hard_copy = ['A'] * (int(background['A']* len (all_seqs)) +3) + ['C'] * (int(background['C']* len (all_seqs))) +  ['G'] * (int(background['G']* len (all_seqs))) + ['T'] * (int(background['T']* len (all_seqs)))
    print len(hard_copy)
    print len (all_seqs)
    for k in range(65, 69):
        soft_copy=[]
        soft_copy.extend(hard_copy)
        shuffle(soft_copy)
        l = 0
        for j in range(0, len(all_seqs)):
            ls=list(all_seqs[j])
            ls[k]= soft_copy[l]
            l = l + 1
            all_seqs[j]="".join(ls)
    
    instances = []
    for seq in all_seqs:
        instances.append(Seq(seq))
	

    m = motifs.create(instances)
    pwm = m.counts.normalize()
    pssm = pwm.log_odds(background)
    peaks_output=[[],[],[],[]]
    
    plt.figure(figsize=(16, 8), dpi=100)
    f,axarr = plt.subplots(4,1)  
 
   
    plasmid_scores=[]
    test_seq = Seq(str(pSC101), IUPAC.unambiguous_dna)
    whole_genome_scores=pssm.calculate(test_seq)
    plasmid_scores.extend(whole_genome_scores)
    outfile= open(path +"List_of_sites_pSC101_calculate.txt".decode("utf-8"), 'w+')
    for i in range(0,len(whole_genome_scores)):
	outfile.write(str(whole_genome_scores[i]) + "\n")
    outfile.close()    
    #ax=axarr[0]
    #ax.hist(whole_genome_scores, label = 'pSC101')
    #ax.set_xlim((-30,30))   
    #ax.legend()    
    test_seq = Seq(str(pSC101), IUPAC.unambiguous_dna).reverse_complement()
    whole_genome_scores=pssm.calculate(test_seq)
    plasmid_scores.extend(whole_genome_scores)
    outfile= open(path +"List_of_sites_pSC101_calculate_rc.txt".decode("utf-8"), 'w+')
    for i in range(0,len(whole_genome_scores)):
	outfile.write(str(whole_genome_scores[i]) + "\n")
    outfile.close()    
    #ax=axarr[1]
    #ax.hist(whole_genome_scores, label = 'pSC101 rc')
    #ax.set_xlim((-30,30))       
    #ax.legend()    
    ax=axarr[0]
    ax.hist(plasmid_scores, color='green', alpha = 0.5)
    ax.set_xlim((-30,30))       
    ax.set_xlabel('Score',size = 9)
    ax.set_ylabel('Number of sites',size = 9)    
    ax.tick_params(axis='both', which='major', labelsize=9)    
    ax.set_title("pSC101",size = 10)
    
    
    plasmid_scores=[]   
    test_seq = Seq(str(pBR322), IUPAC.unambiguous_dna)
    whole_genome_scores=pssm.calculate(test_seq)
    plasmid_scores.extend(whole_genome_scores)
    outfile= open(path +"List_of_sites_pBR322_calculate.txt".decode("utf-8"), 'w+')
    for i in range(0,len(whole_genome_scores)):
	outfile.write(str(whole_genome_scores[i]) + "\n")
    outfile.close() 
    #ax=axarr[2]
    #ax.hist(whole_genome_scores, label = 'pBR322')
    #ax.set_xlim((-30,30))   
    #ax.legend()    
    test_seq = Seq(str(pBR322), IUPAC.unambiguous_dna).reverse_complement()
    whole_genome_scores=pssm.calculate(test_seq)
    plasmid_scores.extend(whole_genome_scores)
    outfile= open(path +"List_of_sites_pBR322_calculate_rc.txt".decode("utf-8"), 'w+')
    for i in range(0,len(whole_genome_scores)):
	outfile.write(str(whole_genome_scores[i]) + "\n")
    outfile.close() 
    #ax=axarr[3]
    #ax.hist(whole_genome_scores, label = 'pBR322 rc')
    #ax.set_xlim((-30,30))    
    #ax.legend()
    ax=axarr[1]
    ax.hist(plasmid_scores, color='orange', alpha = 0.5)
    ax.set_xlim((-30,30))       
    ax.set_xlabel('Score',size = 9)
    ax.set_ylabel('Number of sites',size = 9)    
    ax.tick_params(axis='both', which='major', labelsize=9)    
    ax.set_title("pBR322",size = 10)
    
    plasmid_scores=[]
    test_seq = Seq(str(pUC19), IUPAC.unambiguous_dna)
    whole_genome_scores=pssm.calculate(test_seq)
    plasmid_scores.extend(whole_genome_scores)
    outfile= open(path +"List_of_sites_pUC19_calculate.txt".decode("utf-8"), 'w+')
    for i in range(0,len(whole_genome_scores)):
	outfile.write(str(whole_genome_scores[i]) + "\n")
    outfile.close() 
    #ax=axarr[4]
    #ax.hist(whole_genome_scores, label = 'pUC19')
    #ax.set_xlim((-30,30))   
    #ax.legend()    
    test_seq = Seq(str(pUC19), IUPAC.unambiguous_dna).reverse_complement()
    whole_genome_scores=pssm.calculate(test_seq)   
    plasmid_scores.extend(whole_genome_scores)
    outfile= open(path +"List_of_sites_pUC19_calculate_rc.txt".decode("utf-8"), 'w+')
    for i in range(0,len(whole_genome_scores)):
	outfile.write(str(whole_genome_scores[i]) + "\n")
    outfile.close() 
    #ax=axarr[5]
    #ax.hist(whole_genome_scores, label = 'pUC19 rc')
    #ax.set_xlim((-30,30))    
    #ax.legend()
    ax=axarr[2]
    ax.hist(plasmid_scores, color='blue', alpha = 0.3)
    ax.set_xlim((-30,30))       
    ax.set_xlabel('Score',size = 9)
    ax.set_ylabel('Number of sites',size = 9)    
    ax.tick_params(axis='both', which='major', labelsize=9)    
    ax.set_title("pUC19",size = 10)
    
    plasmid_scores=[]    
    test_seq = Seq(str(pMu), IUPAC.unambiguous_dna)
    whole_genome_scores=pssm.calculate(test_seq)
    plasmid_scores.extend(whole_genome_scores)
    outfile= open(path +"List_of_sites_pMu_calculate.txt".decode("utf-8"), 'w+')
    for i in range(0,len(whole_genome_scores)):
	outfile.write(str(whole_genome_scores[i]) + "\n")
    outfile.close() 
    #ax=axarr[6]
    #ax.hist(whole_genome_scores, label = 'pMu')
    #ax.set_xlim((-30,30))   
    #ax.legend()    
    test_seq = Seq(str(pMu), IUPAC.unambiguous_dna).reverse_complement()
    whole_genome_scores=pssm.calculate(test_seq)   
    plasmid_scores.extend(whole_genome_scores)
    outfile= open(path +"List_of_sites_pMu_calculate_rc.txt".decode("utf-8"), 'w+')
    for i in range(0,len(whole_genome_scores)):
	outfile.write(str(whole_genome_scores[i]) + "\n")
    outfile.close() 
    #ax=axarr[7]
    #ax.hist(whole_genome_scores, label = 'pMu rc')
    #ax.set_xlim((-30,30))    
    #ax.legend()
    ax=axarr[3]
    ax.hist(plasmid_scores, color='red', alpha = 0.5)
    ax.set_xlim((-30,30))       
    ax.set_xlabel('Score',size = 9)
    ax.set_ylabel('Number of sites',size = 9)    
    ax.tick_params(axis='both', which='major', labelsize=9)    
    ax.set_title("Bacteriophage Mu",size = 10)
    
    for i in range(-10,11):
	test_seq=pUC19[(2660 + i  - win_width_l):] + pUC19[0:win_width_r - (len(pUC19)- (2660+i))]
	test_seq = Seq(str(test_seq), IUPAC.unambiguous_dna)
	whole_genome_scores=pssm.calculate(test_seq)   
	#print str(2660 + i) + " ",
	#print whole_genome_scores
    
    plt.tight_layout()
    plt.savefig(path+'motif_score_distribution_plasmids', dpi=100, figsize=(16, 8))
    
        
    
    
def motif_analysis_half_motif():
    peaks_coords=open(path+ "List_of_peaks_merged_2.txt", 'r')
    all_microcin_seqs = []
    all_cfx_seqs = []
    all_seqs = []
    all_seqs_new = []
    peaks_micro_1=[]
    peaks_micro_3=[]
    peaks_cfx_1=[]
    peaks_cfx_3=[]
    background = {'A': 0.245774783354, 'C': 0.2537191331, 'G': 0.254184334046, 'T': 0.246130246797}

    all_unique_peaks=[]
    all_unique_peaks_with_height={}
    for line in peaks_coords:
        peaks = []
        line = line.split('\t')
        sample = str(line[0])
        print sample
        
        for peak in line[1:len(line)]:
            peak = peak.split('-')
            peaks.append(int(peak[0]))
	    if int(peak[0]) not in all_unique_peaks:
		all_unique_peaks.append(int(peak[0]))
		all_unique_peaks_with_height[int(peak[0])]=[]
		all_unique_peaks_with_height[int(peak[0])].append(float(peak[1]))
	    else:
		all_unique_peaks_with_height[int(peak[0])].append(float(peak[1]))	    
    
    common_unique=peaks
    seqs_left=[]   
    seqs_right=[]   
    for i in range(len(common_unique)):
	seq=genomefa[(int(common_unique[i]) - win_width_l + 1):(int(common_unique[i])  -1 )]
	seqs_left.append(seq)
	seq=genomefa[(int(common_unique[i])+4):(int(common_unique[i]) + win_width_r + 1)]
	seqs_right.append(seq)	
	#seq = genomefa[(int(common_unique[i]) - win_width_l + 1):(int(common_unique[i]) + win_width_r + 1)]
	
    all_seqs_left=seqs_left  
    all_seqs_right=seqs_right

    instances_left = []
    instances_right = []
    for seq in all_seqs_left:
        instances_left.append(Seq(str(seq)))
    for seq in all_seqs_right:
	instances_right.append(Seq(str(seq)))
	
    print 'LEFT'
    m = motifs.create(instances_left)
    pwm = m.counts.normalize()
    pssm = pwm.log_odds(background)
    peaks_output=[[],[],[],[]]

    for peak in common_unique:
        test_seq = Seq(str(genomefa[peak - win_width_l + 1 :peak -1 ]), IUPAC.unambiguous_dna)  #
	search_result=pssm.search(test_seq,threshold=-100)
	max_score=-100
	max_pos=0
	b=0
	for position, score in search_result:
	    if position==0:
		max_score=score
		max_pos=position
		b=1
	if b==1:	
	    peaks_output[0].append(peak)
	    peaks_output[1].append(np.mean(all_unique_peaks_with_height[peak]))
	    peaks_output[2].append(max_pos)
	    peaks_output[3].append(float(max_score))	    
    print 'cor coefficient ' + str(pearsonr(peaks_output[1],peaks_output[3]))
          
    
    test_seq = Seq(str(genomefa), IUPAC.unambiguous_dna)
    whole_genome_scores=pssm.calculate(test_seq)
    
    print 'mean height peaks ' + str(np.mean(peaks_output[1]))
    print 'mean motif score peaks ' + str(np.mean(peaks_output[3]))
    print 'mean motif score genome ' + str(np.mean(whole_genome_scores))
    
    plt.figure(figsize=(16, 8), dpi=100)
    f,axarr = plt.subplots(3,1)  
    ax=axarr[0]
    ends=[]
    ends.extend(peaks_output[1])
    ax.hist(ends,label="height")
    ax.legend()
        
    ax=axarr[1]
    ends=[]
    ends.extend(peaks_output[3])
    ax.hist(ends,label="peaks score")
    ax.set_xlim((-30,30))
    ax.legend()
    
    ax=axarr[2]
    ax.hist(whole_genome_scores,label="wholegenome score")
    ax.set_xlim((-30,30))
    ax.legend()
    
    plt.savefig(path+'height_motif_left_score_distributuion_all_peaks_2', dpi=100, figsize=(16, 8))

    outfile= open(path +"List_of_sites_left_all_peaks_2_calculate.txt".decode("utf-8"), 'w+')
    for i in range(0,len(whole_genome_scores)):
	outfile.write(str(whole_genome_scores[i]) + "\n")
    outfile.close()   
    
    test_seq = Seq(str(genomefa), IUPAC.unambiguous_dna).reverse_complement()
    whole_genome_scores=pssm.calculate(test_seq)   
    outfile= open(path +"List_of_sites_left_all_peaks_2_calculate_rc.txt".decode("utf-8"), 'w+')
    for i in range(0,len(whole_genome_scores)):
	outfile.write(str(whole_genome_scores[i]) + "\n")
    outfile.close()    
    
    print 'RIGHT'
    m = motifs.create(instances_right)
    pwm = m.counts.normalize()
    pssm = pwm.log_odds(background)
    peaks_output=[[],[],[],[]]

    for peak in common_unique:
        test_seq = Seq(str(genomefa[peak +4 :peak  + win_width_r + 1]), IUPAC.unambiguous_dna)  #
	search_result=pssm.search(test_seq,threshold=-100)
	max_score=-100
	max_pos=0
	b=0
	for position, score in search_result:
	    if position==0:
		max_score=score
		max_pos=position
		b=1
	if b==1:	
	    peaks_output[0].append(peak)
	    peaks_output[1].append(np.mean(all_unique_peaks_with_height[peak]))
	    peaks_output[2].append(max_pos)
	    peaks_output[3].append(float(max_score))	    
    print 'cor coefficient ' + str(pearsonr(peaks_output[1],peaks_output[3]))
          
    
    test_seq = Seq(str(genomefa), IUPAC.unambiguous_dna)
    whole_genome_scores=pssm.calculate(test_seq)
    
    print 'mean height peaks ' + str(np.mean(peaks_output[1]))
    print 'mean motif score peaks ' + str(np.mean(peaks_output[3]))
    print 'mean motif score genome ' + str(np.mean(whole_genome_scores))
    
    plt.figure(figsize=(16, 8), dpi=100)
    f,axarr = plt.subplots(3,1)  
    ax=axarr[0]
    ends=[]
    ends.extend(peaks_output[1])
    ax.hist(ends,label="height")
    ax.legend()
        
    ax=axarr[1]
    ends=[]
    ends.extend(peaks_output[3])
    ax.hist(ends,label="peaks score")
    ax.set_xlim((-30,30))
    ax.legend()
    
    ax=axarr[2]
    ax.hist(whole_genome_scores,label="wholegenome score")
    ax.set_xlim((-30,30))
    ax.legend()
    
    plt.savefig(path+'height_motif_right_score_distributuion_all_peaks_2', dpi=100, figsize=(16, 8))

    outfile= open(path +"List_of_sites_right_all_peaks_2_calculate.txt".decode("utf-8"), 'w+')
    for i in range(0,len(whole_genome_scores)):
	outfile.write(str(whole_genome_scores[i]) + "\n")
    outfile.close()   
    
    test_seq = Seq(str(genomefa), IUPAC.unambiguous_dna).reverse_complement()
    whole_genome_scores=pssm.calculate(test_seq)   
    outfile= open(path +"List_of_sites_right_all_peaks_2_calculate_rc.txt".decode("utf-8"), 'w+')
    for i in range(0,len(whole_genome_scores)):
	outfile.write(str(whole_genome_scores[i]) + "\n")
    outfile.close()        
    
def filter_sites(threshold):
    sites_coords=open(path+ "List_of_sites_th_0.txt", 'r')
    outfile= open(path + "List_of_sites_th_6_8.txt".decode("utf-8"), 'w+')
    for line in sites_coords:
	line1=line.rstrip().split(' ')
	if float(line1[4]) >=threshold:
	    outfile.write(line)
    
    outfile.close()
    sites_coords.close()
	
def parse_list_of_sites_calculate(genome, infile,outfile):
    #sites_coords=open(path+ "List_of_sites_th_6_8.txt", 'r')
    
    #sites_coords=open(path+ "List_of_sites_th_0_all_peaks_2_calculate.txt", 'r')
    #outfile= open(path + "List_of_sites_th_0_all_peaks_2_calculate_edt.txt".decode("utf-8"), 'w+')
    
    sites_coords=open(path+ infile, 'r')
    outfile= open(path + outfile.decode("utf-8"), 'w+')
    
    print 'len genome ' + str(len(genome))
    #print 'len genomefa ' + str(len(genomefa))
    list_pos={}
    list_lines={}
    #list_lines
    i=0
    j=0
    for line in sites_coords:
	line=line.rstrip()
	line_ar=float(line)
	if line_ar<0 :
	    i+=1
	    continue
	outfile.write(str(i+win_width_l -1 ) + " + " + str(line_ar) + " "+ str(genome[i+win_width_l -1 -1 : i+win_width_l -1])+"\n")
	j+=1
	i+=1	
	
    print j
    sites_coords.close()
    outfile.close()

def parse_list_of_sites_left_calculate(genome, infile,outfile):
    #sites_coords=open(path+ "List_of_sites_th_0_all_peaks_2_calculate.txt", 'r')
    #outfile= open(path + "List_of_sites_th_0_all_peaks_2_calculate_edt.txt".decode("utf-8"), 'w+')
    
    sites_coords=open(path+ infile, 'r')
    outfile= open(path + outfile.decode("utf-8"), 'w+')
    
    print 'len genome ' + str(len(genome))
    list_pos={}
    list_lines={}
    i=0
    j=0
    for line in sites_coords:
	line=line.rstrip()
	line_ar=float(line)
	if line_ar<0 :
	    i+=1
	    continue
	outfile.write(str(i+win_width_l -1 ) + " + " + str(line_ar) + " "+ str(genome[i+win_width_l -1 -1 : i+win_width_l -1])+"\n")
	j+=1
	i+=1	
	
    print j
    sites_coords.close()
    outfile.close()
    
def parse_list_of_sites_right_calculate(genome, infile,outfile):
    #sites_coords=open(path+ "List_of_sites_th_0_all_peaks_2_calculate.txt", 'r')
    #outfile= open(path + "List_of_sites_th_0_all_peaks_2_calculate_edt.txt".decode("utf-8"), 'w+')
    sites_coords=open(path+ infile, 'r')
    outfile= open(path + outfile.decode("utf-8"), 'w+')
    
    print 'len genome ' + str(len(genome))
    list_pos={}
    list_lines={}
    #list_lines
    i=0
    j=0
    for line in sites_coords:
	line=line.rstrip()
	line_ar=float(line)
	if line_ar<0 :
	    i+=1
	    continue
	outfile.write(str(i+-4) + " + " + str(line_ar) + " "+ str(genome[i- 4-1 : i -4])+"\n")
	j+=1
	i+=1	
	
    print j
    sites_coords.close()
    outfile.close()
    
def parse_list_of_sites_calculate_rc(genome, infile,outfile):
    #sites_coords=open(path+ "List_of_sites_th_6_8.txt", 'r')
    #sites_coords=open(path+ "List_of_sites_th_0_all_peaks_2_calculate_rc.txt", 'r')
    #outfile= open(path + "List_of_sites_th_0_all_peaks_2_calculate_rc_edt.txt".decode("utf-8"), 'w+')
    #print 'len genomefa ' + str(len(str(genomefa)))
    #print 'len genomefa ' + str(len(genomefa))
    sites_coords=open(path+ infile, 'r')
    outfile= open(path + outfile.decode("utf-8"), 'w+')

    print 'len genome ' + str(len(genome))    
    
    list_pos={}
    list_lines={}
    #list_lines
    i=0
    
    output=[]
    for line in sites_coords:
	line=line.rstrip()
	line_ar=float(line)
	if line_ar<0 :
	    i+=1
	    continue
	output.append([str(len(genome)-i  -1 - win_width_l -1 -1)   , str(line_ar)] )
	i+=1	
    print len(output)
    for j in range( len(output)-1, -1 ,-1):
	outfile.write(output[j][0] +" - " +output[j][1]+ " "+str(genome[int(output[j][0]) -1 : int(output[j][0])])+"\n")
    sites_coords.close()
    outfile.close()

def parse_list_of_sites():
    #sites_coords=open(path+ "List_of_sites_th_6_8.txt", 'r')
    sites_coords=open(path+ "List_of_sites_th_0.txt", 'r')
    outfile= open(path + "List_of_sites_merged_pos_th_0.txt".decode("utf-8"), 'w+')
    print 'len genomefa ' + str(len(str(genomefa)))
    print 'len genomefa ' + str(len(genomefa))
    list_pos={}
    list_lines={}
    #list_lines

    for line in sites_coords:
	line=line.rstrip()
	line_ar=line.split(' ')
	if len(line_ar)<3:
	    break
	if int(line_ar[1][:-1])<0:
	    pos=len(genomefa)+int(line_ar[1][:-1])+win_width_r
	    strand='-'
	else:
	    pos=int(line_ar[1][:-1])+win_width_l
	    strand='+'
		
	list_lines[pos]=str(pos) + ' ' +strand + ' ' + line_ar[4]
	#print pos
	#print str(pos) + ' ' +strand + ' ' + line_ar[4]
	if pos in list_pos:
	    if list_pos[pos] < float(line_ar[4]):
		list_pos[pos]=float(line_ar[4])	
	else:
	    list_pos[pos]=float(line_ar[4]) 		
	'''
	elif pos+1 in list_pos:
	    if list_pos[pos+1] < float(line_ar[4]):
		list_pos[pos]=float(line_ar[4])  
		del list_pos[pos+1]
	elif pos+2 in list_pos:
	    if list_pos[pos+2] < float(line_ar[4]):
		list_pos[pos]=float(line_ar[4])  
		del list_pos[pos+2]	    
	elif pos+3 in list_pos:
	    if list_pos[pos+3] < float(line_ar[4]):
		list_pos[pos]=float(line_ar[4])  
		del list_pos[pos+3]	    
	elif pos-1 in list_pos:
	    if list_pos[pos-1] < float(line_ar[4]):
		list_pos[pos]=float(line_ar[4])  
		del list_pos[pos-1]	    
	elif pos-2 in list_pos: 
	    if list_pos[pos-2] < float(line_ar[4]):
		list_pos[pos]=float(line_ar[4])  
		del list_pos[pos-2]	    
	elif pos-3 in list_pos:
	    if list_pos[pos-3] < float(line_ar[4]):
		list_pos[pos]=float(line_ar[4])  
		del list_pos[pos-3]
	'''

	    

    print len(list_pos)
    print len(list_lines) 
    for pos in sorted(list_pos.iterkeys()):
    #for pos in list_pos:
	outfile.write(list_lines[pos] + '\n')
    
	
	#if int(line_ar[1][:-1]) <0:	    
	#    outfile.write(line_ar[0] + ' ' + str(len(genomefa)+int(line_ar[1][:-1])) + ' ' + line_ar[1] + ' ' +line_ar[2] +' ' + line_ar[3] + ' ' +line_ar[4]+"\n")
	#else:
	#    outfile.write(line + "\n")
    outfile.close()

    
def site_distribution():
    #peaks_coords=open(path+ "List_of_peaks.txt", 'r')
    sites_coords=open(path+"List_of_sites_merged_pos.txt", 'r')
    fig=plt.figure(figsize=(20,3))  
    f = plt.subplot()  
    seqs=[]
    i=0;
    sites=[]
    for line in sites_coords:	
	line=line.rstrip()
	line=line.split(' ')
	if float(line[2])>=12:
	    pos=int(line[0])
	    sites.append(pos)
	    seq = genomefa[pos - win_width_l + 1:pos + win_width_r + 1]
	    seqs.append(seq)	
    print len(sites)           
    create_motif_plot("motif_from_sites_12.png", seqs)
    
    bins=np.linspace(0, 4647454, 20)
    f.set_xlim(0, 4647454)

    ticks1=[0, 500000, 1000000, 1500000, 2000000, 2500000, 3000000, 3500000, 4000000, 4500000]
    xticknames1=['', '500kb', '1000kb', '1500kb', '2000kb', '2500kb', '3000kb', '3500kb', '4000kb', '4500kb']
    ticks2=[223771, 2729394, 3428115, 3469603, 3600731, 3694454, 4212776]
    xticknames2=['rRNA H', 'rRNA G ', 'rRNA D', 'rRNA B', 'rRNA A', 'rRNA C', 'rRNA E']

    f.set_xticks(ticks1, minor=False)
    f.set_xticklabels(xticknames1, rotation=0, fontsize=10)
    plt.setp(f.set_xticklabels(xticknames1), rotation=0, fontsize=10)

    f.set_xticks(ticks2, minor=True)
    f.set_xticklabels(xticknames2, minor=True)
    plt.setp(f.set_xticklabels(xticknames2, minor=True), rotation=90, fontsize=10)
    f.tick_params(direction='out', pad=18, which='minor')
    f.xaxis.grid(True, which='minor', linewidth=1, linestyle='--', alpha=0.4)
    f.set_title('Sites distribution', fontsize=9)
    
    f.tick_params(axis='both', which='major', labelsize=10)

    xori_ter=0
    yter=1400000
    yori=3711828
    f.hist(sites, bins, alpha=0.3)
    f.plot(yter, xori_ter, 'ro', markersize=5)
    f.plot(yori, xori_ter, 'go', markersize=5)

    plt.tight_layout()
    plt.savefig(path+"Sites_distribution_12.png", dpi=600)  
    plt.close()
    sites_coords.close()
def test_motif_search():
    instances = [Seq("ATGC")]
    m = motifs.create(instances)
    m = motifs.create(instances)
    pwm = m.counts.normalize()
    pssm = pwm.log_odds()    
    test_seq = Seq("GCATCAATGCCCATACGCGTA", IUPAC.unambiguous_dna)
    search_result=pssm.search(test_seq)
    for position, score in search_result:
	print("%d     %5.3f" % (position, score))   
	
    print 'rc'
    test_seq = Seq("GCATCAATGCCCATACGCGTA", IUPAC.unambiguous_dna).reverse_complement()
    search_result=pssm.search(test_seq)
    for position, score in search_result:
	print("%d     %5.3f" % (position, score))    

def calculate_counts_for_TU_plots(sample, peaks,genes):
    quarter=0
    gene_count_1=0
    gene_count_2=0
    gene_count_3=0
    gene_count_4=0
    upstream_count=0
    downstream_count=0
    count_pos_or=0
    count_neg_or=0
    for peak in sorted(peaks.iterkeys()):
	found=0
	beg_prev=0
	end_prev=0
	orient_prev='+'
	check_overlap=0
	for gene in sorted(genes.iterkeys()):
	    beg=gene
	    end=genes[gene][0]
	    orient=genes[gene][1]
	    
	    if check_overlap==1 and orient==orient_prev and beg<=end_prev:
		if peak >=beg and peak <= end:
		    quarter=(end-beg)/4
		    if orient=='+':
			count_pos_or+=1
			if peak>=beg and peak <= beg+quarter:
			    gene_count_1+=1
			    found=1
			    break
			if peak>beg +quarter and peak <= beg+2 * quarter:
			    gene_count_2+=1
			    found=1
			    break
			if peak>beg + 2 * quarter and peak <= beg+3 * quarter:
			    gene_count_3+=1
			    found=1
			    break
			if peak>beg + 3*quarter and peak <= end:
			    gene_count_4+=1
			    found=1
			    break
		    if orient=='-':
			count_neg_or+=1
			if peak>=beg and peak <= beg+quarter:
			    gene_count_4+=1
			    found=1
			    break
			if peak>beg +quarter and peak <= beg+2 * quarter:
			    gene_count_3+=1
			    found=1
			    break
			if peak>beg + 2 * quarter and peak <= beg+3 * quarter:
			    gene_count_2+=1
			    found=1
			    break
			if peak>beg + 3*quarter and peak <= end:
			    gene_count_1+=1
			    found=1
			    break
	    else:		
		if peak >=beg and peak <= end:
		    beg_prev=beg
		    end_prev=end
		    orient_prev=orient		    
		    quarter=(end-beg)/4
		    if orient=='+':
			count_pos_or+=1
			if peak>=beg and peak <= beg+quarter:
			    gene_count_1+=1
			    found=1
			    check_overlap=1
			    continue
			if peak>beg +quarter and peak <= beg+2 * quarter:
			    gene_count_2+=1
			    found=1
			    check_overlap=1
			    continue			    
			if peak>beg + 2 * quarter and peak <= beg+3 * quarter:
			    gene_count_3+=1
			    found=1
			    check_overlap=1
			    continue			    
			if peak>beg + 3*quarter and peak <= end:
			    gene_count_4+=1
			    found=1
			    check_overlap=1
			    continue			    
		    if orient=='-':
			count_neg_or+=1
			if peak>=beg and peak <= beg+quarter:
			    gene_count_4+=1
			    found=1
			    check_overlap=1
			    continue			    
			if peak>beg +quarter and peak <= beg+2 * quarter:
			    gene_count_3+=1
			    found=1
			    check_overlap=1
			    continue			    
			if peak>beg + 2 * quarter and peak <= beg+3 * quarter:
			    gene_count_2+=1
			    found=1
			    check_overlap=1
			    continue			    
			if peak>beg + 3*quarter and peak <= end:
			    gene_count_1+=1
			    found=1
			    check_overlap=1
			    continue			    		    
		if peak > end:
		    beg_prev=beg
		    end_prev=end
		    orient_prev=orient
		    continue
		if peak < beg and peak > end_prev:
		    if peak > end_prev and peak <= (end_prev+ (beg-end_prev)/2):
			if orient_prev=="+":
			    count_pos_or+=1
			    downstream_count+=1
			    found=1
			    break
			if orient_prev=="-":
			    count_neg_or+=1
			    upstream_count+=1
			    found=1
			    break	
		    if peak>(end_prev+ (beg-end_prev)/2) and peak < beg:
			if orient_prev=="+":
			    count_pos_or+=1
			    upstream_count+=1
			    found=1
			    break
			if orient_prev=="-":
			    count_neg_or+=1
			    downstream_count+=1
			    found=1
			    break
	    if found==0:
		print "ACHTUNG!!!"
		return
	    
    sum=upstream_count+downstream_count+gene_count_1+gene_count_2+gene_count_3+gene_count_4
    print sample + "\t" + str(upstream_count) + "\t" + str(gene_count_1)+ "\t" + str(gene_count_2)+ "\t" + str(gene_count_3)+ "\t" + str(gene_count_4)+ "\t" + str(downstream_count)+ "\t" + str(sum) + "-"+str(len(peaks))
    #print count_pos_or 
    #print count_neg_or
    return [upstream_count,gene_count_1,gene_count_2,gene_count_3,gene_count_4,downstream_count]
	

def calculate_counts_for_TU_plots_wide_intergenic_regions(sample, peaks,genes,flank):
    quarter=0
    gene_count_1=0 #1st quarter
    gene_count_2=0 #2nd quarter
    gene_count_3=0 #3rd quarter
    gene_count_4=0 #4th quarter
    upstream_count=0
    downstream_count=0
    count_pos_or=0 
    count_neg_or=0
    for peak in sorted(peaks.iterkeys()): #peaks loop
	found=0
	beg_prev=0
	end_prev=0
	orient_prev='+'
	check_overlap=0
	for gene in sorted(genes.iterkeys()): #genes loop
	    beg=gene #beginnig of the gene
	    end=genes[gene][0] #end of the gene
	    orient=genes[gene][1]     #orientation
	    #wndw= 1 * (end-beg)
	    #wndw=500
	    wndw=flank* (end-beg) #flank
	    if peak >=beg and peak <= end: #peak inside gene
		beg_prev=beg
		end_prev=end
		orient_prev=orient		    
		quarter=(end-beg)/4
		if orient=='+':
		    count_pos_or+=1
		    if peak>=beg and peak <= beg+quarter:
			gene_count_1+=1
			found=1
			check_overlap=1
			continue
		    if peak>beg +quarter and peak <= beg+2 * quarter:
			gene_count_2+=1
			found=1
			check_overlap=1
			continue			    
		    if peak>beg + 2 * quarter and peak <= beg+3 * quarter:
			gene_count_3+=1
			found=1
			check_overlap=1
			continue			    
		    if peak>beg + 3*quarter and peak <= end:
			gene_count_4+=1
			found=1
			check_overlap=1
			continue			    
		if orient=='-':
		    count_neg_or+=1
		    if peak>=beg and peak <= beg+quarter:
			gene_count_4+=1
			found=1
			check_overlap=1
			continue			    
		    if peak>beg +quarter and peak <= beg+2 * quarter:
			gene_count_3+=1
			found=1
			check_overlap=1
			continue			    
		    if peak>beg + 2 * quarter and peak <= beg+3 * quarter:
			gene_count_2+=1
			found=1
			check_overlap=1
			continue			    
		    if peak>beg + 3*quarter and peak <= end:
			gene_count_1+=1
			found=1
			check_overlap=1
			continue			    		    
		
	    if peak > end and peak < end + wndw :
		if orient=="+":
		    count_pos_or+=1
		    downstream_count+=1
		    found=1		
		    continue
		if orient=="-":
		    count_neg_or+=1
		    upstream_count+=1
		    found=1		
		    continue
	    

	    if peak < beg and  peak > beg -wndw: 
		if orient=="+":
		    count_pos_or+=1
		    upstream_count+=1
		    found=1
		    continue
		if orient=="-":
		    count_neg_or+=1
		    downstream_count+=1
		    found=1
		    continue		    
	    if found==1 and peak <= beg - wndw:
		break
	    if found==0:
		#print "ACHTUNG!!!"
		continue
	    
    sum=upstream_count+downstream_count+gene_count_1+gene_count_2+gene_count_3+gene_count_4
    #print sample + "\t" + str(upstream_count) + "\t" + str(gene_count_1)+ "\t" + str(gene_count_2)+ "\t" + str(gene_count_3)+ "\t" + str(gene_count_4)+ "\t" + str(downstream_count)+ "\t" + str(sum) + "-"+str(len(peaks))
    #print count_pos_or 
    #print count_neg_or
    return upstream_count,gene_count_1,gene_count_2,gene_count_3,gene_count_4,downstream_count

def september_get_nearby_genes():
    genes_file="/data/Gyrase/Pipe/Escherichia_coli_K_12_w3110_Mu_genes_annotation.gff"
    peaks_coords=open(path+ "List_of_peaks_MICROCIN_CFX_OXO_2.txt", 'r')
    genes=parse_genes_file(genes_file)
    

    for line in peaks_coords:
        peaks = {}
        line = line.split('\t')
        sample = str(line[0])
        print sample
        heights=[]
        for peak in line[1:len(line)]:
            peak = peak.split('-')
            peaks[int(peak[0])]= float(peak[1])
	    heights.append(float(peak[1]))
	    
	sorted_heights=sorted(heights)
	threshold=sorted_heights[len(sorted_heights)-50]
	peaks_data={}
	j=0;
	for peak in peaks:
	    if peaks[peak]>=threshold:
		
		peaks_data[peak]=[[],[]]
		sorted_genes_keys=sorted(list(genes.keys()))
		for i in range(0,len( sorted(genes.iterkeys()))):
		    if sorted_genes_keys[i]<=peak and genes[sorted_genes_keys[i]][0]>=peak:
		        prev_gene=str(sorted_genes_keys[i-1]) + " " + str(genes[sorted_genes_keys[i-1]][0])+ " " + str(genes[sorted_genes_keys[i-1]][1])+ " " + str(genes[sorted_genes_keys[i-1]][2].split(';')[2])
		        cur_gene=str(sorted_genes_keys[i]) + " " + str(genes[sorted_genes_keys[i]][0])+ " " + str(genes[sorted_genes_keys[i]][1])+ " " + str(genes[sorted_genes_keys[i]][2].split(';')[2])
		        next_gene=str(sorted_genes_keys[i+1]) + " " + str(genes[sorted_genes_keys[i+1]][0])+ " " + str(genes[sorted_genes_keys[i+1]][1])+ " " + str(genes[sorted_genes_keys[i+1]][2].split(';')[2])
		        peaks_data[peak][0].append(prev_gene+ "/" + cur_gene + "/" + next_gene)
		    if i+1 < len(genes) and genes[sorted_genes_keys[i]][0]<=peak and sorted_genes_keys[i+1]>=peak :
			end=genes[sorted_genes_keys[i]][0]
			beg=sorted_genes_keys[i+1]
			if peak < end + (beg-end)/2:
			    cur_gene=str(sorted_genes_keys[i]) + " " + str(genes[sorted_genes_keys[i]][0])+ " " + str(genes[sorted_genes_keys[i]][1])+ " " + str(genes[sorted_genes_keys[i]][2].split(';')[2])
			    next_gene=str(sorted_genes_keys[i+1]) + " " + str(genes[sorted_genes_keys[i+1]][0])+ " " + str(genes[sorted_genes_keys[i+1]][1])+ " " + str(genes[sorted_genes_keys[i+1]][2].split(';')[2])
			    peaks_data[peak][1].append(cur_gene + "/" + next_gene)
			else:
			    cur_gene=str(sorted_genes_keys[i]) + " " + str(genes[sorted_genes_keys[i]][0])+ " " + str(genes[sorted_genes_keys[i]][1])+ " " + str(genes[sorted_genes_keys[i]][2].split(';')[2])
			    next_gene=str(sorted_genes_keys[i+1]) + " " + str(genes[sorted_genes_keys[i+1]][0])+ " " + str(genes[sorted_genes_keys[i+1]][1])+ " " + str(genes[sorted_genes_keys[i+1]][2].split(';')[2])
			    peaks_data[peak][1].append(cur_gene + "/" + next_gene)
		    if sorted_genes_keys[i] > peak:
			break
			    
	for peak in peaks_data:
	    print str(peak) + "\t" + str(peaks[peak])+"\t"+ str(peaks_data[peak][0]) +"\t"+ str(peaks_data[peak][1])
				    
    peaks_coords.close()
    
def TU(filename):
    genes_file="/data/Gyrase/Pipe/Escherichia_coli_K_12_w3110_Mu_genes_annotation.gff"
    operons_file="/data/Gyrase/Pipe/K12_w3110_operons.gff"
    #peaks_coords=open(path+ "List_of_peaks_MICROCIN_CFX_OXO_2.txt", 'r')
    peaks_coords=open(filename, 'r')
    #peaks_coords=open(path+ "List_of_peaks_microcin_1.txt", 'r')
    #sites_coords=open(path + "List_of_sites_merged_pos_th_0.txt", 'r')
    sites_coords=open(path + "List_of_sites_all_peaks_3_calculate_edt.txt", 'r')
    
    genes=parse_genes_file(genes_file)
    operons=parse_operons_file(operons_file)    
    size_of_genes=0
    size_between_genes=0
    
    size_of_operons=0
    size_between_operons=0
    
    end=4647455
    
    client = Client()
    dview = client[:]
    dview.block=False
    view = client.load_balanced_view()
    view.block=True
    #print client.ids
    
    for gene in sorted(genes.iterkeys()):
	size_of_genes=size_of_genes + genes[gene][0]-gene
    size_between_genes=end-size_of_genes
    print 'size of genes ' + str(float(size_of_genes)/(4*end)) + " size intergenic " + str (float(size_between_genes)/(end*2))	
	
    
    for operon in sorted(operons.iterkeys()):
	size_of_operons=size_of_operons + operons[operon][0]-operon
    size_between_operons=end-size_of_operons
    print 'size of operons ' + str(float(size_of_operons)/(4*end)) + " size between operons " + str (float(size_between_operons)/(end*2))	

    
    genes_fraction=[0.0659,0.217,0.217,0.217,0.217,0.0659]
    operons_fraction=[0.2201, 0.1399,  0.1399,  0.1399,  0.1399, 0.2201]	
    

    
    sites = {}
    pks={}
    #print "name\tupstream_count\tgene_count_1\tgene_count_2\tgene_count_3\tgene_count_4\tdownstream_count\tsum_counts-number_of_peaks"        
    for line in sites_coords:	
	line=line.rstrip()
	line = line.split(' ')
	score=float((line[2]))
	if score>=0 and int (line[0]) > 0 and int(line[0])<5000000:
	    sites[int(line[0])]=score
	
    for line in peaks_coords:
	peaks = {}
	line=line.rstrip()
	line = line.split('\t')
	sample = str(line[0])
	#print sample
	for peak in line[1:len(line)]:
	    peak = peak.split('-')
	    if int(peak[0])> 0 and int(peak[0])<5000000:
		peaks[int(peak[0])]=float((peak[1]))
	pks[sample]=peaks
		
    wndws=[0.1, 0.25, 0.5, 1,2,3,4,5,6]
    #wndws=[0.1, 0.25]
    
    op_m_values=[None] * len (wndws)
    gene_m_values=[None] * len (wndws)
    op_c_values=[None] * len (wndws)
    gene_c_values=[None] * len (wndws)
    op_o_values=[None] * len (wndws)
    gene_o_values=[None] * len (wndws)      
    op_a_values=[None] * len (wndws)
    gene_a_values=[None] * len (wndws)  
    
    op_s_values=[None] * len (wndws)
    gene_s_values=[None] * len (wndws)  
    
    hi_g_m=[None] * len (wndws)
    low_g_m=[None] * len (wndws)
    hi_o_m=[None] * len (wndws)
    low_o_m=[None] * len (wndws)  
    
    hi_g_c=[None] * len (wndws)
    low_g_c=[None] * len (wndws)
    hi_o_c=[None] * len (wndws)
    low_o_c=[None] * len (wndws) 
    
    hi_g_o=[None] * len (wndws)
    low_g_o=[None] * len (wndws)
    hi_o_o=[None] * len (wndws)
    low_o_o=[None] * len (wndws)     
    
    hi_g_a=[None] * len (wndws)
    low_g_a=[None] * len (wndws)
    hi_o_a=[None] * len (wndws)
    low_o_a=[None] * len (wndws)     
    
    
    for wndw in wndws: 
	print 'window ' + str(wndw)    
	prefix="_" + str(wndw)	
	for sample in pks:
	    peaks=pks[sample]
	    raw_counts=calculate_counts_for_TU_plots_wide_intergenic_regions(sample + "_peaks_genes", peaks,genes,wndw)
	    #counts = np.divide(raw_counts, genes_fraction)
	    counts=raw_counts
	    sum=counts[0]+counts[1]+counts[2]+counts[3]+counts[4]+counts[5]
	    print sample + "_peaks_genes" + "\t" + str(counts[0]) + "\t" + str(counts[1])+ "\t" + str(counts[2])+ "\t" + str(counts[3])+ "\t" + str(counts[4])+ "\t" + str(counts[5])+ "\t" + str(sum) + "-"+str(len(peaks))
	    TU_plot(sample + "_peaks_genes_wide" +prefix ,counts)	    
	    print 'peaks genes upstream count / downstream count ' + str(float(counts[5])/counts[0])

	    if "MICROCIN" in sample:
		len_micro=len(peaks)
		gene_m_values[wndws.index(wndw)]=float(counts[5])/counts[0]
	    elif "CFX" in sample:
		len_cfx=len(peaks)
		gene_c_values[wndws.index(wndw)]=float(counts[5])/counts[0]
	    elif "OXO" in sample:
		len_oxo=len(peaks)
		gene_o_values[wndws.index(wndw)]=float(counts[5])/counts[0]	  
	    elif "all_peaks":
		gene_a_values[wndws.index(wndw)]=float(counts[5])/counts[0]	  		
		len_all=len(peaks)
		
	    raw_counts=calculate_counts_for_TU_plots_wide_intergenic_regions(sample+ "_peaks_operons", peaks,operons,wndw)
	    #counts = np.divide(raw_counts, operons_fraction)
	    counts=raw_counts
	    sum=counts[0]+counts[1]+counts[2]+counts[3]+counts[4]+counts[5]
	    print sample + "_peaks_operons" + "\t" + str(counts[0]) + "\t" + str(counts[1])+ "\t" + str(counts[2])+ "\t" + str(counts[3])+ "\t" + str(counts[4])+ "\t" + str(counts[5])+ "\t" + str(sum) + "-"+str(len(peaks))
	    TU_plot(sample + "_peaks_operons_wide" +prefix,counts)	    
	    print 'peaks operons upstream count / downstream count ' + str(float(counts[5])/counts[0])
	    
	    if "MICROCIN" in sample:
		op_m_values[wndws.index(wndw)]=float(counts[5])/counts[0]
	    elif "CFX" in sample:
		op_c_values[wndws.index(wndw)]=float(counts[5])/counts[0]
	    elif "OXO" in sample:
		op_o_values[wndws.index(wndw)]=float(counts[5])/counts[0]
	    elif "all_peaks":
		op_a_values[wndws.index(wndw)]=float(counts[5])/counts[0]	  		
		
		
	
	'''    	    
	raw_counts=calculate_counts_for_TU_plots_wide_intergenic_regions("sites_genes", sites,genes,wndw)
	counts = np.divide(raw_counts, genes_fraction)
	counts=raw_counts
	sum=counts[0]+counts[1]+counts[2]+counts[3]+counts[4]+counts[5]
	print "sites_genes" + "\t" + str(counts[0]) + "\t" + str(counts[1])+ "\t" + str(counts[2])+ "\t" + str(counts[3])+ "\t" + str(counts[4])+ "\t" + str(counts[5])+ "\t" + str(sum) + "-"+str(len(sites))
	print 'sites genes upstream count / downstream count ' + str(float(counts[0])/counts[5])
	TU_plot("sites_genes_wide" +prefix,counts)
	gene_s_values[wndws.index(wndw)]=float(counts[5])/counts[0]
	
	raw_counts=calculate_counts_for_TU_plots_wide_intergenic_regions("sites_operons", sites,operons,wndw)
	counts = np.divide(raw_counts, operons_fraction)
	counts=raw_counts
	sum=counts[0]+counts[1]+counts[2]+counts[3]+counts[4]+counts[5]
	print "sites_operons" + "\t" + str(counts[0]) + "\t" + str(counts[1])+ "\t" + str(counts[2])+ "\t" + str(counts[3])+ "\t" + str(counts[4])+ "\t" + str(counts[5])+ "\t" + str(sum) + "-"+str(len(sites))
	print 'sites operons upstream count / downstream count ' + str(float(counts[5])/counts[0])
	TU_plot("sites_operons_wide" +prefix,counts)
	op_s_values[wndws.index(wndw)]=float(counts[5])/counts[0]
	'''  
	
	
	
	data_g_m=[]
	data_o_m=[]
	t = time.time()	
	for j in range (0,50):
	    if j % 5 ==0:
		print str(j*2) + "%",	    
	    results=test(dview,len_micro,genes,operons,wndw)  
	    data_g_m.extend(results[0])
	    data_o_m.extend(results[1])
	t3 = time.time()-t    
	print('%f secs (multicore)' % t3)  	
	print "\nrandom mean upstream count/ downstream count genes " + str(np.mean(data_g_m)) + " sd " + str (np.std(data_g_m))
	print "random mean upstream count/ downstream count operons " + str(np.mean(data_o_m)) + " sd " + str (np.std(data_o_m)) + "\n"

	hi_g_m[wndws.index(wndw)]=np.mean(data_g_m) + 2* np.std(data_g_m)
	low_g_m[wndws.index(wndw)]=np.mean(data_g_m) - 2* np.std(data_g_m)
	hi_o_m[wndws.index(wndw)] = np.mean(data_o_m) +2 * np.std(data_o_m)
	low_o_m[wndws.index(wndw)] = np.mean(data_o_m) -2 * np.std(data_o_m)  	
    
	data_g_c=[]
	data_o_c=[]
	t = time.time()
	for j in range (0,50):
	    if j % 5 ==0:
		print str(j*2) + "%",	    
	    results=test(dview,len_cfx,genes,operons,wndw)  
	    data_g_c.extend(results[0])
	    data_o_c.extend(results[1])
	t3 = time.time()-t    
	print('%f secs (multicore)' % t3)  	
	print "\nrandom mean upstream count/ downstream count genes " + str(np.mean(data_g_c)) + " sd " + str (np.std(data_g_c))
	print "random mean upstream count/ downstream count operons " + str(np.mean(data_o_c)) + " sd " + str (np.std(data_o_c)) + "\n"

	hi_g_c[wndws.index(wndw)]=np.mean(data_g_c) + 2* np.std(data_g_c)
	low_g_c[wndws.index(wndw)]=np.mean(data_g_c) - 2* np.std(data_g_c)
	hi_o_c[wndws.index(wndw)] = np.mean(data_o_c) +2 * np.std(data_o_c)
	low_o_c[wndws.index(wndw)] = np.mean(data_o_c) -2 * np.std(data_o_c) 
	
	data_g_o=[]
	data_o_o=[]
	t = time.time()
	for j in range (0,50):
	    if j % 5 ==0:
		print str(j*2) + "%",	    
	    results=test(dview,len_oxo,genes,operons,wndw)  
	    data_g_o.extend(results[0])
	    data_o_o.extend(results[1])
	t3 = time.time()-t    
	print('%f secs (multicore)' % t3)  	
	print "\nrandom mean upstream count/ downstream count genes " + str(np.mean(data_g_o)) + " sd " + str (np.std(data_g_o))
	print "random mean upstream count/ downstream count operons " + str(np.mean(data_o_o)) + " sd " + str (np.std(data_o_o)) + "\n"

	hi_g_o[wndws.index(wndw)]=np.mean(data_g_o) + 2* np.std(data_g_o)
	low_g_o[wndws.index(wndw)]=np.mean(data_g_o) - 2* np.std(data_g_o)
	hi_o_o[wndws.index(wndw)] = np.mean(data_o_o) +2 * np.std(data_o_o)
	low_o_o[wndws.index(wndw)] = np.mean(data_o_o) -2 * np.std(data_o_o) 	
        
	data_g_a=[]
	data_o_a=[]
	t = time.time()
	for j in range (0,50):
	    if j % 5 ==0:
		print str(j*2) + "%",	    
	    results=test(dview,len_all,genes,operons,wndw)  
	    data_g_a.extend(results[0])
	    data_o_a.extend(results[1])
	t3 = time.time()-t    
	print('%f secs (multicore)' % t3)  	
	print "\nrandom mean upstream count/ downstream count genes " + str(np.mean(data_g_a)) + " sd " + str (np.std(data_g_a))
	print "random mean upstream count/ downstream count operons " + str(np.mean(data_o_a)) + " sd " + str (np.std(data_o_a)) + "\n"

	hi_g_a[wndws.index(wndw)]=np.mean(data_g_a) + 2* np.std(data_g_a)
	low_g_a[wndws.index(wndw)]=np.mean(data_g_a) - 2* np.std(data_g_a)
	hi_o_a[wndws.index(wndw)] = np.mean(data_o_a) +2 * np.std(data_o_a)
	low_o_a[wndws.index(wndw)] = np.mean(data_o_a) -2 * np.std(data_o_a)         
    
    plt.figure(figsize=(16, 8), dpi=100)
    ax = plt.subplot()
    ax.plot(wndws,gene_m_values, label='genes', color='blue')	
    ax.plot(wndws,op_m_values, label='operons', color='red')    
    ax.plot(wndws,low_g_m, '--', label='genes intrerval', color='green')	
    ax.plot(wndws,hi_g_m, '--', color='green')    
    ax.plot(wndws,low_o_m, '--', label='operons interval', color='grey')	
    ax.plot(wndws,hi_o_m, '--',  color='grey')     
    plt.title("distribution  micro" + "")
    ax.set_xticks(wndws)
    ax.legend()
    plt.savefig('distribution_micro.png', dpi=100, figsize=(16, 8))
    plt.close()
    
    
    plt.figure(figsize=(16, 8), dpi=100)
    ax = plt.subplot()  
    ax.plot(wndws,gene_c_values, label='genes', color='blue')   
    ax.plot(wndws,op_c_values, label='operons', color='red')	
    ax.plot(wndws,low_o_c, '--', label='operons interval', color='grey')	
    ax.plot(wndws,hi_o_c, '--', color='grey') 
    ax.plot(wndws,low_g_c, '--', label='genes intrerval', color='green')	
    ax.plot(wndws,hi_g_c, '--', color='green')     
    plt.title("distribution cfx" + "")
    ax.set_xticks(wndws)   
    plt.legend()
    plt.savefig('distribution_cfx.png', dpi=100, figsize=(16, 8))
    plt.close()
    
    plt.figure(figsize=(16, 8), dpi=100)
    ax = plt.subplot()  
    ax.plot(wndws,gene_o_values, label='genes', color='blue')   
    ax.plot(wndws,op_o_values, label='operons', color='red')	
    ax.plot(wndws,low_o_o, '--', label='operons interval', color='grey')	
    ax.plot(wndws,hi_o_o, '--', color='grey') 
    ax.plot(wndws,low_g_o, '--', label='genes intrerval', color='green')	
    ax.plot(wndws,hi_g_o, '--', color='green')     
    plt.title("distribution oxo" + "")
    ax.set_xticks(wndws)   
    plt.legend()
    plt.savefig('distribution_oxo.png', dpi=100, figsize=(16, 8))
    plt.close()    
    
    plt.figure(figsize=(16, 8), dpi=100)
    ax = plt.subplot()  
    ax.plot(wndws,gene_a_values, label='genes', color='blue')   
    ax.plot(wndws,op_a_values, label='operons', color='red')	
    ax.plot(wndws,low_o_a, '--', label='operons interval', color='grey')	
    ax.plot(wndws,hi_o_a, '--', color='grey') 
    ax.plot(wndws,low_g_a, '--', label='genes intrerval', color='green')	
    ax.plot(wndws,hi_g_a, '--', color='green')     
    plt.title("distribution all peaks 3" + "")
    ax.set_xticks(wndws)   
    plt.legend()
    plt.savefig('distribution_all_peaks_3.png', dpi=100, figsize=(16, 8))
    plt.close()        
    
	
    '''
    #
    plt.figure(figsize=(16, 8), dpi=100)
    ax = plt.subplot()    
    plt.plot(wndws,gene_s_values, label='genes', color='blue')
    plt.plot(wndws,op_s_values, label='operons', color='red')
    plt.title("distribution sites" + "")
    ax.set_xticks(wndws)   
    plt.legend()
    plt.savefig('distribution_sites.png', dpi=100, figsize=(16, 8))
    plt.close()    
    '''
    
    sites_coords.close()
    peaks_coords.close()
	

def test ( view,number_of_peaks,genes,operons,wndw):
    num=20
    rand_peaks=[None] *num
    raw_counts_g=[None] *num
    counts_g=[None] *num
    raw_counts_o=[None] *num
    counts_o=[None] *num
    data_g=[None] *num
    data_o=[None] *num
    ar=[None]*20
    ar2=[None]*20
    t = time.time()
    #pr_list = [view.apply_async(calculate_counts_for_TU_plots_wide_intergenic_regions,"random_"+str(i) + "_peaks_genes", rand_peaks[i],genes,wndw) for i in range(num)]
    #view.wait(pr_list)    
    t3 = time.time()-t    
    #print('%f secs (multicore)' % t3)    
    #print pr_list
    
    t = time.time()
    for i in range (0, 20):
	view.targets = [i]
	view.block = False
	rand_peaks[i]=generate_random_peaks(number_of_peaks)
	ar[i] = view.apply(calculate_counts_for_TU_plots_wide_intergenic_regions,"random_"+str(i) + "_peaks_genes", rand_peaks[i],genes,wndw)	
	ar2[i]=view.apply(calculate_counts_for_TU_plots_wide_intergenic_regions,"random_"+str(i)+ "_peaks_operons", rand_peaks[i],operons,wndw)
	#if ar.ready() :
	    #ar.get()	
    view.wait(ar)
    for i in range(20):
	raw_counts_g[i]=ar[i].get()[0]
    view.wait(ar2)
    for i in range(20):
	raw_counts_o[i]=ar2[i].get()[0]    
    
    #print
    t3 = time.time()-t    
    #print('%f secs (multicore)' % t3)  
    

    
    
    #for i in range (0,num):	
	#rand_peaks[i]=generate_random_peaks(len (peaks))	
	#raw_counts_g[i]= view.apply (calculate_counts_for_TU_plots_wide_intergenic_regions,"random_"+str(i) + "_peaks_genes", rand_peaks[i],genes,wndw)
	#raw_counts_o[i]= view.apply(calculate_counts_for_TU_plots_wide_intergenic_regions,"random_"+str(i)+ "_peaks_operons", rand_peaks[i],operons,wndw)
	
    for i in range (0,num):
	counts_g[i]= raw_counts_g[i]
	#print counts_g
	data_g[i]=float(counts_g[i][5])/counts_g[i][0]	
	counts_o[i]=raw_counts_o[i]		
	data_o[i]=float(counts_o[i][5])/counts_o[i][0]
	

    return data_g, data_o
def generate_random_peaks(number):
    peaks={}
    for i in range (0,number):
	peaks[random.randint(0, len(genomefa)-1)]=1
	
    return peaks
def TU_plot(name,counts):
    mean_int=(float)(counts[0]+counts[5])/2
    mean_gene=(float) (counts[1]+counts[2]+counts[3]+counts[4])/4

    plt.figure(figsize=(16, 8), dpi=100)
    plot1 = plt.subplot()
    plot1.set_xlim(-10, 110)
    plot1.set_xticks(np.arange(0, 110, 10))
    xticks = plot1.xaxis.get_major_ticks()
    plot1.set_xlabel('Relative position of gyrase sites, %', size=17)
    plot1.set_ylabel('Number of gyrase sites', size=17)
    #4 bars
    y_pos=np.arange(0,4)    
    plt.bar(y_pos*25, counts[1:5], 25,align='edge',alpha=0.7,color='b', edgecolor = "b" )
    plt.bar(-10, counts[0],10,alpha=0.7,color='red', edgecolor = "b" )
    plt.bar(100, counts[5],10,alpha=0.7,color='g',  edgecolor = "b" )
    
    #statistics 1st bar
    plot2 = plt.subplot()
    interval = scipy.stats.poisson.interval(0.95, mean_int)
    plot2.plot([-7.5, -2.5], [mean_int, mean_int], color='black', lw=0.8) #mean line
    plot2.plot([-10, 0], [interval[0], interval[0]], color='black', lw=0.8) # horizontal line
    plot2.plot([-10, 0], [interval[1], interval[1]], color='black', lw=0.8) # horizontal line
    plot2.plot([-5, -5], [interval[0], interval[1]], color='black', lw=0.8)	# vertical line

    interval = scipy.stats.poisson.interval(0.95, mean_gene)
    plot2.plot([10, 15], [mean_gene, mean_gene], color='black', lw=0.8)
    plot2.plot([10, 15], [interval[0], interval[0]], color='black', lw=0.8)
    plot2.plot([10, 15], [interval[1], interval[1]], color='black', lw=0.8)
    plot2.plot([12.5, 12.5], [interval[0], interval[1]], color='black', lw=0.8)

    plot2.plot([10 + 25, 15 + 25], [mean_gene, mean_gene], color='black', lw=0.8)
    plot2.plot([10 + 25, 15 + 25], [interval[0], interval[0]], color='black', lw=0.8)
    plot2.plot([10 + 25, 15 + 25], [interval[1], interval[1]], color='black', lw=0.8)
    plot2.plot([12.5 + 25, 12.5 + 25], [interval[0], interval[1]], color='black', lw=0.8)

    plot2.plot([10 + 50, 15 + 50], [mean_gene, mean_gene], color='black', lw=0.8)
    plot2.plot([10 + 50, 15 + 50], [interval[0], interval[0]], color='black', lw=0.8)
    plot2.plot([10 + 50, 15 + 50], [interval[1], interval[1]], color='black', lw=0.8)
    plot2.plot([12.5 + 50, 12.5 + 50], [interval[0], interval[1]], color='black', lw=0.8)

    plot2.plot([10 + 75, 15 + 75], [mean_gene, mean_gene], color='black', lw=0.8)
    plot2.plot([10 + 75, 15 + 75], [interval[0], interval[0]], color='black', lw=0.8)
    plot2.plot([10 + 75, 15 + 75], [interval[1], interval[1]], color='black', lw=0.8)
    plot2.plot([12.5 + 75, 12.5 + 75], [interval[0], interval[1]], color='black', lw=0.8)

    interval = scipy.stats.poisson.interval(0.95, mean_int)
    plot2.plot([2.5 + 100, 7.5 + 100], [mean_int, mean_int], color='black', lw=0.8)
    plot2.plot([2.5 + 100, 7.5 + 100], [interval[0], interval[0]], color='black', lw=0.8)
    plot2.plot([2.5 + 100, 7.5 + 100], [interval[1], interval[1]], color='black', lw=0.8)
    plot2.plot([5 + 100, 5 + 100], [interval[0], interval[1]], color='black', lw=0.8)

    #xticks[0].label1.set_visible(False)
    #plt.show()
    plt.title(name)
    plt.savefig(path+name +".png", dpi=100, figsize=(16, 8))
    plt.close()
    
def parse_genes_file(file_name):
    genes_file=open(file_name, 'r')
    genes = {}
    err_count=0
    for line in genes_file:
	if "NC_007779.1_w3110_Mu	ena	transcript" not in line:
	    continue
	line=line.rstrip()
	line = line.split('\t')
	beg = int(line[3])
	if beg not in genes:
	    genes[beg]=[int(line[4]),line[6], line[8]]
	#else:
	    #err_count+=1
	    #print "ACHTUNG 2"
	    #print beg
	    #return
    #
    beg_prev=-1
    end_prev=-1
    or_prev="+"
    mean_overlap=0
    for gene in sorted(genes.iterkeys()):
	if gene<end_prev and or_prev== genes[gene][1]:
	    err_count+=1
	    mean_overlap+=end_prev-gene
	    #print str(beg_prev) + " " + str(end_prev) + " " + str(or_prev) +"\n"
	    #print str(gene) + " " + str(genes[gene][0]) + " " + str(genes[gene][1] + "\n\n")
	beg_prev=gene
	end_prev=genes[gene][0]
	or_prev=genes[gene][1]
    print 'number of overlapping genes ' + str(err_count)
    print 'mean length of overlap ' + str(float(mean_overlap)/err_count)
    genes_file.close()
    return genes

def parse_operons_file(file_name):
    operons_file=open(file_name, 'r')
    operons = {}
    for line in operons_file:
	if "NC_007779.1_w3110_Mu	ena	gene" not in line:
	    continue
	line=line.rstrip()
	line = line.split('\t')
	beg = int(line[3])
	if beg not in operons:
	    operons[beg]=[int(line[4]),line[6]]
	else:
	    print "ACHTUNG 3"
	    return
    operons_file.close()
    return operons

def get_nearby_genes():
    genes_file="/data/Gyrase/Pipe/Escherichia_coli_K_12_w3110_Mu_genes_annotation.gff"
    operons_file="/data/Gyrase/Pipe/K12_w3110_operons_only.gff"
    genes=parse_genes_file(genes_file)
    operons=parse_operons_file(operons_file)    
    
    for gene in sorted(genes.iterkeys()):
	beg=gene
	end=genes[gene][0]
	orient=genes[gene][1]
	inf= genes[gene][2]
	
def intersect_peaks_sites():    
    peaks_coords=open(path+ "List_of_peaks_MICROCIN_CFX.txt", 'r')
    #peaks_coords=open(path+ "List_of_peaks_microcin_1.txt", 'r')
    sites_coords=open(path + "List_of_sites_merged_pos.txt", 'r') 
    for line in peaks_coords:
	peaks = {}
	line=line.rstrip()
	line = line.split('\t')
	sample = str(line[0])
	print sample
	for peak in line[1:len(line)]:
	    peak = peak.split('-')
	    if int(peak[0])> 0 and int(peak[0])<5000000:
		peaks[int(peak[0])]=float((peak[1]))    

		for line in sites_coords:	
		    line=line.rstrip()
		    line = line.split(' ')
		    score=float((line[2]))
		    if score>=0 and int (line[0]) > 0 and int(line[0])<5000000:
			sites[int(line[0])]=score
		    
		for line in peaks_coords:
		    peaks = {}
		    line=line.rstrip()
		    line = line.split('\t')
		    sample = str(line[0])
		    #print sample
		    for peak in line[1:len(line)]:
			peak = peak.split('-')
			if int(peak[0])> 0 and int(peak[0])<5000000:
			    peaks[int(peak[0])]=float((peak[1]))
		    pks[sample]=peaks

def intersect_lists_of_peaks():
    f_mitya=open("/data/Gyrase/Pipe/Fragments_ends/Data_second_mapping/Mitya/_SW_new/List_of_peaks.txt", 'r')
    f_our=open("/data/Gyrase/Pipe/Fragments_ends/Data_second_mapping/case_6_AC_SW_1/List_of_peaks.txt", 'r')

    mitya_peaks={}
    our_peaks={}
    
    for line in f_mitya:
	peaks=[]
	line=line.rstrip()
	
	line=line.split('\t')
	sample= str(line[0].split("_vs_")[0])
	#print sample
	for peak in line[1:len(line)]:
	    peak=peak.split('-')
	    peaks.append(int(peak[0]))
	#print len(peaks)  
	mitya_peaks[sample]=peaks
	
    for line in f_our:
	peaks=[]
	line=line.rstrip()
	line=line.split('\t')
	sample= str(line[0].split("_vs_")[0])
	#print sample
	for peak in line[1:len(line)]:
	    peak=peak.split('-')
	    peaks.append(int(peak[0]))
	#print len(peaks)  
	our_peaks[sample]=peaks    
	
    
    for sample in our_peaks:
	only_our=[]
	only_mitya=[]
	overlap=[]
	our=our_peaks[sample]
	mitya=mitya_peaks[sample]
	overlap_count=0
	for peak in our:
	    if peak in mitya:
		overlap_count+=1
		overlap.append(peak)
	    else:
		only_our.append(peak)
	for peak in mitya:
	    if peak not in our:
		only_mitya.append(peak)
	print sample + "\t" + str(len(our)) + "\t" + str(overlap_count) + "\t" + str (len(mitya))
	create_motif_plot_from_peaks(sample + "_only_our.png",only_our)
	create_motif_plot_from_peaks(sample + "_only_mitya.png",only_mitya)
	create_motif_plot_from_peaks(sample + "_overlap.png",overlap)
	
	
	

def intersect_lists_of_peaks_and_sites():
    #f_peaks=open("/data/Gyrase/Pipe/Fragments_ends/Data_second_mapping/Mitya/_SW_new/List_of_peaks.txt", 'r')
    f_peaks=open(path + "List_of_peaks_score_2.txt", 'r')
    f_sites=open(path+"List_of_sites_th_0_all_peaks_2_calculate_edt.txt", 'r')
    #f_sites=open("/data/Gyrase/Pipe/Fragments_ends/Data_first_mapping/case_6_AC_SW_add_1/List_of_sites_merged_pos.txt", 'r')

    peaks_dic={}
    all_peaks=[]
    all_peaks_height={}
    sites=[]
    for line in f_peaks:
	peaks=[]
	line=line.rstrip()	
	line=line.split('\t')
	sample= str(line[0].split("_vs_")[0])
	#print sample
	for peak in line[1:len(line)]:
	    peak=peak.split(':')
	    peaks.append(int(peak[0]))
	    if peak not in all_peaks:
		all_peaks.append(int(peak[0]))
		all_peaks_height[int(peak[0])]=float(peak[1])
		
	#print len(peaks)  
	peaks_dic[sample]=peaks
	
    for line in f_sites:	
	line=line.rstrip()
	line=line.split(" ")[0]
	#print line
	#line=line.split(' ')
	#sites.append(int(line[0]))   
	sites.append(int(line))   
    
    print "len all peaks " + str(len(all_peaks)) + " len all sites " + str (len(sites))
    overlap_count=0
    
    for peak in all_peaks:
	#if all_peaks.index(peak) % 1000 ==0 :
	    #print all_peaks.index(peak),
	if peak in sites:
	    overlap_count+=1
	    continue
	else:
	    print str(peak) + " "+str(all_peaks_height[peak])
	    continue
	
	for i in range (0,1) :
	    if peak + i in sites :
		overlap_count+=1
		break

    print " peaks " + str (len(all_peaks)) + " overlap " + str(overlap_count) + " sites " + str (len(sites))  
	

def rrna_operons():
    #f_peaks=open("/data/Gyrase/Pipe/Fragments_ends/Data_second_mapping/case_6_AC_SW_1/List_of_peaks.txt", 'r')
    f_peaks=open("/data/Gyrase/Pipe/Fragments_ends/Data_second_mapping/List_of_peaks.txt", 'r')
    #f_sites=open("/data/Gyrase/Pipe/Fragments_ends/Data_first_mapping/case_6_AC_SW_add_1/List_of_sites_th_0_edt.txt", 'r')
    f_sites=open("/data/Gyrase/Pipe/Fragments_ends/Data_first_mapping/case_6_AC_SW_add_1/List_of_sites_merged_pos.txt", 'r')

    fig=plt.figure(figsize=(20,3))  
    f,axarr = plt.subplots(10,1)

    peaks_dic={}
    sites={}
    for line in f_peaks:
	peaks={}
	line=line.rstrip()	
	line=line.split('\t')
	sample= str(line[0].split("_vs_")[0])
	#print sample
	for peak in line[1:len(line)]:
	    peak=peak.split('-')
	    peaks[int(peak[0])]=float(peak[1])  
		
	#print len(peaks)  
	peaks_dic[sample]=peaks
	
    for line in f_sites:	
	line=line.rstrip()
	#line=line.split(":")[0].split(" ")[1]
	#print line
	line=line.split(' ')
	sites[int(line[0])]=float(line[2])
	#sites.append(int(line))   
    peaks_dic["sites"]=sites	
    
    
    A=[3597167, 3602272] 
    B=[3466047, 3471144]
    C=[3690984, 3695995]
    D=[3424645, 3429656]
    E=[4212776, 4217871]
    G=[2725847, 2730930]
    H=[223771, 228875]

    operons=[[3597167, 3602272, '-'], [3466047, 3471144, '-'], [3690984, 3695995, '-'], [3424645, 3429656, '-'], [4212776, 4217871, '+'], [2725847, 2730930, '-'], [223771, 228875, '+']]

    width=5000
    operons_len=[]    
    for j in range(len(operons)):
	operons_len.append(operons[j][1]-operons[j][0])

    print operons_len
    k=0
    for sample in peaks_dic:    
	peaks=peaks_dic[sample]
	operons_stat=[]
	mean=[[],[],[] ]
	for j in range(len(operons)):
	    stats=[0, 0, 0]
	    
	    for peak in peaks:
		if peak>operons[j][0]-width:
		    if peak<operons[j][0]:
			if operons[j][2]=='-':
			    stats[2]=stats[2]+1
			    mean[2].append(peaks[peak])
			elif operons[j][2]=='+':
			    stats[0]=stats[0]+1 
			    mean[0].append(peaks[peak])
			    
		if peak>operons[j][1]:
		    if peak<operons[j][1]+width:
			if operons[j][2]=='-':
			    stats[0]=stats[0]+1
			    mean[0].append(peaks[peak])
			elif operons[j][2]=='+':
			    stats[2]=stats[2]+1 
			    mean[2].append(peaks[peak])
			    
		if peak>=operons[j][0]:
		    if peak<=operons[j][1]:
			stats[1]=stats[1]+1
			mean[1].append(peaks[peak])
	    operons_stat.append(stats)
	    
	print sample
	print operons_stat
	print str(np.mean(mean[0])) + " " +str(np.mean(mean[1])) + " " + str(np.mean(mean[2]))
	
	l=0
	
	if "Un" in sample:
	    continue	
	
	for operon in operons_stat:
	    #print range(l*10,(l+3)*10,10)
	    #print operon
	    ticks2=[10, 50, 90, 130, 170, 210, 250]
	    xticknames2=['rRNA H', 'rRNA G ', 'rRNA D', 'rRNA B', 'rRNA A', 'rRNA C', 'rRNA E']
	
	    
	    axarr[k].set_ylim(0, 50)
	    
	    axarr[k].set_xticks(ticks2)
	    axarr[k].set_xticklabels(xticknames2, fontsize=5)
	    #plt.setp(axarr[k].set_xticklabels(xticknames2, minor=True), rotation=0, fontsize=7)
	    
	   # axarr[k].set_title(sample, fontsize=5)	    
	    axarr[k].text(.5,.8, sample,
	           horizontalalignment='center',
	       transform=axarr[k].transAxes,fontsize=5)	   
	    axarr[k].xaxis.label.set_visible(False)
	    axarr[k].tick_params(axis='both', which='major', labelsize=5)
	    
	    #plot1.set_xticks(np.arange(-120, 112, 10))
	    #plot1.axis([-110, 110, 0.3, 1])	    
	    #axarr[k].tight_layout()
	    axarr[k].bar(range(l*10,(l+3)*10,10), operon, 5,align='edge',alpha=0.7,color='b', edgecolor = "b" )
	    l+=4
	k+=1
	#plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
	f.savefig(path+"Peak_distribution.png", dpi=600) 	
	
	'''
	operons_stat_5000_microcin_IP_Mu_10_mkM_1=[[0, 21], [0, 25], [1, 40], [2, 18], [0, 23], [1, 6], [0, 17]]
	operons_stat_5000_cfx_IP_Mu_10_mkM_1=[[6, 53], [2, 33], [2, 60], [20, 33], [20, 54], [2, 24], [0, 44]]
	operons_stat_5000_microcin_IN_Mu_10_mkM_1=[[0, 0], [1, 1], [2, 0], [2, 1], [2, 1], [1, 0], [0, 1]]
	operons_stat_5000_cfx_IN_Mu_10_mkM_1=[[0, 1], [0, 1], [0, 1], [0, 1], [0, 2], [0, 0], [0, 0]]
	operons_stat_5000_un_IN_Mu_1=[[1, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0]]
	operons_stat_5000_un_IP_Mu_1=[[1, 0], [1, 1], [3, 3], [0, 0], [0, 0], [0, 0], [0, 0]]
	'''
	

	
	'''
	Un_IP_Mu_1=[0, 0, 0, 0, 0, 1, 0]
	Un_IN_Mu_1=[0, 0, 0, 0, 0, 0, 0]
	Cfx_IN_Mu_1=[0, 0, 0, 0, 0, 0, 0]
	Cfx_IP_Mu_10_mkM_1=[5, 6, 4, 5, 7, 4, 5]
	Microcin_IP_Mu_1=[3, 5, 2, 1, 1, 2, 1]
	Microcin_IN_Mu_1=[1, 0, 0, 1, 0, 0, 0]
	operons_length=[5105, 5097, 5011, 5011, 5095, 5083, 5104]
	'''
    f_peaks.close()
    f_sites.close()

def create_summary_table():
    peaks_coords=open(path+ "List_of_peaks.txt", 'r')
    peaks_micro_mu_1={}
    peaks_micro_mu_3={}
    peaks_micro_4={}
    peaks_cfx_mu_1={}
    peaks_cfx_mu_3={}
    peaks_cfx_3={}
    peaks_oxo_mu_1={}
    peaks_oxo_1={}
    peaks_oxo_2={}
    
    for line in peaks_coords:
	peaks = {}
	line = line.split('\t')
	sample = str(line[0])
	#print sample
	for peak in line[1:len(line)]:
	    peak = peak.split('-')
	    peaks[int(peak[0])]=float(peak[1])	    
	if "Microcin_IP_Mu_1" in sample:
	    peaks_micro_mu_1 = peaks
	if "Microcin_IP_Mu_3" in sample:
	    peaks_micro_mu_3 = peaks
	if "Microcin_IP_50_mkM_4" in sample:
	    peaks_micro_4 = peaks	    
	if "Cfx_IP_Mu_10_mkM_1" in sample:
	    peaks_cfx_mu_1 = peaks		
	if "Cfx_IP_Mu_10_mkM_3" in sample:
	    peaks_cfx_mu_3 = peaks
	if "Cfx_IP_10_mkM_3" in sample:
	    peaks_cfx_3 = peaks	
	if "Oxo_IP_Mu_120_mkM_1" in sample:
	    peaks_oxo_mu_1 = peaks
	if "Oxo_IP_120_mkM_1" in sample:
	    peaks_oxo_1 = peaks
	if "Oxo_IP_120_mkM_2" in sample:
	    peaks_oxo_2 = peaks	
    peaks_coords.close()
    #peaks_coords=open(path+ "List_of_peaks_score_1_motif_2.txt", 'r')
    peaks_coords=open(path+ "List_of_peaks_all_motif_3.txt", 'r')
    all_peaks_score_1={}
    for line in peaks_coords:    
	    line=line.rstrip()
	    line = line.split('\t')
	    sample = str(line[0])
	    for peak in line[1:len(line)]:
		    peak = peak.split(':')    
		    all_peaks_score_1[int(peak[0])]=float((peak[1]))    
    peaks_coords.close()
    
    outfile=open(path + "Summary_table_edt.txt","w+")
    outfile.write("peak" + "\t" +"score" + "\t"+ "Cfx_mu_1" + "\t" +"Cfx_mu_3" + "\t" + "Cfx_3" + "\t" + "Micro_mu_1" + "\t" + "Micro_mu_3" + "\t" + "Micro_4"+ "\t" )
    outfile.write("Oxo_mu_1" + "\t" + "Oxo_1" + "\t" + "Oxo_2\n")
    na="0"
    for peak in sorted(all_peaks_score_1.iterkeys()):
	outfile.write(str(peak) + "\t" +str(all_peaks_score_1[peak] ))
	if peak in peaks_cfx_mu_1:
	    outfile.write("\t" +str(peaks_cfx_mu_1[peak] ))
	else:
	    outfile.write("\t" +na )
	if peak in peaks_cfx_mu_3:
	    outfile.write("\t" +str(peaks_cfx_mu_3[peak] ))
	else:
	    outfile.write("\t" +na )
	if peak in peaks_cfx_3:
	    outfile.write("\t" +str(peaks_cfx_3[peak] ))
	else:
	    outfile.write("\t"+na  )
	if peak in peaks_micro_mu_1:
	    outfile.write("\t" +str(peaks_micro_mu_1[peak] ))
	else:
	    outfile.write("\t" +na)
	if peak in peaks_micro_mu_3:
	    outfile.write("\t" +str(peaks_micro_mu_3[peak] ))
	else:
	    outfile.write("\t" +na )
	if peak in peaks_micro_4:
	    outfile.write("\t" +str(peaks_micro_4[peak] ))
	else:
	    outfile.write("\t"+na  )
	if peak in peaks_oxo_mu_1:
	    outfile.write("\t" +str(peaks_oxo_mu_1[peak] ))
	else:
	    outfile.write("\t" +na )
	if peak in peaks_oxo_1:
	    outfile.write("\t" +str(peaks_oxo_1[peak] ))
	else:
	    outfile.write("\t" +na )
	if peak in peaks_oxo_2:
	    outfile.write("\t" +str(peaks_oxo_2[peak]) + "\n")
	else:
	    outfile.write("\t"+na+ "\n" )
    
    outfile.close()
	
	
    
    
    
    
def calculate_correlation():
    
    outfile=open(path + "Summary_table.txt","r")
    
    ar=[]
    for i in range (0,9):
	ar.append([])
    
    print ar	
    for line in outfile:    
	line = line.rstrip().split("\t")
	if "peak" not in line and len(line)==11 and '' not in line:		
	    #print line
	    #print line[i+2]		
	    for i in range(0,9):
		ar[i].append(float(line[i+2]))
    
    for i in range(0,9):
	for j in range(0,9):
	    print str(pearsonr(ar[i],ar[j])[0]),
	print "\n"
    
    for i in range(0,9):
	for j in range(0,9):
	    print str(pearsonr(ar[i],ar[j])[1]),
	print "\n"
    outfile.close()    
def retrieve_seqs():
    outfile=open(path + "Summary_table_edt.txt","r")    
    peaks={}
    for line in outfile:    
	if 'peak' in line:
	    continue
	
	line = line.rstrip().split("\t")
	sc=[None]* (len(line)-1)
	#print line
	for i in range(0,len(line)-1):
	    #print line[i+1],
	    sc[i]=float(line[i+1])
		
	peaks[int(line[0])]=sc
    outfile.close()  
    
    outfile=open(path + "List_of_sites_all_peaks_2_calculate_edt.txt","r")    
    sites={}
    for line in outfile:
	line = line.rstrip().split(" ")
	sc=[None]* (len(line)-1)
	sc[0]=line[1]
	sc[1]=float(line[2])

	sites[int(line[0])]=sc
    outfile.close()     
    
    print 'high height peaks'
    th=50
    for peak in peaks:
	if peaks[peak][1]>th or peaks[peak][2]>th or peaks[peak][3]>th: 
	    if peaks[peak][4]>th or peaks[peak][5]>th or peaks[peak][6]>th:
		if peaks[peak][7]>th or peaks[peak][8]>th or peaks[peak][9]>th:
		    print str(peak) + " " + str(peaks[peak][0])
    print 'strong motif weak peak'
    th=10
    for peak in peaks:
	if peaks[peak][0]>10:
	    if peaks[peak][1]<th and peaks[peak][2]<th and peaks[peak][3]<th: 
		if peaks[peak][4]<th and peaks[peak][5]<th and peaks[peak][6]<th:
		    if peaks[peak][7]<th and peaks[peak][8]<th and peaks[peak][9]<th:
			print str(peak) + " " + str(peaks[peak][0])
			
    print 'weak motif strong peak'
    th=50
    for peak in peaks:
	if peaks[peak][0]<5:
	    if peaks[peak][1]>th or peaks[peak][2]>th or peaks[peak][3]>th: 
		if peaks[peak][4]>th or peaks[peak][5]>th or peaks[peak][6]>th:
		    if peaks[peak][7]>th or peaks[peak][8]>th or peaks[peak][9]>th:
			print str(peak) + " " + str(peaks[peak][0])			
    
    print 'strong sites'
    th=13
    for site in sites:
	if sites[site][1]>th:
	    print str(site )+ str(sites[site])
	    
    print 'very weak sites'
    th=-22	    
    for site in sites:
	if sites[site][1]<th:
	    print str(site )+ str(sites[site])    
	    
    
    print 'strong sites RC'
    outfile=open(path + "List_of_sites_all_peaks_2_calculate_rc_edt.txt","r")    
    sites={}
    for line in outfile:
	line = line.rstrip().split(" ")
	sc=[None]* (len(line)-1)
	sc[0]=line[1]
	sc[1]=float(line[2])

	sites[int(line[0])]=sc
    outfile.close()    
    th=12
    for site in sites:
	if sites[site][1]>th:
	    print str(site )+ str(sites[site])   
	    
	    
    print 'strong sites pSC101'
    outfile=open(path + "List_of_sites_pSC101_calculate_edt.txt","r")    
    sites={}
    for line in outfile:
	line = line.rstrip().split(" ")
	sc=[None]* (len(line)-1)
	sc[0]=line[1]
	sc[1]=float(line[2])

	sites[int(line[0])]=sc
    outfile.close()    
    th=9
    for site in sites:
	if sites[site][1]>th:
	    print str(site )+ str(sites[site])  
	    
    print 'strong sites pSC101 rc'
    outfile=open(path + "List_of_sites_pSC101_calculate_rc_edt.txt","r")    
    sites={}
    for line in outfile:
	line = line.rstrip().split(" ")
	sc=[None]* (len(line)-1)
	sc[0]=line[1]
	sc[1]=float(line[2])

	sites[int(line[0])]=sc
    outfile.close()    
    th=9
    for site in sites:
	if sites[site][1]>th:
	    print str(site )+ str(sites[site])     
    
    print 'strong sites pBR322'
    outfile=open(path + "List_of_sites_pBR322_calculate_edt.txt","r")    
    sites={}
    for line in outfile:
	line = line.rstrip().split(" ")
	sc=[None]* (len(line)-1)
	sc[0]=line[1]
	sc[1]=float(line[2])

	sites[int(line[0])]=sc
    outfile.close()    
    th=10
    for site in sites:
	if sites[site][1]>th:
	    print str(site )+ str(sites[site])      
	    
	    
    print 'strong sites pBR322 rc'
    outfile=open(path + "List_of_sites_pBR322_calculate_rc_edt.txt","r")    
    sites={}
    for line in outfile:
	line = line.rstrip().split(" ")
	sc=[None]* (len(line)-1)
	sc[0]=line[1]
	sc[1]=float(line[2])

	sites[int(line[0])]=sc
    outfile.close()    
    th=10
    for site in sites:
	if sites[site][1]>th:
	    print str(site )+ str(sites[site])    

    print 'strong sites pUC19'
    outfile=open(path + "List_of_sites_pUC19_calculate_edt.txt","r")    
    sites={}
    for line in outfile:
	line = line.rstrip().split(" ")
	sc=[None]* (len(line)-1)
	sc[0]=line[1]
	sc[1]=float(line[2])

	sites[int(line[0])]=sc
    outfile.close()    
    th=5
    for site in sites:
	if sites[site][1]>th:
	    print str(site )+ str(sites[site])      
		
		
    print 'strong sites pUC19 rc'
    outfile=open(path + "List_of_sites_pUC19_calculate_rc_edt.txt","r")    
    sites={}
    for line in outfile:
	line = line.rstrip().split(" ")
	sc=[None]* (len(line)-1)
	sc[0]=line[1]
	sc[1]=float(line[2])

	sites[int(line[0])]=sc
    outfile.close()    
    th=5
    for site in sites:
	if sites[site][1]>th:
	    print str(site )+ str(sites[site])    
    
def create_score_file_for_peaks_different_motif():
    #peaks_coords=open(path+ "List_of_peaks_merged_1.txt", 'r')
    #peaks_coords_motif=open(path+ "List_of_peaks_merged_2.txt", 'r')
    peaks_coords=open(path+ "List_of_peaks.txt", 'r')
    peaks_coords_motif=open(path+ "List_of_peaks_merged_3.txt", 'r')    
    outfile= open(path + "List_of_peaks_all_motif_3.txt".decode("utf-8"), 'w+')
    all_microcin_seqs = []
    all_cfx_seqs = []
    all_seqs = []
    all_seqs_new = []
    background = {'A': 0.245774783354, 'C': 0.2537191331, 'G': 0.254184334046, 'T': 0.246130246797}
    
    all_unique_peaks=[]
    all_unique_peaks_with_height={}
    
    for line in peaks_coords_motif:
        peaks = []
        line = line.rstrip().split('\t')
        sample = str(line[0])
        print sample
        for peak in line[1:len(line)]:
            peak = peak.split('-')
            peaks.append(int(peak[0]))
            
        seqs = []
        for i in range(len(peaks)):
            seq = genomefa[(int(peaks[i]) - win_width_l + 1):(int(peaks[i]) + win_width_r + 1)]
            seqs.append(seq)
        print 'len(seqs)=' + str(len(seqs))

        if len(seqs) == 0:
            continue
            #return

	
	all_seqs=seqs
	
	hard_copy = ['A'] * (int(background['A']* len (all_seqs)) +3) + ['C'] * (int(background['C']* len (all_seqs))) +  ['G'] * (int(background['G']* len (all_seqs))) + ['T'] * (int(background['T']* len (all_seqs)))
	
	if len(peaks)==2631 or len(peaks)==479:
	    hard_copy = ['A'] * (int(background['A']* len (all_seqs)) +3) + ['C'] * (int(background['C']* len (all_seqs))) +  ['G'] * (int(background['G']* len (all_seqs))) + ['T'] * (int(background['T']* len (all_seqs)))
	if len(peaks)==5396:
	    hard_copy = ['A'] * (int(background['A']* len (all_seqs)) +2) + ['C'] * (int(background['C']* len (all_seqs))) +  ['G'] * (int(background['G']* len (all_seqs))) + ['T'] * (int(background['T']* len (all_seqs)))
	if len(peaks)==12565:
	    hard_copy = ['A'] * (int(background['A']* len (all_seqs)) +5) + ['C'] * (int(background['C']* len (all_seqs))) +  ['G'] * (int(background['G']* len (all_seqs))) + ['T'] * (int(background['T']* len (all_seqs)))

	sm = (int(background['A']* len (all_seqs)) +5) + (int(background['C']* len (all_seqs))) +  (int(background['G']* len (all_seqs))) +  (int(background['T']* len (all_seqs)))
	
	print str(len (peaks)) + " " + str(sm)
	for k in range(65, 69):
	    soft_copy=[]
	    soft_copy.extend(hard_copy)
	    shuffle(soft_copy)
	    l = 0
	    for j in range(0, len(all_seqs)):
		ls=list(all_seqs[j])
		ls[k]= soft_copy[l]
		l = l + 1
		all_seqs[j]="".join(ls)
    
	
	instances = []
	for seq in seqs:
	    instances.append(Seq(seq))
	    
    
	m = motifs.create(instances)
	pwm = m.counts.normalize()
	pssm = pwm.log_odds(background)
	
   
    for line in peaks_coords:
        peaks = []
        line = line.rstrip().split('\t')
        sample = str(line[0])
        print sample
        for peak in line[1:len(line)]:
            peak = peak.split('-')
            peaks.append(int(peak[0]))

	peaks_output=[[],[],[],[]]
    
	#print "Peak Av_height  Position : score "
    
	#for peak in all_unique_peaks:
	for peak in peaks:
	    test_seq = Seq(str(genomefa[peak - win_width_l + 1 :peak + win_width_r + 1 ]), IUPAC.unambiguous_dna)  #
	    search_result=pssm.search(test_seq,threshold=-100)
	    max_score=-100
	    max_pos=0
	    b=0
	    for position, score in search_result:
		if position==0:
		#if score>max_score:
		    max_score=score
		    max_pos=position
		    b=1
	    if b==1:	
		peaks_output[0].append(peak)
		peaks_output[1].append(0)
		peaks_output[2].append(max_pos)
		peaks_output[3].append(float(max_score))	    
		#print("%d  %5.3f  %d  %5.3f" % (peak, np.mean(all_unique_peaks_with_height[peak]), max_pos, max_score))    
	outfile.write(str(sample))
	for i in range(0, len(peaks_output[1])):
	    outfile.write("\t" + str (peaks_output[0][i]) + ":" +str (peaks_output[3][i]))
	outfile.write("\n")
	#print 'cor coefficient ' + str(pearsonr(peaks_output[1],peaks_output[3])) 
    outfile.close()
    peaks_coords_motif.close()
    peaks_coords.close()
    
    
def create_score_file_for_peaks(prefix):
    peaks_coords=open(path+ "List_of_peaks_merged_" +prefix + ".txt", 'r')
    outfile= open(path + "List_of_peaks_score_"+prefix + ".txt".decode("utf-8"), 'w+')
    all_microcin_seqs = []
    all_cfx_seqs = []
    all_seqs = []
    all_seqs_new = []
    background = {'A': 0.245774783354, 'C': 0.2537191331, 'G': 0.254184334046, 'T': 0.246130246797}

    all_unique_peaks=[]
    all_unique_peaks_with_height={}
    for line in peaks_coords:
        peaks = []
        line = line.rstrip().split('\t')
        sample = str(line[0])
        print sample
        for peak in line[1:len(line)]:
            peak = peak.split('-')
            peaks.append(int(peak[0]))
            
        seqs = []
        for i in range(len(peaks)):
            seq = genomefa[(int(peaks[i]) - win_width_l + 1):(int(peaks[i]) + win_width_r + 1)]
            seqs.append(seq)
        print 'len(seqs)=' + str(len(seqs))

        if len(seqs) == 0:
            continue
            #return

	all_seqs=seqs
	
	hard_copy = ['A'] * (int(background['A']* len (all_seqs)) +3) + ['C'] * (int(background['C']* len (all_seqs))) +  ['G'] * (int(background['G']* len (all_seqs))) + ['T'] * (int(background['T']* len (all_seqs)))
	
	if len(peaks)==2631 or len(peaks)==479:
	    hard_copy = ['A'] * (int(background['A']* len (all_seqs)) +3) + ['C'] * (int(background['C']* len (all_seqs))) +  ['G'] * (int(background['G']* len (all_seqs))) + ['T'] * (int(background['T']* len (all_seqs)))
	if len(peaks)==5396:
	    hard_copy = ['A'] * (int(background['A']* len (all_seqs)) +2) + ['C'] * (int(background['C']* len (all_seqs))) +  ['G'] * (int(background['G']* len (all_seqs))) + ['T'] * (int(background['T']* len (all_seqs)))
	if len(peaks)==12565:
	    hard_copy = ['A'] * (int(background['A']* len (all_seqs)) +5) + ['C'] * (int(background['C']* len (all_seqs))) +  ['G'] * (int(background['G']* len (all_seqs))) + ['T'] * (int(background['T']* len (all_seqs)))

	sm = (int(background['A']* len (all_seqs)) +3) + (int(background['C']* len (all_seqs))) +  (int(background['G']* len (all_seqs))) +  (int(background['T']* len (all_seqs)))
	
	print str(len (peaks)) + " " + str(sm)
	for k in range(65, 69):
	    soft_copy=[]
	    soft_copy.extend(hard_copy)
	    shuffle(soft_copy)
	    l = 0
	    for j in range(0, len(all_seqs)):
		ls=list(all_seqs[j])
		ls[k]= soft_copy[l]
		l = l + 1
		all_seqs[j]="".join(ls)
    
	
	instances = []
	for seq in seqs:
	    instances.append(Seq(seq))
	    
    
	m = motifs.create(instances)
	pwm = m.counts.normalize()
	pssm = pwm.log_odds(background)
	peaks_output=[[],[],[],[]]
    
	#print "Peak Av_height  Position : score "
	
	#for peak in all_unique_peaks:
	for peak in peaks:
	    test_seq = Seq(str(genomefa[peak - win_width_l + 1 :peak + win_width_r + 1 ]), IUPAC.unambiguous_dna)  #
	    search_result=pssm.search(test_seq,threshold=-100)
	    max_score=-100
	    max_pos=0
	    b=0
	    for position, score in search_result:
		if position==0:
		#if score>max_score:
		    max_score=score
		    max_pos=position
		    b=1
	    if b==1:	
		peaks_output[0].append(peak)
		peaks_output[1].append(0)
		peaks_output[2].append(max_pos)
		peaks_output[3].append(float(max_score))	    
	outfile.write(str(sample))
	for i in range(0, len(peaks_output[1])):
	    outfile.write("\t" + str (peaks_output[0][i]) + ":" +str (peaks_output[3][i]))
	outfile.write("\n")
    outfile.close()

def ROC_curves():
        peaks_coords=open(path+ "List_of_peaks_score1.txt", 'r')
	sites_coords=open(path + "List_of_sites_merged_pos_th_0.txt", 'r')
	outfile= open(path + "tpr_fpr1.txt".decode("utf-8"), 'w+')
	prec=0.1
	max_score=16.0
	sites={}
	peaks={}
    
	TPR=[]
	FPR=[]
	scr=[]
	not_peaks=0
    
	for line in sites_coords:	
	    line=line.rstrip()
	    line = line.split(' ')
	    score=float((line[2]))
	    if score>=0 and int (line[0]) > 0 and int(line[0])<5000000:
		sites[int(line[0])]=score
    
    
	for line in peaks_coords:
	    line=line.rstrip()
	    line = line.split('\t')
	    sample = str(line[0])
	    for peak in line[1:len(line)]:
		pk = peak.split(':')
		#print pk
		peak=int(pk[0])
		if str(pk[1])[0]=='$': 
		    height=list(pk[1])
		    height[0]='-'
		    height=float("".join(height))
		    #print height
		else:
		    height=float((pk[1]))
		if peak in peaks:
		    if peaks[peak]<height:
			peaks[peak]=height
		else:
		    peaks[peak]=height
    
	for	site in sites:
	    if site not in peaks:
		not_peaks+=1
    
	#print np.arange(0.0,max_score,prec)
	for score in np.arange(0.0,max_score,prec):
	    print score
	    tp=0
	    fn=0
	    tn=0
	    fp=0
	    for peak in peaks:
		if peaks[peak]>=score:
		    if score>=15:
			print peak
		    tp+=1
	    for site in sites:
		if site not in peaks and sites[site]>=score:
		    fp+=1
	    scr.append(score)
	    TPR.append(float(tp)/ len(peaks))
	    FPR.append(float(fp)/not_peaks)
    
    
    
    
	tpr_string="TPR"
	fpr_string="FPR"
	scr_string="score"
	for i in range (0, len(TPR)):
	    tpr_string=tpr_string + "\t" + str(TPR[i])
	    fpr_string=fpr_string + "\t" + str(FPR[i])
	    scr_string=scr_string + "\t" + str(scr[i])
	tpr_string=tpr_string + "\n" 
	fpr_string=fpr_string + "\n" 
	scr_string=scr_string + "\n" 
    
	outfile.write(scr_string)
	outfile.write(tpr_string)
	outfile.write(fpr_string)
    
	outfile.close()
    
	plt.figure(figsize=(16, 8), dpi=100)
	plt.plot(FPR, TPR)
	plt.title("ROC")
	#plt.show()
	plt.savefig(path+"roc.png", dpi=100, figsize=(16, 8))
	plt.close() 
    
	plt.figure(figsize=(16, 8), dpi=100)
	plt.plot(scr, TPR,label='TPR', color='blue')
    
	plt.plot(scr, FPR, label = 'FPR', color='red')
    
	plt.plot([0, 16], [0.05,0.05], '--')
	plt.title("FPR TPR")
	plt.legend()
	#plt.show()
	plt.savefig(path+"fpr_tpr1.png", dpi=100, figsize=(16, 8))
	plt.close()     
    
    

def new_files_averaged_lines():
    files=[ "Oxo_IP_Mu_120_mkM_1_ends.wig", "Oxo_IP_120_mkM_1_ends.wig", "Oxo_IP_120_mkM_2_ends.wig", 
	   "Microcin_IP_50_mkM_4_ends.wig", "Un_IN_Mu_1_new_ends.wig",   "Un_IP_Mu_1_new_ends.wig", "Un_IN_1_ends.wig", "Un_IN_2_ends.wig"]
    mask_array=[]
    for k in range(0,4647999,1):
	if (k >= d1_b and k <= d1_e):
	    mask_array.append(True)
	elif (k >= d2_b and k <= d2_e):
	    mask_array.append(True)
	elif (k >= d3_b and k <= d3_e):
	    mask_array.append(True)
	else:
	    mask_array.append(False)    
    for fl in files:
	print fl
	pars_com = Parsers()
	wigin = open(path+fl, 'r')
	ends=pars_com.Wig_pars_return_normalized(wigin) 
	control=[]
	#create_gc_plot_fitting(fl, peaks_coord, ends, control)    
	
	
	mc = np.ma.masked_array(return_averaged_line(ends), mask=mask_array)
	
    
	plt.figure(figsize=(16, 8), dpi=100)
	plt.plot(xcoord,mc, '.', label='initial averaged 200kb', color='blue')
	
    
	plt.title(fl[:-4] + "_200kb")
	plt.legend()
	plt.savefig(path + fl[:-4] + '_200kb.png', dpi=100, figsize=(16, 8))
	plt.close()

def pitctures_for_presentation():
    
    
    ap=venn_diagram()
    all_peaks1=ap[0]
    all_peaks2=ap[1]
    all_peaks3=ap[2]
    
    

    
    
    
    create_score_file_for_peaks("1")
    create_score_file_for_peaks("2")
    #ROC_curves()
    
    
    print "cor 1"
    list_of_peaks_merged= open(path + "List_of_peaks_merged_1.txt".decode("utf-8"), 'r')
    pk_heights={}
    for line in list_of_peaks_merged:    
	    line=line.rstrip()
	    line = line.split('\t')
	    sample = str(line[0])
	    for peak in line[1:len(line)]:
		    pk = peak.split('-')    
		    pk_heights[int(pk[0])]=float((pk[1])) 
    list_of_peaks_scores= open(path + "List_of_peaks_score_1.txt".decode("utf-8"), 'r')
    pk_scores={}
    for line in list_of_peaks_scores:    
	    line=line.rstrip()
	    line = line.split('\t')
	    sample = str(line[0])
	    for peak in line[1:len(line)]:
		    pk = peak.split(':')    
		    pk_scores[int(pk[0])]=float((pk[1])) 
    heights=[]
    scores=[]    
    for pk in pk_scores.iterkeys():
	    heights.append(pk_heights[pk])
	    scores.append(pk_scores[pk])    
    print "len cor 1 " + str(len(scores)    )
    print 'cor coefficient ' + str(pearsonr(heights,scores))    
    list_of_peaks_merged.close()
    
    
    print "cor 2"
    list_of_peaks_merged= open(path + "List_of_peaks_merged_2.txt".decode("utf-8"), 'r')
    pk_heights={}
    for line in list_of_peaks_merged:    
	    line=line.rstrip()
	    line = line.split('\t')
	    sample = str(line[0])
	    for peak in line[1:len(line)]:
		    pk = peak.split('-')    
		    pk_heights[int(pk[0])]=float((pk[1]))    
    list_of_peaks_scores= open(path + "List_of_peaks_score_2.txt".decode("utf-8"), 'r')
    pk_scores={}
    for line in list_of_peaks_scores:    
	    line=line.rstrip()
	    line = line.split('\t')
	    sample = str(line[0])
	    for peak in line[1:len(line)]:
		    pk = peak.split(':')    
		    pk_scores[int(pk[0])]=float((pk[1])) 
    
    heights=[]
    scores=[]    
    for pk in pk_scores.iterkeys():
	    heights.append(pk_heights[pk])
	    scores.append(pk_scores[pk])
    print "len cor 2 " + str(len(scores))
    print 'cor coefficient ' + str(pearsonr(heights,scores))
    
    list_of_peaks_merged.close()
    
    '''
    print len(all_peaks1)
    print len(all_peaks2)
    #all_peaks2={}
    #for peak in all_peaks:
    #    all_peaks2[peak]=1
    genes_file="/data/Gyrase/Pipe/Escherichia_coli_K_12_w3110_Mu_genes_annotation.gff"
    operons_file="/data/Gyrase/Pipe/K12_w3110_operons_only.gff"
    genes=parse_genes_file(genes_file)
    operons=parse_operons_file(operons_file) 
    wndws=[0.1, 0.25, 0.5, 1,2,3,4,5,6]
    sample="all_peaks2"
    op_s_values=[None] * len (wndws)
    gene_s_values=[None] * len (wndws)
    for wndw in wndws: 
	print 'window ' + str(wndw)    
	prefix="_" + str(wndw)
	raw_counts=calculate_counts_for_TU_plots_wide_intergenic_regions(sample + "_peaks_genes", all_peaks2,genes,wndw)
	counts=raw_counts
	TU_plot(sample + "_peaks_genes_wide" +prefix ,counts)
	gene_s_values[wndws.index(wndw)]=float(counts[5])/counts[0]
    
	raw_counts=calculate_counts_for_TU_plots_wide_intergenic_regions(sample+ "_peaks_operons", all_peaks2,operons,wndw)
	counts=raw_counts
	TU_plot(sample + "_peaks_operons_wide" +prefix,counts)	 
	op_s_values[wndws.index(wndw)]=float(counts[5])/counts[0]
    
    plt.figure(figsize=(16, 8), dpi=100)
    ax = plt.subplot()
    
    plt.plot(wndws,gene_s_values, label='genes', color='blue')
    plt.plot(wndws,op_s_values, label='operons', color='red')
    plt.title("distribution all peaks2" + "")
    ax.set_xticks(wndws)   
    plt.legend()
    plt.savefig(path+'distribution_al_peaks2.png', dpi=100, figsize=(16, 8))
    plt.close() 
    '''    

def combine_bed_graph():
    f1=open(u"/data/Gyrase/Pipe/Fragments_ends/Transcription/GSM2516397_RNAseq_3end_Cyto_1.bedGraph", 'r')
    f2=open(u"/data/Gyrase/Pipe/Fragments_ends/Transcription/GSM2516400_RNAseq_3end_Cyto_2.bedGraph", 'r')
    f3=open(u"/data/Gyrase/Pipe/Fragments_ends/Transcription/GSM2516401_RNAseq_3end_Nucleoid_1.bedGraph", 'r')
    f4=open(u"/data/Gyrase/Pipe/Fragments_ends/Transcription/GSM2516402_RNAseq_3end_Nucleoid_2.bedGraph", 'r')
    
    o1=open(u"/data/Gyrase/Pipe/Fragments_ends/Transcription/Cyto_1.wig", 'w+')
    o2=open(u"/data/Gyrase/Pipe/Fragments_ends/Transcription/Cyto_2.wig", 'w+')
    o11=open(u"/data/Gyrase/Pipe/Fragments_ends/Transcription/Cyto_1.bedGraph", 'w+')
    o22=open(u"/data/Gyrase/Pipe/Fragments_ends/Transcription/Cyto_2.bedGraph", 'w+')    
    o3=open(u"/data/Gyrase/Pipe/Fragments_ends/Transcription/Nucleoid_1.wig", 'w+')
    o4=open(u"/data/Gyrase/Pipe/Fragments_ends/Transcription/Nucleoid_2.wig", 'w+')
    
    o5=open(u"/data/Gyrase/Pipe/Fragments_ends/Transcription/merged.wig", 'w+')
    
    o11.write("track type=bedGraph name=\"TopHat - read coverage\"\n")
    o22.write("track type=bedGraph name=\"TopHat - read coverage\"\n")
    
    o1.write("track type=wiggle_0 name=\"Cyto_1\" description=\"decription\" autoScale=off viewLimits=0.0:25.0\n")
    o1.write("fixedStep chrom=U00096 start=1 step=1\n")
    o2.write("track type=wiggle_0 name=\"Cyto_2\" description=\"decription\" autoScale=off viewLimits=0.0:25.0\n")
    o2.write("fixedStep chrom=U00096 start=1 step=1\n")
    o3.write("track type=wiggle_0 name=\"Nucleoid_1\" description=\"decription\" autoScale=off viewLimits=0.0:25.0\n")
    o3.write("fixedStep chrom=U00096 start=1 step=1\n")
    o4.write("track type=wiggle_0 name=\"Nucleoid_2\" description=\"decription\" autoScale=off viewLimits=0.0:25.0\n")
    o4.write("fixedStep chrom=U00096 start=1 step=1\n")    
    
    o5.write("track type=wiggle_0 name=\"merged\" description=\"decription\" autoScale=off viewLimits=0.0:25.0\n")
    o5.write("fixedStep chrom=U00096 start=1 step=1\n")     
    
    total=[[],[],[],[]]
    prev=0
    for line in f1:
	if line.startswith("track"):
	    continue
	line=line.rstrip().split("\t")
	
	for i in range(int(line[1]), int(line[2]) ):
	    o1.write(line[3] + "\n")
	    total[0].append(int(line[3]))
	for i in range(int(line[1]), int(line[2])):
	    o11.write("U00096" + "\t" + str(i) + "\t" +str (i+1) +"\t" + line[3] + "\n")
    print line[1] + "\t" + line[2] + "\t" + line[3]
    prev=0
    for line in f2:
	if line.startswith("track"):
	    continue
	line=line.rstrip().split("\t")
	if prev != int(line[1]):
	    print "WARNING2"
	    print line[1]
	    print prev
	    break	
	for i in range(int(line[1]), int(line[2]) ):
	    o2.write(line[3] + "\n")
	    total[1].append(int(line[3]))
	for i in range(int(line[1]), int(line[2]) ):
	    o22.write("U00096" + "\t" + str(i) + "\t" +str (i+1) +"\t" + line[3] + "\n")	
	prev=int(line[2])
    print line[1] + "\t" + line[2] + "\t" + line[3]
    prev=0
    for line in f3:
	if line.startswith("track"):
	    continue
	line=line.rstrip().split("\t")
	for i in range(int(line[1]), int(line[2]) ):
	    o3.write(line[3] + "\n")
	    total[2].append(int(line[3]))
    print line[1] + "\t" + line[2] + "\t" + line[3]	    
    prev=0
    for line in f4:
	if line.startswith("track"):
	    continue
	line=line.rstrip().split("\t")
	for i in range(int(line[1]), int(line[2]) ):
	    o4.write(line[3] + "\n")	
	    total[3].append(int(line[3]))
    print line[1] + "\t" + line[2] + "\t" + line[3]
    print len(total[0])
    print len(total[1])
    print len(total[2])
    print len(total[3])
    max_len=max(len(total[0]),len(total[1]),len(total[2]),len(total[3]),len(genomefa_2))    
    for i in range(0, max_len):
	sm=0
	if i>=len(total[0]):
	    sm=sm + 0
	else:
	    sm=sm+total[0][i]
	if i>=len(total[1]):
	    sm=sm + 0
	else:
	    sm=sm+total[1][i]	
	if i>=len(total[2]):
	    sm=sm + 0
	else:
	    sm=sm+total[2][i]	
	if i>=len(total[3]):
	    sm=sm + 0
	else:
	    sm=sm+total[3][i]	
	o5.write(str(sm) + "\n")
    
    
    f1.close()
    f2.close()
    f3.close()
    f4.close()
    o1.close()
    o2.close()
    o3.close()
    o4.close()
    o5.close()
def check_transcription_level():
    #print genome2.features[:10]
    wigfile=open(path + "Transcription/merged.wig", 'r')
    gfffile=open("/data/Gyrase/Pipe/Escherichia_coli_K_12_w3110_Mu_genes_annotation.gff", 'r')
    o1=open(path + "/Transcription/transcription.txt", 'w+')
    o2=open(path + "/Transcription/gene_expression.txt", 'w+')
    wig=[]
    gff=[]
    gb={}
    gb_syn={}
    gb_w3110_syn={}
    notfound=0
    found=0
    found_flag=0
    for line in wigfile:
	if line.startswith("track") or line.startswith("fixedStep") :
	    continue
	#line=line.rstrip().split("\t")   
	wig.append(int(line))
    
    for feat in genome2.features:
	if feat.type=='gene':
	    #print feat
	    #print feat.location.start
	    #print np.mean(wig[feat.location.start: feat.location.end])
	    o1.write (str(feat.qualifiers['gene'])[2:-2] + "\t" + str(feat.location.start) + "\t" + str(feat.location.end)  + "\t" + str(feat.location.strand)+ "\t" + str(np.mean(wig[feat.location.start: feat.location.end])) +"\n")
	    gb[str(feat.qualifiers['gene'])[2:-2]]=str(feat.location.start) + "\t" + str(feat.location.end)  + "\t" + str(feat.location.strand)+ "\t" + str(np.mean(wig[feat.location.start: feat.location.end]))
	    gb_syn[str(feat.qualifiers['locus_tag'])[2:-2]]=str(feat.qualifiers['gene'])[2:-2] + "\t" +str(feat.qualifiers['gene_synonym'])+ "\t"+ str(feat.location.start) + "\t" + str(feat.location.end)  + "\t" + str(feat.location.strand)+ "\t" + str(np.mean(wig[feat.location.start: feat.location.end]))
	    
    for feat in genome_w3110.features:
	if feat.type=='CDS':
	    #print feat
	    #print feat.location.start
	    #print np.mean(wig[feat.location.start: feat.location.end])
	    #print feat
	    #print str(feat.qualifiers['note'])
	    note=str(feat.qualifiers['note']).split(";")
	    for i in range(0,len(note)):
		if ":" in note[i]:
		    note=note[i]
		    break
	    #print note

	    note=note.split(":")[2]
	    if note.endswith("']"):
		note=note[:-2]
	    #print note
	    gb_w3110_syn[note]=str(feat.qualifiers['gene'])[2:-2] + "\t" +str(feat.location.start) + "\t" + str(feat.location.end)  + "\t" + str(feat.location.strand)

    strand_mismatch=0
    for tag in gb_w3110_syn:
	genename=gb_w3110_syn[tag].split("\t")[0]
	found_flag=0
	start=gb_w3110_syn[tag].split("\t")[1]
	end=gb_w3110_syn[tag].split("\t")[2]
	strand=gb_w3110_syn[tag].split("\t")[3]	
	if tag in gb_syn:
	    found+=1
	    found_flag=1
	    strand2=gb_syn[tag].split("\t")[4]
	    expression=gb_syn[tag].split("\t")[5] 
	    start2=gb_syn[tag].split("\t")[2]
	    end2=gb_syn[tag].split("\t")[3]
	    o2.write(genename + "\t"  + tag + "\t"+start+ "\t"+ end + "\t"+strand + "\t"+expression+ "\n")
	    if strand!=strand2:
		strand_mismatch+=1
		print genename + "\t"  +start+ "\t"+ end + "\t"+strand  + "\t"+start2+ "\t"+ end2 + "\t"+strand2
	else: 	    
	    for tag2 in gb_syn:
		if genename in gb_syn[tag2]:
		    found+=1
		    found_flag=1
		    strand2=gb_syn[tag2].split("\t")[4]
		    expression=gb_syn[tag2].split("\t")[5]
		    start2=gb_syn[tag2].split("\t")[2]
		    end2=gb_syn[tag2].split("\t")[3]		    
		    o2.write(genename + "\t"  + tag + "\t"+start+ "\t"+ end + "\t"+strand + "\t"+expression +"\n")
		    if strand!=strand2:
			strand_mismatch+=1
			print genename + "\t"  +start+ "\t"+ end + "\t"+strand  + "\t"+start2+ "\t"+ end2 + "\t"+strand2
		    break
	    if found_flag==0:
		notfound+=1
		#print genename
	    
    for line in gfffile:
	found_flag=0
	if not line.startswith("NC_007779.1_w3110_Mu	ena	gene"):
	    continue
	line=line.rstrip().split("\t")
	name=line[8].split(';')[1].split("=")[1]
	gff.append(name)

    print notfound
    print found
    print strand_mismatch
    gfffile.close()
    wigfile.close()
    o1.close()
    o2.close()
    
    

def get_genes_for_transcription():
    st=open(path + "Summary_table_edt.txt".decode("utf-8"), 'r')
    ge = open(path + "gene_expression.txt".decode("utf-8"), 'r')
    
    pks = {}
    pks_expression={}
    
    cfx={}
    micro={}
    oxo={}
    micro_cfx={}
    micro_cfx_oxo={}
    
    for line in st:
	if line.startswith("peak"):
	    continue
	line1 = line.rstrip().split('\t')
	pks[int(line1[0])]=line.rstrip()
	tmp_m=0
	tmp_c=0
	tmp_o=0
	
	height=[]
	
	if line1[2]!=0:
	    tmp_c+=1
	if line1[3]!=0:
	    tmp_c+=1	
	if line1[4]!=0:
	    tmp_c+=1	    
	
	if line1[5]!=0:
	    tmp_m+=1	
	if line1[6]!=0:
	    tmp_m+=1
	if line1[7]!=0:
	    tmp_m+=1	    
	    
	if line1[8]!=0:
	    tmp_o+=1
	if line1[9]!=0:
	    tmp_o+=1
	if line1[10]!=0:
	    tmp_o+=1
	
	if tmp_c>=2:
	    cfx[int(line1[0])]=[]
	    height.append(float(line1[2]))
	    height.append(float(line1[3]))
	    height.append(float(line1[4]))
	    cfx[int(line1[0])][0]=[np.mean(height)]
	    cfx[int(line1[0])][0]=[float(line1[1])]
	 
	if tmp_m>=2:
	    micro[int(line1[0])]=[]
	    height.append(float(line1[5]))
	    height.append(float(line1[6]))
	    height.append(float(line1[7]))
	    micro[int(line1[0])][0]=[np.mean(height)]
	    micro[int(line1[0])][0]=[float(line1[1])]
	    
	if tmp_o>=2:
	    oxo[int(line1[0])]=[]
	    height.append(float(line1[8]))
	    height.append(float(line1[9]))
	    height.append(float(line1[10]))
	    oxo[int(line1[0])][0]=[np.mean(height)]
	    oxo[int(line1[0])][0]=[float(line1[1])]
	    
	for peak in micro:
	    if peak not in micro_cfx :
		micro_cfx[peak]=micro[peak]
	for peak in cfx:
	    if peak not in micro_cfx :
		micro_cfx[peak]=cfx[peak]
	    else:
		height=micro_cfx[peak][0]
		height=(1.0 *height + cfx[peak][0])/2
		micro_cfx[peak][0]=height
	
	micro_cfx_oxo=micro_cfx
	
	for peak in oxo:
	    if peak not in micro_cfx_oxo :
		micro_cfx_oxo[peak]=oxo[peak]
	    else:
		height=micro_cfx_oxo[peak][0] *2
		height=(1.0 *height + oxo[peak][0])/2
		micro_cfx_oxo[peak][0]=height	    
	
	
    print str(len(micro)) + ' ' + str(len(cfx)) + ' ' + str(len(oxo)) 
    
    genes={}
    expression=[]
    for line in ge:
	line1 = line.rstrip().split('\t')
	genes[line1[0]]=line.rstrip()
	expression.append(float(line1[5]))

    print sorted(expression)
    threshold=sorted(expression)[len (expression) - 100]
    print threshold

    wnd = 1000
    genes_peaks={}
    genes_hg = {}
    genes_ms = {}
    
    for gene in genes:

	start = int(genes[gene].split('\t')[2])
	end = int(genes[gene].split('\t')[3])
	strand = int(genes[gene].split('\t')[4])
	expression = float(genes[gene].split('\t')[5])
	if expression>= threshold:
	    genes_peaks[gene] = [[], [], []]
	    genes_hg [gene] = [[], [], []]
	    genes_ms[gene] = [[], [], []]
	    for peak in pks:
		if genes[gene].split("\t")[5] >= threshold:
		    if peak >= start and peak <= end:
			genes_peaks[gene][1].append(peak)
			genes_ms[gene][1].append(float(pks[peak].split('\t')[1]))

			hg = 0
			for i in range(1, 10):
			    hg = hg + float(pks[peak].split('\t')[i + 1])
			genes_hg[gene][1].append(float(hg / 9))

		    if strand == 1:
			if peak < start and peak >= start - wnd:
			    genes_peaks[gene][0].append(peak)
			    genes_ms[gene][0].append(float(pks[peak].split('\t')[1]))
			    hg = 0
			    for i in range(1, 10):
				hg = hg + float(pks[peak].split('\t')[i + 1])
			    genes_hg[gene][0].append(float(hg / 9))
			if peak > end and peak <= end + wnd:
			    genes_peaks[gene][2].append(peak)
			    genes_ms[gene][2].append(float(pks[peak].split('\t')[1]))
			    hg = 0
			    for i in range(1, 10):
				hg = hg + float(pks[peak].split('\t')[i + 1])
			    genes_hg[gene][2].append(float(hg / 9))

		    if strand == -1:
			if peak < start and peak >= start - wnd:
			    genes_peaks[gene][2].append(peak)
			    genes_ms[gene][2].append(float(pks[peak].split('\t')[1]))
			    hg = 0
			    for i in range(1, 10):
				hg = hg + float(pks[peak].split('\t')[i + 1])
			    genes_hg[gene][2].append(float(hg / 9))
			if peak > end and peak <= end + wnd:
			    genes_peaks[gene][0].append(peak)
			    genes_ms[gene][0].append(float(pks[peak].split('\t')[1]))
			    hg = 0
			    for i in range(1, 10):
				hg = hg + float(pks[peak].split('\t')[i + 1])
			    genes_hg[gene][0].append(float(hg / 9))
    sum =[0,0,0,0,0,0,0,0,0]
    for gene in genes_peaks:
	print gene + "\t" + str(len(genes_peaks[gene][0]))+ "\t" + str(len(genes_peaks[gene][1]))+ "\t" + str(len(genes_peaks[gene][2])),
	print  "\t" + str(np.mean(genes_ms[gene][0])) + "\t" + str(np.mean(genes_ms[gene][1])) + "\t" + str(np.mean(genes_ms[gene][2])),
	print  "\t" + str(np.mean(genes_hg[gene][0])) + "\t" + str(np.mean(genes_hg[gene][1])) + "\t" + str(np.mean(genes_hg[gene][2]))

	for i in range(0,3):
	    sum[i]=sum[i]+ len(genes_peaks[gene][i])
	for i in range(3,6):
	    if np.mean(genes_ms[gene][i-3]) >=0 :
		sum[i]=sum[i]+ np.mean(genes_ms[gene][i-3])
	for i in range(6,9):
	    if np.mean(genes_hg[gene][i-6]) >=0 :
		sum[i]=sum[i]+ np.mean(genes_hg[gene][i-6])

    print sum
    for peak in pks:
	pks_expression[peak]=[[],[],[]]
	for gene in genes:
	    if genes[gene].split("\t")[5] >= threshold:
		start=int(genes[gene].split('\t')[2])
		end = int( genes[gene].split('\t')[3])
		strand = int( genes[gene].split('\t')[4])
		expression=float(genes[gene].split('\t')[5])
		if peak>=start and peak <= end:
		    pks_expression[peak][1].append(genes[gene])
		if strand==1:
		    if peak<start and peak>=start - wnd:
			pks_expression[peak][0].append(genes[gene])
		    if peak>end  and peak<=end +  wnd:
			pks_expression[peak][2].append(genes[gene])
		if strand==-1:
		    if peak<start and peak>=start - wnd:
			pks_expression[peak][2].append(genes[gene])
		    if peak>end  and peak<=end +  wnd:
			pks_expression[peak][0].append(genes[gene])
    motif_score=[]
    motif_score_in=[]
    motif_score_up = []
    motif_score_down = []
    height=[[],[],[],[],[],[],[],[],[]]
    height_mean=[]
    expr_in=[]
    expr_up=[]
    expr_down=[]
    #for peak in pks_expression:
	#print pks_expression[peak]

    for peak in pks:
	#print peak
	motif_score.append(float(pks[peak].split('\t')[1]))
	hg=0
	for i in range(1,10):
	    height[i-1].append(float(pks[peak].split('\t')[i+1]))
	    hg=hg+float(pks[peak].split('\t')[i+1])
	height_mean.append(float(hg/9))

	expr=[]
	for line in pks_expression[peak][1]:
	    #print line
	    line=line.split("\t")
	    expr.append(float(line[5]))
	if len(expr)>0:
	    expr_in.append(np.mean(expr))
	    motif_score_in.append(float(pks[peak].split('\t')[1]))

	expr=[]
	for line in pks_expression[peak][0]:
	    line=line.split("\t")
	    expr.append(float(line[5]))
	if len(expr) > 0:
	    expr_up.append(np.mean(expr))
	    motif_score_up.append(float(pks[peak].split('\t')[1]))

	expr=[]
	for line in pks_expression[peak][2]:
	    line=line.split("\t")
	    expr.append(float(line[5]))
	if len(expr) > 0:
	    expr_down.append(np.mean(expr))
	    motif_score_down.append(float(pks[peak].split('\t')[1]))

    print pearsonr(motif_score_in,expr_in)
    print pearsonr(motif_score_up, expr_up)
    print pearsonr(motif_score_down, expr_down)
    print pearsonr(height_mean, expr_down)

def september_get_genes_for_transcription():
    st=open(path + "Summary_table_edt.txt".decode("utf-8"), 'r')
    #ge = open(path + "gene_expression_mu_correction.txt".decode("utf-8"), 'r')
    ge = open(path + "operon_expression.txt".decode("utf-8"), 'r')
    
    pks = {}
    pks_expression={}
    
    cfx={}
    micro={}
    oxo={}
    micro_cfx={}
    micro_cfx_oxo={}
    
    for line in st:
	if line.startswith("peak"):
	    continue
	line1 = line.rstrip().split('\t')
	tmp_m=0
	tmp_c=0
	tmp_o=0
	
	height=[]
	
	if line1[2]!='0':
	    tmp_c+=1
	if line1[3]!='0':
	    tmp_c+=1	
	if line1[4]!='0':
	    tmp_c+=1	    
	
	if line1[5]!='0':
	    tmp_m+=1	
	if line1[6]!='0':
	    tmp_m+=1
	if line1[7]!='0':
	    tmp_m+=1	    
	    
	if line1[8]!='0':
	    tmp_o+=1
	if line1[9]!='0':
	    tmp_o+=1
	if line1[10]!='0':
	    tmp_o+=1
	
	height=[]
	if tmp_c>=2:
	    cfx[int(line1[0])]=[0,0]
	    height.append(float(line1[2]))
	    height.append(float(line1[3]))
	    height.append(float(line1[4]))
	    cfx[int(line1[0])][0]=np.mean(height)
	    cfx[int(line1[0])][1]=float(line1[1])
	
	height=[] 
	if tmp_m>=2:
	    micro[int(line1[0])]=[0,0]
	    height.append(float(line1[5]))
	    height.append(float(line1[6]))
	    height.append(float(line1[7]))
	    micro[int(line1[0])][0]=np.mean(height)
	    micro[int(line1[0])][1]=float(line1[1])
	    #print int(line1[0]) , 
	    #print height ,
	    #print float(line1[1])
	
	height=[]    
	if tmp_o>=2:
	    oxo[int(line1[0])]=[0,0]
	    height.append(float(line1[8]))
	    height.append(float(line1[9]))
	    height.append(float(line1[10]))
	    oxo[int(line1[0])][0]=np.mean(height)
	    oxo[int(line1[0])][1]=float(line1[1])
	    
    for peak in micro:
	if peak not in micro_cfx :
	    micro_cfx[peak]=micro[peak]
    for peak in cfx:
	if peak not in micro_cfx :
	    micro_cfx[peak]=cfx[peak]
	else:
	    height=micro_cfx[peak][0]
	    #print height
	    height=(1.0 *height + cfx[peak][0])/2
	    #height=max(height,cfx[peak][0])
	    micro_cfx[peak][0]=height
    
    micro_cfx_oxo=micro_cfx
    
    for peak in oxo:
	if peak not in micro_cfx_oxo :
	    micro_cfx_oxo[peak]=oxo[peak]
	else:
	    height=micro_cfx_oxo[peak][0] *2
	    height=(1.0 *height + oxo[peak][0])/2
	    #height=micro_cfx_oxo[peak][0] 
	    #height=max(height,oxo[peak][0])	    
	    micro_cfx_oxo[peak][0]=height	    
	
    print 'PART 1'
    print str(len(micro)) + ' ' + str(len(cfx)) + ' ' + str(len(oxo)) 
    
    genes={}
    expression=[]
    for line in ge:
	line1 = line.rstrip().split('\t')
	name=line1[0]
	start = int(line1[2])
	end = int(line1[3])
	strand = int(line1[4])
	expr = float(line1[5])
	if name in genes:
	    name=name + "_"+str(start)+ "_"+str(end)
	genes[name]=line.rstrip()	
	#genes[line1[0]]=line.rstrip()
	expression.append(float(line1[5]))

    #print sorted(expression)
    threshold=sorted(expression)[len (expression) - 100]
    print threshold
    #threshold=0
    threshold=sorted(expression)[100]
    print threshold
    wnd = 1000
    genes_peaks={}
    genes_hg = {}
    genes_ms = {}
    
    peaks={}
    peaks["micro"]=micro
    peaks["cfx"]=cfx
    peaks["oxo"]=oxo
    peaks["micro_cfx"]=micro_cfx
    peaks["micro_cfx_oxo"]=micro_cfx_oxo
    
    j=0
    for sample in sorted(list(peaks.keys())):
	print "======================= " + sample + " ======================="
    
	for gene in genes:
	    start = int(genes[gene].split('\t')[2])
	    end = int(genes[gene].split('\t')[3])
	    strand = int(genes[gene].split('\t')[4])
	    expression = float(genes[gene].split('\t')[5])
	    if j==-1:
		print gene ,
		print start ,
		print end ,
		print strand ,
		print expression 
	    if expression<= threshold:
		genes_peaks[gene] = [[], [], []]
		genes_hg [gene] = [[], [], []]
		genes_ms[gene] = [[], [], []]
		for peak in peaks[sample]:
		    if peak >= start and peak <= end:
			genes_peaks[gene][1].append(peak)
			genes_ms[gene][1].append(peaks[sample][peak][1])    
			genes_hg[gene][1].append(peaks[sample][peak][0])

		    if strand == 1:
			if peak < start and peak >= start - wnd:
			    genes_peaks[gene][0].append(peak)
			    genes_ms[gene][0].append(peaks[sample][peak][1])
			    genes_hg[gene][0].append(peaks[sample][peak][0])
			    
			if peak > end and peak <= end + wnd:
			    genes_peaks[gene][2].append(peak)
			    genes_ms[gene][2].append(peaks[sample][peak][1])
			    genes_hg[gene][2].append(peaks[sample][peak][0])
			    

		    if strand == -1:
			if peak < start and peak >= start - wnd:
			    genes_peaks[gene][2].append(peak)
			    genes_ms[gene][2].append(peaks[sample][peak][1])
			    genes_hg[gene][2].append(peaks[sample][peak][0])
			    
			if peak > end and peak <= end + wnd:
			    genes_peaks[gene][0].append(peak)
			    genes_ms[gene][0].append(peaks[sample][peak][1])
			    genes_hg[gene][0].append(peaks[sample][peak][0])
				
	
	
	j=0
	for gene in genes_peaks:
	    start = int(genes[gene].split('\t')[2])
	    end = int(genes[gene].split('\t')[3])
	    strand = int(genes[gene].split('\t')[4])
	    expression = float(genes[gene].split('\t')[5])	
	    
	    if  len (genes_peaks[gene][0]) + len (genes_peaks[gene][1]) +  len (genes_peaks[gene][2])>0 and j==-1:
		print gene ,
		print start ,
		print end ,
		print strand ,
		print expression 
	    
	    if len (genes_peaks[gene][0])>0 and j==-1 :		
		print 'UP' ,
		for i in range (0, len (genes_peaks[gene][0])):
		    print genes_peaks[gene][0][i] ,
		    print genes_hg[gene][0][i] ,
		    print genes_ms[gene][0][i]
		    #motif_score_up.append(genes_ms[gene][0][i])
		    #height_up.append(genes_hg[gene][0][i])	
		    
	    if len (genes_peaks[gene][1])>0 and j==-1:		
		print 'IN' ,
		for i in range (0, len (genes_peaks[gene][1])):
		    print genes_peaks[gene][1][i] ,
		    print genes_hg[gene][1][i] ,
		    print genes_ms[gene][1][i]
		    #motif_score_up.append(genes_ms[gene][0][i])
		    #height_up.append(genes_hg[gene][0][i])	
		    
	    if len (genes_peaks[gene][2])>0 and j==-1:		
		print 'DOWN' ,
		for i in range (0, len (genes_peaks[gene][2])):
		    print genes_peaks[gene][2][i] ,
		    print genes_hg[gene][2][i] ,
		    print genes_ms[gene][2][i] 
		    #motif_score_up.append(genes_ms[gene][0][i])
		    #height_up.append(genes_hg[gene][0][i])	
	    if  len (genes_peaks[gene][0]) + len (genes_peaks[gene][1]) +  len (genes_peaks[gene][2])>0 and j==-1:
		print '\n'
	    
	    
	j+=1
	motif_score_in=[]
	motif_score_up=[]
	motif_score_down=[]
	height_in=[]
	height_up=[]
	height_down=[]
	av_motif_score_in=[]
	av_motif_score_up=[]
	av_motif_score_down=[]
	av_height_in=[]
	av_height_up=[]
	av_height_down=[]
	cum_motif_score_in=[]
	cum_motif_score_up=[]
	cum_motif_score_down=[]
	cum_height_in=[]
	cum_height_up=[]
	cum_height_down=[]
	max_motif_score_in=[]
	max_motif_score_up=[]
	max_motif_score_down=[]
	max_height_in=[]
	max_height_up=[]
	max_height_down=[]	
	list_expression_in=[]	
	list_expression_up=[]	
	list_expression_down=[]
	single_expression_in=[]	
	single_expression_up=[]	
	single_expression_down=[]	
	number_of_peaks_in=[]
	number_of_peaks_up=[]
	number_of_peaks_down=[]
	
	for gene in genes_peaks:
	    start = int(genes[gene].split('\t')[2])
	    end = int(genes[gene].split('\t')[3])
	    strand = int(genes[gene].split('\t')[4])
	    expression = float(genes[gene].split('\t')[5])	
	    
	    
	    for i in range (0, len (genes_peaks[gene][0])):
		list_expression_up.append(expression)
		motif_score_up.append(genes_ms[gene][0][i])
		height_up.append(genes_hg[gene][0][i])
		
	    for i in range (0, len (genes_peaks[gene][1])):
		list_expression_in.append(expression)
		motif_score_in.append(genes_ms[gene][1][i])
		height_in.append(genes_hg[gene][1][i])
		
	    for i in range (0, len (genes_peaks[gene][2])):
		list_expression_down.append(expression)
		motif_score_down.append(genes_ms[gene][2][i])
		height_down.append(genes_hg[gene][2][i])		    
	    
	    if (len(genes_peaks[gene][0])>0):
		single_expression_up.append(expression)
		av_motif_score_up.append(np.mean(genes_ms[gene][0]))
		av_height_up.append(np.mean(genes_hg[gene][0]))
		cum_motif_score_up.append(np.sum(genes_ms[gene][0]))
		cum_height_up.append(np.sum(genes_hg[gene][0]))
		max_motif_score_up.append(max(genes_ms[gene][0]))
		max_height_up.append(max(genes_hg[gene][0]))		
		number_of_peaks_up.append(len(genes_peaks[gene][0]))
	    
	    if (len(genes_peaks[gene][1])>0):
		single_expression_in.append(expression)
		av_motif_score_in.append(np.mean(genes_ms[gene][1]))
		av_height_in.append(np.mean(genes_hg[gene][1]))
		cum_motif_score_in.append(np.sum(genes_ms[gene][1]))
		cum_height_in.append(np.sum(genes_hg[gene][1]))
		max_motif_score_in.append(max(genes_ms[gene][1]))
		max_height_in.append(max(genes_hg[gene][1]))			
		number_of_peaks_in.append(len(genes_peaks[gene][1]))
	    
	    if (len(genes_peaks[gene][2])>0):
		single_expression_down.append(expression)
		av_motif_score_down.append(np.mean(genes_ms[gene][2]))
		av_height_down.append(np.mean(genes_hg[gene][2]))
		cum_motif_score_down.append(np.sum(genes_ms[gene][2]))
		cum_height_down.append(np.sum(genes_hg[gene][2]))
		max_motif_score_down.append(max(genes_ms[gene][2]))
		max_height_down.append(max(genes_hg[gene][2]))			
		number_of_peaks_down.append(len(genes_peaks[gene][2]))	    
	    
	print sample
	print 'Motif score'
	print len(motif_score_up) ,
	print len(motif_score_in) ,
	print len(motif_score_down)
	
	print ' List MS up ', 
	print pearsonr(motif_score_up,list_expression_up)   
	print ' List MS in ',
	print pearsonr(motif_score_in,list_expression_in) 
	print ' List MS down ',
	print pearsonr(motif_score_down,list_expression_down) 
	print ' AV MS up ',
	print pearsonr(av_motif_score_up,single_expression_up) 
	print ' AV MS in ',
	print pearsonr(av_motif_score_in,single_expression_in)  
	print ' AV MS down ',
	print pearsonr(av_motif_score_down,single_expression_down)  
	print ' Cum MS up ',
	print pearsonr(cum_motif_score_up,single_expression_up) 
	print ' Cum MS in ',
	print pearsonr(cum_motif_score_in,single_expression_in)  
	print ' Cum MS down ',
	print pearsonr(cum_motif_score_down,single_expression_down) 
	print ' Max MS up ',
	print pearsonr(max_motif_score_up,single_expression_up) 
	print ' Max MS in ',
	print pearsonr(max_motif_score_in,single_expression_in)  
	print ' Max MS down ',
	print pearsonr(max_motif_score_down,single_expression_down)  	
	
	
	
	
	print 'Height'
	print len(height_up) ,
	print len(height_in) ,
	print len(height_down) 	
	print ' List height up ',
	print pearsonr(height_up,list_expression_up)   
	print ' List height in ',
	print pearsonr(height_in,list_expression_in) 
	print ' List height down ',
	print pearsonr(height_down,list_expression_down) 
	print ' AV height up ',
	print pearsonr(av_height_up,single_expression_up) 
	print ' AV height in ',
	print pearsonr(av_height_in,single_expression_in)  
	print ' AV height down ',
	print pearsonr(av_height_down,single_expression_down)  
	print ' Cum height up ',
	print pearsonr(cum_height_up,single_expression_up) 
	print ' Cum height in ',
	print pearsonr(cum_height_in,single_expression_in)  
	print ' Cum height down ',
	print pearsonr(cum_height_down,single_expression_down)  		
	print ' Max height up ',
	print pearsonr(max_height_up,single_expression_up) 
	print ' max height in ',
	print pearsonr(max_height_in,single_expression_in)  
	print ' max height down ',
	print pearsonr(max_height_down,single_expression_down) 	
	print 'Number of peaks'

	print pearsonr(number_of_peaks_up,single_expression_up) 
	print pearsonr(number_of_peaks_in,single_expression_in)  
	print pearsonr(number_of_peaks_down,single_expression_down)  	
    
    

def september_transcription_of_operons():
    operons_file="/data/Gyrase/Pipe/K12_w3110_operons.gff"
    genes_file = open(path + "gene_expression_mu_correction.txt".decode("utf-8"), 'r')
    outfile=open(path+"operon_expression.txt", 'w+')
    operons=parse_operons_file(operons_file)  
    operon_expression={}
    genes={}
    for line in genes_file:
	line1 = line.rstrip().split('\t')
	name=line1[0]
	start = int(line1[2])
	end = int(line1[3])
	strand = int(line1[4])
	expression = float(line1[5])
	if name in genes:
	    name=name + "_"+str(start)+ "_"+str(end)
	genes[name]=[start,end,strand,expression]
	
    for operon in sorted(operons.iterkeys()):
	operon_expression[operon]=[]
	for gene in genes:
	    if genes[gene][0]>=operon-1 and genes[gene][0]<=operons[operon][0]:
		operon_expression[operon].append(genes[gene][3])

    i=0
    for operon in sorted(operons.iterkeys()):
	if len (operon_expression[operon])>=1:
	    i+=1
	    if operons[operon][1]=="+":
		strand="+1"
	    if operons[operon][1]=="-":
		strand="-1"	    
	    outfile.write ("id_" + str(i) + "\t" + "id_" + str(i) + "\t" + str(operon) + "\t" + str(operons[operon][0])+ "\t" + strand + "\t" + str(np.mean(operon_expression[operon])) + "\n" )
    outfile.close()
    genes_file.close()
    
def september_check_closely_distributed_peaks():
    peaks_coords=open(path+ "List_of_peaks_MICROCIN_CFX_OXO_2.txt", 'r')
    st=open(path + "Summary_table_edt.txt".decode("utf-8"), 'r')
    
    scores={}
    for line in st:
	if line.startswith("peak"):
	    continue
	line1 = line.rstrip().split('\t')	
	scores[int(line1[0])]=float(line1[1])	

	    
    for line in peaks_coords:
        peaks = []
        line = line.split('\t')
        sample = str(line[0])
        print sample
        heights=[]
        for peak in line[1:len(line)]:
	    peak = peak.split('-')
	    peaks.append(int(peak[0]))
	    heights.append(float(peak[1]))
	    
	
	print len(peaks)
	
	for i in range(0,len(peaks)):
	    seed=peaks[i]
	    seeds=[]
	    seeds.append(str(peaks[i]) + "-" + str(heights[i]) + "-" + str(scores[peaks[i]]))
	    j=i
	    while True:
		j+=1
		if j<len(peaks) and peaks[j]-peaks[i]<=100:
		    seeds.append(str(peaks[j]) + "-" + str(heights[j]) + "-" +str(scores[peaks[j]]))
		else:
		    break
	    if len(seeds)>=5:
		print seeds
		
    peaks_coords.close()
    st.close()
		
#case_2_FE_FE()
#case_3_FE_SW()

#outfile_list_of_peaks= open(path + "List_of_peaks.txt".decode("utf-8"), 'w+')

#case_6_AC_SW()
#peaks_distribution()
#peak_analysis()
#atgc_content()
#motif_analysis()
#parse_list_of_sites()
#test_motif_search()
#site_distribution()
#peaks_analysis2()
#TU()

#outfile_list_of_peaks= open(path + "List_of_peaks.txt".decode("utf-8"), 'w+')
#case_7_Mitya()
#peaks_distribution()
#rrna_operons()

#venn_diagram()
#pitctures_for_presentation()
#TU()
#create_score_file_for_peaks_different_motif()
#create_summary_table()
#retrive_seqs()
#peak_analysis()
#motif_analysis()
#create_score_file_for_peaks()
#test_motif_search()
#parse_list_of_sites_calculate(genomefa,"List_of_sites_th_0_all_peaks_2_calculate.txt","List_of_sites_all_peaks_2_calculate_edt.txt")
#parse_list_of_sites_calculate_rc(genomefa,"List_of_sites_th_0_all_peaks_2_calculate_rc.txt","List_of_sites_all_peaks_2_calculate_rc_edt.txt")

#motif_analysis_plasmid()

#parse_list_of_sites_calculate(pSC101,"List_of_sites_pSC101_calculate.txt","List_of_sites_pSC101_calculate_edt.txt")
#parse_list_of_sites_calculate_rc(pSC101,"List_of_sites_pSC101_calculate_rc.txt","List_of_sites_pSC101_calculate_rc_edt.txt")

#parse_list_of_sites_calculate(pBR322,"List_of_sites_pBR322_calculate.txt","List_of_sites_pBR322_calculate_edt.txt")
#parse_list_of_sites_calculate_rc(pBR322,"List_of_sites_pBR322_calculate_rc.txt","List_of_sites_pBR322_calculate_rc_edt.txt")

#parse_list_of_sites_calculate(pUC19,"List_of_sites_pUC19_calculate.txt","List_of_sites_pUC19_calculate_edt.txt")
#parse_list_of_sites_calculate_rc(pUC19,"List_of_sites_pUC19_calculate_rc.txt","List_of_sites_pUC19_calculate_rc_edt.txt")

#retrieve_seqs()
#motif_analysis_half_motif()
#parse_list_of_sites_left_calculate(genomefa,"List_of_sites_left_all_peaks_2_calculate.txt","List_of_sites_left_all_peaks_2_calculate_edt.txt")
#parse_list_of_sites_right_calculate(genomefa,"List_of_sites_right_all_peaks_2_calculate.txt","List_of_sites_right_all_peaks_2_calculate_edt.txt")

#motif_analysis_random()
#intersect_lists_of_peaks_and_sites()

#filter_sites(6.88)
#parse_list_of_sites()
#site_distribution()

#peaks_analysis2()

#TU()
#combine_bed_graph()
#check_transcription_level()

#intersect_lists_of_peaks()
#intersect_lists_of_peaks_and_sites()
#rrna_operons()
#create_file()

#parse_list_of_sites()
#ROC_curves()

#new_files()

#September 2017
#venn_diagram()
#september_create_motif_pictures_and_files()
#motif_analysis()
#motif_analysis_plasmid()

'''
parse_list_of_sites_calculate(genomefa,"List_of_sites_th_0_all_peaks_3_calculate.txt","List_of_sites_all_peaks_3_calculate_edt.txt")
parse_list_of_sites_calculate_rc(genomefa,"List_of_sites_th_0_all_peaks_3_calculate_rc.txt","List_of_sites_all_peaks_3_calculate_rc_edt.txt")

parse_list_of_sites_calculate(pSC101,"List_of_sites_pSC101_calculate.txt","List_of_sites_pSC101_calculate_edt.txt")
parse_list_of_sites_calculate_rc(pSC101,"List_of_sites_pSC101_calculate_rc.txt","List_of_sites_pSC101_calculate_rc_edt.txt")

parse_list_of_sites_calculate(pBR322,"List_of_sites_pBR322_calculate.txt","List_of_sites_pBR322_calculate_edt.txt")
parse_list_of_sites_calculate_rc(pBR322,"List_of_sites_pBR322_calculate_rc.txt","List_of_sites_pBR322_calculate_rc_edt.txt")

parse_list_of_sites_calculate(pUC19,"List_of_sites_pUC19_calculate.txt","List_of_sites_pUC19_calculate_edt.txt")
parse_list_of_sites_calculate_rc(pUC19,"List_of_sites_pUC19_calculate_rc.txt","List_of_sites_pUC19_calculate_rc_edt.txt")

parse_list_of_sites_calculate(pMu,"List_of_sites_pMu_calculate.txt","List_of_sites_pMu_calculate_edt.txt")
parse_list_of_sites_calculate_rc(pMu,"List_of_sites_pMu_calculate_rc.txt","List_of_sites_pMu_calculate_rc_edt.txt")
'''
#september_calculate_correlation()
#september_get_nearby_genes()
#TU(path+ "List_of_peaks_MICROCIN_CFX_OXO_2.txt")
#TU(path+ "List_of_peaks_merged_3.txt")

#create_score_file_for_peaks("3")
#create_score_file_for_peaks_different_motif()
#create_summary_table()

#september_transcription_of_operons()
#september_get_genes_for_transcription()
#september_check_closely_distributed_peaks()
case_7_Mitya()