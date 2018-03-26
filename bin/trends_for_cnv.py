from __future__ import division
import glob
import collections
import matplotlib.pyplot as plt
import statistics
from math import log,ceil
from sympy import *
import argparse
from Bio import SeqIO
from subprocess import Popen
import os.path

def main():

	parser = argparse.ArgumentParser(description='Tool Help')
	
	parser.add_argument("-gff","--annotation-file", help="GFF file")
	parser.add_argument("-r","--reference", help="Genome reference .fasta file")
	parser.add_argument("-i1","--bam1", help="Case .bam file")
	parser.add_argument("-i2","--bam2", help="Control .bam file")
	parser.add_argument("-w","--window-factor-size", help="Window factor size")

	args = parser.parse_args()
	bam_files = [args.bam1,args.bam2]
	
	genome_reference_fasta = args.reference
	window = args.window_factor_size
	annotation_file = args.annotation_file

	result = make_depth_files(window,bam_files,genome_reference_fasta)

	chroms = result[0]
	windows = result[1]

	if os.path.exists("cnv.txt"):
		new_file = open("cnv.txt","w")
		new_file.close()
	
	for i in range(len(chroms)):
		cnv(chroms[i],windows[i],annotation_file,bam_files)
	

def make_depth_files(window,bams_files,genome_reference_fasta):
	chroms = {}
	list_chroms = []
	list_window = []
	for seq_record in SeqIO.parse(genome_reference_fasta,"fasta"):
		chroms[str(seq_record.id)] = len(str(seq_record.seq))
		

	sorted_chroms = collections.OrderedDict(sorted(chroms.items()))
	for c,l in sorted_chroms.items():
		list_chroms.append(c)
		list_window.append(ceil(l/float(window)))

	bams = bams_files
	sorted_bams = []

	for bam in bams:
		name = bam[:bam.find(".bam")]
		sorted_bams.append(name+"_sorted.bam")
		
		pipe = Popen("samtools sort -O bam "+bam+" > "+name+"_sorted.bam",shell=True)
		pipe.wait()	
		pipe = Popen("samtools index "+name+"_sorted.bam",shell=True)
		pipe.wait()	



	for bam in sorted_bams:
		name = bam[:bam.find("_")]
		for c,l in sorted_chroms.items():
			pipe = Popen("samtools depth -r "+c+":0-"+str(l)+" "+bam+" > "+name+"_"+c+".cov",shell=True)
			pipe.wait()
	return list_chroms,list_window

def search_cov_files(path):
	
	return glob.glob(path)

def calculate_median_chrom(cov_files,window):

	median = []
	tam = len(cov_files)
	for i in range(tam):
		f = open(cov_files[i])

		count = 0
		depth = []
		
		for line in f:
			d = line.split()
			prof = int(d[2])
			if count < window:
				depth.append(prof)
				count += 1
			else:
				depth.append(prof)
				median.append(statistics.median(depth))
				depth = []
				count = 0
		f.close()

	return median

def calculate_median_segment(cov_files,window,median_chrom):

	tam = len(cov_files)
	cnv_list = []
	pos_list = []
	for i in range(tam):
		f = open(cov_files[i])

		count = 0
		depth = []

		for line in f:
			d = line.split()
			prof = int(d[2])
			pos = int(d[1])
			if count == 0:
				pos_list.append(pos)

			if count < window:
				depth.append(prof)
				count += 1
			else:
				depth.append(prof)
				median_segment = statistics.median(depth)
				depth = []
				count = 0
						
				cnv_value =  (median_segment/float(median_chrom))
				cnv_list.append(cnv_value)
		f.close()
	return cnv_list

def solve_trend(iterations,cnv_case,cnv_control,window,k,annotation_file):
	cnv_file = open("cnv.txt","a+")

	x = []
	y = []
	gene_x = []
	gene_y = []
	a = Symbol('a')
	b = Symbol('b')
	
	for i in range(iterations):
		
		sol = solve( [( ( (a*cnv_control[i]-b*cnv_case[i])/(a*cnv_control[i]) )*100) ],[a,b])
		try:
			u = float(str(sol[a])[:str(sol[a]).find("*") ])
			y.append(log(u,10))
			x.append((i+1)*window )
			gff_file = open(annotation_file)
			if log(u,10) >= 0.1 or log(u,10) <= -0.1:
				for gff_line in gff_file:
					s = gff_line.split()
					if len(s) == 9 and s[0] == k and s[2]=="gene":
						start = int(s[3])
						end = int(s[4])
						gene_id = s[8][s[8].find("=")+1:s[8].find(";")]
						descr = s[8].split(";")[2]
						descr = descr[descr.find("=")+1:]
						descr = descr.replace("+"," ")
						descr = descr.replace("%2C",",")
						descr = descr.replace("%2C",",")
						descr = descr.replace("%28","(")
						descr = descr.replace("%29",")")
						desr = descr.replace("%2B","+")
						descr = descr.replace("%2F","/")
						descr = descr.replace("%27","'")
						if s[6] == "+":
							if (i+1)*window >= start and (i+1)*window <= end:
								cnv_file.write(k+"\t"+gene_id+"\t"+descr+"\t+\t"+str((i+1)*window)+"\t"+str((i+2)*window)+"\t"+str(log(u,10))+"\n")
								cnv_file.flush()
								gene_x.append((i+1)*window)
								gene_y.append(log(u,10))
						else:
							if (i+1)*window >= start and (i+1)*window <= end:
								cnv_file.write(k+"\t"+gene_id+"\t"+descr+"\t-\t"+str((i+1)*window)+"\t"+str((i+2)*window)+"\t"+str(log(u,10))+"\n")
								cnv_file.flush()

								gene_x.append((i+1)*window)
								gene_y.append(log(u,10))


		except Exception,e:
			print e
			pass
	cnv_file.close()
	return gene_x,gene_y,x,y

def cnv(k,w,annotation_file,bams):
	
	fig, ax = plt.subplots()

	'''
	chromosome median
	'''
	
	
	window = w
	
	prefix = bams[1]
	prefix = prefix[:prefix.find(".bam")]
	files_case = search_cov_files(prefix+"_"+k+".cov")
	
	median_case = calculate_median_chrom(files_case,window)
	median_chrom_case = statistics.median(median_case)
	
	prefix = bams[0]
	prefix = prefix[:prefix.find(".bam")]
	files_control = search_cov_files(prefix+"_"+k+".cov")
	
	median_control = calculate_median_chrom(files_control,window)
	median_chrom_control = statistics.median(median_control)


	cnv_control = []
	cnv_case = []

	chroms = [k]

	'''
	segment median
	'''

	#for c in chroms:
	prefix = bams[1]
	prefix = prefix[:prefix.find(".bam")]

	files_case = search_cov_files(prefix+"_"+k+".cov")
	prefix = bams[0]
	prefix = prefix[:prefix.find(".bam")]
	files_control = search_cov_files(prefix+"_"+k+".cov")
	
		

	tam = 1
	cnv_case = calculate_median_segment(files_case,window,median_chrom_case)
	cnv_control = calculate_median_segment(files_control,window,median_chrom_control)




	x = []
	y = []



	iterations = 0
	

	if len(cnv_control) > len(cnv_case):
		iterations = len(cnv_case)
	else:
		iterations = len(cnv_control)

	for i in range(iterations):
		x.append((i+1)*window )
		y.append(cnv_control[i])
	plt.plot(x,y,label="control")
	
	x = []
	y = []

	for i in range(iterations):
		x.append((i+1)*window )
		y.append(cnv_case[i])
	plt.plot(x,y,label="case")


	trend_analisys = solve_trend(iterations,cnv_case,cnv_control,window,k,annotation_file)
	gene_x = trend_analisys[0]
	gene_y = trend_analisys[1]
	x = trend_analisys[2]
	y = trend_analisys[3]


	ax.plot(gene_x,gene_y, 'o',color="black",label=k+" genes")

	plt.plot(x,y,label="Trends")
	plt.setp(ax.get_xticklabels(), rotation=25, fontsize = 6, horizontalalignment='right')
	plt.xlabel('Position')
	plt.ylabel('Log(Solution)')
	plt.title('CNV ('+k+")")
	legend = ax.legend(loc='center left', shadow=True,bbox_to_anchor=(1, 0.5))

	plt.savefig("CNV_analisys_"+k+".png", bbox_inches='tight')




if __name__ == "__main__":
	main()
