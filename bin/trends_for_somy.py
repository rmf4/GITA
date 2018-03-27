
import glob
from math import sqrt
import collections
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import statistics
from Bio import SeqIO
import argparse
import subprocess
from subprocess import Popen

def main():

	parser = argparse.ArgumentParser(description='Tool Help')
	
	
	parser.add_argument("-r","--reference", help="Genome reference .fasta file")
	parser.add_argument("-i1","--bam1", help="Case .bam file")
	parser.add_argument("-i2","--bam2", help="Control .bam file")
	parser.add_argument("-w","--window-size", help="Window factor size")
	parser.add_argument("-p","--ploidy", help="Ploidy")

	args = parser.parse_args()
	bam_files = [args.bam1,args.bam2]
	
	genome_reference_fasta = args.reference
	window = int(args.window_size)

	result = make_depth_files(window,bam_files,genome_reference_fasta)

	chroms = result


	ploidy = int(args.ploidy)

	somy(chroms,bam_files,ploidy,window)

	
	
def make_depth_files(window,bams_files,genome_reference_fasta):
	chroms = {}
	list_chroms = []
	for seq_record in SeqIO.parse(genome_reference_fasta,"fasta"):
		chroms[str(seq_record.id)] = len(str(seq_record.seq))
		

	sorted_chroms = collections.OrderedDict(sorted(chroms.items()))
	for c,l in sorted_chroms.items():
		list_chroms.append(c)


	bams = bams_files
	sorted_bams = []

	for bam in bams:
		name = bam[:bam.find(".bam")]
		sorted_bams.append(name+"_sorted.bam")
		pipe = Popen("samtools",stdout=subprocess.PIPE,stderr=subprocess.PIPE,stdin=subprocess.PIPE,shell=True)
		pipe.wait()
		version = pipe.communicate()[1].split("\n")[2]
		if version == 'Version: 0.1.19-96b5f2294a':
			pipe = Popen("samtools sort "+bam+" "+name+"_sorted",shell=True)
			pipe.wait()
		else:
			pipe = Popen("samtools sort -O bam "+bam+" > "+name+"_sorted.bam",shell=True)
			pipe.wait()

		pipe = Popen("samtools index "+name+"_sorted.bam",shell=True)
		pipe.wait()	

	for bam in sorted_bams:
		name = bam[:bam.find("_")]
		for c,l in sorted_chroms.items():
			pipe = Popen("samtools depth -r "+c+":0-"+str(l)+" "+bam+" > "+name+"_"+c+".cov",shell=True)
			pipe.wait()
	return list_chroms

def search_cov_files(path):
	return glob.glob(path)

def calculate_median_genome(cov_files,window):
	tam = len(cov_files)
	median = []
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

def list_chroms(genome_reference_fasta):
	chroms = {}
	list_chroms = []
	for seq_record in SeqIO.parse(genome_reference_fasta,"fasta"):
		chroms[str(seq_record.id)] = len(str(seq_record.seq))
		

	sorted_chroms = collections.OrderedDict(sorted(chroms.items()))
	for c,l in sorted_chroms.items():
		list_chroms.append(c)
	return list_chroms

def somy(chroms,bams,ploidy,window):

	case_name = bams[0][bams[0].rfind("/")+1:bams[0].find(".bam")]
	control_name = bams[1][bams[1].rfind("/")+1:bams[1].find(".bam")]

	somy_table = open("somy.txt","w")
	somy_table.write("\t"+case_name+"\t"+control_name+"\nSomy Value\t")
	somy_table.flush()
	
	'''
	genome median
	'''

	
	prefix = bams[1]
	prefix = prefix[:prefix.find(".bam")]
	files_case = search_cov_files(prefix+"*.cov")

	median_case = calculate_median_genome(files_case,window)
	median_genome_case = statistics.median(median_case)


	prefix = bams[0]
	prefix = prefix[:prefix.find(".bam")]

	files_control = search_cov_files(prefix+"*.cov")


	median_control = calculate_median_genome(files_control,window)
	median_genome_control = statistics.median(median_control)



	control = []
	case = []

	
	'''
	chrom median
	'''
	for c in chroms:

		prefix = bams[1]
		prefix = prefix[:prefix.find(".bam")]

		files_case = search_cov_files(prefix+"_"+c+".cov")
		median_case = calculate_median_chrom(files_case,window) 
		median_chrom_case = statistics.median(median_case)
		
		somy_case = ploidy * (float(median_chrom_case)/median_genome_case)
		somy_table.write(str(somy_case)+"\t")
		somy_table.flush()
		case.append(somy_case)

		prefix = bams[0]
		prefix = prefix[:prefix.find(".bam")]
		files_control = search_cov_files(prefix+"_"+c+".cov")
		median_control = calculate_median_chrom(files_control,window) 
		median_chrom_control = statistics.median(median_control)		

		somy_control = ploidy * (float(median_chrom_control)/median_genome_control)
		somy_table.write(str(somy_control)+"\n")
		somy_table.flush()
		control.append(somy_control)



	fig, ax = plt.subplots()
	n = len(chroms)
	ind = np.arange(n)
	width = 0.35 


	p1 = ax.bar(ind , control,   width, color='red',label=control_name)
	p2 = ax.bar(ind + width, case, width, color='blue',label=case_name)

	label = ()

	for f in range(1,n+1):
		label += (f,)

	plt.xticks(ind+width/2)
	ax.set_xticklabels(label, fontsize=6)

	plt.xlabel('Chromosomes')
	plt.ylabel('Somy')
	plt.title('Somy Analysis')
	legend = ax.legend(loc='center left', shadow=True,bbox_to_anchor=(1, 0.5))
	plt.savefig("Somy_analisys.png", bbox_inches='tight')


if __name__ == "__main__":
	main()