#!/usr/bin/python
import sys, getopt
import time
import fileinput
import re
import pprint
import operator
from collections import defaultdict
from tempfile import NamedTemporaryFile
from subprocess import call
from os import remove
from Queue import Queue
from threading import Thread
import numpy as np
import logging
import Levenshtein
#import networkx as nx

MIN_CLUSTER = 5
MIN_HOMOPOL_CLUSTER = 3
MIN_HOMOPOL_LENGTH = 8
MAX_GROUPS = 5

def main(argv):
	
#	try:
#		opts, args = getopt.getopt(argv,"h1:2:",["list1=","list2="])	
#	except getopt.GetoptError:
#		print 'test.py -1 <file1> -2 <file2>'
#		sys.exit(2)
#	for opt, arg in opts:
#		if opt == '-h':
#			print 'test.py -1 <file1> -2 <file2>'
#			sys.exit()
#		elif opt in ("-1", "--list1"):
#			file1 = arg
#			name1 = file1.rsplit(".",1)[0]
#		elif opt in ("-2", "--list2"):
#			file2 = arg
#			name2 = file2.rsplit(".",1)[0]
#			
#D08NBACXX110923:4:2305:3011:115838      147     1       11963   0       101M    =       11725   -339    GCCTTTGAGAGGTCACAGGGTCTTGATGCTGTGGTCTTCATCTGCAGGTGTCTGACTTCCAGCAACTGCTGGCCTGTGCCAGGGTGCAAGCTGAGCACTGG   DCAA:CDCDDDCDDCDDBDCCCCCEDFDDHGHHFIIIIIHEFDCCDGGIIIJIHEIHEIGJIGGFHEGHCFGGIGHGBIGGGGHHGIIDGF>BEDDAD@?8
#C04N3ACXX110925:1:2304:13680:152165     147     1       11963   0       101M    =       11720   -344    GCCTTTCAGAGGTCACAGGGTCTTGATGCTGTGGTCTTCATCTGCAGGTGTCTGCCTTCCAGCAACTGCTGGCCTGTGCCAGGGTGCAAGCTGAGCACTGG   #########DCC>553;?;;(;;33:CABB?;.)7=.=)3G@78C3HFB=(?B0(*?0*?BED?0):D@GB=ECA3FC<:DBA+<A+<++3??BA2++:11
#D08NBACXX110923:4:2308:4472:182599      147     1       11964   0       101M    =       11730   -335    CCTTTGAGAGGTCACAGGGTCTTGATGCTGTGGTCTTCATCTGCAGGTGTCTGACTTCCAGCAACTGCTGGCCTTTGCCAGGGTGCAAGCTGAGCACTGGA   <CDCCDDCCCDDDDDDBC?CCCEC@FEEBHC>HE>CGGJIIGBIGHGIIHGFFFCIIHGEGH@JIIIIIGEGFA2HCIJIGFJFJGJIHHDHHFDFFFCCC
	NUM_THREADS = 2

	#clusters_queue = Queue()
	#for i in range(NUM_THREADS):
	#	worker = Thread(target=worker_cluster, args=(i, clusters_queue,))
	#	worker.setDaemon(True)
	#	worker.start()
	logging.basicConfig(level=logging.DEBUG,format='[%(levelname)s] (%(threadName)-10s) %(message)s',)


	MIN_MISSING_ADJ = 5
	MAX_DIST_READS = 100
	sorted = 1;
	last_line = [];
	cluster = defaultdict(list)
	cluster_hsupport = defaultdict(int)
	line_count = 0

	for line in fileinput.input():
		line_count = line_count+1
		tokens = line.split()

		if ( line_count % 10000000 == 0 ):
			sys.stderr.write(str(line_count)+" "+time.strftime("%Y-%m-%d %H:%M")+" "+tokens[2]+" "+tokens[3]+" "+'\n')


		if ( last_line ):
			#Is this sam file ordered?
			if ( int(last_line[3]) > int(tokens[3]) and last_line[2] == tokens[2] ):
				exit(1)

			if ( tokens[5] != "*" ):
				matches = []
				num = ""
				length_matches = 0
				for char in tokens[5]:
					if ( char.isdigit() ):
						num += char
					else:
						matches.append([int(num),ord(char)])
						num = ""
						length_matches += 1
			
				#Is this a partial alignment?
				first_char = matches[0][1]
				last_char = matches[-1][1]

				if ( length_matches >= 2 and ( first_char == 83 or last_char == 83 or first_char == 72 or last_char == 72 ) ):
					start = 0
					count = 0
					coord_count = 0

					for segment in matches:
						type = segment[1]
						val = segment[0]
			
						start = int(tokens[3])+int(coord_count)
						if ( type == 77 or type == 73): # type == M or I
							count += val
							coord_count += val
						elif ( type == 68 ): # type == D
							coord_count += val
						elif ( type == 83 ): # type == S
							if ( count == 0 ):
								seq = revcomp(tokens[9][count:count+val])
								qual = tokens[10][count:count+val][::-1]
								extremity = "down"
								cluster[start].append(tokens+[extremity]+[matches]+[seq]+[qual]);
							else:
								seq = tokens[9][count:count+val]
								qual = tokens[10][count:count+val]
								extremity = "up"
								cluster[start].append(tokens+[extremity]+[matches]+[seq]+[qual]);
							count += val
						elif ( type == 72 ): # type == H
							cluster_hsupport[start]+=1
					
					#Cluster reads by breakpoint coordinate (start);
					cremove_keys = []
					#logging.debug(str(len(cluster)))
					for ckey in cluster:
						#If no more reads can be added to this cluster, proccess it.
						if ( ckey+MAX_DIST_READS < int(tokens[3]) ):
							#if ( len(cluster[ckey]) >= MIN_CLUSTER or len(cluster[ckey]) >= MIN_HOMOPOL_CLUSTER):
							softclips = len(cluster[ckey])
							hardclips = cluster_hsupport[ckey]

							#logging.debug(str(ckey)+" "+str(softclips)+" "+str(hardclips))
							#if ( softclips + hardclips >= MIN_CLUSTER and softclips >= 3):
							if ( softclips >= MIN_CLUSTER ):
							#	clusters_queue.put((cluster[ckey],ckey,hardclips))	
								ProcessCluster(cluster[ckey],ckey,hardclips)
							cremove_keys.append(ckey)
						
					for ckey in cremove_keys:
						del cluster[ckey]
						del cluster_hsupport[ckey]

					##########
					#
					# DEBUG: Print clusters
					#
					##########
					#print "-"+tokens[3]
					#for ckey in cluster:
					#	logging.debug("-"+str(ckey)+":")
					#	for read in cluster[ckey]:
					#		logging.debug("-"+str(read))
					#print "--"
		last_line = tokens
		
	for ckey in cluster:
		softclips = len(cluster[ckey])
		hardclips = cluster_hsupport[ckey]

		if ( softclips >= MIN_CLUSTER ):
			ProcessCluster(cluster[ckey],ckey,hardclips)
	
def ProcessCluster(cluster,start,hsupport):
	#print cluster 
	seqs2 = []
	quals2 = []
	align_quals = []
	extremities = []
	ids = []
	count = 1;

	for alignment in ( cluster ):
		#Trim bases with low quality
		s = 0
		for c in alignment[-1]:
			if ( ord(c) <= 40 ):
				s=s+1
			else:
				break
		seq = alignment[-2][s:]
		qual = alignment[-1][s:]
		e = 0
		for c in alignment[-1][::-1]:
			if ( ord(c) <= 40 ):
				e=e+1
			else:
				break

		if ( e != 0 ):
			seqs2.append(seq[:-e]);
			quals2.append(qual[:-e]);
		else:
			seqs2.append(seq);
			quals2.append(qual);

		align_quals.append(alignment[4])
		extremities.append(alignment[-4])
		ids.append(alignment[0])
		count+=1

	#logging.debug(seq2)
	clusters = f_clust(seqs2,quals2,align_quals,extremities)
	#clusters = []
				

	nclusters = len(clusters)
	if ( nclusters > 0 ):
		for c_id,c_seq,c_sup,c_aqual,c_ext,c_poly,c_polymax in clusters:
			if ( len(c_seq[0]) >= 10 ):
				if ( c_poly*1.0/len(c_seq[0]) <= 0.5 ):
					print cluster[0][2]+" "+str(start)+" "+str(c_id)+" "+str(c_sup+int(hsupport/nclusters))+" "+c_seq[0]+" "+c_seq[1]+" "+"_".join(c_ext)+" "+"_".join(c_aqual)
				else:
					print cluster[0][2]+" "+str(start)+" "+str(c_id)+" "+str(c_sup+int(hsupport/nclusters))+" PolyT("+str(c_polymax)+")"+c_seq[0]+" "+c_seq[1]+" "+"_".join(c_ext)+" "+"_".join(c_aqual)
	#print "----"
	#else:
		#try homopol cluster
	#	homopol_cluster = f_homopol(seqs2,quals2,align_quals,extremities)
	#	if ( len(homopol_cluster) > 0 ):
	#		for c_id,c_seq,c_sup,c_aqual,c_ext in homopol_cluster:
	#			print cluster[0][2]+" "+str(start)+" "+str(c_id)+" "+str(c_sup)+" "+c_seq[0]+" "+c_seq[1]+" "+"_".join(c_ext)+" "+"_".join(c_aqual)

def LevenshteinDistance_old(s, t):
	# degenerate cases
	if (s == t):
		return 0;
	if ( len(s) == 0 ):
		return len(t)
	if ( len(t) == 0 ):
		return len(s);

	# create two work vectors of integer distances
	v0 = range(len(t) + 1);
	v1 = [0]*(len(t) + 1);
	return 0

	# initialize v0 (the previous row of distances)
	# this row is A[0][i]: edit distance for an empty s
	# the distance is just the number of characters to delete from t
	for i in range(len(s)):
		# calculate v1 (current row distances) from the previous row v0

		# first element of v1 is A[i+1][0]
		#   edit distance is delete (i+1) chars from s to match empty t
		v1[0] = i + 1;

		# use formula to fill in the rest of the row
		for j in range(len(t)):
			cost = 0
			if (s[i] == t[j]):
				cost =  0
			else:
				cost = 1;
			v1[j + 1] = min(v1[j] + 1, v0[j + 1] + 1, v0[j] + cost);

		# copy v1 (current row) to v0 (previous row) for next iteration
		for j in range(len(v0)):
			v0[j] = v1[j];

	return v1[len(t)];

def LevenshteinDistance_faster1(s1, s2):
    if len(s1) > len(s2):
        s1, s2 = s2, s1

    distances = range(len(s1) + 1)
    for i2, c2 in enumerate(s2):
        distances_ = [i2+1]
        for i1, c1 in enumerate(s1):
            if c1 == c2:
                distances_.append(distances[i1])
            else:
                distances_.append(1 + min((distances[i1], distances[i1 + 1], distances_[-1])))
        distances = distances_
    return distances[-1]

def LevenshteinDistance(source, target):
	if len(source) < len(target):
		return LevenshteinDistance(target, source)

    # So now we have len(source) >= len(target).
	#if len(target) == 0:
	#	return len(source)

    # We call tuple() to force strings to be used as sequences
    # ('c', 'a', 't', 's') - numpy uses them as values by default.
	source = np.array(tuple(source))
	target = np.array(tuple(target))

    # We use a dynamic programming algorithm, but with the
    # added optimization that we only need the last two rows
    # of the matrix.
	previous_row = np.arange(target.size + 1)
	for s in source:
        # Insertion (target grows longer than source):
		current_row = previous_row + 1

        # Substitution or matching:
        # Target and source items are aligned, and either
        # are different (cost of 1), or are the same (cost of 0).
		current_row[1:] = np.minimum(current_row[1:],np.add(previous_row[:-1], target != s))

        # Deletion (target grows shorter than source):
		current_row[1:] = np.minimum(current_row[1:],current_row[0:-1] + 1)

		previous_row = current_row

	return previous_row[-1]

def homo_start(sequence):
	count_A = 0
	count_T = 0

	count = 0
	first = 0
	second = 0
	for char in sequence:
		if ( char == 'A' ):
			count_A += 1
			if ( first ):
				first = 0
		else:
			if ( first == 1 ):
				second = 1
			else:
				first = 1
			if (second == 1 ):
				if ( count >= 2 ):
					break
		count += 1

	count = 0
	first = 0
	second = 0
	for char in sequence:
		if ( char == 'T'):
			count_T += 1
			if ( first ):
				first = 0
		else:
			if ( first == 1 ):
				second = 1
			else:
				first = 1
			if (second == 1 ):
				if ( count >= 2 ):
					break
		count += 1

	return max(count_A,count_T)

def f_homopol(sequences,qualities,align_quals,extremities):
	homopol_groups = [];
	homopol_qual_groups = [];
	homopol_aqual_groups = [];
	homopol_ext_groups = [];
	support = 0
	sort_seq = sorted(enumerate(sequences),key=lambda x:len(x[1]), reverse=True)

	for index, seq in sort_seq:
		qual = qualities[index]
		a_qual = align_quals[index]
		ext = extremities[index]
		if ( homo_start(seq) >= MIN_HOMOPOL_LENGTH ):
			homopol_groups.append(seq)
			homopol_qual_groups.append(qual)
			homopol_aqual_groups.append(a_qual)
			homopol_ext_groups.append(ext)
			support += 1
	return_list = []
	#print "Searching for Homopols"


	if ( len(homopol_groups) >= MIN_HOMOPOL_CLUSTER ):
		return_list.append([0,consensus(homopol_groups,homopol_qual_groups),support,homopol_aqual_groups,homopol_ext_groups])
	return return_list

def f_clust(sequences,qualities,align_quals,extremities):
	count_seq = 0
	sort_seq = sorted(enumerate(sequences),key=lambda x:len(x[1]), reverse=True)
	groups = defaultdict(list);
	groups_qual = defaultdict(list);
	groups_aqual = defaultdict(list);
	groups_ext = defaultdict(list);
	groups_poly = defaultdict(int);
	groups_polymax = defaultdict(int);

	#notinc = range(len(sequences))
	first = 1

	for index, seq in sort_seq:
		qual = qualities[index]
		a_qual = align_quals[index]
		ext = extremities[index]
		added = False
		min_dist = 9999
		min_idx = -1

		#Remove PolyT from the beggining of the sequence
		max_poly = 0
		poly_count = 0
		poly_true = 0
		for c in seq:
			if ( c == 'T' ):
				poly_count=poly_count+1
			else:
				break
		if ( poly_count >= 3 ):
			seq=seq[poly_count:]
			poly_true = 1
		
	

		#Create first group
		if ( first == 1 ):
			groups[0].append(seq)
			groups_qual[0].append(qual)
			groups_aqual[0].append(a_qual)
			groups_ext[0].append(ext)
			groups_poly[0]=int(groups_poly[0])+poly_true
			if ( poly_true and groups_polymax[0] < poly_count ):
				groups_polymax[0]=poly_count
			first = 0
		else:
			if ( len(groups) >= MAX_GROUPS ):
				break
			for gp,gplist in groups.iteritems():
				#Make both strings the same length
				seq1 = gplist[0][:len(seq)]
				#dist = LevenshteinDistance(seq1,seq)
				dist = Levenshtein.distance(seq1,seq)
				#logging.debug(seq+" "+seq1+" "+str(dist))
				#logging.debug(groups)
				#dist = 0
				if ( dist < min_dist and ( dist < 3 or dist < round(len(seq)*0.25) ) ):
					min_dist = dist
					min_idx = gp
					added = True
				if ( dist == 0 ):
					break
			
			if (not added):
				#create group
				idx = len(groups)
				groups[idx].append(seq)
				groups_qual[idx].append(qual)
				groups_aqual[idx].append(a_qual)
				groups_ext[idx].append(ext)
				groups_poly[idx]=groups_poly[idx]+poly_true
				if ( poly_true and groups_polymax[idx] < poly_count ):
					groups_polymax[idx]=poly_count
			else:
				#Append to an existing group
				groups[min_idx].append(seq)
				groups_qual[min_idx].append(qual)
				groups_aqual[min_idx].append(a_qual)
				groups_ext[min_idx].append(ext)
				groups_poly[min_idx]=groups_poly[min_idx]+poly_true
				if ( poly_true and groups_polymax[min_idx] < poly_count ):
					groups_polymax[min_idx]=poly_count

#	print "--- Std cluster ---"
#	print sequences
#	print qualities
#	print "--"
#	print groups
#	print groups_qual
#	print "----"

	index = 0
	return_list = []
	for gp,gplist in groups.iteritems():
		if ( len (gplist) >= MIN_CLUSTER ):
			return_list.append([index,consensus(gplist,groups_qual[gp]),len(gplist),groups_aqual[gp],groups_ext[gp],groups_poly[gp],groups_polymax[gp]])
			index += 1

	return return_list


def consensus(sequences,quals):
#	print "------ Consensus"
#	print sequences
#	print quals
#	print "------"
	consensus = ""
	consensus_quals = ""
	len_sequences = len(sequences)

	for nt_i in range(len(sequences[0])):
		nts = defaultdict(int);
		nts_count = defaultdict(int);
		idx = 0
		for seq in sequences:
			if ( nt_i < len(seq) and seq[nt_i] != 'N' ):
				nts[seq[nt_i]] += ord(quals[idx][nt_i]); 
				nts_count[seq[nt_i]] += 1;
			idx+=1

		biggest_c = -1
		biggest_n = ""
		biggest_q = ""

		for key,sum in nts.iteritems():
			if ( sum > biggest_c ):
				biggest_c = sum
				biggest_n = key
				biggest_q = chr(int(sum/nts_count[key]))

		consensus += biggest_n
		consensus_quals += biggest_q

	return [consensus,consensus_quals]

def revcomp(sequence):
	#Store original sequence case
	case = [0]*len(sequence);
	count = 0
	for c in sequence:
		case[count] = c.isupper()
		count+=1
	case = case[::-1]	

	#Complement bases
	seq_dict = {'A':'T','T':'A','G':'C','C':'G','N':'N'};
	tseq = sequence.upper()
	rev = "".join([seq_dict[base] for base in reversed(tseq)])

	#Update original lower/upper case
	count = 0
	for i in case:
		if ( not i ):
			rev[count]=rev[count].lower()
		count += 1;
	return rev	

if __name__ == "__main__":
	main(sys.argv[1:])

