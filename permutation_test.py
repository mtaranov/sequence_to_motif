#script that utilizes MOODS python library to scan 1kb sequences from hg19 with HOCOMOCO motifs
import argparse
import pysam
import os 
import numpy as np 
import MOODS.parsers
import MOODS.tools
import MOODS.scan
from itertools import izip
import random
import pickle 
import pdb 

def parse_args():
    parser=argparse.ArgumentParser(description='permutation test for a specified motif')
    parser.add_argument('--pwm',help='pwm_file_for_motif')
    parser.add_argument('--reference',help='reference fasta file')
    parser.add_argument('--out_prefix',help='output file prefix')
    parser.add_argument('--chrom_sizes',help='file containing sizes of chromosomes')
    parser.add_argument('--background_freqs',help='file containing the background allele frequencies')
    parser.add_argument('--p_val',type=float,help='p-value threshold for inferring motif presence in a sequence',default=0.001)
    parser.add_argument('--bin_size',type=int,help='bin sizes for scanning reference genome',default=1000)
    parser.add_argument('--slide',type=int,help='number of bases to slide each window',default=500)
    parser.add_argument('--num_tests',type=int,help='number of permutation tests',default=100)
    parser.add_argument('--chrom',help='chromosome',default=None)
    parser.add_argument('--use_pickle',help='provide a pre-generated pickle file to use in computing FDR, OTHER INPUTS WILL BE IGNORED!',default=None) 
    return parser.parse_args()

#helper function to convert numpy array to tuple (MOODS needs tuples for scanning)
def totuple(a):
    try:
        return tuple(totuple(i) for i in a)
    except TypeError:
        return a
    
def parse_chrom_sizes(chrom_sizes_file):
    chrom_sizes=dict()
    with open(chrom_sizes_file) as f:
        data=f.read().strip().split('\n')
        for line in data:
            tokens=line.split('\t')
            chrom_sizes[tokens[0]]=int(tokens[1])
    return chrom_sizes 

def parse_background_freqs(background_freqs_file):
    data=open(background_freqs_file,'r').read().strip().split('\n')
    bg_freqs=[0,0,0,0]
    for line in data:
        tokens=line.split('\t')
        if tokens[0]=="A":
            bg_freqs[0]=float(tokens[1])
        elif tokens[0]=="C":
            bg_freqs[1]=float(tokens[1])
        elif tokens[0]=="G":
            bg_freqs[2]=float(tokens[1])
        elif tokens[0]=="T":
            bg_freqs[3]=float(tokens[1])
    print str(bg_freqs)
    return tuple(bg_freqs) 

#perfom initial genome scan prior to permutation 
def initial_scan(chrom_sizes,args,reference,scanner,threshold):
    initial_collisions=0 
    potential_hits=dict()
    if args.chrom!=None:
        chroms=[args.chrom]
    else: 
        chroms=chrom_sizes.keys()
    for chrom in chroms:
        print("scanning:"+str(chrom))
        pos_start=0
        pos_end=pos_start+args.bin_size
        while pos_end < chrom_sizes[chrom]:            
            if pos_start%1000000==0: 
                print("pos_start:"+str(pos_start)+"/"+str(chrom_sizes[chrom])) 
            #get the next genome bin to scan 
            seq=reference.fetch(chrom,pos_start,pos_end)
            #scan!
            results=scanner.scan(seq)[0] #only 1 motif in our test 
            scores=[r.score for r in results]
            offsets=[r.pos for r in results]
            #check if each hit is greater than the acceptance threshold
            for hit_index in range(len(scores)):
                if scores[hit_index]> threshold:
                    #accept!
                    #get 10-bp window for motif position
                    offset_floor=10*(offsets[hit_index]/10)
                    potential_hits[tuple([chrom,pos_start,pos_end])]=dict()
                    if offset_floor in potential_hits[tuple([chrom,pos_start,pos_end])]:
                        #we have a hash collision!
                        print("Hash collision on initial scan!!")
                        initial_collisions+=1
                    potential_hits[tuple([chrom,pos_start,pos_end])][offset_floor]=dict()
                    potential_hits[tuple([chrom,pos_start,pos_end])][offset_floor]['offset']=[offsets[hit_index]]
                    potential_hits[tuple([chrom,pos_start,pos_end])][offset_floor]['score']=[scores[hit_index]]
            #update by the sliding window!
            pos_start+=args.slide
            pos_end+=args.slide
    print('completed initial scan!')
    print(str(len(potential_hits.keys())))
    return potential_hits,initial_collisions

#perform a permutation test on specified genome region 
def permutation_test(region,bins_in_region,reference,scanner,args):
    updated_bins=bins_in_region
    collisions=0 
    seq=reference.fetch(region[0],region[1],region[2])
    seq_len=len(seq)
    for test_id in range(args.num_tests):
        #permute the sequence
        shuffled = ''.join(random.sample(seq, seq_len))
        results=scanner.scan(shuffled)[0]
        scores=[r.score for r in results]
        offsets=[r.pos for r in results]
        offset_bins=[10*(offset/10) for offset in offsets]
        offset_bin_dict=dict()
        for i in range(len(offset_bins)):
            if offset_bins[i] in offset_bin_dict:
                print("Collision in permutation test!")
                #pdb.set_trace()
                collisions+=1 
            offset_bin_dict[offset_bins[i]]=[offsets[i],scores[i]]
        for bin_index in bins_in_region.keys():
            if bin_index not in offset_bin_dict:
                updated_bins[bin_index]['offset'].append(None)
                updated_bins[bin_index]['score'].append(0)
            else:
                updated_bins[bin_index]['offset'].append(offset_bin_dict[bin_index][0])
                updated_bins[bin_index]['score'].append(offset_bin_dict[bin_index][1]) 
    return updated_bins,collisions

#calculates false discovery rate from permutation test 
def calculate_fdr(potential_hits,args):
    bed_entries=[] 
    for region in potential_hits:
        chrom=region[0]
        start_pos=region[1] 
        for cur_bin in potential_hits[region]:
            #get the average offset
            offset_positions=[p for p in potential_hits[region][cur_bin]['offset'] if p is not None]
            mean_offset=sum(offset_positions)/float(len(offset_positions))
            adjusted_pos=start_pos+int(round(mean_offset))
            #compute fdr
            true_score=potential_hits[region][cur_bin]['score'][0]
            permuted_scores=np.array(potential_hits[region][cur_bin]['score'][1::])
            fdr=np.sum(permuted_scores>=true_score)/float(len(permuted_scores))
            bed_entries.append(chrom+'\t'+str(adjusted_pos)+'\t'+str(fdr))
    bed_file='\n'.join(bed_entries)
    #write the output bed file
    outf=open(args.out_prefix,'w')
    outf.write(bed_file)
    
            
            

def main():
    args=parse_args()
    if args.use_pickle!=None:
        potential_hits=pickle.load(open(args.use_pickle,'rb'))
        calculate_fdr(potential_hits,args)
        exit()

    chrom_sizes=parse_chrom_sizes(args.chrom_sizes)
    print("parsed chromosome sizes")

    bg = parse_background_freqs(args.background_freqs)
    print("got background frequencies")

    #prep the reference fasta
    reference=pysam.FastaFile(args.reference)

    #get the matrix files
    matrix_file_name=args.pwm
    matrix=totuple(np.transpose(np.loadtxt(matrix_file_name,skiprows=1)))
    motif_name=matrix_file_name.split('.')[0]
    threshold=MOODS.tools.threshold_from_p(matrix,bg,args.p_val)
    #create the moods scanner
    scanner = MOODS.scan.Scanner(7)
    scanner.set_motifs([matrix], bg, [threshold],)
    print('read in motif matrix and threshold')
    #perform preliminary scan to identify all regions of potential motif hits
    potential_hits,initial_collisions= initial_scan(chrom_sizes,args,reference,scanner,threshold)
    print("completed initial genome scan, there were: "+str(initial_collisions) + " collisions!") 
    #perform permutation test by looking only at positions where hits were observed
    counter=0
    permutation_test_collisions=0 
    for region in potential_hits:
        print(str(counter))
        potential_hits[region],collisions=permutation_test(region,potential_hits[region],reference,scanner,args)
        permutation_test_collisions+=collisions
        counter+=1 
    print("completed permutation tests, there were: "+str(permutation_test_collisions)+" collisions!")
    #pickle the results in case they must be examined later 
    print("pickling!")
    pickle.dump(potential_hits, open(args.out_prefix, 'wb'))
    #calculate FDR and generate bed file of motif hits to FDR values
    calculate_fdr(potential_hits,args) 

if __name__=='__main__':
    main()
