#script that utilizes MOODS python library to scan 1kb sequences from hg19 with HOCOMOCO motifs
import argparse
import pysam
import os 
import numpy as np 
import MOODS.parsers
import MOODS.tools
import MOODS.scan
from itertools import izip 
import pdb 

def parse_args():
    parser=argparse.ArgumentParser(description='provide a director of PWM motif files and source fasta file as well the interval sizes to scan')
    parser.add_argument('--pwm_dir',help='directory of PWM files to process')
    parser.add_argument('--reference',help='reference fasta file')
    parser.add_argument('--out_prefix',help='output file prefix')
    parser.add_argument('--chrom_sizes',help='file containing sizes of chromosomes')
    parser.add_argument('--background_freqs',help='file containing the background allele frequencies')
    parser.add_argument('--p_val',type=float,help='p-value threshold for inferring motif presence in a sequence',default=0.001)
    parser.add_argument('--bin_size',type=int,help='bin sizes for scanning reference genome',default=1000)
    parser.add_argument('--num_hits_per_motif',type=int,help="maximum instances of a single motif in a sequence that will be recorded",default=1)
    parser.add_argument('--chrom',help="optionally, specify the chromosomes of interest for the parser to focus on, comma-separated list of chromosomes",default=None)
    parser.add_argument('--positions_bed',help="input bed file containing sequence positions to scan",default=None)
    parser.add_argument('--binarize',help="remove this flag if you don't want to binarize the scores",action='store_true') 
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

def calculate_bg_from_scratch(reference,positions_bed_file):
    bg_freqs=[0,0,0,0]
    numbases=0
    positions_bed=open(positions_bed_file,'r').read().strip().split('\n')
    positions_bed=[entry.split('\t') for entry in positions_bed]
    for region in positions_bed:
        seq=(reference.fetch(region[0],int(region[1]),int(region[2]))).lower()
        numA=seq.count('a')
        numC=seq.count('c')
        numG=seq.count('g') 
        numT=seq.count('t') 
        numbases=numbases+numA+numC+numG+numT
        bg_freqs[0]+=numA
        bg_freqs[1]+=numC
        bg_freqs[2]+=numG
        bg_freqs[3]+=numT
    bg_freqs=[float(i)/numbases for i in bg_freqs]
    print(str(bg_freqs))
    return bg_freqs

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

def global_scan(args,chrom_sizes,num_motifs,scanner,thresholds,reference,motif_names,output_dir):         
    #iterate through each chromosome in the reference sequence to scan for motifs
    if args.chrom==None:
        chroms=chrom_sizes.keys()
    else:
        chroms=args.chrom.split(',')
    for chrom in chroms:
        print("scanning:"+str(chrom))
        #pre-allocate the output numpy array with zeros
        num_sequence_bins=chrom_sizes[chrom]/args.bin_size
        chrom_motif_mat=np.zeros((num_sequence_bins,args.num_hits_per_motif*num_motifs))
        chrom_pos_mat=np.zeros((num_sequence_bins,2))
        
        pos_start=0
        pos_end=pos_start+args.bin_size
        bin_index=0 
        while pos_end < chrom_sizes[chrom]:
            if pos_start%1000000==0: 
                print("pos_start:"+str(pos_start)+"/"+str(chrom_sizes[chrom])) 
            #get the next genome bin to scan 
            seq=reference.fetch(chrom,pos_start,pos_end)
            #scan!
            results=scanner.scan(seq)
            #threshold the motif scores to binary values, either take the max score or the top 3 scores to threshold 
            for motif_index in range(num_motifs):
                #get the top n scores for each motif
                results_cur_motif=[r.score for r in results[motif_index]]
                results_cur_motif.sort(reverse=True)
                #pad to the desired length
                results_cur_motif+=[0]*(args.num_hits_per_motif-len(results_cur_motif))
                #truncate to the desired length
                results_cur_motif=results_cur_motif[0:args.num_hits_per_motif]
                if args.binarize==False:
                    chrom_motif_mat[bin_index,motif_index*args.num_hits_per_motif:(motif_index+1)*args.num_hits_per_motif]=results_cur_motif
                else:
                    chrom_motif_mat[bin_index,motif_index*args.num_hits_per_motif:(motif_index+1)*args.num_hits_per_motif]=[int(m>thresholds[motif_index])for m in results_cur_motif]
            chrom_pos_mat[bin_index][0]=pos_start
            chrom_pos_mat[bin_index][1]=pos_end                                 
            #update indices 
            pos_start=pos_end
            pos_end=pos_start+args.bin_size
            bin_index+=1

        #save output numpy pickles for the chromosome
        np.save("/".join([output_dir,".".join([chrom,"mat"])]),chrom_motif_mat)
        np.save("/".join([output_dir,".".join([chrom,"pos"])]),chrom_pos_mat)
        outf_names=open(output_dir+'/motif_names.txt','w')
        outf_names.write('\n'.join(motif_names))
        print("finished processing chromosome:"+str(chrom))
        
def scan_specified_positions(args,positions,num_motifs,scanner,thresholds,reference,motif_names,output_dir):
    num_sequence_bins=len(positions)
    motif_mat=np.zeros((num_sequence_bins,args.num_hits_per_motif*num_motifs))
    pos_mat=np.empty((num_sequence_bins,3),dtype=object)
    bin_index=0 
    for position in positions:
        seq=reference.fetch(position[0],int(position[1]),int(position[2]))
        #scan!
        if bin_index%10000==0:
            print(str(bin_index))
        results=scanner.scan(seq)
        for motif_index in range(num_motifs):
            #get the top n scores for each motif
            results_cur_motif=[r.score for r in results[motif_index]]
            results_cur_motif.sort(reverse=True)
            #pad to the desired length
            results_cur_motif+=[0]*(args.num_hits_per_motif-len(results_cur_motif))
            #truncate to the desired length
            results_cur_motif=results_cur_motif[0:args.num_hits_per_motif]
            if args.binarize==False:
                motif_mat[bin_index,motif_index*args.num_hits_per_motif:(motif_index+1)*args.num_hits_per_motif]=results_cur_motif
            else:
                motif_mat[bin_index,motif_index*args.num_hits_per_motif:(motif_index+1)*args.num_hits_per_motif]=[int(m>thresholds[motif_index])for m in results_cur_motif]
        pos_mat[bin_index][0]=position[0]
        pos_mat[bin_index][1]=position[1]
        pos_mat[bin_index][2]=position[2]
        bin_index+=1
    #save the output numpy pickle
    #save output numpy pickles for the chromosome
    np.save("/".join([output_dir,"mat"]),motif_mat)
    np.save("/".join([output_dir,"pos"]),pos_mat)
    outf_names=open(output_dir+'/motif_names.txt','w')
    outf_names.write('\n'.join(motif_names))
        

def main():
    args=parse_args()
    chrom_sizes=parse_chrom_sizes(args.chrom_sizes)
    print("parsed chromosome sizes")
    
    #prep the reference fasta
    reference=pysam.FastaFile(args.reference)


    if args.background_freqs!=None:
    	bg = parse_background_freqs(args.background_freqs)
    else:
        print("No background_freqs argument provided, computing frequencies from positions_bed argument") 
        bg = calculate_bg_from_scratch(reference,args.positions_bed)
    print("got background frequencies")

    #get the matrix files
    matrix_file_names=['/'.join([args.pwm_dir,i]) for i in os.listdir(args.pwm_dir) if i.endswith('.pwm')]

    #parse the pwm matrices 
    matrices=[totuple(np.transpose(np.loadtxt(f,skiprows=1))) for f in matrix_file_names]
    motif_names=[f.split('.')[0] for f in matrix_file_names]
    num_motifs=len(matrices)
    
    #get the p-value cutoff thresholds 
    thresholds=[MOODS.tools.threshold_from_p(m,bg,args.p_val) for m in matrices]

    #create the moods scanner
    scanner = MOODS.scan.Scanner(7)
    scanner.set_motifs(matrices, bg, thresholds, )
    print('read in motif matrices and thresholds')

    #create the output s_m, s_3m, s_tf directories
    output_dir=args.out_prefix+"_"+str(args.num_hits_per_motif)+"_"+str(args.p_val)
    try:
        os.makedirs(output_dir)
    except:
        print("Directory already exists:"+output_dir) 
        
    #decide whether we are scanning globally or at specific positions from a bed file 
    if args.positions_bed==None:
        global_scan(args,chrom_sizes,num_motifs,scanner,thresholds,reference,motif_names,output_dir)
    else:
        positions=[i.split('\t') for i in open(args.positions_bed,'r').read().strip().split('\n')]
        scan_specified_positions(args,positions,num_motifs,scanner,thresholds,reference,motif_names,output_dir)


if __name__=='__main__':
    main()
