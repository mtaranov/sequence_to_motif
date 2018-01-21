#Convert matrix mapping sequences to motifs to a matrix of sequences to transcription factors
import argparse
import numpy as np
import pdb
import os 

def parse_args():
    parser=argparse.ArgumentParser(description='provide a directory of sequence-to-motif matrices and associated annotation files, as obtained by running "scan_motifs.py", and this script will map the motifs to transcription factors')
    parser.add_argument('--input_dir',help="directory containing the output from scan_motifs.py script")
    parser.add_argument('--num_hits_per_motif',type=int,help="maximum instances of a single motif in a sequence",default=1)
    parser.add_argument('--motif_names_file',help="file containing the list of motif names along the columns of the SxM matrix", default="motif_names.txt")
    return parser.parse_args()


def get_tf_to_motif_dict(input_dir,motif_file_name):
    motifs=open("/".join([input_dir,motif_file_name]),'r').read().strip().split('\n')
    tf_to_motifs=dict()
    for i in range(len(motifs)):
        motif =motifs[i]
        if motif not in tf_to_motifs:
            tf_to_motifs[motif]=[i]
        else:
            tf_to_motifs[motif].append(i)
    return tf_to_motifs


def main():
    args=parse_args()
    tf_to_motifs=get_tf_to_motif_dict(args.input_dir,args.motif_names_file)
    matrix_file_names=['/'.join([args.input_dir,i]) for i in os.listdir(args.input_dir) if i.endswith('.mat.npy')]
    num_tfs=len(tf_to_motifs.keys())
    tfs=tf_to_motifs.keys()
    #write the tf list to file
    outf_tf_names=open(args.input_dir+'/tf_names.txt','w')
    outf_tf_names.write('\n'.join(tfs))
    #read in the motif SxM matrix
    for matrix_file_name in matrix_file_names:
        matrix=np.load(matrix_file_name)
        num_seqs=matrix.shape[0] 
        #create a collapsed SxTF matrix
        tf_mat=np.zeros((num_seqs,num_tfs*args.num_hits_per_motif))
        for tf_index in range(len(tfs)):
            tf=tfs[tf_index]
            motif_columns=[matrix[:,i*args.num_hits_per_motif:(i+1)*args.num_hits_per_motif] for i in tf_to_motifs[tf]]
            pairwise_max=np.max(motif_columns,axis=0)
            tf_mat[:,tf_index*args.num_hits_per_motif:(tf_index+1)*args.num_hits_per_motif]=pairwise_max
        #save to output file!
        current_file=matrix_file_name.replace('mat','tf')
        np.save(current_file,tf_mat)
        print("processed:"+str(matrix_file_name))
        
if __name__=="__main__":
    main() 
