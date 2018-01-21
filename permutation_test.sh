#for chrom in `seq 2 12`# X Y
#do 
#    python permutation_test.py --pwm pwm/ATF3_HUMAN.H10MO.A.pwm --reference hg19/hg19.genome.fa --out_prefix ATF3_p5e-5_permutation_test_$chrom.dat --chrom_sizes hg19/hg19.chrom.sizes --background_freqs hg19_background_freqs.txt --p_val 0.00005 --bin_size 10000 --slide 500 --num_tests 100 --chrom chr$chrom &
#done
python permutation_test.py --pwm pwm/ATF3_HUMAN.H10MO.A.pwm --reference hg19/hg19.genome.fa --out_prefix ATF3_p5e-5_permutation_test_$chrom.dat --chrom_sizes hg19/hg19.chrom.sizes --background_freqs hg19_background_freqs.txt --p_val 0.00005 --bin_size 10000 --slide 500 --num_tests 100 --chrom chr$chrom --use_pickle ATF3_p5e-5_permutation_test_Y.dat
