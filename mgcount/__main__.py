## APR-2021  Andrea Hita Ardiaca

## -- Modules import
import argparse
import os
import re
import pandas as pd
import numpy as np
import multiprocessing as mp
from tempfile import TemporaryDirectory
from gtfparse import read_gtf

## -- MGcount source code
from mgcount import hierarchassign as hra
from mgcount import multigraph as mg
from mgcount import expreport as count

def main():

    ##--------------- I/O files
    parser = argparse.ArgumentParser(
        description = """MGcount RNA-seq quantification pipeline \n""")
    
    parser.add_argument('--bam_infiles',
                        '-i',
                        type=str,
                        help = 'Alignment input file names \n',
                        required = True)
 
    parser.add_argument('--sample_id',
                        type=str,
                        help = 'SampleID input file names \n',
                        required = False, 
                        default = '')
    
    parser.add_argument('--outdir',
                        '-o',
                        type=str,
                        help='Output directory path \n',
                        required = True)
    
    parser.add_argument('--gtf',
                        type=  str,
                        help = 'Annotations file name \n',
                        required = True)
    
    parser.add_argument('-T',
                        '--n_cores',
                        type = int,
                        help = 'Number of cores for parallelization \n',
                        required = False,
                        default = 1)

    parser.add_argument('--subs',
                        type = int,
                        help = '''Optional subsapling number of alignments 
                                  to train the multigraph network \n''',
                        required = False,
                        default = 0)
                            
    parser.add_argument('-p',
                        '--paired_flag',
                        action = 'store_true',
                        help = 'Paired end flag \n')

    parser.add_argument('-s',
                        '--strand_option',
                        type = str,
                        help = '''Options available are 0: unstranded, 
                                  1: forwardstranded (default)
                                  or 2:reverse_stranded \n''',
                        choices = ['0','1','2'],
                        required = False,
                        default = '1')
    
    parser.add_argument('-th_low',
                        type = float,
                        help = '''Low minimal threshold for feature-to-feature
                                   multi-mapping fraction\n''',
                        required = False,
                        default = 0.01)    
                        
    parser.add_argument('-th_high',
                        type = float,
                        help = '''High minimal threshold for feature-to-feature
                                   multi-mapping fraction\n''',
                        required = False,
                        default = 0.75)    
    
    parser.add_argument('-seed',
                        type = float,
                        help = '''Optional fixed seed for random numbers generation 
                                  during communities detection\n''',
                        required = False,
                        default = None) 

    # parser.add_argument('--feat_r',
    #                     type = str,
    #                     help = '''GTF feature type for rRNA reads
    #                               to be assigned to\n''',
    #                     required = False,
    #                     default = 'exon')
    
    # parser.add_argument('--feat_output_r',
    #                     type = str,
    #                     help = '''GTF field name for which to summarize
    #                               expresion of rRNA assigned reads\n''',
    #                     required = False,
    #                     default = 'gene_name')
    
    # parser.add_argument('--feat_biotype_r',
    #                     type = str,
    #                     help = '''GTF field name defining biotype for
    #                               rRNA features\n''',
    #                     required = False,
    #                     default = 'gene_biotype')

    # parser.add_argument('--min_overlap_r',
    #                     type = str,
    #                     help = '''Minimal feature-alignment overlapping fraction
    #                               for assignation (default 0.75) for long round \n''',
    #                     required = False,
    #                     default = '1')
    
    # parser.add_argument('--ml_flag_r',
    #                     type = int,
    #                     help = '''Multi-loci graph detection
    #                               based groups frag for rRNA round\n''',
    #                     required = False,
    #                     default = 0)
    
    ## ---- Long round config.
    parser.add_argument('--feat_output_long',
                        type = str,
                        help = '''GTF field name for which to summarize 
                                  expresion of longRNA assigned reads\n''',
                        required = False,
                        default = 'gene_name')

    parser.add_argument('--feat_biotype_long',
                        type = str,
                        help = '''GTF field name defining biotype for
                                  longRNA features\n''',
                        required = False,
                        default = 'gene_biotype')

    parser.add_argument('--min_overlap_long',
                        type = str,
                        help = '''Minimal feature-alignment overlapping fraction
                                  for assignation (default 0.75) for long round \n''',
                        required = False,
                        default = '1')

    parser.add_argument('--ml_flag_long',
                        type = int,
                        help = '''Multi-loci graph detection based 
                                  groups flag for long round\n''',
                        required = False,
                        default = 1)

    
    ## ---- Small round config.
    parser.add_argument('--feat_small',
                        type = str,
                        help = '''GTF feature type for smallRNA reads 
                                  to be assigned to\n''',
                        required = False,
                        default = 'transcript')

    parser.add_argument('--feat_output_small',
                        type = str,
                        help = '''GTF field name for which to summarize 
                                  counts of smallRNA assigned reads\n''',
                        required = False,
                        default = 'transcript_name')

    parser.add_argument('--feat_biotype_small',
                        type = str,
                        help = '''GTF field name defining biotype for
                                  smallRNA features\n''',
                        required = False,
                        default = 'transcript_biotype')

    parser.add_argument('--min_overlap_small',
                        type = str,
                        help = '''Minimal feature-alignment overlapping fraction
                                  for assignation (default 1) for long round \n''',
                        required = False,
                        default = '1')

    parser.add_argument('--ml_flag_small',
                        type = int,
                        help = '''Multi-loci graph detection 
                                  based groups flag for small round\n''',
                        required = False,
                        default = 1)

    ## ---- featureCounts software path                    
    parser.add_argument('--featureCounts_path',
                        type = str,
                        help = ''' Path to featureCounts executable file 
                                   (by default, the software is looked for 
                                    on /user/bin/featureCounts) \n''',
                        required = False,
                        default = '/usr/bin/featureCounts')
    
    parser.add_argument('--btyperounds_filename',
                        help = '''Optional .csv file with gene_biotype 
                                  assignation_round custom defined pairs \n''',
                        required = False,
                        default = None)

    
    ##--------------- Parse arguments
    args = parser.parse_args()

    ## Required arguments
    outdir = os.path.abspath(args.outdir) + '/'
    gtf_filename = os.path.abspath(args.gtf)
    n_cores = args.n_cores

    ## Params
    end = ('Paired' if args.paired_flag else 'Single')
    strand = args.strand_option
    n_sub = args.subs
    th_high = args.th_high
    th_low = args.th_low
    seed = args.seed
    
    ## Additional optinal arguments
    fcpath = args.featureCounts_path
    btyperounds_filename = args.btyperounds_filename

    ## Get input files
    infiles = [line.rstrip() for line in open(os.path.abspath(args.bam_infiles))]
    if '' in infiles: infiles.remove('')
    
    ## Select default crounds file if non-user-input defined
    if btyperounds_filename == None:
        ##btype_crounds = config.get_btypes_by_crounds() for major compatibility
        ##btype_crounds = pd.read_csv(utils.resource_path('btypes_crounds.csv')) for compiled version
        btype_crounds = pd.read_csv(os.path.dirname(__file__) + '/data/btypes_crounds.csv')
    else: btype_crounds = pd.read_csv(btyperounds_filename)

    ## Counting rounds configuration
    crounds = pd.DataFrame({
        'r':[##'r',
             'small',
             'long_exon',
             'long_intron'],
        
        'annot':[##'r',
                 'small',
                 'long',
                 'long'],
        
        'feat':[##args.feat_r,
                args.feat_small,
                'exon',
                'gene'],
        
        'attr':[##args.feat_output_r,
                args.feat_output_small,
                args.feat_output_long,
                args.feat_output_long],
        
        'btype':[##args.biotype_r,
                 args.feat_biotype_small,
                 args.feat_biotype_long,
                 args.feat_biotype_long],
        
        'ml' : [##args.ml_flag_r,
                args.ml_flag_small,
                args.ml_flag_long,
                'Derived' if args.ml_flag_long else False],

        'minov' : [##args.minov_r,
                args.min_overlap_small,
                args.min_overlap_long,
                args.min_overlap_long],
        
        'suf' : [##'',
                 '',
                 '-exon',
                 '-intron']})
    
    ## Get gene biotypes
    gtf = read_gtf(gtf_filename)

    ## Create directory for temporary files
    if not os.path.exists(outdir): os.mkdir(outdir)
    tmpdir = TemporaryDirectory(prefix = outdir + '.mg_')
    tmppath = tmpdir.name + '/'
    
    ## keep all tmp files
    ##tmppath = outdir + 'tmp/'
    ##if not os.path.exists(tmppath): os.mkdir(tmppath)

    ## Set method to start parallel child processes to "spawn"
    mp.set_start_method("spawn")
    
    
    ## RUN COUNTING PIPELINE
    ## ---------------------------------------------------------------------------
    
    ## ----- Hierarchical read. assignation loop by sample
    crounds = hra.read_assignation_loop(infiles, outdir, fcpath, tmppath, gtf_filename,
                                        crounds, btype_crounds, n_cores, end, strand)

    ## ----- Multi-loci group extraction
    ml = mg.MG(infiles, outdir, tmppath, gtf, crounds, btype_crounds, n_cores,
               n_sub, th_low, th_high, seed)

    ## ----- Expression. level summarization
    counts = count.extract_count_matrix(infiles, outdir, tmppath, crounds, n_cores, ml)

    ## ---------------------------------------------------------------------------
    
    ##--------------- Save count matrix
    ## Save only nonempty features
    ## Will label samples by using filename fragment that is distinct in all samples
    nonempty_idx = np.where(counts.sum(axis=1)!= 0)[0].tolist()
    if len(args.sample_id) > 0:
        sid = np.array(pd.read_csv(args.sample_id, header = None)[0])
        counts.columns = sid

    counts.iloc[nonempty_idx].to_csv(
        outdir + 'count_matrix.csv', header = True, index = True)
    
    ## Save feature metadata table
    feats_metadata = count.extract_feat_metadata(outdir, tmppath, gtf, crounds, infiles[0], ml)
    feats_metadata.reindex(counts.iloc[nonempty_idx].index).to_csv(
        outdir + 'feature_metadata.csv')

    print('......')
    print('Voilà.')
    print('Olé!')
    return 0

if __name__ == "__main__":
    main()

