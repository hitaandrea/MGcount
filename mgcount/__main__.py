## APR-2021  Andrea Hita Ardiaca

## -- Modules import
import argparse
import os
import re

import pandas as pd
import numpy as np
import multiprocessing as mp

from warnings import simplefilter
from tempfile import TemporaryDirectory
from gtfparse import read_gtf
from shutil import rmtree
from scipy import sparse
from scipy.io import mmwrite

## -- MGcount source code
from mgcount import hierarchassign as hra
from mgcount import multigraph as mg
from mgcount import expreport as count

def main():

    parser = argparse.ArgumentParser(
        description = """MGcount RNA-seq quantification pipeline v1.0.2 \n""")
        
    ##--------------- Required I/O
    parser_required = parser.add_argument_group('required arguments')
    
    parser_required.add_argument('--bam_infiles',
                                 '-i',
                                 type = os.path.abspath,
                                 help = 'Alignment input file names \n',
                                 required = True)
 
    parser_required.add_argument('--outdir',
                                 '-o',
                                 type = os.path.abspath,
                                 help ='Output directory path \n',
                                 required = True)
    
    parser_required.add_argument('--gtf',
                                 type = os.path.abspath,
                                 help = 'Annotations file name \n',
                                 required = True)
    
    ##--------------- Optional arguments
    parser.add_argument('--sample_id',
                        type=str,
                        help = 'Optional sampleID names \n',
                        required = False, 
                        default = '')
    
    parser.add_argument('-T',
                        '--n_cores',
                        type = int,
                        help = 'Number of cores for parallelization \n',
                        required = False,
                        default = 1)
                            
    parser.add_argument('-p',
                        '--paired_flag',
                        action = 'store_true',
                        help = 'Paired end flag \n')

    parser.add_argument('-m',
                        '--mtx',
                        action = 'store_true',
                        help = 'Output count matrix in sparse format \n')
    
    parser.add_argument('-s',
                        '--strand_option',
                        type = str,
                        help = '''Options available are 0: unstranded, 
                                  1: forwardstranded (default)
                                  or 2:reverse_stranded \n''',
                        choices = ['0','1','2'],
                        required = False,
                        default = '1')
    
    parser.add_argument('--th_low',
                        type = float,
                        help = '''Low minimal threshold for feature-to-feature
                                   multi-mapping fraction\n''',
                        required = False,
                        default = 0.01)    
                        
    parser.add_argument('--th_high',
                        type = float,
                        help = '''High minimal threshold for feature-to-feature
                                   multi-mapping fraction\n''',
                        required = False,
                        default = 0.75)    
    
    parser.add_argument('--seed',
                        type = float,
                        help = '''Optional fixed seed for random numbers generation 
                                  during communities detection\n''',
                        required = False,
                        default = None)
    
    ## ---- Long round config.
    parser.add_argument('--feature_output_long',
                        type = str,
                        help = '''GTF field name for which to summarize 
                                  expresion of longRNA assigned reads\n''',
                        required = False,
                        default = 'gene_name')
    
    parser.add_argument('--feature_biotype_long',
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
    parser.add_argument('--feature_small',
                        type = str,
                        help = '''GTF feature type for smallRNA reads 
                                  to be assigned to\n''',
                        required = False,
                        default = 'transcript')
    
    parser.add_argument('--feature_output_small',
                        type = str,
                        help = '''GTF field name for which to summarize 
                                  counts of smallRNA assigned reads\n''',
                        required = False,
                        default = 'transcript_name')
    
    parser.add_argument('--feature_biotype_small',
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
    
    ## Suppress future warnings
    simplefilter(action='ignore', category=FutureWarning)

    ## Set method to start parallel child processes to "spawn"
    mp.set_start_method("spawn")
    
    ##--------------- Parse arguments
    args = parser.parse_args()
    
    ## Required arguments
    outdir = args.outdir
    gtf_filename = os.path.abspath(args.gtf)
    n_cores = args.n_cores

    ## Params
    end = ('Paired' if args.paired_flag else 'Single')
    mtxF = args.mtx
    strand = args.strand_option
    th_high = args.th_high
    th_low = args.th_low
    seed = args.seed
    
    ## Additional optinal arguments
    fcpath = args.featureCounts_path
    btyperounds_filename = args.btyperounds_filename

    ## Get input files
    infiles = [os.path.abspath(line.rstrip()) for line in open(args.bam_infiles)]
    if '' in infiles: infiles.remove('')
    
    ## Select default crounds file if non-user-input defined
    if btyperounds_filename == None:
        ##btype_crounds = config.get_btypes_by_crounds() for major compatibility
        ##btype_crounds = pd.read_csv(utils.resource_path('btypes_crounds.csv')) for compiled version
        btype_crounds = pd.read_csv(os.path.dirname(__file__) + '/data/btypes_crounds.csv')
    else: btype_crounds = pd.read_csv(btyperounds_filename)
    
    ## Counting rounds configuration
    crounds = pd.DataFrame({
        'r':['small',
             'long_exon',
             'long_intron'],
        
        'annot':['small',
                 'long',
                 'long'],
        
        'feat':[args.feature_small,
                'exon',
                'gene'],
        
        'attr':[args.feature_output_small,
                args.feature_output_long,
                args.feature_output_long],
        
        'btype':[args.feature_biotype_small,
                 args.feature_biotype_long,
                 args.feature_biotype_long],
        
        'ml' : [args.ml_flag_small,
                args.ml_flag_long,
                'Derived' if args.ml_flag_long else False],
        
        'minov' : [args.min_overlap_small,
                args.min_overlap_long,
                args.min_overlap_long],
        
        'suf' : ['',
                 '-exon',
                 '-intron']})
    
    ## Get gene biotypes
    gtf = read_gtf(gtf_filename)

    ## Create directory for temporary files
    if not os.path.exists(outdir): os.mkdir(outdir)
    tmpdir = TemporaryDirectory(prefix = os.path.join(outdir, '.mg_'))
    tmppath = os.path.abspath(tmpdir.name)
    
    ## keep all tmp files
    ##tmppath = os.path.join(outdir,'tmp')
    ##if not os.path.exists(tmppath): os.mkdir(tmppath)
    
    try: 
        ## RUN COUNTING PIPELINE
        ## ---------------------------------------------------------------------------
        
        ## ----- Hierarchical read. assignation loop by sample
        crounds = hra.read_assignation_loop(infiles, fcpath, tmppath, gtf_filename,
                                            crounds, btype_crounds, n_cores, end, strand)
        
        ## ----- Multi-loci group extraction
        ml = mg.MG(infiles, outdir, tmppath, gtf, crounds, btype_crounds, n_cores,
                   th_low, th_high, seed) 
        
        ## ----- Expression. level summarization
        counts = count.extract_count_matrix(infiles, tmppath, crounds, n_cores, ml)
        
        ## ---------------------------------------------------------------------------
        
        ##--------------- Save count matrix

        ## Save only nonempty features
        nonempty_idx = np.where(counts.sum(axis=1)!= 0)[0].tolist()

        ## Will label samples by using filename fragment that is distinct in all samples
        if len(args.sample_id) > 0:
            sid = np.array(pd.read_csv(args.sample_id, header = None)[0])
            counts.columns = sid

        ## Save count matrix
        if mtxF:
             mmwrite(
                 os.path.join(outdir,'count_matrix.mtx'),
                 sparse.csr_matrix(counts.iloc[nonempty_idx].values))
        else:
            counts.iloc[nonempty_idx].to_csv(
                os.path.join(outdir,'count_matrix.csv'), header = True, index = True)
            
        ## Save feature metadata table
        feats_metadata = count.extract_feat_metadata(tmppath, gtf, crounds, infiles[0], ml)
        feats_metadata.reindex(counts.iloc[nonempty_idx].index).to_csv(
            os.path.join(outdir,'feature_metadata.csv'))
        
    finally:
        if os.path.exists(tmppath):
            rmtree(tmppath)
        
    print('......')
    print('Voilà.')
    print('Olé!')
    return 0

if __name__ == "__main__":
    main()

