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
from mgcount import utils 

def main():

    parser = argparse.ArgumentParser(
        description = """MGcount RNA-seq quantification pipeline v1.0.2 \n""")
        
    ##--------------- Required I/O
    parser_required = parser.add_argument_group('required arguments')
    
    parser_required.add_argument('--bam_infiles',
                                 '-i',
                                 type = os.path.abspath,
                                 help = 'Alignment (.bam) input file names\n',
                                 required = True)
 
    parser_required.add_argument('--outdir',
                                 '-o',
                                 type = os.path.abspath,
                                 help ='Output directory path \n',
                                 required = True)
    
    parser_required.add_argument('--gtf',
                                 type = os.path.abspath,
                                 help = 'Annotations (.gtf) file name \n',
                                 required = True)
    
    ##--------------- Optional arguments
    parser.add_argument('--sample_names',
                        type= os.path.abspath,
                        help = '''Optional sample names to be used as matrix column ids
                                (alternatively, "bam_infiles" file names will be used)"\n''',
                        required = False, 
                        default = None)
    
    parser.add_argument('-T',
                        '--n_cores',
                        type = int,
                        help = 'Number of cores for parallelization (default 1) \n',
                        required = False,
                        default = 1)
                            
    parser.add_argument('-p',
                        '--paired_flag',
                        action = 'store_true',
                        help = 'Paired end flag (default 1)\n')

    parser.add_argument('-m',
                        '--mtx',
                        action = 'store_true',
                        help = 'Output count matrix in sparse format (default False)\n')
    
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
                                   multi-mapping fraction (default 0.01)\n''',
                        required = False,
                        default = 0.01)    
                        
    parser.add_argument('--th_high',
                        type = float,
                        help = '''High minimal threshold for feature-to-feature
                                   multi-mapping fraction (default 0.75)\n''',
                        required = False,
                        default = 0.75)    
    
    parser.add_argument('--seed',
                        type = float,
                        help = '''Optional fixed seed for random numbers generation 
                                  during communities detection (default random integer)\n''',
                        required = False,
                        default = None)
    
    ## ---- Long round config.
    parser.add_argument('--feature_output_long',
                        type = str,
                        help = '''GTF field name for which to summarize 
                                  expresion of longRNA assigned reads (default: gene_name)\n''',
                        required = False,
                        default = 'gene_name')
    
    parser.add_argument('--feature_biotype_long',
                        type = str,
                        help = '''GTF field name defining biotype for
                                  longRNA features (default: gene_biotype)\n''',
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
                                  groups flag for long round.
                                  Options are 1: enable (default) and 0: disable\n''',
                        required = False,
                        default = 1)
    
    ## ---- Small round config.
    parser.add_argument('--feature_small',
                        type = str,
                        help = '''GTF feature type for smallRNA reads 
                                  to be assigned to (default: transcript)\n''',
                        required = False,
                        default = 'transcript')
    
    parser.add_argument('--feature_output_small',
                        type = str,
                        help = '''GTF field name for which to summarize 
                                  counts of smallRNA assigned reads
                                  (default: transcript_name)\n''',
                        required = False,
                        default = 'transcript_name')
    
    parser.add_argument('--feature_biotype_small',
                        type = str,
                        help = '''GTF field name defining biotype for
                                  smallRNA features
                                  (default: transcript_biotype)\n''',
                        required = False,
                        default = 'transcript_biotype')
    
    parser.add_argument('--min_overlap_small',
                        type = str,
                        help = '''Minimal feature-alignment overlapping fraction
                                  for assignation (default 1) for small round \n''',
                        required = False,
                        default = '1')

    parser.add_argument('--ml_flag_small',
                        type = int,
                        help = '''Multi-loci graph detection 
                                  based groups flag for small round.
                                  Options are 1: enable (default) and 0: disable\n''',
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
                                  assignation_round custom defined pairs
                                  (alternatively, in-built btyperounds will be used)\n''',
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

    ## Params
    end = ('Paired' if args.paired_flag else 'Single')
    mtxF = args.mtx
    strand = args.strand_option
    th_high = args.th_high
    th_low = args.th_low
    seed = args.seed
    n_cores = args.n_cores

    ## Additional optinal arguments
    fcpath = args.featureCounts_path
    btyperounds_filename = args.btyperounds_filename

    
    ## --------------- Read inputs
    
    ## Get input files
    infiles = pd.Series([os.path.abspath(line.rstrip()) for line in open(args.bam_infiles)])
    if '' in infiles: infiles.remove('')

    ## Get annotations file
    gtf = read_gtf(gtf_filename)
    
    ## If not parse as argument, define sample_names as filenames, with full-paths if duplicates 
    infiles.index = [re.sub('/','_',re.sub('.bam','',infile)) for infile in infiles]

    tvar = [re.sub('.bam','', utils.path_leaf(infile)) for infile in infiles]
    if len(tvar) == len(set(tvar)):
        infiles.index = tvar
    del tvar
    
    if args.sample_names is not None:
        infiles.index = [os.path.abspath(line.rstrip()) for line in open(args.sample_names)]

    ## Select default crounds file if non parsed as argument
    ##btype_crounds = config.get_btypes_by_crounds() for major compatibility
    ##btype_crounds = pd.read_csv(utils.resource_path('btypes_crounds.csv')) for compiled version
    if btyperounds_filename is None:
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
    
    ## Create directory for temporary files
    if not os.path.exists(outdir): os.mkdir(outdir)
    tmpdir = TemporaryDirectory(prefix = os.path.join(outdir, '.mg_'))
    tmppath = os.path.abspath(tmpdir.name)
    
    ## keep all tmp files
    ## tmppath = os.path.join(outdir,'tmp')
    ## if not os.path.exists(tmppath): os.mkdir(tmppath)
    
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
        counts, features = count.extract_count_matrix(infiles, tmppath, crounds,
                                                      n_cores, ml, mtxF)
        ## ---------------------------------------------------------------------------

        
        ## --------------- Write outputs
                
        ## Save count matrix
        if mtxF:
             mmwrite(os.path.join(outdir,'count_matrix.mtx'), counts)
        else:
            counts.columns = infiles.index
            counts.to_csv(os.path.join(outdir,'count_matrix.csv'), header=True, index=True)

        ## Save sample metadata table
        infiles.to_csv(os.path.join(outdir,'sample_metadata.csv'), index=True)
            
        ## Save feature metadata table
        feats_metadata = count.extract_feat_metadata(tmppath, gtf, crounds, infiles.index[0], ml)
        feats_metadata.reindex(features).to_csv(os.path.join(outdir,'feature_metadata.csv'))
        
    finally:
        if os.path.exists(tmppath):
            rmtree(tmppath)
        
    print('......')
    print('Voilà.')
    print('Olé!')
    return 0

if __name__ == "__main__":
    main()

