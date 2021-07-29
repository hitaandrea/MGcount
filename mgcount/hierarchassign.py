## APR-2021  Andrea Hita Ardiaca

## -- Modules import
import os
import re
import sys
import pandas as pd
import numpy as np
import multiprocessing as mp
import pysam 

from shutil import copyfile
from functools import partial

## MGcount source code
from mgcount import utils

def read_assignation_loop(infiles, outdir, fcpath, tmppath, gtf_filename,
                          crounds, btype_crounds, n_cores, end, strand):

    print("--------------------------------------------------------")
    print("Assigning alignments to features (1/3)")
    print("--------------------------------------------------------")
    
    ## Split annotation files by cround
    for i in crounds['annot'].unique():

        btypes = btype_crounds['biotype'][btype_crounds['assignation_round']==i]
        tag = crounds['btype'][crounds['annot'] == i].iloc[0]
        pattern = re.compile('|'.join([tag+' "'+bt+ '"' for bt in btypes]))
        gtf_subset = open(tmppath + i + '_annot.gtf','w')
        
        if len(btypes) != 0:
            gtf_subset.writelines([line for line in open(gtf_filename) if pattern.search(line)])
        
        if os.stat(tmppath + i + '_annot.gtf').st_size == 0:
            crounds = crounds.drop(crounds[crounds['annot'] == i].index).reset_index(drop=True)
            
    ## Hiearhihcal assignation parallelized by sample
    n_cores = min(mp.cpu_count(), n_cores)
    func_arh = partial(assign_rna_reads_hierarchically,
                       outdir = outdir,
                       fcpath = fcpath,
                       tmppath = tmppath,
                       crounds = crounds,
                       end = end,
                       strand = strand)
                       
    with mp.Pool(processes = n_cores) as pool:
        result_list = pool.map(func_arh, infiles)

    return crounds



def assign_rna_reads_hierarchically(infile, outdir, fcpath, tmppath, crounds, end, strand):

    ##sn = re.sub('/','_',re.sub('.bam','',infile))
    sn = re.sub('.bam','', utils.path_leaf(infile))

    ## Copy alignment file
    copyfile(infile, tmppath + sn + '_alignment.bam')
    
    for i in range(0, crounds.shape[0]):
        
        ## Counting round
        r = crounds['r'][i]
        print(r + ' feat. assignation round for: ' + sn)
        
        ## Call feature counts
        call_featurecounts(outdir, sn, fcpath, tmppath, end, strand, crounds.iloc[i])
        
        if not os.path.isfile(tmppath + sn + '_alignment.bam.featureCounts'):
            sys.exit('''Assignation failed. Please check your gtf or/and the 
                        arguments configuration provided for "feat", "attr" & "biotype''')

        ## Filter reads with at least 1 assigned alignment from input .bam file
        os.system('grep "Assigned" ' + tmppath + sn +
                  '_alignment.bam.featureCounts ' +
                  '| cut -f 1 >> ' + tmppath + sn + '_fids.txt')
        utils.remove_duplicated_lines(tmppath + sn + '_fids.txt',
                                      tmppath + sn + '_fids2.txt')     
        filter_bamreads(tmppath + sn + '_alignment.bam',
                        tmppath + sn + '_alignment2.bam',
                        tmppath + sn + '_fids2.txt')
        copyfile(tmppath + sn + '_alignment2.bam',
                 tmppath + sn + '_alignment.bam')
        os.remove(tmppath + sn + '_alignment2.bam')
        
        ## Export alignment-feature assignation records
        os.system('grep "Assigned" ' +
                  tmppath + sn + '_alignment.bam.featureCounts >> ' +
                  tmppath + sn + '_fc_' + r)
        ##copyfile(tmppath + sn + '_alignment.bam.featureCounts',
        ##         tmppath + sn + '_fc2_' + r)
        os.remove(tmppath + sn + '_alignment.bam.featureCounts')
    
    ## Remove temporary files
    os.remove(tmppath + sn + '_fids.txt')
    os.remove(tmppath + sn + '_fids2.txt')
    
    return 0


def call_featurecounts(outdir, sn, fcpath, tmppath, end, strand, cround):

    ## Change temporarily the path 
    cwd = os.getcwd()
    os.chdir(tmppath)

    ## Run program
    os.system(fcpath +

              ## ---- alignment input file
              ' -o ' + tmppath +
              sn + '_counts_' + cround['r']  + '.csv' +

              ## ---- annotations
              ' -a ' + tmppath + cround['annot'] + '_annot.gtf' +
               
              ' -F ' + 'GTF' +
              ' -t '+ cround['feat'] +
              ' -g ' + cround['attr'] +
              ##'-A '
               
              ## ---- level of summarization
              ## -f if not specify, count at atttribute level
               
              ## ---- overlap between reads and features
              ' -O' +
              ##' --minOverlap 15' +
              ' --fracOverlap ' + cround['minov'] +
              ## '--fracOverlapFeature +
              ' --largestOverlap' +
              ## '--nonOverlap '
              ## '--nonOverlapFeature ' 
              ## '--readExtension5 '
              ## '--readExtension3 '
              ## '--read2pos '
               
              ## ---- multi-mapping reads
              ' -M' +
              ' --fraction' +
               
              ## ---- read filtering
              ## '-Q '
              ## '--splitOnly '
              ## '--nonSplitOnly '
              ## '--primary  '
              ## '--ignoreDup '
               
              ## ---- strandness
              ' -s ' + strand +
               
              ## ---- exon-exon junctions
              ## '-J '
              ## '-G '
               
              ## ---- parameters specific to paired end reads
              (' -p ' if end == 'Paired' else '') +
              ## '-B '
              ## '-P '
              ## '-d '
              ## '-D '
              ## '-C '
              ## '--donotsort '
               
              ## number of CPU threads
              ##'-T '
              ## read groups
              ## '-byReadGroup'
               
              ## assignment results for each read
              ' -R CORE ' +
              ## '--RPath '
               
              ## miscellaneous
              ## '--tmpDir '
              ## '--maxMOp '
              ## '--verbose '
              ## '-v '
              tmppath + sn + '_alignment.bam ' +
              '>> /dev/null 2>&1')
        
    ## Re-change path to original
    os.chdir(cwd)
    

def filter_bamreads(bam_infile, bam_outfile, reads_list_infile):

    # Get selected read nemaes
    fq = open(reads_list_infile).readlines()
    fq = [x.strip() for x in fq]
    fq = set(fq)

    # Suppress missing index warning when reading file
    save = pysam.set_verbosity(0) 
    infile = pysam.AlignmentFile(bam_infile)
    pysam.set_verbosity(save)

    outfile = pysam.AlignmentFile(bam_outfile, template= infile, mode= 'wb')

    # Write selected reads in output alignment file
    for aln in infile.fetch(until_eof = True):
        if aln.query_name not in fq:
            outfile.write(aln)

    ## Exit
    infile.close()
    outfile.close()

    
