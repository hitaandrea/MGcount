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

def read_assignation_loop(infiles, fcpath, tmppath, gtf_filename,
                          crounds, btype_crounds, n_cores, end, strand):

    print("--------------------------------------------------------")
    print("Assigning alignments to features (1/3)")
    print("--------------------------------------------------------")
    
    ## Split annotation files by cround
    for i in crounds['annot'].unique():

        btypes = btype_crounds['biotype'][btype_crounds['assignation_round']==i]
        tag = crounds['btype'][crounds['annot'] == i].iloc[0]
        pattern = re.compile('|'.join([tag+' "'+bt+ '"' for bt in btypes]))
        gtf_subset = open(os.path.join(tmppath, i + '_annot.gtf'),'w')
        
        if len(btypes) != 0:
            gtf_subset.writelines([line for line in open(gtf_filename) if pattern.search(line)])
        
        if os.stat(os.path.join(tmppath, i + '_annot.gtf')).st_size == 0:
            crounds = crounds.drop(crounds[crounds['annot'] == i].index).reset_index(drop=True)
    
    ## Hiearhihcal assignation parallelized by sample
    n_cores = min(mp.cpu_count(), n_cores)
    func_arh = partial(assign_rna_reads_hierarchically,
                       fcpath = fcpath,
                       tmppath = tmppath,
                       crounds = crounds,
                       end = end,
                       strand = strand)
    
    with mp.Pool(processes = n_cores) as pool:
        result_list = pool.starmap(func_arh, zip(infiles.tolist(), infiles.index.tolist()))

    return crounds



def assign_rna_reads_hierarchically(infile, sn, fcpath, tmppath, crounds, end, strand):

    ## Copy alignment file
    copyfile(infile, os.path.join(tmppath, sn + '_alignment.bam'))
    
    for i in range(0, crounds.shape[0]):
        
        ## Counting round
        r = crounds['r'][i]
        print(r + ' assignation round for: ' + sn)
        
        ## Call feature counts
        call_featurecounts(sn, fcpath, tmppath, end, strand, crounds.iloc[i])
    
        if not os.path.isfile(os.path.join(tmppath, sn + '_alignment.bam.featureCounts')):
            sys.exit('READ ASSIGNATION FAILED. Please check that:\n'
                     '1. The gtf contains the required non-empty fields as defined by the '
                     'assignation arguments "feature", "feature_output" & "feature_biotype'
                     '2. featureCounts is available on the path specified by '
                     '--featureCounts_path argument (default: /user/bin/featureCounts)\n'
                     '3. Input bam files are not empty or corrupted')
            
        ## Filter reads with at least 1 assigned alignment from input .bam file
        extract_assigned_readIds(os.path.join(tmppath, sn + '_alignment.bam.featureCounts'),
                                 os.path.join(tmppath, sn + '_fc_' + r),
                                 os.path.join(tmppath, sn + '_fids.txt'))
        filter_bamreads(os.path.join(tmppath, sn + '_alignment.bam'),
                        os.path.join(tmppath, sn + '_alignment2.bam'),
                        os.path.join(tmppath, sn + '_fids.txt')) ##'_fids2.txt')
        copyfile(os.path.join(tmppath, sn + '_alignment2.bam'),
                 os.path.join(tmppath, sn + '_alignment.bam'))
        
        ## Remove temporary files
        os.remove(os.path.join(tmppath, sn + '_alignment2.bam'))
        os.remove(os.path.join(tmppath, sn + '_alignment.bam.featureCounts'))
    
    ## Remove temporary files
    os.remove(os.path.join(tmppath, sn + '_fids.txt'))
    
    return 0


def call_featurecounts(sn, fcpath, tmppath, end, strand, cround):

    ## Change temporarily the path 
    cwd = os.getcwd()
    os.chdir(tmppath)

    ## Run program
    os.system(fcpath +

              ## ---- alignment input file
              ' -o ' + '"' +
              os.path.join(tmppath, sn + '_counts_' + cround['r']  + '.csv') +  '" ' +

              ## ---- annotations
              ' -a ' + '"' +
              os.path.join(tmppath, cround['annot'] + '_annot.gtf') +  '" ' +
               
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
              '"' + os.path.join(tmppath, sn + '_alignment.bam') + '" ' +
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

    
def extract_assigned_readIds(fc_infile, fc_outfile, reads_list_outfile):
    
    ## Open file connections
    infile = open(fc_infile, 'r')
    outfile = open(fc_outfile, 'w')
    rids_outfile = open(reads_list_outfile, 'w')

    ## Initialize read_ids seen set
    rids_seen = set()

    ## Iterate over fc file            
    for fc_line in infile.readlines():
        aline = re.findall(r'Assigned', fc_line)
        if aline:
            outfile.write(fc_line)
            rid = fc_line.split('\t')[0]
            if rid not in rids_seen:
                rids_outfile.write(rid+'\n')
                rids_seen.add(rid)

    ## Exit
    infile.close()
    outfile.close()
    rids_outfile.close()
           
