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
        result_list = pool.map(
           func_arh, [x for x in zip(infiles.tolist(), infiles.index.tolist())])
        
    return crounds



def assign_rna_reads_hierarchically(infile, fcpath, tmppath, crounds, end, strand):

    bf, sn = infile
    bam_infile = os.path.join(tmppath, sn + '.bam')
    rids_file = os.path.join(tmppath, sn + '_fids.txt')
    ## print('--> Hierarchical assignations for ' + sn)

    ## Copy alignment file
    copyfile(bf, bam_infile)

    ## Assignation LOOP:
    for i in range(0, crounds.shape[0]):
        
        ## Counting round
        r = crounds['r'][i]
        bam_outfile = os.path.join(tmppath, sn + '_' + r + '.bam')

        fc_infile = os.path.join(tmppath, bam_infile.rsplit('/')[-1] + '.featureCounts')
        fc_outfile = os.path.join(tmppath, sn + '_fc_' + r)
        print('--> ' + r + ' assignation round for: ' + sn)
        
        ## Call feature counts
        call_featurecounts(bam_infile, sn, fcpath, tmppath, end, strand, crounds.iloc[i])
        if not os.path.isfile(fc_infile):
            sys.exit('READ ASSIGNATION FAILED for' + sn + '. Please check that:\n'
                     '1. The gtf contains the required non-empty fields as defined by the '
                     'assignation arguments "feature", "feature_output" & "feature_biotype'
                     '2. featureCounts is available on the path specified by '
                     '--featureCounts_path argument (default: /user/bin/featureCounts)\n'
                     '3. Input bam files are not empty or corrupted')
        else:
            ## Filter reads with at least 1 assigned alignment from input .bam file
            extract_assigned_readIds(fc_infile, fc_outfile, rids_file)
            filter_bamreads(bam_infile, bam_outfile, rids_file)

            ## Remove raw fc, update bam, and proceed with next iteration
            os.remove(fc_infile);
            bam_infile = bam_outfile
    
    ## Remove temporary files
    os.remove(os.path.join(tmppath, sn + '_fids.txt'))
    
    return 0


def call_featurecounts(bam_infile, sn, fcpath, tmppath, end, strand, cround):

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
              '"' + bam_infile + '" ' +
              '>> /dev/null 2>&1')
        
    ## Re-change path to original
    os.chdir(cwd)
    

def filter_bamreads(bam_infile, bam_outfile, rids_infile):

    # Get selected read nemaes
    fq = open(rids_infile).readlines()
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

    
def extract_assigned_readIds(fc_infile, fc_outfile, rids_outfile):
    
    ## Open file connections
    infile = open(fc_infile, 'r')
    outfile = open(fc_outfile, 'w')
    rids_outfile = open(rids_outfile, 'w')
    
    ## Iterate over fc file            
    for fc_line in infile.readlines():
        aline = re.findall(r'Assigned', fc_line)
        if aline:
            outfile.write(fc_line)
            rids_outfile.write(fc_line.split('\t')[0]+'\n')           
    ## Exit
    infile.close()
    outfile.close()
           
