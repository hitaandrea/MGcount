## APR-2021  Andrea Hita Ardiaca

## -- Modules import
import os
import re
import pandas as pd
import numpy as np
import multiprocessing as mp
import itertools as it
import random as rnd
import scipy.sparse
from igraph import Graph

from scipy.io import mmwrite, mmread
from copy import deepcopy
from functools import partial

## MGcount source code
from mgcount import utils

def MG(infiles, outdir, tmppath, gtf, crounds, btype_crounds, n_cores, 
       th_low, th_high, seed): 

    print("--------------------------------------------------------")
    print("Extracting multi-loci groups (2/3)")
    print("--------------------------------------------------------")
    
    out_ml = dict()
    samples = [re.sub('.bam','', utils.path_leaf(infile)) for infile in infiles]

    ## Setup multicore processing
    n_cores = min(mp.cpu_count(), n_cores)
    for i in np.where(crounds['ml']==True)[0]:

        ## ---------------- Get non-empty features for mg generation
        cround = crounds.iloc[i,]
        feat_btypes = pd.merge(
            gtf[[cround['attr'], cround['btype']]][gtf['feature']==cround['feat']],
            btype_crounds, left_on = cround['btype'], right_on = 'biotype')
            
        if cround['btype'] != 'biotype':
            feat_btypes.drop(['biotype'], axis='columns', inplace=True)
        
        ## Get sample file_names
        fc_infiles = [tmppath + sn + '_fc_' + cround['r'] for sn in samples]

        ## Get non-empty features (reduces mg matrix size and speeds-up computation)
        counts = np.zeros(pd.read_table(tmppath+samples[0]+'_counts_'+ cround['r']+'.csv', 
                           sep = '\t', skiprows = 1, usecols = [6]).shape[0])
        for sn in samples:
            counts = counts + pd.read_table(tmppath+sn+'_counts_'+ cround['r']+'.csv', 
                                             sep = '\t', skiprows = 1, usecols = [6]).iloc[:,0].to_numpy() 
        
        feats = pd.read_table(tmppath+samples[0]+'_counts_'+cround['r']+'.csv',
                               sep = '\t', skiprows = 1, usecols = [0]).iloc[
                                   np.where(counts!= 0)[0].tolist(),0].reset_index(
                                       drop = True)    
                                       
        ## Get gene biotypes
        feats = pd.DataFrame({cround['attr']:feats}).merge(
            feat_btypes[feat_btypes['assignation_round'] ==
                        cround['annot']].drop_duplicates().reset_index(drop = True),
            on=cround['attr'], how='left')
        feats['dupli_entries'] = feats[cround['attr']].duplicated(keep=False)
        feats = feats.drop_duplicates().reset_index(drop = True)
        feats.loc[feats[cround['attr']].duplicated(keep=False), cround['btype']] = 'Hybrid'
        feats = feats.drop_duplicates().reset_index(drop = True)
        feats.to_csv(tmppath + 'feats_'+ cround['r'] + '.csv')
        
        ## ---------------- Extract adjacency matrix
        print("Building " + cround['r'] + " multi-mappers graph" )

        ## Parallelize graph matrix generation by sample, asynchronously
        func_eam = partial(extract_adjacency_matrix,
                           feat_list = feats[cround['attr']].tolist())

        with mp.Pool(processes = n_cores) as pool:
            mgm_list = pool.map(func_eam, fc_infiles)

        ## Add matrices
        mgm = mgm_list[0]
        mmwrite(tmppath+samples[0]+'mgm_'+cround['r']+ '.mtx', mgm_list[0])
        for k in range(1,len(mgm_list)):
            mmwrite(tmppath+samples[k]+'mgm_'+cround['r']+ '.mtx', mgm_list[k])
            mgm = mgm + mgm_list[k]

        ## Store output matrix
        mmwrite(outdir + 'multigraph_matrix_'+ cround['r'] +'.mtx', mgm)
        ##mgm = mmread(outdir + 'multigraph_matrix_'+ cround['r'] +'.mtx')
        
        ## ---------------- Build graph and extract multiloci groups
        print("Detecting " + cround['r'] + " communities from multi-mappers graph")
        ml = pd.DataFrame()
        
        ## LOOP by biotype if sRNA
        btypeloop = feats[cround['btype']].unique() if cround['r'] == 'small' else [''] 

        for btype in btypeloop:

            ## Get biotype adjacency matrix subset
            if btype != '':
                fids = np.where(feats[cround['btype']] == btype)[0]
                inM = mgm.tocsr()[fids,:][:,fids]; inFeats = feats.iloc[fids.tolist()]
            else:
                inM = mgm.tocsr(); inFeats = feats

            if inFeats.shape[0] == 0:
                continue
                                    
            ## Get multi-loci groups
            df = get_graph_loci_groups(inM, inFeats, cround, th_low, th_high, seed)

            ## Append output
            btype1 = cround['annot'] if btype == '' else btype
            df['community_id'] = [btype1+'-cl-'+str(cid) if cid is not None else None
                                 for cid in df['community_id']]
            df['feature'] = [df['community_name'].iloc[i] if df['community_flag'].iloc[i]
                               else df[cround['attr']].iloc[i]
                               for i in range(df.shape[0])]
            ml = ml.append(df)
            print('   --> ' + btype1 + ' multi-loci groups successfully extracted.')

        ## Store multi-loci groups list
        ml = ml.set_index(cround['attr'])
        ml = ml.reindex(index = feats[cround['attr']].tolist())
        ml = ml.reset_index()
        ml.to_csv(outdir+'multigraph_communities_'+cround['r']+'.csv', index = False)
        out_ml[cround['annot']] = ml
        
    return out_ml


def extract_adjacency_matrix(fc_infile, feat_list):
                         
    ## -------- Import feature counts output data
    N = len(feat_list)
    M = np.zeros(shape = (len(feat_list),len(feat_list)), dtype = np.uint32)
    vdict = dict(zip(feat_list,range(0,len(feat_list))))

    if os.stat(fc_infile).st_size == 0:
        M_sparse = scipy.sparse.csc_matrix(M)
        return M_sparse

    os.system('sort -o ' + fc_infile + ' ' + fc_infile)
    fc_in = open(fc_infile, 'r')
    currentRead = fc_in.readline().rstrip().split('\t')
    connected_vertices = [int(vdict[vname]) for vname in currentRead[3].split(',')]
    i = 0; nreads = 0
    maxIter = sum(1 for line in open(fc_infile))-1

    ## -------- LOOP: Adjacency matrix generation
    while i <= maxIter:
        newRead = fc_in.readline().rstrip().split('\t')

        if newRead[0] != currentRead[0]: 
            edges = list(it.product(connected_vertices, repeat = 2))
            connected_vertices = []
            currentRead = newRead
            for idx in edges:
                M[idx] = M[idx] + 1
            nreads+=1

        if i == maxIter:
            break
        
        for vname in newRead[3].split(','):
            connected_vertices.append(int(vdict[vname]))
        i+=1

    print('   --> Multi-mappers adjacency matrix successfully extracted from ' +
          str(nreads) + ' reads and ' +
          str(i-1) + ' alignments!')
    M_sparse = scipy.sparse.csc_matrix(M)
    M_sparse.dtype = np.int32
    return M_sparse


def get_graph_loci_groups(inM, inFeats, cround, th_low, th_high, seed):
    
    ## ---- Scale adj. matrix
    D = deepcopy(np.diag(inM.todense()))
    D[np.where(D == 0)] = 1
    sM = inM.multiply(scipy.sparse.csr_matrix(1/D)).todense()
    np.fill_diagonal(sM,0)

    ## If no multimappers, no graph is built, just exit
    if np.all(sM == 0):
       return pd.DataFrame({cround['attr']:inFeats[cround['attr']],
                            cround['btype']:inFeats[cround['btype']],
                            'community_flag':False,'community_id':None,
                            'community_name':None,'community_biotype':None})
    
    ## ---- Generate multi-mappers graph
    sources, targets = sM.nonzero()
    g = Graph(n = len(D), edges = zip(sources, targets), directed=True,
              edge_attrs={'weight': np.array(sM[sources, targets])[0]},
              vertex_attrs={'weight': np.log10(D), 'naln': D})
    g.es[np.where(np.array(g.es['weight']) > 1)[0].tolist()]['weight'] = 1.0

    ## ---- Multi-loci groups automatic recognition
    
    ## Filter weak edges
    gcf = deepcopy(g)

    ## Infomap algorithm community detection
    if cround['annot'] == 'long':
        gcf.delete_edges(np.where(np.array(g.es['weight']) < th_high)[0].tolist())
    elif np.median(g.degree()) < 50:
        th_arr = np.array([th_low, th_high])
        gcf.delete_edges(np.where(np.array(g.es['weight']) < th_arr[(np.abs(
            th_arr - np.median(g.es['weight']))).argmin()])[0].tolist())

    ## Set seed before runing infomap algorithm
    if seed is not None:
        rnd.seed(seed)
    
    ## Run infomap (Rosvall et al. 2008)
    g_communities = gcf.community_infomap(trials = 100)
    g.vs['cl'] = g_communities.membership
    
    ## Initialize ml. output structure
    out = pd.DataFrame({cround['attr'] :inFeats[cround['attr']],
                        cround['btype']:inFeats[cround['btype']],
                        'mb':g.vs['cl'],
                        'naln':g.vs['naln']})
    out = out.merge(out.rename(columns = {'naln':'naln_community'}
        ).groupby(['mb']).naln_community.sum().reset_index(), on = 'mb')

    out['community_flag'] =  (out['mb'].isin(out['mb'].loc[out.duplicated('mb')])) & (out['naln_community'] > 50)
    out['community_id']=None; out['community_name']=None
    out['community_name_dedup']=None; out['community_biotype'] = None

    ## Extract multi-loci communities list
    cl_names = out[cround['attr']].loc[out['community_flag'] == False].tolist()
    dups = []
    
    ## ---- Multi-loci groups ouput generation
    for cl in out['mb'].loc[out['community_flag']].unique():

        ## Get gene ids in multi-loci gropu
        cl_id = cl
        ##cl_members = np.where(np.array(g_communities.membership) == cl)[0]    
        cl_members = np.where(np.array(g.vs['cl']) == cl)[0]    
        cl_feats = inFeats.iloc[cl_members]
        cl_feats = cl_feats.assign(nalign = D[cl_members].tolist())

        ## Define multi-loci community name
        tmp = define_multiloci_featname_paste(cl_feats, cround)
        cl_name = tmp['cn']; cl_btype = tmp['btype']
        
        ## Track possible name duplicates
        cl_names.append(cl_name)
        s = sum([i == cl_name for i in cl_names])
        cl_name_suf = cl_name + '_' + str(s)
        if s > 1:
            dups.append(cl_name)
            
        ## Add new multi-loci group to output
        out.loc[out.mb == cl, 'community_id'] = cl_id
        out.loc[out.mb == cl, 'community_name'] = cl_name
        out.loc[out.mb == cl, 'community_name_dedup'] = cl_name_suf
        out.loc[out.mb == cl, 'community_biotype'] = cl_btype

    ## Deduplicate duplicated multi-loci group names
    out.loc[out['community_name'].isin(dups), 'community_name'] = out.loc[
        out['community_name'].isin(dups), 'community_name_dedup']
    del(out['community_name_dedup'])
    del(out['mb'])
    
    ## Return multi-loci groups data.frame
    return out


def define_multiloci_featname_topcommon(cl_feats, cround):

    topF = cl_feats.iloc[np.argmax(np.array(cl_feats['nalign']))]
    cn = [utils.lcs(topF[cround['attr']], g) for g in cl_feats[cround['attr']]]
    cn = [i.pop() for i in cn if len(i) > 0 ]

    if (len(cn) > 1) & (any([len(i)>3 for i in cn])):
        cn = [i for i in cn if len(re.sub('-','',i)) > 3]

    if len(cn) > 1:
        cn = utils.most_frequent(cn)
    else:
        cn = cn[0]
        
    return {'cn':cn.rstrip('-'), 'btype':topF[cround['btype']]}


def define_multiloci_featname_paste(cl_feats, cround):

    ## Identify if all loci are pseudogenes
    anyNonP = any("pseudogene" not in s for s in cl_feats[cround['btype']].tolist())

    ## If not, exclude pseudogenes from multi-loci feature name
    if anyNonP:
        cl_feats = cl_feats.loc[["pseudogene" not in s
                                 for s in cl_feats[cround['btype']].tolist()]]

    ## Get top. aligned loci within community as community biotype
    topF = cl_feats.iloc[np.argmax(np.array(cl_feats['nalign']))]
    bt = topF[cround['btype']]

    ## Collapse feat.names from individual loci to generate group name
    tocoll = cl_feats.loc[cl_feats[cround['btype']] == bt]
    tocoll = tocoll.sort_values('nalign', ascending = False)[cround['attr']].tolist()
    
    if len(tocoll) > 8:
        cn = '_'.join(['_'.join(tocoll[0:8]),'etc'])
    else:
        cn = '_'.join(tocoll)
    
    return {'cn':cn, 'btype': bt}
