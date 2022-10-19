## APR-2021  Andrea Hita Ardiaca

## -- Modules import
import pandas as pd
import multiprocessing as mp
import numpy as np
import os as os
import re as re
from scipy import sparse
from functools import partial

def extract_count_matrix(infiles, tmppath, crounds, n_cores, ml=None, mtxF=False):

    print("--------------------------------------------------------")
    print("Computing expression matrix (3/3)")
    print("--------------------------------------------------------")

    ## Get features
    features = pd.Series(get_counts(sn=infiles.index[0], tmppath=tmppath, crounds=crounds, ml=ml, mtxF=False).index)
    
    ## Set number of cores
    n_cores = min(mp.cpu_count(), n_cores)

    ## Build count_matrix by sample
    if mtxF is False:
        with mp.Pool(processes = n_cores) as pool:
            func_gc = partial(get_counts,
                              tmppath = tmppath,
                              crounds = crounds,
                              ml = ml,
                              mtxF = mtxF)
            counts = pd.concat(pool.map(func_gc, infiles.index), axis = 1)
    else:
        with mp.Pool(processes = n_cores) as pool:
            func_gc = partial(get_counts,
                              tmppath = tmppath,
                              crounds = crounds,
                              ml = ml,
                              mtxF = mtxF)
            counts = sparse.csr_matrix(sparse.hstack(pool.map(func_gc, infiles.index)))

    ## Get only non-empty features
    nonempty_idx = np.where(counts.sum(axis=1)!= 0)[0].tolist()
    features = features.iloc[nonempty_idx]
    features.index = features
    
    if mtxF is False:
        counts = counts.iloc[nonempty_idx]
    else:
        counts = counts[nonempty_idx,:]
    
    ## Extract feat. metadata table
    return [counts, features]
        

def get_counts(sn, tmppath, crounds, ml, mtxF):

    print(sn)

    ## By counting round, retrieve output and combine
    out = pd.DataFrame()
    for i in range(0, crounds.shape[0]):

        cround = crounds.iloc[i]

        ## Import fractionated counts
        df = pd.read_table(
            os.path.join(tmppath, sn + '_counts_' + cround['r']  + '.csv'),
            sep = '\t', skiprows = 1, usecols = [0,6], index_col = 0)
        df.columns = ['counts']

        ## If small non-coding RNA, group by multi-loci recognized groups
        if cround['ml'] != 0:
            df = df.merge(ml[cround['annot']][['feature',cround['attr']]],
              left_index = True, right_on = cround['attr'])
        else:
            df['feature'] = df.index 
            df = df[['feature','counts']]

        df = df.groupby(['feature'], as_index=False).agg({'counts':sum})

        ## Append counting round
        df['feature'] = df['feature'] + cround['suf']
        out = out.append(df)
            
    out = out.set_index('feature')

    if mtxF is False:
        return out
    else:
        return sparse.csr_matrix(out)


def extract_feat_metadata(tmppath, gtf, crounds, onesn, ml):
    
    ## Initialize data.frame
    colnames = ['feature','assignation_round','annotations_subset',
                'feature_type','feature_output','feature_biotype','community_flag']
    feats_metadata = pd.DataFrame(columns = colnames)

    for i in range(0, crounds.shape[0]):
        cround = crounds.iloc[i]
        counts = pd.read_table(
            os.path.join(tmppath, onesn + '_counts_' + cround['r']  + '.csv'),
            sep = '\t', skiprows = 1, usecols=[0], index_col=0)
        
        ## ---- Fill feats-metadata frame
        if cround['ml']:
            countsml = counts.merge(
                ml[cround['annot']][['feature',cround['attr']]],
                left_index = True,
                right_on = cround['attr'])
            countsml = countsml.drop_duplicates('feature').reset_index(drop=True)

            df = pd.DataFrame({'feature':countsml['feature']})
            df['feature_output'] = cround['attr']
            df['assignation_round'] = cround['r']
            df['annotations_subset'] = cround['annot']
            df['feature_type'] = cround['feat']
            
            df = df.merge(ml[cround['annot']][[
                'feature','community_biotype','community_flag',
                cround['attr'],cround['btype']]], 
                          on = 'feature')

            df['feature_biotype'] = [df['community_biotype'].iloc[i] if df['community_flag'].iloc[i]
                         else df[cround['btype']].iloc[i] for i in range(df.shape[0])]

            df = df.drop(columns = [cround['btype'], 'community_biotype'])
            df = df.drop_duplicates('feature').reset_index(drop = True)
            df = df[colnames]    
        else:

            df = pd.DataFrame({'feature':counts.index.tolist()})
            df['feature_output'] = cround['attr']
            df['assignation_round'] = cround['r']
            df['annotations_subset'] = cround['annot']
            df['feature_type'] = cround['feat']

            feat_btypes = df.merge(
                gtf[[cround['attr'],cround['btype']]][
                    gtf['feature']==cround['feat']].drop_duplicates(#cround['attr']
                    ).reset_index(drop=True),
                left_on = 'feature', right_on = cround['attr'], how = 'left')

            feat_btypes['dupli_entries'] = feat_btypes[cround['attr']].duplicated(keep=False)
            feat_btypes = feat_btypes.drop_duplicates().reset_index(drop = True)

            feat_btypes.loc[feat_btypes['feature'].duplicated(keep=False), cround['btype']] = 'Hybrid'

            feat_btypes = feat_btypes.drop_duplicates().reset_index(drop = True)
            feat_btypes = feat_btypes.set_index('feature')
            feat_btypes = feat_btypes.reindex(df['feature'])

            df['feature_biotype'] = feat_btypes[cround['btype']].tolist()
            
        df['feature'] = df['feature']+cround['suf']
        feats_metadata = feats_metadata.append(df)
                            
    ## Save output
    feats_metadata = feats_metadata.set_index('feature')
    return feats_metadata
 

