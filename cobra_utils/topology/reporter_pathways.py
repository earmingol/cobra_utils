# -*- coding: utf-8 -*-

from __future__ import absolute_import

import numpy as np
import pandas as pd
import scipy.stats as stats

from sklearn.utils import resample
from cobra_utils import query
from _collections import defaultdict


def reporter_pathways(model, p_val_df, pathways=None, rxn_pathways_association=None, verbose=True):
    '''
    This function computes an aggregate p-value for each pathway (SubSystem in the metabolic reconstruction) based on the
    network topology of the metabolic reconstruction. It takes the p-value for differential expression of each gene and
    compute the aggregate p-value for the neighbor reactions of a given pathway.

    Parameters
    ----------
    model : cobra.core.Model.Model
        A cobra model.

    p_val_df : pandas.DataFrame
        A dataframe with gene names as index. It have to contains the p-values for the differential expression
        of the respective indexing genes.

    pathways : array-like
        An array or list containing pathway names (str) to be considered.

    rxn_pathways_association : dict
        A dictionary where the keys are the pathways and the values a list of reactions (RxnIDs) that belong to those
        pathways.

    verbose : boolean, True by default.
        A variable to enable or disable the printings of this function.

    Returns
    -------
    path_p_values : pandas.DataFrame
        A dataframe reporting the respective p-values for the pathways that had associated genes containing a p-value
        in p_val_matrix. Additionally, the corrected, mean and std Z values as well as gene number for the given pathway
        are reported in each case.
    '''
    if verbose:
        print('Running reporter pathways analysis')
    # Drop nan genes
    df = p_val_df.dropna(how='all', axis=0)

    # Evaluate information of dataframe
    if 'value' not in list(df.columns):
        cols = list(df.columns)
        df.rename(columns={cols[0] : 'value'}, inplace=True)

    if not isinstance(df.index, str):
        df.index = df.index.map(str)

    # Get gene Z scores
    gene_Z_scores = pd.DataFrame(stats.norm.ppf(df.values) * -1.0, index=df.index, columns=['value'])

    # Convert inf values to numerical values
    gene_Z_scores = gene_Z_scores.replace(np.inf, 15.0)
    gene_Z_scores = gene_Z_scores.replace(-np.inf, -15.0)
    gene_Z_scores = gene_Z_scores.dropna()

    # Genes - Rxn - SubSystems info
    if rxn_pathways_association is None:
        rxn_info = query.rxn_info_from_genes(model=model,
                                             genes=list(df.index),
                                             verbose=verbose)
    else:
        records = []
        for key, val in rxn_pathways_association.items():
            for reaction in val:
                rxn = model.reactions.get_by_id(reaction)
                if len(rxn.genes) != 0:
                    for gene in rxn.genes:
                        records.append((rxn.id, str(gene.id), key))
        rxn_info = pd.DataFrame.from_records(records, columns=['RxnID', 'GeneID', 'SubSystem'])

    if pathways is not None:
        rxn_info = rxn_info.loc[rxn_info.SubSystem.isin(pathways)]

    rxn_info = rxn_info[['GeneID', 'SubSystem']]
    rxn_info = rxn_info.loc[rxn_info.GeneID != '']
    rxn_info = rxn_info.loc[rxn_info.SubSystem != '']
    rxn_info.drop_duplicates(inplace=True)

    unique_pathways = rxn_info.SubSystem.unique()

    # For each pathway calculate the aggregate Z-score and keep track of the number of neighbouring genes
    Z_scores = np.empty((len(unique_pathways), 4))
    Z_scores[:] = np.nan
    Z_scores = pd.DataFrame(Z_scores, index=unique_pathways, columns=['Z-score', 'Mean-Z', 'Std-Z', 'Genes-Number'])

    for path in unique_pathways:
        path_genes = rxn_info.loc[rxn_info.SubSystem == path]['GeneID'].unique().tolist()
        path_genes = list(set(path_genes).intersection(set(gene_Z_scores.index)))

        if len(path_genes) > 0:
            Z_scores.loc[path, 'Z-score'] = np.nansum(gene_Z_scores.loc[path_genes]['value'].values) / np.sqrt(len(path_genes))
            Z_scores.loc[path, 'Mean-Z'] = np.nanmean(gene_Z_scores.loc[path_genes]['value'].values)
            Z_scores.loc[path, 'Std-Z'] = np.nanstd(gene_Z_scores.loc[path_genes]['value'].values)
            Z_scores.loc[path, 'Genes-Number'] = len(path_genes)

    # Remove the metabolites which have no Z-scores
    Z_scores = Z_scores.loc[~Z_scores['Z-score'].isna()]

    # Correct for background by calculating the mean Z-score for random sets of the same size as the ones that
    # were found for the pathways
    for i, size in enumerate(Z_scores['Genes-Number'].unique()):
        size = int(size)
        # Sample 100000 sets for each size. Sample with replacement
        n_samples = 100000

        random_Z_set = np.empty((n_samples, size))

        for j in range(size):
            random_Z_set[:, j] = resample(gene_Z_scores.values, n_samples=n_samples).flatten()

        bg_Z = np.nansum(random_Z_set, axis=1) / np.sqrt(size)
        mean_bg_Z = np.nanmean(bg_Z)
        std_bg_Z = np.nanstd(bg_Z)

        Z_scores.loc[Z_scores['Genes-Number'] == size, 'Z-score'] = (Z_scores.loc[Z_scores['Genes-Number'] == size, 'Z-score'].values - mean_bg_Z) / std_bg_Z

    # Calculate p-values
    path_p_values = Z_scores['Z-score'].apply(lambda x: 1.0 - stats.norm.cdf(x)).to_frame()
    path_p_values.rename(columns={'Z-score': 'p-value'}, inplace=True)

    # Report results
    path_p_values['corrected Z'] = Z_scores['Z-score'].values
    path_p_values['mean Z'] = Z_scores['Mean-Z'].values
    path_p_values['std Z'] = Z_scores['Std-Z'].values
    path_p_values['gene number'] = Z_scores['Genes-Number'].values

    #Sort p-values from smallest value.
    path_p_values.sort_values(by='p-value', ascending=True, inplace=True)
    return path_p_values





