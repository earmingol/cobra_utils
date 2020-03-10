# -*- coding: utf-8 -*-

from __future__ import absolute_import

import numpy as np
import pandas as pd
import scipy.stats as stats

from sklearn.utils import resample
from cobra_utils import query


def reporter_metabolites(model, p_val_df, genes=None, verbose=True):
    '''
    This function computes an aggregate p-value for each metabolite based on the network topology of the metabolic
    reconstruction. It takes the p-value for differential expression of each gene and compute the aggregate p-value
    for the neighbor reactions of a given metabolite.

    More information on:
    https://www.pnas.org/cgi/doi/10.1073/pnas.0406811102

    This code was adapted from RAVEN 2.0 code available on:
    https://github.com/SysBioChalmers/RAVEN/blob/master/core/reporterMetabolites.m

    Parameters
    ----------
    model : cobra.core.Model.Model
        A cobra model.

    p_val_df : pandas.DataFrame
        A dataframe with gene names as index. It have to contains the p-values for the differential expression
        of the respective indexing genes.

    genes : array-like
        An array or list containing gene names (str) to be considered.

    verbose : boolean, True by default.
        A variable to enable or disable the printings of this function.

    Returns
    -------
    met_p_values : pandas.DataFrame
        A dataframe reporting the respective p-values for the metabolites that had associated genes containing a p-value
        in p_val_matrix. Additionally, the corrected, mean and std Z values as well as gene number for the given metabolite
        are reported in each case.
    '''
    if verbose:
        print('Running reporter metabolites analysis')
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

    # Mets - Genes info
    met_info = query.met_info_from_model(model=model,
                                         verbose=verbose)
    met_info = met_info.loc[met_info.GeneID.isin(list(gene_Z_scores.index))]

    if genes is not None:
        met_info = met_info.loc[met_info.GeneID.isin(genes)]

    unique_mets = met_info.MetID.unique()

    met_info = met_info[['MetID', 'GeneID']]
    met_info = met_info.loc[met_info.GeneID != '']
    met_info.drop_duplicates(inplace=True)

    # For each metabolite calculate the aggregate Z-score and keep track of the number of neighbouring genes
    Z_scores = np.empty((len(unique_mets), 4))
    Z_scores[:] = np.nan
    Z_scores = pd.DataFrame(Z_scores, index=unique_mets, columns=['Z-score', 'Mean-Z', 'Std-Z', 'Genes-Number'])

    for met in unique_mets:
        met_genes = met_info.loc[met_info.MetID == met]['GeneID'].unique()
        met_genes = list(set(met_genes).intersection(set(gene_Z_scores.index)))

        if len(met_genes) > 0:
            Z_scores.loc[met, 'Z-score'] = np.nansum(gene_Z_scores.loc[met_genes]['value'].values) / np.sqrt(len(met_genes))
            Z_scores.loc[met, 'Mean-Z'] = np.nanmean(gene_Z_scores.loc[met_genes]['value'].values)
            Z_scores.loc[met, 'Std-Z'] = np.nanstd(gene_Z_scores.loc[met_genes]['value'].values)
            Z_scores.loc[met, 'Genes-Number'] = len(met_genes)

    # Remove the metabolites which have no Z-scores
    Z_scores = Z_scores.loc[~Z_scores['Z-score'].isna()]

    # Correct for background by calculating the mean Z-score for random sets of the same size as the ones that
    # were found for the metabolites
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
    met_p_values = Z_scores['Z-score'].apply(lambda x: 1.0 - stats.norm.cdf(x)).to_frame()
    met_p_values.rename(columns={'Z-score': 'p-value'}, inplace=True)

    # Report results
    met_p_values['corrected Z'] = Z_scores['Z-score'].values
    met_p_values['mean Z'] = Z_scores['Mean-Z'].values
    met_p_values['std Z'] = Z_scores['Std-Z'].values
    met_p_values['gene number'] = Z_scores['Genes-Number'].values

    #Sort p-values from smallest value.
    met_p_values.sort_values(by='p-value', ascending=True, inplace=True)
    return met_p_values





