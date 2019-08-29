# -*- coding: utf-8 -*-

from __future__ import absolute_import


def get_rxn_ids(model):
    '''
    This function returns a list of IDs of all reactions in the model.

    Parameters
    ----------
    model : cobra.core.Model.Model
        A cobra model.

    Returns
    -------
    rxns : list
        A list containing all IDs of rxns in the model.
    '''
    rxns = []
    for reaction in model.reactions:
        rxns.append(reaction.id)
    rxns = list(set(rxns))
    return rxns


def get_gene_ids(model):
    '''
    This function returns a list of IDs of all genes in the model.

    Parameters
    ----------
    model : cobra.core.Model.Model
        A cobra model.

    Returns
    -------
    genes : list
        A list containing all IDs of genes in the model.
    '''
    genes = []
    for gene in model.genes:
        genes.append(gene.id)
    genes = list(set(genes))
    return genes


def get_met_ids(model):
    '''
    This function returns a list of IDs of all metabolites in the model.

    Parameters
    ----------
    model : cobra.core.Model.Model
        A cobra model.

    Returns
    -------
    mets : list
        A list containing all IDs of metabolites in the model.
    '''
    mets = []
    for met in model.metabolites:
        mets.append(met.id)
    mets = list(set(mets))
    return mets