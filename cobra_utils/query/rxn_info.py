# -*- coding: utf-8 -*-

from __future__ import absolute_import

import pandas as pd
from cobra_utils.query import get_ids

import warnings


def rxn_info_from_metabolites(model, metabolites, verbose=True):
    '''
    This function looks for all the reactions where the metabolites in the list participate. Also, it retrieves the genes
    associated to those reactions.

    Parameters
    ----------
    model : cobra.core.Model.Model
        A cobra model.

    metabolites : array-like
        An iterable object containing a list of metabolite ids present in the model.

    verbose : boolean, True by default.
        A variable to enable or disable the printings of this function.

    Returns
    -------
    rxn_gene_association : pandas.DataFrame
        A pandas dataframe containing the information retrieved. The columns are :
        'MetID', 'MetName', 'RxnID', 'RxnName', 'GeneID', 'Subsystem', 'RxnFormula'
    '''
    if verbose:
        print('Using list of metabolites to get reactions where they participate. Also, getting genes of those reactions.')

    rxn_gene_association = []
    mets_ = set(metabolites).intersection(set(get_ids.get_met_ids(model)))
    if verbose:
        excluded = set(metabolites) - mets_
        warnings.warn('{} are not in the model'.format(excluded))

    for metabolite in mets_:
        met = model.metabolites.get_by_id(metabolite)
        for rxn in met.reactions:
            if len(rxn.genes) != 0:
                for gene in rxn.genes:
                    rxn_gene_association.append(
                        (met.id, met.name, rxn.id, rxn.name, str(gene.id), rxn.subsystem, rxn.reaction))
            else:
                rxn_gene_association.append(
                    (met.id, met.name, rxn.id, rxn.name, '', rxn.subsystem, rxn.reaction))

    labels = ['MetID', 'MetName', 'RxnID', 'RxnName', 'GeneID', 'Subsystem', 'RxnFormula']
    rxn_gene_association = pd.DataFrame.from_records(rxn_gene_association, columns=labels)
    if verbose:
        print('Information correctly obtained.')
    return rxn_gene_association


def rxn_info_from_reactions(model, reactions, verbose=True):
    '''
    This function looks for all the reactions and genes that are associated from a list of reactions ids.

    Parameters
    ----------
    model : cobra.core.Model.Model
        A cobra model.

    reactions : array-like
        An iterable object containing a list of reaction ids present in the model.

    verbose : boolean, True by default.
        A variable to enable or disable the printings of this function.

    Returns
    -------
    rxn_gene_association : pandas.DataFrame
        A pandas dataframe containing the information retrieved. The columns are :
        'RxnID', 'RxnName', 'GeneID', 'SubSystem', 'RxnFormula'
    '''
    if verbose:
        print('Using list of reactions to get their information and genes associated.')

    rxn_gene_association = []
    rxns_ = set(reactions).intersection(set(get_ids.get_rxn_ids(model)))
    if verbose:
        excluded = set(reactions) - rxns_
        warnings.warn('{} are not in the model'.format(excluded))

    for reaction in rxns_:
        rxn = model.reactions.get_by_id(reaction)
        if len(rxn.genes) != 0:
            for gene in rxn.genes:
                rxn_gene_association.append((rxn.id, rxn.name, str(gene.id), rxn.subsystem, rxn.reaction))
        else:
            rxn_gene_association.append((rxn.id, rxn.name, '', rxn.subsystem, rxn.reaction))
    labels = ['RxnID', 'RxnName', 'GeneID', 'SubSystem', 'RxnFormula']
    rxn_gene_association = pd.DataFrame.from_records(rxn_gene_association, columns=labels)
    if verbose:
        print('Information correctly obtained.')
    return rxn_gene_association


def rxn_info_from_genes(model, genes, verbose=True):
    '''
    This function looks for all the reactions and genes that are associated from a list of gene ids.

    Parameters
    ----------
    model : cobra.core.Model.Model
        A cobra model.

    genes : array-like
        An iterable object containing a list of gene ids present in the model.

    verbose : boolean, True by default.
        A variable to enable or disable the printings of this function.

    Returns
    -------
    rxn_gene_association : pandas.DataFrame
        A pandas dataframe containing the information retrieved. The columns are :
        'GeneID', 'RxnID', 'RxnName', 'SubSystem', 'RxnFormula'
    '''
    if verbose:
        print('Using list of genes to get the reactions associated and their information.')

    rxn_gene_association = []
    genes_ = set(genes).intersection(set(get_ids.get_gene_ids(model)))
    if verbose:
        excluded = set(genes) - genes_
        warnings.warn('{} are not in the model'.format(excluded))

    for gene in genes_:
        g = model.genes.get_by_id(gene)
        for rxn in g.reactions:
                rxn_gene_association.append((str(g.id), rxn.id, rxn.name, rxn.subsystem, rxn.reaction))
    labels = ['GeneID', 'RxnID', 'RxnName', 'SubSystem', 'RxnFormula']
    rxn_gene_association = pd.DataFrame.from_records(rxn_gene_association, columns=labels)
    if verbose:
        print('Information correctly obtained.')
    return rxn_gene_association


def rxn_info_from_model(model, verbose=True):
    '''
    This function looks for all the reactions in the model and returns their respective information.

    Parameters
    ----------
    model : cobra.core.Model.Model
        A cobra model.

    verbose : boolean, True by default.
        A variable to enable or disable the printings of this function.

    Returns
    -------
    rxn_gene_association : pandas.DataFrame
        A pandas dataframe containing the information retrieved. The columns are :
        'RxnID', 'RxnName', 'GeneID', 'SubSystem', 'RxnFormula'
    '''
    if verbose:
        print('Getting information for all reactions in the model.')

    rxn_gene_association = []
    for rxn in model.reactions:
        if len(rxn.genes) != 0:
            for gene in rxn.genes:
                rxn_gene_association.append((rxn.id, rxn.name, str(gene.id), rxn.subsystem, rxn.reaction))
        else:
            rxn_gene_association.append((rxn.id, rxn.name, '', rxn.subsystem, rxn.reaction))
    labels = ['RxnID', 'RxnName', 'GeneID', 'SubSystem', 'RxnFormula']
    rxn_gene_association = pd.DataFrame.from_records(rxn_gene_association, columns=labels)
    if verbose:
        print('Information correctly obtained.')
    return rxn_gene_association