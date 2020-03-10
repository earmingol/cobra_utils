# -*- coding: utf-8 -*-

from __future__ import absolute_import

import pandas as pd
from cobra_utils.query import get_ids

import warnings


def met_info_from_metabolites(model, metabolites, verbose=True):
    '''
    This function looks for all the metabolites in a list and find their reaction association. Also, it retrieves the genes
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
    met_rxn_gene_association : pandas.DataFrame
        A pandas dataframe containing the information retrieved. The columns are :
        'MetID', 'MetName', 'RxnID', 'RxnName', 'GeneID', 'Subsystem', 'RxnFormula'
    '''
    if verbose:
        print('Using list of metabolites to get information where they participate. Also, getting associated reactions and genes.')

    met_rxn_gene_association = []
    mets_ = set(metabolites).intersection(set(get_ids.get_met_ids(model)))
    if verbose:
        excluded = set(metabolites) - mets_
        warnings.warn('{} are not in the model'.format(excluded))

    for metabolite in mets_:
        met = model.metabolites.get_by_id(metabolite)
        for rxn in met.reactions:
            if len(rxn.genes) != 0:
                for gene in rxn.genes:
                    met_rxn_gene_association.append((met.id, met.name, rxn.id, rxn.name, str(gene.id), rxn.subsystem, rxn.reaction))
            else:
                met_rxn_gene_association.append((met.id, met.name, rxn.id, rxn.name, '', rxn.subsystem, rxn.reaction))

    labels = ['MetID', 'MetName', 'RxnID', 'RxnName', 'GeneID', 'Subsystem', 'RxnFormula']
    met_rxn_gene_association = pd.DataFrame.from_records(met_rxn_gene_association, columns=labels)
    if verbose:
        print('Information correctly obtained.')
    return met_rxn_gene_association


def met_info_from_reactions(model, reactions, verbose=True):
    '''
    This function looks for all the metabolites involved in reactions that are in a list.

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
    met_rxn_gene_association : pandas.DataFrame
        A pandas dataframe containing the information retrieved. The columns are :
        'RxnID', 'RxnName', 'MetID', 'MetName', 'GeneID', 'Subsystem', 'RxnFormula'
    '''
    if verbose:
        print('Using list of reactions to get information about metabolites and genes associated.')

    met_rxn_gene_association = []
    rxns_ = set(reactions).intersection(set(get_ids.get_rxn_ids(model)))
    if verbose:
        excluded = set(reactions) - rxns_
        warnings.warn('{} are not in the model'.format(excluded))

    for reaction in rxns_:
        rxn = model.reactions.get_by_id(reaction)
        for met in rxn.metabolites:
            if len(rxn.genes) != 0:
                for gene in rxn.genes:
                    met_rxn_gene_association.append((rxn.id, rxn.name, met.id, met.name, str(gene.id), rxn.subsystem, rxn.reaction))
            else:
                met_rxn_gene_association.append((rxn.id, rxn.name, met.id, met.name,  '', rxn.subsystem, rxn.reaction))
    labels = ['RxnID', 'RxnName', 'MetID', 'MetName', 'GeneID', 'Subsystem', 'RxnFormula']
    met_rxn_gene_association = pd.DataFrame.from_records(met_rxn_gene_association, columns=labels)
    if verbose:
        print('Information correctly obtained.')
    return met_rxn_gene_association


def met_info_from_genes(model, genes, verbose=True):
    '''
    This function looks for all the metabolites involved in reactions that are associated to a list of gene ids.

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
    met_rxn_gene_association : pandas.DataFrame
        A pandas dataframe containing the information retrieved. The columns are :
        'GeneID', 'MetID', 'MetName', 'RxnID', 'RxnName', 'SubSystem', 'RxnFormula'
    '''
    if verbose:
        print('Using list of genes to get the metabolites associated and their information.')

    met_rxn_gene_association = []
    genes_ = set(genes).intersection(set(get_ids.get_gene_ids(model)))
    if verbose:
        excluded = set(genes) - genes_
        warnings.warn('{} are not in the model'.format(excluded))

    for gene in genes_:
        g = model.genes.get_by_id(gene)
        for rxn in g.reactions:
            for met in rxn.metabolites:
                met_rxn_gene_association.append((str(g.id), met.id, met.name, rxn.id, rxn.name, rxn.subsystem, rxn.reaction))
    labels = ['GeneID', 'MetID', 'MetName', 'RxnID', 'RxnName', 'SubSystem', 'RxnFormula']
    met_rxn_gene_association = pd.DataFrame.from_records(met_rxn_gene_association, columns=labels)
    if verbose:
        print('Information correctly obtained.')
    return met_rxn_gene_association


def met_info_from_model(model, verbose=True):
    '''
    This function looks for all the metabolites in the model and returns their respective information.

    Parameters
    ----------
    model : cobra.core.Model.Model
        A cobra model.

    verbose : boolean, True by default.
        A variable to enable or disable the printings of this function.

    Returns
    -------
    met_rxn_gene_association : pandas.DataFrame
        A pandas dataframe containing the information retrieved. The columns are :
        'MetID', 'MetName', 'GeneID', 'RxnID', 'RxnName', 'SubSystem', 'RxnFormula'
    '''
    if verbose:
        print('Getting information for all metabolites in the model.')

    met_association = []
    for met in model.metabolites:
        for rxn in met.reactions:
            if len(rxn.genes) != 0:
                for gene in rxn.genes:
                    met_association.append((met.id, met.name, str(gene.id), rxn.id, rxn.name, rxn.subsystem, rxn.reaction))
            else:
                met_association.append((met.id, met.name,  '', rxn.id, rxn.name, rxn.subsystem, rxn.reaction))
    labels = ['MetID', 'MetName', 'GeneID', 'RxnID', 'RxnName', 'SubSystem', 'RxnFormula']
    met_association = pd.DataFrame.from_records(met_association, columns=labels)
    if verbose:
        print('Information correctly obtained.')
    return met_association