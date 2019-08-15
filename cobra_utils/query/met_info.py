# -*- coding: utf-8 -*-

from __future__ import absolute_import

import pandas as pd

def met_info_from_model(model, verbose=True):
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
    met_rxn_gene_association : pandas.DataFrame
        A pandas dataframe containing the information retrieved. The columns are :
        'MetID', GeneID', 'RxnID', 'RxnName', 'SubSystem', 'RxnFormula'
    '''
    if verbose:
        print('Getting information for all metabolites in the model.')

    met_association = []
    for met in model.metabolites:
        for rxn in met.reactions:
            if len(rxn.genes) != 0:
                for gene in rxn.genes:
                    met_association.append((str(met.id), str(gene.id), rxn.id, rxn.name, rxn.subsystem, rxn.reaction))
            else:
                met_association.append((str(met.id), '', rxn.id, rxn.name, rxn.subsystem, rxn.reaction))
    labels = ['MetID', 'GeneID', 'RxnID', 'RxnName', 'SubSystem', 'RxnFormula']
    met_association = pd.DataFrame.from_records(met_association, columns=labels)
    if verbose:
        print('Information correctly obtained.')
    return met_association