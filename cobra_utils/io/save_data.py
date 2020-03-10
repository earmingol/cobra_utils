# -*- coding: utf-8 -*-

from __future__ import absolute_import

import cobra


def save_model(filename, format='matlab', verbose=True, **kwargs):
    '''
    This function saves a metabolic reconstruction into a given format.

    Parameters
    ----------
    filename : str
        Filename of the model to save. It is preferable to use absolute path.

    format : str, 'matlab' by default.
        Format of the file containing the model. Options to use:
        'json' for .json file
        'matlab' for .mat file
        'sbml' for .xml file
        'yaml' for  .yaml or .yml file

    verbose : boolean, True by default
        A variable to enable or disable the printings of this function.
    '''
    if verbose:
        print('Saving genome-scale model')
    try:
        if format == 'json':
            model = cobra.io.save_json_model(filename, **kwargs)
        elif format == 'matlab':
            model = cobra.io.save_matlab_model(filename, **kwargs)
        elif format == 'sbml':
            model = cobra.io.write_sbml_model(filename, **kwargs)
        elif format == 'yaml':
            model = cobra.io.save_yaml_model(filename, **kwargs)
        else:
            raise NotImplementedError("Format {} not implemented. Specify a correct format for the model".format(format))
    except:
        raise ImportError("The file has an incorrect format or does not match with implemented formats")
    if verbose:
        print('Model correctly saved.')
    return model