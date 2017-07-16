import numpy as np
import pandas as pd
import scipy.optimize as opt
path = "C:/Users/RJ/Desktop/testdir/test.xlsx"


def import_data(path):
    """
    ONLY WORKS IF DATA IS IN PERFECT ORDER. WILL WRITE UPDATED FUNCTION IN THE FUTURE TO HANDLE ACTUAL DATA LAYOUT.

    Parameters
    ----------
    path: location of file

    Returns
    -------
    imported data from excel file without labels

    """
    excel = pd.read_excel(path, header=None)
    return excel


def labels(data):
    row = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    column = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12']
    row = np.array([row[i] for i in range(data.shape[0])])
    column = np.array([column[i] for i in range(data.shape[1])])
    return row, column


def add_sample_name(data, orientation=0):
    """
    Parameters
    ----------
    data: imported data
    orientation: same convention as axis

    Returns
    -------
    Adds sample name replicates depending on orientation

    """
    sample = ['S1', 'S1', 'S2', 'S2', 'S3', 'S3', 'S4', 'S4', 'S5', 'S5', 'S6', 'S6']
    rowIndex = pd.MultiIndex.from_arrays(
        [np.array([sample[i] for i in range(data.shape[0])]), labels(data)[0]])
    colIndex = pd.MultiIndex.from_arrays(
        [np.array([sample[i] for i in range(data.shape[1])]), labels(data)[1]])
    if orientation == 0:
        df = pd.DataFrame(data=data.values, index=rowIndex, columns=labels(data)[1])
        return df
    if orientation == 1:
        df = pd.DataFrame(data=data.values, index=labels(data)[0], columns=colIndex)
        return df


def get_replicates(raw_data, orientation=0):
    """
    Parameters
    ----------
    raw_data: raw data from file
    orientation: same convention as axis function. used to decide the axis of replicates to be returned

    Returns
    -------
    Replicate separation depending on plate orientation

    Warning
    -------
    You're not supposed to call list/tuple with []. idk how to get rid of it or what to change it to

    """

    if np.shape(raw_data)[0] % 2 != 0 or np.shape(raw_data)[1] % 2 != 0:
        return 'uneven replicates! cannot split evenly'
    if orientation == 0:
        return np.split(raw_data, np.shape(raw_data)[0]/2, axis=0)
    if orientation == 1:
        return np.split(raw_data, np.shape(raw_data)[1]/2, axis=1)


def calculate_cv(raw_data, orientation=0):
    """
    Parameters
    ----------
    raw_data: raw data from file
    orientation: same convention as axis function. used to determine axis calculations

    Returns
    -------

    """
    if orientation == 0:
        average = np.array([np.mean(x, axis=0) for x in get_replicates(raw_data, orientation=0)])
        std = np.array([np.std(x, axis=0, ddof=1) for x in get_replicates(raw_data, orientation=0)])
        with np.errstate(divide='ignore', invalid='ignore'):
            cv = pd.DataFrame(np.array(np.nan_to_num(std/average)) * 100)
            return cv
    if orientation == 1:
        average = np.array([np.mean(x, axis=1) for x in get_replicates(raw_data, orientation=1)])
        std = np.array([np.std(x, axis=1, ddof=1) for x in get_replicates(raw_data, orientation=1)])
        with np.errstate(divide='ignore', invalid='ignore'):
            cv = pd.DataFrame(np.array(np.nan_to_num(std/average)) * 100)
            return cv


def get_concentrations(starting_concentration, dilution_ratio, n_dilutions, graph_type='inhibition'):
    """
    Parameters
    ----------
    starting_concentration: desired starting concentration for assay
    dilution_ratio: desired dilution ratio for serial dilutions
    n_dilutions: number of dilutions
    graph_type: inhibition or drc

    Returns
    -------
    x-axis concentrations for assays in non log

    """
    if graph_type == 'inhibition':
        concentrations = np.array(starting_concentration * np.power(dilution_ratio, range(n_dilutions)))
        return concentrations
    if graph_type == 'drc':
        concentrations = starting_concentration / np.power(dilution_ratio, range(n_dilutions))
        return concentrations


def log_dilution(concentrations):
    """
    Parameters
    ----------
    concentrations: concentrations

    Returns
    -------
    log(concentrations)

    """
    return np.log10(concentrations)


def inhibition(concentration, bottom, top, logIC50, hill_slope):
    """
    Parameters
    ----------
    concentration: list of concentrations in log
    bottom: bottom of curve
    top: top of curve
    logIC50: LogIC50
    hill_slope: HillSlope

    Returns
    -------
    Response

    """
    return bottom + (top - bottom)/(1+np.power(10, ((logIC50-concentration) * hill_slope)))


def inhibition_coefficients(data, concentrations):
    """
    Parameters
    ----------
    data: labeled data
    concentrations: concentrations in log

    Returns
    -------
    coefficients for each replicate

    """
    coefficient_storage = []
    for i, j in data.groupby(level=0):
        replicate_average = np.mean(j, axis=0)
        coefficients, d = opt.curve_fit(inhibition, concentrations, replicate_average)
        curve_coefficients = dict(zip(['top', 'bottom', 'logIC50', 'hill_slope'], coefficients))
        coefficient_storage.append(curve_coefficients)
    coefficient_storage = pd.DataFrame(data=coefficient_storage)
    return coefficient_storage

