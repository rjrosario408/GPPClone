import numpy as np
import pandas as pd
import scipy.optimize as opt
import seaborn as sns
location = "C:/Users/RJ/Desktop/testdir/lmao_test.xlsx"


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


def labels(imported_data):
    """
    Parameters
    ----------
    imported_data: takes in imported data

    Returns
    -------
    plate labels for rows and columns

    """
    row = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    column = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12']
    row = np.array([row[i] for i in range(imported_data.shape[0])])
    column = np.array([column[i] for i in range(imported_data.shape[1])])
    return row, column


def add_sample_name(imported_data, orientation=0):
    """
    Parameters
    ----------
    imported_data: imported data
    orientation: same convention as axis

    Returns
    -------
    Adds sample name replicates depending on orientation

    """
    sample = ['S1', 'S1', 'S2', 'S2', 'S3', 'S3', 'S4', 'S4', 'S5', 'S5', 'S6', 'S6']
    rowIndex = pd.MultiIndex.from_arrays(
        [np.array([sample[i] for i in range(imported_data.shape[0])]), labels(imported_data)[0]])
    colIndex = pd.MultiIndex.from_arrays(
        [np.array([sample[i] for i in range(imported_data.shape[1])]), labels(imported_data)[1]])
    if orientation == 0:
        df = pd.DataFrame(data=imported_data.values, index=rowIndex, columns=labels(imported_data)[1])
        return df
    if orientation == 1:
        df = pd.DataFrame(data=imported_data.values, index=labels(imported_data)[0], columns=colIndex)
        return df


def calculate_cv(data):
    """
    Parameters
    ----------
    data: labeled data

    Returns
    -------
    CV for each replicate

    """
    if not isinstance(data.index, pd.core.index.MultiIndex):
        data = data.T
    average = np.array([np.mean(j, axis=0) for i, j in data.groupby(level=0)])
    std = np.array([np.std(j, axis=0, ddof=1) for i, j in data.groupby(level=0)])
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
        concentrations = pd.DataFrame(starting_concentration * np.power(dilution_ratio, range(n_dilutions)))
        return concentrations
    if graph_type == 'drc':
        concentrations = pd.DataFrame(starting_concentration * np.power(dilution_ratio, range(n_dilutions)))
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


def transform_y(data):
    """
    Parameters
    ----------
    data: named data

    Returns
    -------
    100 - feature scaled *100 value of response

    """
    if not isinstance(data.index, pd.core.index.MultiIndex):
        data = data.T
    y_transform = []
    for i, j in data.groupby(level=0):
        length = j.shape[1]-1
        x_min = np.average([j.iloc[0, 0], j.iloc[1, 0]])
        x_max = np.average([j.iloc[0, length], j.iloc[1, length]])
        y_transform.append(np.average(j.apply(lambda x: 100-((x-x_min)/(x_max-x_min)) * 100), axis=0)[1:length])
    return pd.DataFrame(data=y_transform, columns=[i for i in range(1, len(y_transform[1])+1)])


def inhibition_coefficients(response, concentrations):
    """
    Parameters
    ----------
    response: transformed y response
    concentrations: concentrations in log

    Returns
    -------
    coefficients for each replicate

    """
    coefficient_storage = []
    concentrations = concentrations.iloc[:, 0]
    response = np.array(response)
    for i in response:
        coefficients, d = opt.curve_fit(inhibition, concentrations, i[:])
        curve_coefficients = dict(zip(['top', 'bottom', 'logIC50', 'hill_slope'], coefficients))
        coefficient_storage.append(curve_coefficients)
    coefficient_storage = pd.DataFrame(data=coefficient_storage,
                                       index=[('S'+str(i+1)) for i in range(len(coefficient_storage))]).round(decimals=3)
    return coefficient_storage


def graph(concentrations, response, cv, fit):
    """
    Parameters
    ----------
    concentrations: concentrations in non log
    response: transformed y data
    cv: right now using CV as error calculation as a place holder. will change in the future for actual error calc
    fit: coefficients from inhibitions coefficients

    Returns
    -------
    inhibition vs log dilution graph and a table with all the fit coefficients

    """
    concentrations = concentrations.iloc[:, 0]
    cv = cv.iloc[:, 1:cv.shape[1]-1]
    sns.plt.axes(xscale='log')
    sns.plt.xlabel("Log Dilutions")
    sns.plt.ylabel("%inhibition")
    sns.plt.title("%inhibition vs Log Dilutions")
    for i in range(response.shape[0]):
        sns.plt.errorbar(concentrations, response.values[i], yerr=cv.values[i], fmt='-o', label=('S'+str(i+1)))
    sns.plt.legend()
    sns.plt.table(cellText=fit.values, colWidths=[0.25] * len(fit.columns), rowLabels=fit.index, colLabels=fit.columns,
                  cellLoc='center', rowLoc='center', loc='bottom')
