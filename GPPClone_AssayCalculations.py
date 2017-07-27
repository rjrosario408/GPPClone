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
    row = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
    column = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12']
    row = np.array([row[i] for i in range(imported_data.shape[0])])
    column = np.array([column[i] for i in range(imported_data.shape[1])])
    row_index = pd.MultiIndex.from_arrays(
        [np.array([sample[i] for i in range(imported_data.shape[0])]), row])
    col_index = pd.MultiIndex.from_arrays(
        [np.array([sample[i] for i in range(imported_data.shape[1])]), column])
    if orientation == 0:
        df = pd.DataFrame(data=imported_data.values, index=row_index, columns=column)
        return df
    if orientation == 1:
        df = pd.DataFrame(data=imported_data.values, index=row, columns=col_index)
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
    with np.errstate(divide='ignore', invalid='ignore'):
        cv = data.std(axis=0, level=0, ddof=1)/data.mean(axis=0, level=0) * 100
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
        concentrations = pd.DataFrame(starting_concentration * np.power(dilution_ratio, range(n_dilutions)),
                                      index=[i for i in range(1, n_dilutions+1)])
        return concentrations
    if graph_type == 'drc':
        concentrations = pd.DataFrame(starting_concentration / np.power(dilution_ratio, range(n_dilutions)),
                                      index=[i for i in range(1, n_dilutions+1)])
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


def inhibition(concentration, bottom, top, log_ic50, hill_slope):
    """
    Parameters
    ----------
    concentration: list of concentrations in log
    bottom: bottom of curve
    top: top of curve
    log_ic50: LogIC50
    hill_slope: HillSlope

    Returns
    -------
    Inhibition Response

    """
    return bottom + (top - bottom)/(1+np.power(10, ((log_ic50-concentration) * hill_slope)))


def drc(concentration, bottom, top, log_ec50, hill_slope):
    """
    Parameters
    ----------
    concentration: list of concentrations in log
    bottom: bottom of curve
    top: top of curve
    log_ec50: LogIC50
    hill_slope: HillSlope

    Returns
    -------
    Dose Response

    """
    return bottom + (top - bottom)/(1+np.power(10, ((log_ec50-concentration) * hill_slope)))


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
        x_min = j.loc[:, '1'].mean()
        x_max = j.loc[:, '12'].mean()
        y_transform.append(j.loc[:, '2':'11'].apply
                           (lambda x: 100-((x-x_min)/(x_max-x_min)) * 100).mean(axis=0, level=0))
    return pd.concat(y_transform)


def error_transform_y(data):
    """
    Parameters
    ----------
    data: named data

    Returns
    -------
    STD between transformed replicates

    """
    if not isinstance(data.index, pd.core.index.MultiIndex):
        data = data.T
    y_transform_error = []
    for i, j in data.groupby(level=0):
        x_min = j.loc[:, '1'].mean()
        x_max = j.loc[:, '12'].mean()
        y_transform_error.append(j.loc[:, '2':'11'].apply
                                 (lambda x: 100-((x-x_min)/(x_max-x_min)) * 100).std(axis=0, level=0, ddof=1))
    return pd.concat(y_transform_error)


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
    for i, j in response.groupby(level=0):
        coefficients, d = opt.curve_fit(inhibition, concentrations, j.iloc[0, :])
        curve_coefficients = dict(zip(['top', 'bottom', 'logIC50', 'hill_slope'], coefficients))
        coefficient_storage.append(curve_coefficients)
    coefficient_storage = pd.DataFrame(
        data=coefficient_storage, index=[('S'+str(i+1)) for i in range(len(coefficient_storage))]).round(decimals=3)
    return coefficient_storage


def normalize_y_drc(data):
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
        x_min = j.loc[:, '1'].mean()
        x_max = j.loc[:, '12'].mean()
        y_transform.append(j.loc[:, '1':'11'].apply
                           (lambda x: 100 - ((x-x_min)/(x_max-x_min)) * 100).mean(axis=0, level=0))
    return pd.concat(y_transform)


def error_normalize_ydrc(data):
    """
    Parameters
    ----------
    data: named data

    Returns
    -------
    STD between transformed replicates

    """
    if not isinstance(data.index, pd.core.index.MultiIndex):
        data = data.T
    y_transform_error = []
    for i, j in data.groupby(level=0):
        x_min = j.loc[:, '1'].mean()
        x_max = j.loc[:, '12'].mean()
        y_transform_error.append(j.loc[:, '1':'11'].apply
                                 (lambda x: 100-((x-x_min)/(x_max-x_min)) * 100).std(axis=0, level=0, ddof=1))
    return pd.concat(y_transform_error)


def drc_coefficients(response, concentrations):
    """
    Parameters
    ----------
    response: transformed y response
    concentrations: concentrations in log

    Returns
    -------
    drc coefficients for each replicate

    """
    coefficient_storage = []
    concentrations = concentrations.iloc[:, 0]
    for i, j in response.groupby(level=0):
        coefficients, d = opt.curve_fit(drc, concentrations, j.iloc[0, :])
        curve_coefficients = dict(zip(['bottom', 'top', 'logEC50', 'hill_slope'], coefficients))
        coefficient_storage.append(curve_coefficients)
    coefficient_storage = pd.DataFrame(
        data=coefficient_storage, index=[('S'+str(i+1)) for i in range(len(coefficient_storage))]).round(decimals=3)
    return coefficient_storage


def inhibition_graph(data, concentrations):
    """
    Parameters
    ----------
    data: named data
    concentrations: desired x values in non log.

    Returns
    -------
    inhibition response vs. log dilutions and table with calculated fit parameter

    """
    log_x = log_dilution(concentrations)
    concentrations = concentrations.iloc[:, 0]
    y = transform_y(data)
    error = error_transform_y(data)
    fit = inhibition_coefficients(y, log_x)

    sns.plt.axes(xscale='log')
    sns.plt.xlabel("Log Dilutions")
    sns.plt.ylabel("%inhibition")
    sns.plt.title("%inhibition vs Log Dilutions")

    for i in range(y.shape[0]):
        sns.plt.errorbar(concentrations, y.values[i], yerr=error.values[i], fmt='-o', label=('S'+str(i+1)))
    sns.plt.legend()
    sns.plt.table(cellText=fit.values, colWidths=[0.25] * len(fit.columns), rowLabels=fit.index, colLabels=fit.columns,
                  cellLoc='center', rowLoc='center', loc='bottom')


def drc_graph(data, concentrations):
    """
    Parameters
    ----------
    data: named data
    concentrations: desired x values in non log.

    Returns
    -------
    drc response vs. log dilutions and table with calculated fit parameter

    """
    log_x = log_dilution(concentrations)
    concentrations = concentrations.iloc[:, 0]
    y = normalize_y_drc(data)
    error = error_normalize_ydrc(data)
    fit = drc_coefficients(y, log_x)

    sns.plt.axes(xscale='log')
    sns.plt.xlabel("Log Dose")
    sns.plt.ylabel("Response")
    sns.plt.title("Dose Response Curve")

    for i in range(y.shape[0]):
        sns.plt.errorbar(concentrations, y.values[i], yerr=error.values[i], fmt='-o', label=('S'+str(i+1)))
    sns.plt.legend()
    sns.plt.table(cellText=fit.values, colWidths=[0.25] * len(fit.columns), rowLabels=fit.index, colLabels=fit.columns,
                  cellLoc='center', rowLoc='center', loc='bottom')
