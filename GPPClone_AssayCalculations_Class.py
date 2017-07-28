import numpy as np
import pandas as pd
import scipy.optimize as opt
import seaborn as sns
location = "C:/Users/RJ/Desktop/testdir/test20170720/Raw.xlsx"


class Nab (object):
    """
    Take in drc or inhibition data and return graph
    """

    def __init__(self, path, orientation=0):
        self.path = path
        self.orientation = orientation
        self.data = self.add_sample_name()

    def add_sample_name(self):
        excel = pd.read_excel(self.path, header=None)
        sample = ['S1', 'S1', 'S2', 'S2', 'S3', 'S3', 'S4', 'S4', 'S5', 'S5', 'S6', 'S6']
        row = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
        column = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12']
        row = np.array([row[i] for i in range(excel.shape[0])])
        column = np.array([column[i] for i in range(excel.shape[1])])
        row_index = pd.MultiIndex.from_arrays(
            [np.array([sample[i] for i in range(excel.shape[0])]), row])

        col_index = pd.MultiIndex.from_arrays(
            [np.array([sample[i] for i in range(excel.shape[1])]), column])
        if self.orientation == 0:
            df = pd.DataFrame(data=excel.values, index=row_index, columns=column)
            return df
        if self.orientation == 1:
            df = pd.DataFrame(data=excel.values, index=row, columns=col_index)
            return df

    def calculate_cv(self):
        if not isinstance(self.data.index, pd.core.index.MultiIndex):
            self.data = self.data.T

        with np.errstate(divide='ignore', invalid='ignore'):
            cv = self.data.std(axis=0, level=0, ddof=1) / self.data.mean(axis=0, level=0) * 100
            return cv

    @staticmethod
    def get_concentrations(starting_concentration, dilution_ratio, n_dilutions, graph_type='inhibition'):
        if graph_type == 'inhibition':
            concentrations = pd.DataFrame(starting_concentration * np.power(dilution_ratio, range(n_dilutions)),
                                          index=[i for i in range(1, n_dilutions + 1)])
            return concentrations
        if graph_type == 'drc':
            concentrations = pd.DataFrame(starting_concentration / np.power(dilution_ratio, range(n_dilutions)),
                                          index=[i for i in range(1, n_dilutions + 1)])
            return concentrations

    @staticmethod
    def log_dilution(concentrations):
        return np.log10(concentrations)

    @staticmethod
    def inhibition(concentration, bottom, top, log_ic50, hill_slope):
        return bottom + (top - bottom) / (1 + np.power(10, ((log_ic50 - concentration) * hill_slope)))

    def transform_y(self):
        if not isinstance(self.data.index, pd.core.index.MultiIndex):
            self.data = self.data.T
        y_transform = []
        for i, j in self.data.groupby(level=0):
            x_min = j.loc[:, '1'].mean()
            x_max = j.loc[:, '12'].mean()
            y_transform.append(j.loc[:, '2':'11'].apply
                               (lambda x: 100 - ((x - x_min) / (x_max - x_min)) * 100).mean(axis=0, level=0))
        return pd.concat(y_transform)

    def error_transform_y(self):
        if not isinstance(self.data.index, pd.core.index.MultiIndex):
            self.data = self.data.T
        y_transform_error = []
        for i, j in self.data.groupby(level=0):
            x_min = j.loc[:, '1'].mean()
            x_max = j.loc[:, '12'].mean()
            y_transform_error.append(j.loc[:, '2':'11'].apply
                                     (lambda x: 100 - ((x - x_min) / (x_max - x_min)) * 100).std(axis=0, level=0,
                                                                                                 ddof=1))
        return pd.concat(y_transform_error)

    def inhibition_coefficients(self):
        coefficient_storage = []
        concentrations = Nab.log_dilution(Nab.get_concentrations(1000, 2, 10)).iloc[:, 0]
        for i, j in self.transform_y().groupby(level=0):
            coefficients, d = opt.curve_fit(Nab.inhibition, concentrations, j.iloc[0, :])
            curve_coefficients = dict(zip(['top', 'bottom', 'logIC50', 'hill_slope'], coefficients))
            coefficient_storage.append(curve_coefficients)
        coefficient_storage = pd.DataFrame(
            data=coefficient_storage, index=[('S' + str(i + 1)) for i in range(len(coefficient_storage))]).round(
            decimals=3)
        return coefficient_storage

    def __repr__(self):
        return "<%s>" % self.data

