import GPPClone_AssayCalculations as GPPCA
import pandas as pd
import pandas.util.testing as pdt


def test_calculate_cv_20170720():
    raw = "C:/Users/RJ/Desktop/testdir/test20170720/Raw.xlsx"
    calculate_cv_path = "C:/Users/RJ/Desktop/testdir/test20170720/calculate_cv.xlsx"
    calculate_cv_code = GPPCA.calculate_cv((GPPCA.add_sample_name(GPPCA.import_data(raw))))
    calculate_cv_gpp = pd.read_excel(calculate_cv_path)
    calculate_cv_gpp.columns = calculate_cv_gpp.columns.astype(str)
    pdt.assert_almost_equal(calculate_cv_code, calculate_cv_gpp)


def test_log_dilution_20170720():
    log_dilution_path = "C:/Users/RJ/Desktop/testdir/test20170720/log_dilution.xlsx"
    log_dilution_code = GPPCA.log_dilution(GPPCA.get_concentrations(1000, 2, 10))
    log_dilution_gpp = pd.read_excel(log_dilution_path)
    pdt.assert_almost_equal(log_dilution_code, log_dilution_gpp)


def test_calculate_y_transform_20170720():
    raw = "C:/Users/RJ/Desktop/testdir/test20170720/Raw.xlsx"
    transform_y_path = "C:/Users/RJ/Desktop/testdir/test20170720/y_transform.xlsx"
    transform_y_code = GPPCA.transform_y(GPPCA.add_sample_name(GPPCA.import_data(raw)))
    transform_y_gpp = pd.read_excel(transform_y_path)
    transform_y_gpp.columns = transform_y_gpp.columns.astype(str)
    pdt.assert_almost_equal(transform_y_code, transform_y_gpp)


def test_inhibition_coefficients():
    raw = "C:/Users/RJ/Desktop/testdir/test20170720/Raw.xlsx"
    inhibition_coefficients_path = "C:/Users/RJ/Desktop/testdir/test20170720/inhibition_coefficients.xlsx"
    response = GPPCA.transform_y(GPPCA.add_sample_name(GPPCA.import_data(raw)))
    logx = GPPCA.log_dilution(GPPCA.get_concentrations(1000, 2, 10))
    inhibition_coefficients_code = GPPCA.inhibition_coefficients(response, logx)
    inhibition_coefficients_gpp = pd.read_excel(inhibition_coefficients_path)
    pdt.assert_frame_equal(inhibition_coefficients_code, inhibition_coefficients_gpp)


def test_normalize_y_drc_20170726():
    raw = "C:/Users/RJ/Desktop/testdir/test20170726/Raw.xlsx"
    normalize_y_path = "C:/Users/RJ/Desktop/testdir/test20170726/normalize_y_drc.xlsx"
    normalize_y_code = GPPCA.normalize_y_drc(GPPCA.add_sample_name(GPPCA.import_data(raw)))
    normalize_y_gpp = pd.read_excel(normalize_y_path)
    normalize_y_gpp.columns = normalize_y_gpp.columns.astype(str)
    pdt.assert_almost_equal(normalize_y_code, normalize_y_gpp)


def test_drc_coefficients():
    raw = "C:/Users/RJ/Desktop/testdir/test20170726/Raw.xlsx"
    drc_coefficients_path = "C:/Users/RJ/Desktop/testdir/test20170726/drc_coefficients.xlsx"
    response = GPPCA.normalize_y_drc(GPPCA.add_sample_name(GPPCA.import_data(raw)))
    logx = GPPCA.log_dilution(GPPCA.get_concentrations(2000000, 2, 11, graph_type='drc'))
    drc_coefficients_code = GPPCA.drc_coefficients(response, logx)
    drc_coefficients_gpp = pd.read_excel(drc_coefficients_path)
    pdt.assert_frame_equal(drc_coefficients_code, drc_coefficients_gpp)

