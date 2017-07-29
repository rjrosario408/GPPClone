import GPPClone_AssayCalculations_Class as GPPCA
import pandas as pd
import pandas.util.testing as pdt


def test_calculate_cv_20170720():
    raw = GPPCA.Inhibition("C:/Users/RJ/Desktop/testdir/test20170720/Raw.xlsx")
    calculate_cv_path = "C:/Users/RJ/Desktop/testdir/test20170720/calculate_cv.xlsx"
    calculate_cv_code = raw.calculate_cv()
    calculate_cv_gpp = pd.read_excel(calculate_cv_path)
    calculate_cv_gpp.columns = calculate_cv_gpp.columns.astype(str)
    pdt.assert_almost_equal(calculate_cv_code, calculate_cv_gpp)


def test_log_dilution_20170720():
    log_dilution_path = "C:/Users/RJ/Desktop/testdir/test20170720/log_dilution.xlsx"
    log_dilution_code = GPPCA.Inhibition.log_dilution(GPPCA.Inhibition.get_concentrations(1000, 2, 10, 'inhibition'))
    log_dilution_gpp = pd.read_excel(log_dilution_path)
    pdt.assert_almost_equal(log_dilution_code, log_dilution_gpp)


def test_calculate_y_transform_20170720():
    raw = GPPCA.Inhibition("C:/Users/RJ/Desktop/testdir/test20170720/Raw.xlsx")
    transform_y_path = "C:/Users/RJ/Desktop/testdir/test20170720/y_transform.xlsx"
    transform_y_code = raw.transform_y()
    transform_y_gpp = pd.read_excel(transform_y_path)
    transform_y_gpp.columns = transform_y_gpp.columns.astype(str)
    pdt.assert_almost_equal(transform_y_code, transform_y_gpp)


def test_inhibition_coefficients():
    raw = GPPCA.Inhibition("C:/Users/RJ/Desktop/testdir/test20170720/Raw.xlsx")
    inhibition_coefficients_path = "C:/Users/RJ/Desktop/testdir/test20170720/inhibition_coefficients.xlsx"
    inhibition_coefficients_code = raw.coefficients()
    inhibition_coefficients_gpp = pd.read_excel(inhibition_coefficients_path)
    pdt.assert_frame_equal(inhibition_coefficients_code, inhibition_coefficients_gpp)


def test_normalize_y_drc_20170726():
    raw = GPPCA.DoseResponse("C:/Users/RJ/Desktop/testdir/test20170726/Raw.xlsx")
    transform_y_path = "C:/Users/RJ/Desktop/testdir/test20170726/normalize_y_drc.xlsx"
    transform_y_code = raw.transform_y()
    transform_y_gpp = pd.read_excel(transform_y_path)
    transform_y_gpp.columns = transform_y_gpp.columns.astype(str)
    pdt.assert_almost_equal(transform_y_code, transform_y_gpp)


def test_drc_coefficients():
    raw = GPPCA.DoseResponse("C:/Users/RJ/Desktop/testdir/test20170720/Raw.xlsx")
    inhibition_coefficients_path = "C:/Users/RJ/Desktop/testdir/test20170726/drc_coefficients.xlsx"
    inhibition_coefficients_code = raw.coefficients()
    inhibition_coefficients_gpp = pd.read_excel(inhibition_coefficients_path)
    pdt.assert_frame_equal(inhibition_coefficients_code, inhibition_coefficients_gpp)

