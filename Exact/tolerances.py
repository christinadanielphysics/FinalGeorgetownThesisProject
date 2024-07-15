from numpy import finfo

TOLERANCE_computational_zero_for_denominator_approximate = 1e-5 # for the filter_lesser and the filter_greater functions in the Denominator.py file
# high accuracy with 1e-7 for 4-site (this tolerance is for excluding small spectral function weights)
TOLERANCE_combine = 1e-5 # binning in the frequency domain, basically, in the combine.py file and in the filter_lesser and filter_greater functions in the Denominator.py file

MAXITER_brent = 10000
TOLERANCE_computational_zero_for_denominator = 4*finfo(float).eps
TOLERANCE_brent_for_polynomial = 4*finfo(float).eps
TOLERANCE_brent_for_derivative = 4*finfo(float).eps


