import numpy as np
from math import log,sqrt,exp,isnan

#####################Utils####################

def remove_none(array):
	"""Removes the None values.

	Identifies the None values in an array and returns the array without them.
	Only 1-dimensional arrays are supported

	Parameters
	----------
	array: array-like
		the array from which to remove the None values
	"""
	mask=np.array([a is None for a in array])
	return np.array(array)[mask==False]

def is_boolean(s):
	""" Checks whether string s is equal to 0 or 1"""
	#TODO: also chekc for True and False and checks that doesn't break anything in upload
	return s =="0" or s == "1"
	
def is_float(s):
	""" Checks whether string s is a number"""
	try:
		float(s)
		return True
	except ValueError:
		return False
	
def moving_average(sequence, n_rec=1):

	""" Smoothes a sequence.

	Smoothes a sequence using the moving average method.
	New x(t) is the average of x(t-1), x(t), and x(t+1).
	x(0) is 2*x(0)+x(1)

	Parameters
	----------
	sequence : array-like
		the sequence to smooth.
	n_rec : int
		how many iterations of the moving average should be applied to the sequence.
	"""

	if(len(sequence)<=1):
		return sequence
	
	if n_rec<=0:
		print("warning: moving average: n_rec should be positive")
		return sequence
		
	seq=np.zeros(np.shape(sequence))
	seq[0]=sequence[0]
	seq[-1]=sequence[-1]
	seq[0:-1] = seq[0:-1] + sequence[1:]
	seq = seq + sequence
	seq[1:] = seq[1:] + sequence[:-1]
	seq=seq/3
	
	if n_rec==1:
		return seq
	else:
		return moving_average(seq, n_rec=n_rec-1)

# interpolate a model of the time serie with an exponential
# returns the exponential of the sampled model
def exp_interpolation(timeseries, get_model=False):
	"""Exponential interpolation of a time-series
	
	Fit an exponential curve to a time-series.
	For a time-series with values v and time t,
	fits a curve of the form v = a x exp( b x t)
	by minimizing the mean-squared error.

	Formulaes for a and b were taken from:
	http://mathworld.wolfram.com/LeastSquaresFittingExponential.html

	Parameters
	----------
	timeseries : array-like, array-like
		The 2 arrays should have the same size.
		The first array contains the timestamps.
		The second array contains the values of the series.

	get_model : boolean, default: False
		Whether to return the coefficients a and b.

	Returns
	-------
	time,series(, a, b) : array-like, array-like (, float, float)
		Returns a tuple of 2 or 4 elements.
		The first is an array with the timestamps (as input).
		The second is an array with the corresponding interpolated values.
		If get_model is True, the third and forth elements are the a and b coefficients.
	
	"""
	time,series=timeseries
	logseries = [log(s) for s in series]
	
	# old way:
	# model = linear_model.LinearRegression(fit_intercept=True)
	# model.fit([[t] for t in time],logseries)
	# loginterpolation=model.predict([[t] for t in time])
	# interpolation = [exp(l) for l in loginterpolation]
	
	a = ( np.sum([x*x*y for x,y in zip(time,series)])*np.sum([y*logy for y,logy in zip(series,logseries)]) - np.sum([x*y for x,y in zip(time,series)])*np.sum([x*y*logy for x,y,logy in zip(time,series,logseries)]) ) / ( np.sum([y for y in series])*np.sum([x*x*y for x,y in zip(time,series)]) - pow(np.sum([x*y for x,y in zip(time,series)]),2) )
		
	b = ( np.sum([y for y in series])*np.sum([x*y*logy for x,y,logy in zip(time,series,logseries)]) - np.sum([x*y for x,y in zip(time,series)])*np.sum([y*logy for y,logy in zip(series,logseries)]) ) / ( np.sum([y for y in series])*np.sum([x*x*y for x,y in zip(time,series)]) - pow(np.sum([x*y for x,y in zip(time,series)]),2) )
	
	A=exp(a)
	B=b
	
	interpolation=[A*exp(B*t) for t in time]
	if get_model:
		return time,interpolation,A,B
	else:
		return time,interpolation







