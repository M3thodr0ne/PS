import warnings

warnings.filterwarnings('error')
try:
	warnings.warn(Warning())
except Warning:
	print 'Warning was raised as an exception!'
	
bo = true
print(bo)