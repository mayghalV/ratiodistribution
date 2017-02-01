# ratiodistribution
Statisical implementations for the ratio distribution including error bounds estimation

##Motivation
A regession has predicted a value for variable X of mu_x with error bounds (often at 95% confidence) e_x.  
Similarly, another independent regression has predicted a value for variable Y  of mu_y with error e_y.

ie:

X = mu_x +- e_x
Y = mu_y +- e_y

##Questions
*If we define Z as the ratio of X and Y, ie Z:=X/Y, what is the distribution of this ratio?
*What are the error bounds of Z with a% confidence, ie e_z? (0%<a<100%)

##Credit and more information
This package implements the formulas from Hinkley, D. V. (December 1969). "On the Ratio of Two Correlated Normal Random Variables". Biometrika. 56 (3): 635–639.
Also see https://en.wikipedia.org/wiki/Ratio_distribution#Gaussian_ratio_distribution for more information.

##Methods:
###ratiodistribution.sdFromError(error, confidence)
Returns standard deviation of the errors given the error value and its confidence.
####Parameters:
error: e_x from above (float)
confidence: % confidence of error bound (float)
####Returns:
Standard deviation of the distribution of errors (float)


###ratiodistribution.pdf(z, mu_x, sd_x, mu_y, sd_y)
Probability density function for ratio distribution
####Parameters:
z: Value of z (float)  
mu_x: predicted value of X (float)  
sd_x: standard deviation of the distribution of errors for X (float)  
mu_y: predicted value of Y (float)  
sd_y: standard deviation of the distribution of errors for Y (float)  
####Returns:
Probability density function (float)

###ratiodistribution.cdf(z, mu_x, sd_x, mu_y, sd_y)  
Cumulative distribution function for ratio distribution
####Parameters:
z: Value of z (float)  
mu_x: predicted value of X (float)  
sd_x: standard deviation of the distribution of errors for X (float)  
mu_y: predicted value of Y (float)  
sd_y: standard deviation of the distribution of errors for Y (float)  
####Returns:
Cumulative distribution function (float)

###ratiodistribution.invCdf(cdf, mu_x, sd_x, mu_y, sd_y)
Inverse cumulative distribution function for ratio distribution
####Parameters:
cdf: Cumulative distribution (float)  
mu_x: predicted value of X (float)  
sd_x: standard deviation of the distribution of errors for X (float)  
mu_y: predicted value of Y (float)  
sd_y: standard deviation of the distribution of errors for Y (float)  
####Returns:
z such that P(Z<=z)=cdf (float)  

###ratiodistribution.error(mu_x, sd_x, mu_y, sd_y, confidence=0.95)
Calculates the error bounds of Z with confidence *confidence*. Note given the nature of the ratio distribution the size of the errors might not be symmetric around the mean.
####Parameters:
mu_x: predicted value of X (float)  
e_x: error in estimation of X at confidence level *confidence*(float)  
mu_y: predicted value of Y (float)  
e_y: error in estimation of Y at confidence level *confidence* (float)  
confidence: confidence interval of e_x, e_y and new bounds  
####Returns:
[lowerValue, upperValue] at confidence level *confidence* (list of 2 floats)  


###ratiodistribution.plotPdf(mu_x, e_x, mu_y, e_y, confidence=0.95)  
Plots shape of pdf of Z
####Parameters:
mu_x: predicted value of X (float)  
e_x: error in estimation of X at confidence level *confidence*(float)  
mu_y: predicted value of Y (float)  
e_y: error in estimation of Y at confidence level *confidence* (float)  
confidence: confidence interval of e_x, e_y  
####Returns:
None