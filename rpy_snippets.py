'''
This file is meant to test that the installation of rpy went without a hitch

Some things I test
Plotting
'''

#built in objects in R
import rpy2.robjects as robjects

#the main interpreter
from rpy2.robjects import r as R

'''
************
* PLOTTING *
************
'''
#can aliases for R functions!
rplot = R("plot")

#If you don't make aliases... kind of disgustting
R("jpeg(\"test.jpeg\")")

x = [1, 2, 3, 4]
y = [1, 2, 4, 8]
rplot(x, y, type = "b", main = "Test", xlab = "x", ylab = "y")

'''
*************
* CONSTANTS *
*************
'''
#Yes you can print them!
print(R('pi'))

'''
***********
* IMPORTS *
***********
'''
#You can also import packages from R and use them!