
#########################################
#
# David L Gibbs
# draft of probabiltiy model for single sample GSA.
#
# There are three hypothesis.
# 1.) The genes in the gene set are disregulated DOWN
# 2.) The genes in the gene set are not disregulated
# 3.) The genes in the gene set are disregulated UP.

# We set a prior for each. To start, this can be uniform, or we
# can use an informative prior taken from the GTEx project.

# The prior should state: the probability that the genes
# (in the set) are disregulated.  To do that, we will
# only use the middle distribution. If they are evenly
# spaced out in the GTEx data, then there *IS* a chance
# that they could be disregulated.  If they are already
# bunched up on one side, then the probability is low.
# They are already far from "normal" in the middle.
# What to do if they are bunched low in GTEx and high
# in the observed measurements??

# I guess the prior needs to be HYPOTHESIS specific.
# for H_low, the prior needs to check _where_ they are
# located... if they are in the middle, or up high, then
# they have the potential for DISREGULATION DOWN!!

# If the genes are bunched low... they do not have the
# potential for DISREGULATION LOW, but they _do_ have the
# potential for DISREGULATION HIGH!!

# Beta function parameters will probably need to be adjusted.

#

from scipy.stats import beta
import numpy as np

hLow = beta(0.5,2)
hMid = beta(2,2)
hHigh = beta(2,0.5)

# getting probabilties #
# hLow.sf(0.4)
# 1 - hMid.cdf(0.4) - hMid.sf(0.6)
# hHigh.cdf(0.4)

########################################

def pmodel(dat1, priors, mode):

    from scipy.stats import beta
    import numpy as np

    hLow = beta(0.5,2)
    hMid = beta(2,2)
    hHigh = beta(2,0.5)

    p1 = priors[0] # h1
    p2 = priors[1] # h2
    p3 = priors[2] # h3

    for di in dat1:
        p1 = p1 * hLow.sf(di)
        p2 = p2 * min(hMid.sf(di), hMid.cdf(di)) * 2.0
        p3 = p3 * hHigh.cdf(di)

    ptot = p1+p2+p3
    p1 /= ptot
    p2 /= ptot
    p3 /= ptot
    #print(dat1)
    print("low  " + str(p1) + "mid  " + str(p2) + "high " + str(p3))

    if mode == 'high':
        return(p3)
    elif mode == 'low':
        return(p1)
    else:
        return(p2)
'''
########################################
dat2 = [-0.3, 0.4, 0.44, 0.33]
p1 = 0.333 # h1
p2 = 0.333 # h2
p3 = 0.333 # h3

print(pmodel(dat2, [p1, p2, p3], 'mid'))

for di in dat2:
    p1 = p1 * hLow.sf(di)
    p2 = p2 * min(hMid.sf(di), hMid.cdf(di)) * 2.0
    p3 = p3 * hHigh.cdf(di)

ptot = p1+p2+p3
p1 /= ptot
p2 /= ptot
p3 /= ptot
print("low  " + str(p1) + "\nmid  " + str(p2) + "\nhigh " + str(p3))
print()

########################################
dat3 = [0.88, 0.89, 0.9]

p1 = 0.333 # h1
p2 = 0.333 # h2
p3 = 0.333 # h3

for di in dat3:
    p1 = p1 * hLow.sf(di)
    p2 = p2 * min(hMid.sf(di), hMid.cdf(di)) * 2.0
    p3 = p3 * hHigh.cdf(di)

ptot = p1+p2+p3
p1 /= ptot
p2 /= ptot
p3 /= ptot
print("low  " + str(p1) + "\nmid  " + str(p2) + "\nhigh " + str(p3))
print()

########################################
dat4 = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

p1 = 0.333 # h1
p2 = 0.333 # h2
p3 = 0.333 # h3

for di in dat4:
    p1 = p1 * hLow.sf(di)
    p2 = p2 * min(hMid.sf(di), hMid.cdf(di)) * 2.0
    p3 = p3 * hHigh.cdf(di)

ptot = p1+p2+p3
p1 /= ptot
p2 /= ptot
p3 /= ptot
print("low  " + str(p1) + "\nmid  " + str(p2) + "\nhigh " + str(p3))
print()

########################################

dat4 = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

p1 = 0.333 # h1
p2 = 0.333 # h2
p3 = 0.333 # h3

for di in dat4:
    p1 = p1 * hLow.sf(di)
    p2 = p2 * min(hMid.sf(di), hMid.cdf(di)) * 2.0
    p3 = p3 * hHigh.cdf(di)

ptot = p1+p2+p3
p1 /= ptot
p2 /= ptot
p3 /= ptot
print("low  " + str(p1) + "\nmid  " + str(p2) + "\nhigh " + str(p3))
print()

########################################
# even a strong prior, in the face of convincing data,
# cannot move the posterior.

dat4 = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

p1 = 0.05 # h1 -- low
p2 = 0.05 # h2 -- middle
p3 = 0.9 # h3 -- high

for di in dat4:
    p1 = p1 * hLow.sf(di)
    p2 = p2 * min(hMid.sf(di), hMid.cdf(di)) * 2.0
    p3 = p3 * hHigh.cdf(di)

ptot = p1+p2+p3
p1 /= ptot
p2 /= ptot
p3 /= ptot
print("low  " + str(p1) + "\nmid  " + str(p2) + "\nhigh " + str(p3))
print()
#############
'''
