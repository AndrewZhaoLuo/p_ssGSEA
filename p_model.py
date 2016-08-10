
'''
David L Gibbs \\

Probabiltiy model for single sample GSA. \\

There are three hypothesis. \\

1.) The genes in the gene set are disregulated DOWN \\

2.) The genes in the gene set are not disregulated \\

3.) The genes in the gene set are disregulated UP. \\

We set a prior for each. To start, this can be uniform, or we \\
can use an informative prior taken from the GTEx project. \\

The prior should state: the probability that the genes \\
(in the set) are disregulated.  To do that, we will \\
only use the middle distribution. If they are evenly \\
spaced out in the GTEx data, then there *IS* a chance \\
that they could be disregulated.  If they are already \\
bunched up on one side, then the probability is low. \\
They are already far from "normal" in the middle. \\
What to do if they are bunched low in GTEx and high \\
in the observed measurements?? \\

I guess the prior needs to be HYPOTHESIS specific. \\
for H_low, the prior needs to check _where_ they are \\
located... if they are in the middle, or up high, then \\
they have the potential for DISREGULATION DOWN!! \\

If the genes are bunched low... they do not have the \\
potential for DISREGULATION LOW, but they _do_ have the \\
potential for DISREGULATION HIGH!! \\

Future work: \\

Main issue with current algorithm is the normalization I am doing, which guarantees values \\
fall between 0 and 1 to fit the current distributions. The issue arises since if I get a 0 \\
or 1; some extreme expression value especially in series (ie I evaluate several values close to 0) \\
this invariably leads to probabilities becoming very skewed. Evaluating several 1's will lead
to high's prior to be 1 and mid and low to 0 due to rounding. As a result, If I have a sequence
of expressions 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0. I would say the gene set is underexpressed
but because of those first few 1's, the probability of underexpressed is 0. I suggest
maybe creating distributions of mid, low, and high for every gene. \\

-Andrew
'''

from scipy.stats import beta

def pmodel(dat1, priors, mode):
    """
    Given normalized expression data, returns the probability that the genes are high, low, or mid expressed

    :param dat1: a list of gene expression values normalized to (0, 1).
    :type dat1: list

    :param priors: a 3 element list with the prior probability that the genes are under, mid, or overexpressed. The first element should be the proability of low expression, second for mid expression, and third for last expression.
    :type priors: list

    :param mode: a string either high, low, or mid representing the probability that should be returned in the end
    :type mode: float

    :return: a probability measuring certainty the expression in the given set are expressed according to the pattern given by mode
    """

    hLow = beta(0.5,2)
    hMid = beta(2,2)
    hHigh = beta(2,0.5)

    p1 = priors[0]
    p2 = priors[1]
    p3 = priors[2]

    for di in dat1:
        if di == 0:
            di = 0.00001
        if di == 1:
            di = 0.99999
        #print(p1,'\t', p2,'\t', p3)
        p1 = p1 * (hLow.sf(di))
        p2 = p2 * (min(hMid.sf(di), hMid.cdf(di)) * 2.0)
        p3 = p3 * (hHigh.cdf(di))
        
        ptot = p1 + p2 + p3
        if ptot == 0:
            print(p1, '\t', p2, '\t', p3)
            print(dat1)
            exit()
        p1 /= ptot
        p2 /= ptot
        p3 /= ptot

    if mode == 'high':
        return(p3)
    elif mode == 'low':
        return(p1)
    return(p2)