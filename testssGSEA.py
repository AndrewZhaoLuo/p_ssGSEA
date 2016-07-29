from cache_codec import *

dump_ssGSEA_scores("BC")

scores = load_ssGSEA_scores("BC")["FARMER_BREAST_CANCER_CLUSTER_6"]

#v = scores.values()
#div = abs(max(v) - min(v))

#print(div)
print(max(scores.values()))
print(min(scores.values()))
print(scores[38])
