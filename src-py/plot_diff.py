import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import seaborn as sns

with open("/mnt/scratch2/avi/meta-map/kraken/counts/hi.txt") as f:
    hi = pd.read_table(f, header=None).set_index(0)

with open("/mnt/scratch2/avi/meta-map/kraken/counts/lo.txt") as f:
    lo = pd.read_table(f, header=None).set_index(0)

with open("/mnt/scratch2/avi/meta-map/kraken/counts/tr.txt") as f:
    tr = pd.read_table(f, header=None, sep=" ").set_index(0)

with open("/mnt/scratch2/avi/meta-map/kraken/counts/krk.txt") as f:
    kr = pd.read_table(f, header=None).set_index(0)

ct = pd.concat([tr, lo, hi, kr], axis=1)
ct.columns = ["truth", 'MM', 'Unique', "kraken"]

print (ct.corr(method="spearman"))

mm_ards = np.abs(ct["truth"] - ct["MM"]) / ct["truth"]
kr_ards = np.abs(ct["truth"] - ct["kraken"]) / ct["truth"]
un_ards = np.abs(ct["truth"] - ct["Unique"]) / ct["truth"]

#sns.distplot(mm_ards,  hist=False, label="MM")
#sns.distplot(un_ards,  hist=False, label="Unique")
#sns.distplot(kr_ards,  hist=False, label="Unique")

sns.regplot(x="truth", y="MM", data=ct, label="MM")
sns.regplot(x="truth", y="Unique", data=ct, label="Unique", scatter_kws={'alpha':0.1})
sns.regplot(x="truth", y="kraken", data=ct, label="Kraken", scatter_kws={'alpha':0.3})
sns.regplot(x="truth", y="truth", data=ct, label="truth", scatter_kws={'alpha':0.3})
plt.legend()
plt.show()
plt.savefig("un.pdf")


print "MARD: MM"
print np.mean( np.abs(ct["truth"] - ct["MM"]) / ct["truth"] )
print "MARD: Unique"
print np.mean( np.abs(ct["truth"] - ct["Unique"]) / ct["truth"] )
print "MARD: Kraken"
print np.mean( np.abs(ct["truth"] - ct["kraken"]) / ct["truth"] )
