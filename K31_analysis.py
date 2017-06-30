import sys
sys.path.append("/scicore/home/neher/neher/HIV/hivwholeseq")
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
plt.ion()

from hivwholeseq.patients.patients import load_patients, Patient

def dt(d1,d2):
	return (d1.toordinal()-d2.toordinal())/365.25


if __name__=="__main__":
	pats = load_patients(csv=True)
	fmt = "%d/%m/%Y"
	res = []
	cutoff = 0.01
	region = 'gag'
	for pcode, pat in pats.iterrows():
		try:
			EDI = datetime.strptime(pat["infect date best"], fmt)
			P = Patient(pat)
			aft = P.get_allele_frequency_trajectories(region)[0]
			for si, (scode, sample) in enumerate(P.samples.iterrows()):
				try:
					date = datetime.strptime(sample["date"], fmt)
					af = aft[si]
					mask = (1.0 - af.max(axis=0))>cutoff
					div = ((1.0 - (af**2).sum(axis=0))*mask).mean()
					print(EDI, date, div)
					res.append((dt(date, EDI), div))
				except:
					print(scode, "didn't work")

		except:
			print("skipping patient ", pcode)

	res = np.array(res)
	plt.plot(res[:,0], res[:,1], 'o')
