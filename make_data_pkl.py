from typing import NamedTuple
import copy, sys
import pickle
from my_bp import RejectedSample, format_poly, get_rejected_sample

with open("./testdata10M.dat") as f:
    s1 = [format_poly(f.readline()) for _ in range(4)]
    sk = f.readline()

    rejected_samples = [[],[],[],[]]
    for _ in range(999999999):
        try:
            sample = get_rejected_sample(f)
            rejected_samples[sample.poly_idx].append(sample)
        except IndexError:
            break

with open("rejected_samples_0.pkl", "wb") as f:
    pickle.dump({"rejected_samples": rejected_samples[0]},f)
with open("rejected_samples_1.pkl", "wb") as f:
    pickle.dump({"rejected_samples": rejected_samples[1]},f)
with open("rejected_samples_2.pkl", "wb") as f:
    pickle.dump({"rejected_samples": rejected_samples[2]},f)
with open("rejected_samples_3.pkl", "wb") as f:
    pickle.dump({"rejected_samples": rejected_samples[3]},f)

print(f'{len(rejected_samples[0])} rejected samples for poly_idx 0')
print(f'{len(rejected_samples[1])} rejected samples for poly_idx 1')
print(f'{len(rejected_samples[2])} rejected samples for poly_idx 2')
print(f'{len(rejected_samples[3])} rejected samples for poly_idx 3')
