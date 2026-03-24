from math import dist
from typing import NamedTuple
import copy, sys
import pickle
import numpy as np
import matplotlib.pyplot as plt

from python.inequalities import (
    extract_inequality_coefficients,
    get_ineqsign,
    Inequality,
)
import python.version as version
from python.propagation_data import (
    LWEInstance,
    PropagationData,
    PropagationDataStep,
)
from python.helpers import (
    IneqType,
    print_v
)
from python.run import (
    create_graph_inequalities,
    propagate
)
from python.solve import solve
from recover import process_propagation_data
from dsa_dist import get_dist_z

beta = 78
gamma1 = 2**17
q = 2**23 - 2**13 + 1
NO_X = True

class RejectedSample(NamedTuple):
    coeff: int
    coeff_idx: int
    poly_idx: int
    c: list[int]
    x: list[list[int]]

def format_poly(poly) -> list[int]:
    return [int(x) for x in poly.split(",")[:-1]]

def rot(poly: list[int]) -> list[int]:
    tmp = poly[255]
    poly[1:256] = poly[0:255]
    poly[0] = -tmp
    return poly

def xtimes(poly: list[int], times: int) -> list[int]:
    t_poly = copy.deepcopy(poly)
    for _ in range(times):
        t_poly = rot(t_poly)
    return t_poly

def inner_prod(poly1: list[int], poly2: list[int]) -> int:
    return sum(x * y for x, y in zip(poly1, poly2))

def poly_mul(poly1: list[int], poly2: list[int]) -> list[int]:
    res = [0] * 256
    for i in range(256):
        res[i] = inner_prod(poly1[::-1], xtimes(poly2, i))
    return res[::-1]

def get_rejected_sample(f) -> RejectedSample:
    coeff = int(f.readline().split(":")[1])
    coeff_idx = int(f.readline().split(":")[1])
    poly_idx = int(f.readline().split(":")[1])
    c = format_poly(f.readline().split(":")[1])
    x = None
    if not NO_X:
        f.readline()
        x = [format_poly(f.readline()) for _ in range(4)]
    return RejectedSample(coeff, coeff_idx, poly_idx, c, x)

def test_ineq(sample: RejectedSample, s1: list[list[int]]) -> bool:
    # for i in range(4):
    #     x = poly_mul(s1[i], sample.c)
    #     if sample.x[i] != x:
    #         print(f"Not match {i}")
    
    poly_idx = sample.poly_idx
    x = poly_mul(s1[poly_idx], sample.c)
    lb = -beta
    ub = beta

    if sample.coeff > 0:
        lb = sample.coeff - gamma1
    else:
        ub = sample.coeff + gamma1 - 1

    print(f"{lb} <= {x[sample.coeff_idx]} <= {ub}")
    xi = inner_prod(s1[sample.poly_idx][::-1], xtimes(sample.c, 255-sample.coeff_idx))
    if xi != x[sample.coeff_idx]:
        print(f"Inner prod {xi} does not match poly mul {x[sample.coeff_idx]}")

    return lb <= x[sample.coeff_idx] <= ub


def main():
    STEPS = 20
    STEP_SIZE = 1
    USE_BEST_STEP = True
    SCA_OBS = True
    PERFECT_INEQ = False

    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <n_rej|--pck>")
        sys.exit(1)

    arg = sys.argv[1]
    try:
        n_rej = int(arg)
        is_numeric_arg = True
    except ValueError:
        is_numeric_arg = False

    attacked_s1_idx = 0
    if is_numeric_arg:
        with open("./testdata1M_noX.dat") as f:
            s1 = [format_poly(f.readline()) for _ in range(4)]
            sk = f.readline()
            rejected_samples = [get_rejected_sample(f) for _ in range(n_rej)]

        inequalities = []
        list_obs = []
        dist_rej_z = get_dist_z()[-2*beta:-1]

        for sample in rejected_samples:
            # if not test_ineq(sample, s1):
            #     print("Inequality not satisfied for sample:", sample)

            if attacked_s1_idx == sample.poly_idx:
                lb = -beta
                ub = beta
                ci = xtimes(sample.c, 255-sample.coeff_idx)

                if SCA_OBS:
                    observed = (((256*q + sample.coeff) % (256*q))  - (gamma1 - beta)) % 256
                    if PERFECT_INEQ:
                        start = 60
                        stop = 80
                        if (start <= observed <= stop) and sample.coeff > 0:
                            lb = observed - beta
                            inequalities.append(Inequality(ci, IneqType.GE, lb, True, 1))
                            #inequalities.append(Inequality(ci, IneqType.LE, beta, True, 1))
                        elif (156 - start >= observed >= 156 - stop) and sample.coeff < 0:
                            ub = observed - beta - 1
                            inequalities.append(Inequality(ci, IneqType.LE, ub, True, 1))
                            #inequalities.append(Inequality(ci, IneqType.GE, -beta, True, 1))

                        else:
                            continue
                    else:
                        if (58 <= observed <= 98):
                            prob_pos = dist_rej_z[observed + 15]
                            prob_neg = dist_rej_z[156 - observed - 15]
                            denom = prob_pos + prob_neg
                            prob_pos = prob_pos / denom
                            prob_neg = prob_neg / denom
                            inequalities.append(Inequality(ci, IneqType.GE, observed - beta, "UNK", prob_pos))
                                #inequalities.append(Inequality(ci, IneqType.LE, observed - beta - 1, "UNK", prob_neg))
                        else:
                            continue

                    # if sample.coeff < 0:
                    #     list_obs.append(observed)
                else:
                    if sample.coeff > 0:
                        lb = sample.coeff - gamma1
                        inequalities.append(Inequality(ci, IneqType.GE, lb, True, 1))
                    else:
                        ub = sample.coeff + gamma1 - 1
                        inequalities.append(Inequality(ci, IneqType.LE, ub, True, 1))

                
        
        # plt.hist(list_obs, bins=100)
        # plt.savefig("obs_hist.png")
        
        propagation_data = PropagationData.new(
            s1[attacked_s1_idx][::-1], inequalities, 0, 0, 0
        )

        with open("bp_state.pkl", "wb") as f:
            pickle.dump({"s1": s1,"propagation_data": propagation_data,"inequalities": inequalities},f)


    elif arg == "--pkl":
        with open("bp_state.pkl", "rb") as f:
            state = pickle.load(f)
            s1 = state["s1"]
            propagation_data = state["propagation_data"]
            inequalities = state["inequalities"]


    ### run_with_inequality
    key_priori_dist = {i: np.float64(1/5) for i in range(-2, 3)}
    # key_priori_dist = {-2: np.float64(0.3), -1: np.float64(0.15), 0: np.float64(0.1), 1: np.float64(0.15), 2: np.float64(0.3)}
    g = create_graph_inequalities(
        inequalities,
        key_priori_dist
    )
        
    success_bp = propagate(
        s1[attacked_s1_idx][::-1],
        g,
        STEPS,
        STEP_SIZE,
        propagation_data=propagation_data,
    )
    print("BP alone " + "succeeded" if success_bp else "did not succeed" + ".")
    step_idx = (
        max(
            propagation_data.steps,
            key=lambda x: propagation_data.steps[x].correct_coefficients,
        )
        if USE_BEST_STEP
        else -1
    )
    step_list = sorted(propagation_data.steps.items(), key=lambda x: x[0])
    pos = 0
    for idx, _ in reversed(step_list):
        if idx == step_idx:
            break
        pos += 1
    print_v(
        f"Using step {step_idx} which is on position {pos}, has {propagation_data.steps[step_idx].correct_coefficients} correct coefficients and {propagation_data.steps[step_idx].recovered_coefficients} recovered coefficients."
    )
    propagation_data.recovered_coefficients = propagation_data.steps[
        step_idx
    ].recovered_coefficients
    # success = solve(
    #     propagation_data,
    #     block_size=block_size,
    #     run_reduction=run_reduction,
    #     max_beta=max_beta,
    #     perform=run_reduction,
    #     add_fplll=add_fplll,
    #     step=step_idx,
    #     max_enum=max_enum,
    #     step_rank=pos,
    # )
    process_propagation_data(
        propagation_data,
        False,
        True,
        False,
        False,
        "pdf",
        None,
        50,
        "archive_type",
        "data",
        None,
    )
if __name__ == "__main__":
    main()
