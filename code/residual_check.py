"""Numerical check of the tree ODE and of the expansion of P(tau_o <= 1/2).

h' = lam (1+h)(1+h-t)/(1-2t+h), h(0)=0, on [0,1/2];  H = h(1/2).
P(root frozen by 1/2) = 1 - exp(lam)/(1+H)^3.
R = 4H/lam - log(1/H) - 1 is the error in H = (lam/4)[log(1/H)+1+o(1)];
K = P/(3 lam/4) - (L - log L) should tend to log 4 - 1/3,
and K - (log 4 - 1/3) = R - log(4H/(lam L)) + O(lam L^2).
"""
import json
from mpmath import mp, mpf, odefun, log, exp

mp.dps = 40
target = log(4) - mpf(1) / 3
out = []
for e in [3, 4, 6, 8, 10, 12]:
    lam = mpf(10) ** (-e)
    f = odefun(lambda t, h: lam * (1 + h) * (1 + h - t) / (1 - 2 * t + h), 0, mpf(0))
    H = f(mpf(1) / 2)
    P = 1 - exp(lam) / (1 + H) ** 3
    L = log(1 / lam)
    R = 4 * H / lam - log(1 / H) - 1
    K = P / (3 * lam / 4) - (L - log(L))
    g = log(4 * H / (lam * L))
    row = {"lam": f"1e-{e}", "H": float(H), "P": float(P), "K": float(K),
           "K_minus_target": float(K - target), "R": float(R), "log_4H_over_lamL": float(g),
           "R_minus_log": float(R - g), "L_times_R": float(L * R), "target": float(target)}
    print(row, flush=True)
    out.append(row)
json.dump(out, open(__file__.replace(".py", "_output.json"), "w"), indent=1)
