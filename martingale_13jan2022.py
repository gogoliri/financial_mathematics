import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt

m = 100000
dt = 1 / 365
T = 1
t = 0
n = 365  # number of day
S0 = 1
K = 1.1
r = 0.01
sigma = 0.2

ST = np.zeros(m)
STs = np.zeros(m)

for i in range(m):
    eps = np.random.normal(0, 1)
    ST[i] = S0 * np.exp((r - 0.5 * sigma ** 2) * T + sigma * eps * np.sqrt(T))
    STs[i] = S0 * np.exp((r + sigma ** 2 - 0.5 * sigma ** 2) * T + sigma * eps * np.sqrt(T))


def blscall(St, K, r, sigma, T, t):
    stau = np.sqrt(T - t)
    d1 = (np.log(St / K) + (r + 0.5 * (sigma ** 2)) * (T - t)) / (sigma * stau)
    d2 = d1 - sigma * stau
    bsc = St * norm.cdf(d1) - K * np.exp(-r * (T - t)) * norm.cdf(d2)
    return bsc


Q = np.mean(ST > K)
Qs = np.mean(STs > K)

C_martingale = S0 * Qs - K * np.exp(-r * T) * Q
C_bls = blscall(S0, K, r, sigma, T, t)
print(f"Call option pricing from Black-Scholes: {C_bls}")
print(f"Call option pricing from Martingale: {C_martingale}")

# ex2

C_MC = np.exp(-r * T) * np.mean(np.maximum(ST-K,0))
C_MC_s = S0 * np.mean((np.maximum(STs-K, 0) / STs))


print(f"Call option pricing using risk neutral: {C_MC}")
print(f"Call option pricing using alternative choice numeraire: {C_MC_s}")
