import numpy as np
import math


def diff(alpha, derivative=1):
    n = len(alpha)

    b = np.zeros(n)
    b[derivative] = math.factorial(derivative)
    A = np.power.outer(alpha, np.arange(n)).T
    beta = np.linalg.solve(A, b)

    alpha_n = A[n - 1, :] * alpha  # alpha**n
    coef = np.inner(alpha_n, beta)
    EG = np.dot(A, beta) - b

    eG = sum(abs(EG))

    eM = np.inner(abs(alpha_n), abs(beta))

    if abs(coef) <= (eG * eM):
        coef = np.inner(alpha_n * alpha, beta)
        factor = math.factorial(n + 1) / coef
        order = n - 1
    else:
        factor = math.factorial(n) / coef
        order = n - 2

    return beta, factor, order
