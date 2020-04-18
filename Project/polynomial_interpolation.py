import numpy as np

# Returns Lagrande coef. for polynome centered at x = 0, evaluated at x = 0 + eval
def lagrange_coef(coefs, eval=.5):
    array = np.empty(len(coefs))
    for i, coef in enumerate(coefs):
        other = coefs[:i] + coefs[i+1:]
        other = np.array(other)
        num = np.product(eval - other)
        den = np.product(coef - other)
        array[i] = num / den
    return array

if __name__ == '__main__':
    print("Coeff. 3rd order polynomial at left wall:\n\t",  lagrange_coef([0, +1/2, +3/2, +5/2],        eval=-.5))
    print("Coeff. 3rd order polynomial at right wall:\n\t", lagrange_coef([0, -1/2, -3/2, -5/2],        eval=+.5))
    print("Coeff. 4th order polynomial at left wall:\n\t",  lagrange_coef([0, +1/2, +3/2, +5/2, +7/2],  eval=-.5))
    print("Coeff. 4th order polynomial at right wall:\n\t", lagrange_coef([0, -1/2, -3/2, -5/2, -7/2],  eval=+.5))