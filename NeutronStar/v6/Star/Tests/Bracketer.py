import numpy as np 

def func(val, u):
    return(u**(3.0/2)-val)

def D_func(u):
    return (3.0/2.0)*u**(1.0/2.0)


def bracketer(f, init_val):

    ulast = init_val
    i = 0

    while True:

        i = i+1
        unext = ulast - func(f,ulast)/D_func(ulast)
        if abs(unext-ulast ) <0.00001:
            return unext
        else:
            ulast = unext

if __name__ == "__main__":
    val = 1000
    result = bracketer(val, 2)
    print(result)
