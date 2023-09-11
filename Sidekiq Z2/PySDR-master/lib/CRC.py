"""
Cyclic Redundancy Check functions

"""

import numpy as np

# Returns XOR of 'a' and 'b'
# (both of same length)
def xor(a, b):
 
    # initialize result
    result = np.array([])
 
    # Traverse all bits, if bits are
    # same, then XOR is 0, else 1
    for i in range(1, len(b)):
        if a[i] == b[i]:
            result = np.append(result, [0])
        else:
            result = np.append(result, [1])
 
    return result
 
 
# Performs Modulo-2 division
def mod2div(dividend, divisor):
 
    # Number of bits to be XORed at a time.
    pick = np.size(divisor)
 
    # Slicing the dividend to appropriate
    # length for particular step
    tmp = dividend[0: pick]
 
    while pick < np.size(dividend):
 
        if tmp[0] == 1:
 
            # replace the dividend by the result
            # of XOR and pull 1 bit down
            tmp = np.append(xor(divisor, tmp), dividend[pick])
 
        else:   # If leftmost bit is '0'
            # If the leftmost bit of the dividend (or the
            # part used in each step) is 0, the step cannot
            # use the regular divisor; we need to use an
            # all-0s divisor.
            tmp = np.append(xor(np.ones(pick), tmp), dividend[pick])

        # increment pick to move further
        pick += 1
 
    # For the last n bits, we have to carry it out
    # normally as increased value of pick will cause
    # Index Out of Bounds.
    if tmp[0] == 1:
        tmp = xor(divisor, tmp)
    else:
        tmp = xor(np.zeros(pick), tmp)
 
    checkword = tmp
    return checkword
 
# Function used at the sender side to encode
# data by appending remainder of modular division
# at the end of data.
 
 
def encodeData(data, key):
 
    l_key = np.size(key)
 
    # Appends n-1 zeroes at end of data
    appended_data = np.append(data, np.zeros(l_key-1))
    remainder = mod2div(appended_data, key)
 
    # Append remainder in the original data
    codeword = np.append(data, remainder)
    return codeword
 
def CRCcheck(codeword, key):

    checkword = mod2div(codeword, key)
    
    if np.all(checkword == 0):
        error = 0
    else:
        error = 1

    return error


# Driver code
data = np.array([1,0,0,1,0,0])
key = np.array([1,1,0,1])

codeword = encodeData(data, key)
error = CRCcheck(codeword, key)
# print(error)
