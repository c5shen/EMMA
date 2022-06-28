"""
Created on Nov 14, 2012

@author: Siavash Mirarab
"""

# 1.19.2022 - Copied from SEPP/UPP by Chengze


def gcd(a, b):
    """Return greatest common divisor using Euclid's Algorithm."""
    while b:
        a, b = b, a % b
    return a


def lcm(a, b):
    """Return lowest common multiple."""
    return a * b // gcd(a, b)
