import random
from math import gcd

def is_coprime_with_list(num, num_list):
    """
    Check if `num` is coprime with all numbers in `num_list`.
    """
    return all(gcd(num, other) == 1 for other in num_list)

def generate_pairwise_coprime(count, start=2, end=50):
    """
    Generate a list of `count` pairwise coprime integers.
    
    Parameters:
        count (int): Number of coprime integers to generate.
        start (int): Minimum value for the integers.
        end (int): Maximum value for the integers.
    
    Returns:
        list: A list of pairwise coprime integers.
    """
    coprime_list = []
    
    while len(coprime_list) < count:
        num = random.randint(start, end)
        if is_coprime_with_list(num, coprime_list):
            coprime_list.append(num)
    
    return coprime_list