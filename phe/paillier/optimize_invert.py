def gcd_extended(a, b):
    if a == 0:
        return 0, 1, b
    x1, y1, g = gcd_extended(b%a, a)

    x = y1 - (b//a) * x1;
    y = x1;

    return x, y, g

def extended_euclidian_inverse(a, b):
    x, y, g = gcd_extended(a, b)
    if g != 1:
        raise ValueError('%d has no inverse mod %d' % (a, b))
    else:
        d = (x % b + b) % b;
    return d

def invert(a, b):
    """
    The multiplicitive inverse of a in the integers modulo b.

    :return int: x, where a * x == 1 mod b
    """
    return extended_euclidian_inverse(a, b)
