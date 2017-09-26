import StringIO

import numpy as np

def parse_data_str(s, shape):
    return np.loadtxt(StringIO.StringIO(s)).reshape(*shape).astype(int)

DATA1 = dict(
    mask = parse_data_str("""
    0 0 0 0
    0 0 0 0
    0 0 0 0
    0 0 0 0
    0 0 0 0

    0 0 0 0
    0 0 0 0
    0 1 0 0
    0 0 0 0
    0 0 0 0

    0 0 0 0
    0 0 0 0
    0 0 0 0
    0 0 0 0
    0 0 0 0
    """, shape=(3, 5, 4)),
    N3=1,
    N2=6,
    N1=6,
)

DATA2 = dict(
    mask = parse_data_str("""
    0 0 0 0
    0 0 0 0
    0 0 0 0
    0 0 0 0
    0 0 0 0

    0 0 0 0
    0 0 0 0
    0 1 1 0
    0 0 0 0
    0 0 0 0

    0 0 0 0
    0 0 0 0
    0 0 0 0
    0 0 0 0
    0 0 0 0
    """, shape=(3, 5, 4)),
    N3=2,
    N2=5*2,
    N1=12,
)

DATA4 = dict(
    mask = parse_data_str("""
    0 0 0 0
    0 0 0 0
    0 0 0 0
    0 0 0 0
    0 0 0 0

    0 0 0 0
    0 1 1 0
    0 1 1 0
    0 0 0 0
    0 0 0 0

    0 0 0 0
    0 0 0 0
    0 0 0 0
    0 0 0 0
    0 0 0 0
    """, shape=(3, 5, 4)),
    N3=4,
    N2=2*4+2*4,
    N1=2*4+4*2,
)

DATA3 = dict(
    mask = parse_data_str("""
    0 0 0 0
    0 0 0 0
    0 0 0 0
    0 0 0 0
    0 0 0 0

    0 0 0 0
    0 1 0 0
    0 1 1 0
    0 1 0 0
    0 0 0 0

    0 0 0 0
    0 0 0 0
    0 0 0 0
    0 0 0 0
    0 0 0 0
    """, shape=(3, 5, 4)),
    N3=4,
    N2=2*4 + 3*3 + 1
)
