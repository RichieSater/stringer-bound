"""Finite Bernstein-cover witnesses for the unique-top ``L=4`` family.

The verifier derives every polynomial and checks every local coefficient in
exact rational arithmetic.  The data below specify only prefix-free dyadic
partitions and, where needed, nonnegative integer multipliers.
"""

Z_LE_A_FINAL_LEAVES = (
    ('0000', 'T', None),
    ('00010', 'T', None),
    ('00011', 'M4', (43931, 1)),
    ('0010', 'T', None),
    ('00110', 'T', None),
    ('00111', 'M4', (490678, 1)),
    ('01000', 'M4', (14783, 1)),
    ('0100100', 'M4', (27645, 1)),
    ('0100101', 'M4', (5053830, 1)),
    ('0100110', 'T', None),
    ('0100111', 'M4', (64274491, 1)),
    ('01010', 'M4', (5692665, 1)),
    ('01011', 'G4', None),
    ('01100', 'T', None),
    ('0110100', 'T', None),
    ('0110101', 'G4', None),
    ('011011', 'T', None),
    ('0111', 'G4', None),
    ('10000', 'M4', (2389, 1)),
    ('10001', 'G4', None),
    ('10010', 'M4', (9509, 1)),
    ('10011', 'G4', None),
    ('10100', 'M4', (3637, 1)),
    ('10101', 'G4', None),
    ('10110', 'M4', (15232, 1)),
    ('10111', 'G4', None),
    ('11', 'G4', None),
)

B_GE_V_FINAL_LEAVES = (
    ('00', 'T', None),
    ('0100', 'G4', None),
    ('01010', 'G4', None),
    ('01011', 'T', None),
    ('01100', 'G4', None),
    ('01101', 'T', None),
    ('01110', 'G4', None),
    ('01111', 'T', None),
    ('1', 'T', None),
)

MIDDLE_LOW_Z_LEAVES = (
    ('00', 'T', None),
    ('010', 'G4', None),
    ('0110000', 'T', None),
    ('0110001', 'G4', None),
    ('011001', 'T', None),
    ('0110100', 'T', None),
    ('0110101', 'G4', None),
    ('011011', 'T', None),
    ('0111', 'G4', None),
    ('1', 'T', None),
)

LOWER_LOW_Z_LEAVES = (
    ('000', 'T', None),
    ('001000', 'G4', None),
    ('0010010', 'T', None),
    ('0010011', 'G4', None),
    ('00101', 'G4', None),
    ('00110', 'T', None),
    ('00111000', 'G4', None),
    ('00111001', 'T', None),
    ('0011101', 'G4', None),
    ('001111', 'T', None),
    ('010', 'T', None),
    ('01100', 'T', None),
    ('011010', 'G4', None),
    ('011011', 'T', None),
    ('0111', 'T', None),
    ('1000', 'G4', None),
    ('100100', 'G4', None),
    ('1001010', 'T', None),
    ('1001011', 'G4', None),
    ('10011', 'G4', None),
    ('101', 'G4', None),
    ('1100', 'G4', None),
    ('1101', 'T', None),
    ('1110', 'G4', None),
    ('11110', 'T', None),
    ('111110', 'G4', None),
    ('111111', 'T', None),
)

LOWER_MID_CORNER_LEAVES = (
    ('00', 'T', None),
    ('01', 'G4', None),
    ('1', 'T', None),
)
