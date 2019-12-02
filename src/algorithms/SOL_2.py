def DIST_NAIF(x, y):
    return DIST_NAIF_REC(x, y, 0, 0, 0, float('inf'))

def DIST_NAIF_REC(x, y, i, j, c, dist):
    """
    x: Dna
    y: Dna
    i: indice dans [0..|x|]
    j: indice dans [0..|y|]
    dist: coup du meilleur alignement connu avant cet appel

    return
    dist meilleur coup connu apre cet appel
    """

    if i == len(x.seq) and j == len(x.seq):
        if c < dist:
            dist = c
    else:
        if i<len(x.seq) and j<len(y.seq):
            dist = DIST_NAIF_REC(x, y, i+1, j+1, c+Dna.c_atomique(x.seq[i+1], y.seq[j+1]), dist)
        if i<len(x.seq):
            dist = DIST_NAIF_REC(x, y, i+1, j, c+Dna.c_ins, dist)
        if j<len(y.seq):
            dist = DIST_NAIF_REC(x, y, i, j+1, c+Dna.c_del, dist)

    return dist
