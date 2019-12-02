

class DnaParser:
    def __init__(self, filename=""):

        self.x = Dna();
        self.y = Dna();
        self.filename = filename

    def parse(self):

        if self.filename != "":

            with open(self.filename, "r") as file:
                n = int(file.readline())
                m = int(file.readline())

                x = file.readline()
                x = x.strip().split(' ')

                y = file.readline()
                y = y.strip().split(' ')


                self.x = Dna(x)
                self.y = Dna(y)

        return self.x, self.y



class Dna:
    '''
    cette class permet de faire des operation sur une sequence d'adn
    '''
    cins = 2
    cdel = 2
    csub_concord = 3
    csub_nonconcord = 4

    def __init__(self, seq=[]):

        self.seq = list(seq)
        self.l = len(seq)

    def __str__(self):
        return "|".join(self.seq)

    @staticmethod
    def alignement(a1, a2):
        '''
        a1: tuple (xb, yb) d'adn
        a2: tuple (x, y) d'adn


        return true | false , si a1 est un alignement de a2
        '''
        xbar, ybar = a1
        x, y = a2

        if not  (Dna.pi(xbar) == x.seq and  Dna.pi(ybar) == y.seq and len(xbar.seq) == len(ybar.seq)):
            #print("pi xbar = ", Dna.pi(xbar))
            return False

        for i in range(len(xbar.seq)):
            if xbar.seq[i] == ybar.seq[i] == "_":
                return False
        return True

    @staticmethod
    def pi(xbar):
        return list(filter(lambda e: e != "_", xbar.seq))
        #return xbar.seq.replace("_", "")
    @staticmethod
    def c(xbar, ybar):
        coup = 0
        for i in range(len(xbar.seq)):
            coup += Dna.c_atomique(xbar.seq[i], ybar.seq[i])
        return coup

    @staticmethod
    def c_atomique(xbar_i, ybar_i):
        if xbar_i == "_":
            return Dna.cins
        if ybar_i == "_":
            return Dna.cdel

        if xbar_i == ybar_i:
            return 0
        if ((xbar_i, ybar_i) == ("A", "T")) or \
           ((xbar_i, ybar_i) == ("T", "A")) or \
           ((xbar_i, ybar_i) == ("G", "C")) or \
           ((xbar_i, ybar_i) == ("C", "G")):
           return Dna.csub_concord

        return Dna.csub_nonconcord

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
    #print("i=",i, "j=",j)
    if i == len(x.seq) and j == len(y.seq):
        if c < dist:
            dist = c
            print(dist)

    else:
        if i<len(x.seq) and j<len(y.seq):
            dist = DIST_NAIF_REC(x, y, i+1, j+1, c+Dna.c_atomique(x.seq[i], y.seq[j]), dist)
            #print(dist)
        if i<len(x.seq):
            dist = DIST_NAIF_REC(x, y, i+1, j, c+Dna.cins, dist)
            #print(dist)
        if j<len(y.seq):
            dist = DIST_NAIF_REC(x, y, i, j+1, c+Dna.cdel, dist)
            #print(dist)
    return dist

def printD(x, y, D):

    print(" ", y)
    for i in range(len(x.seq)):
        print(x.seq[i], end="|")
        for j in range(len(y.seq)):
            print(D[i][j], end="|")
        print()

def SOL_1(x, y, D):
    xbar_seq = []
    ybar_seq = []

    i = len(x.seq)-1
    j = len(y.seq)-1
    currentD = D[i][j]
    while i > 0 and j > 0:
        #print("i=", i, ", j=", j)
        if(D[i-1][j] + Dna.cins == currentD):
            currentD = D[i-1][j]
            ybar_seq.append("_")
            xbar_seq.append(x.seq[i])
            i = i-1


            #insetion


        elif D[i][j-1] + Dna.cdel == currentD:
            currentD = D[i][j-1]
            #delete
            xbar_seq.append("_")
            ybar_seq.append(y.seq[j])

            j = j-1


        elif D[i-1][j-1] + Dna.c_atomique(x.seq[i], y.seq[j]) == currentD:
            currentD = D[i-1][j-1]

            xbar_seq.append(x.seq[i])
            ybar_seq.append(y.seq[j])

            i = i-1
            j = j-1

        # case 0
    xbar_seq.append(x.seq[i])
    ybar_seq.append(y.seq[j])


    return Dna(list(reversed(xbar_seq))), Dna(list(reversed(ybar_seq)))
def DIST_1(x, y):
    """
    x, y adn

    return tableau D contenant toute les valeur de D(i, j)
    """
    D = [[Dna.cdel*j for j in range(len(y.seq))]]

    #contruction de la matrice
    for i in range(1, len(x.seq)):
        D.append([Dna.cins*i]+[None]*(len(y.seq)-1))




    #remplissage
    for j in range(1, len(y.seq)):
        for i in range(1, len(x.seq)):

            D[i][j] = min(D[i-1][j] + Dna.cins,\
                          D[i][j-1] + Dna.cdel,\
                          D[i-1][j-1] + Dna.c_atomique(x.seq[i], y.seq[j]))

    return D

if __name__ == "__main__":

    parser = DnaParser("../data/Inst_0000010_44.adn")



    x, y = parser.parse()
    print(x)
    print(y)
    import time
    t1 = time.time()
    d = DIST_1(x, y)
    printD(x, y, d)
    xbar, ybar = SOL_1(x, y, d)
    print(xbar)
    print(ybar)
    t2 = time.time() - t1
    print("Execution time = ", t2)
