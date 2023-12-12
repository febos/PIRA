import numpy as np

def Nussinov(seq, minloop = 3):

    seq = seq.upper()
    N   = len(seq)

    bps = {"GU":1,"UG":1,
           "AU":2,"UA":2,
           "GC":3,"CG":3}

    D = np.zeros((N,N))
    K = {}

    for h in range(N):
        for i in range(N-h):

            j = i+h

            if j - i - 1 < minloop:

                D[i,j] = 0

            else:

                bestk, bestscorek = -1, -1

                for k in range(i, j-minloop):

                    if seq[k] + seq[j] in bps:
                        scorek = D[i, k - 1] + D[k + 1, j - 1] + bps[seq[k] + seq[j]]
                        if scorek > bestscorek:
                            bestk, bestscorek = k, scorek

                if bestscorek >= D[i, j - 1]:

                    K[(i, j)] = bestk
                    D[i, j]   = bestscorek

                else:
                    D[i, j] = D[i, j - 1]

    print(D)
    print(K)
                            

if __name__ == "__main__":

    #seq = "ACGUACGCUAGCUGCUCGAUCGUCGAUCGAUCGACGCUAGCGCGUCGGGU"
    seq = "GGUCCAC"

    dbn = Nussinov(seq)

    print(seq)
    print(dbn)
