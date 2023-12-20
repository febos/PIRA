import numpy as np

def Nussinov(seq, minloop = 3):

    seq = seq.upper()
    N   = len(seq)    # sequence length

    # base pair scores
    bps = {"GU":1,"UG":1,
           "AU":2,"UA":2,
           "GC":3,"CG":3}

    # score matrix
    D = np.zeros((N,N))
    # base pair dictionary for traceback
    K = {}

    # going diagonale by diagonale,
    # skipping first [minloop] diagonals
    # to avoid too short hairpins
    for h in range(minloop + 1, N):
        # row index
        for i in range(N - h):
            # column index
            j = i + h

            bestk, bestscorek = -1, -1
            # searching for the optimal base pair
            for k in range(i, j - minloop):

                if seq[k] + seq[j] in bps:
                    scorek = D[i, k - 1] + D[k + 1, j - 1] + bps[seq[k] + seq[j]]
                    if scorek > bestscorek:
                        bestk, bestscorek = k, scorek
            # see if no base pair is not better
            if bestscorek >= D[i, j - 1]:

                K[(i, j)] = bestk
                D[i, j]   = bestscorek
            # if it is
            else:
                D[i, j] = D[i, j - 1]

    # traceback
    print(D)
    print(K)
                            

if __name__ == "__main__":

    #seq = "ACGUACGCUAGCUGCUCGAUCGUCGAUCGAUCGACGCUAGCGCGUCGGGU"
    seq = "GGUCCAC"

    dbn = Nussinov(seq)

    print(seq)
    print(dbn)
