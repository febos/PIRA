import numpy as np

### https://doi.org/10.1073/pnas.77.11.6309

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
    dbn = ['.' for ch in seq] # dbn template

    queue = {(0, N - 1), } # start from the last cell

    while queue:

        newq = set()

        # for each cell in the queue
        for i,j in queue:
            # if we have (k, j) base pair
            if (i,j) in K:
                k = K[(i,j)]
                # if we have more base pairs before k  
                if D[(i, k - 1)] > 0:
                    newq.add((i, k - 1))
                # if we have more base pairs inside [k, j]
                if D[(k + 1, j - 1)] > 0:
                    newq.add((k + 1, j - 1))
                # add (k, j) bp to dbn
                dbn[k] = '('
                dbn[j] = ')'
            # if j is unpaired
            else:
                # if we have more base pairs inside [i, j - 1]
                if D[(i, j - 1)] > 0:
                    newq.add((i, j - 1))
        # update queue
        queue = newq
    
    return ''.join(dbn)
                            

if __name__ == "__main__":

    seq = "ACGUACGCUAGCUGCUCGAUCGUCGAUCGAUCGACGCUAGCGCGUCGGGU"
    #seq = "GGUUCCAC"

    dbn = Nussinov(seq)

    print(seq)
    print(dbn)
