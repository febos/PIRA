import numpy as np

### https://doi.org/10.1073/pnas.77.11.6309

# base pair scores
bps = {"GU":-1,"UG":-1,
       "AU":-2,"UA":-2,
       "GC":-3,"CG":-3}


def BackTrack(begin, end, K, minloop, partial = False):

    queue = {(begin, end), } # start from the last cell
    basepairs = []

    while queue:

        newq = set()

        # for each cell in the queue
        for i,j in queue:
            # if we have (k, j) base pair
            if (i,j) in K:
                k = K[(i,j)]
                # if we may have more base pairs before k  
                if (k - 1) - i > minloop:
                    newq.add((i, k - 1))
                # if we may have more base pairs inside [k, j]
                if (j - 1) - (k + 1) > minloop and not partial:
                    newq.add((k + 1, j - 1))
                # add (k, j) bp
                basepairs.append((k,j))
            # if j is unpaired
            else:
                # if we may have more base pairs inside [i, j - 1]
                if (j - 1) - i > minloop:
                    newq.add((i, j - 1))
        # update queue
        queue = newq

    return sorted(basepairs)


def Ekj(seq,k,j,inners):

    return bps[seq[k]+seq[j]]


def Nussinov(seq, minloop = 3):

    seq = seq.upper()
    N   = len(seq)    # sequence length

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

            bestk, bestscorek = -1, 10**9
            # searching for the optimal base pair
            for k in range(i, j - minloop):

                if seq[k] + seq[j] in bps:
                    inners = BackTrack(k + 1, j - 1, K, minloop, partial = True)
                    scorek = D[i, k - 1] + D[k + 1, j - 1] + Ekj(seq,k,j,inners)
                    if scorek < bestscorek:
                        bestk, bestscorek = k, scorek
            # see if no base pair is not better
            if bestscorek <= D[i, j - 1]:

                K[(i, j)] = bestk
                D[i, j]   = bestscorek
            # if it is
            else:
                D[i, j] = D[i, j - 1]

    # traceback
    dbn = ['.' for ch in seq] # dbn template
    predictedbps = BackTrack(0, N - 1, K, minloop)
    for v,w in predictedbps:
        dbn[v] = '('
        dbn[w] = ')'
    
    return D[0, N - 1], ''.join(dbn)
                            

if __name__ == "__main__":

    seq = 'GGGCGGCUAGCUCAGCGGAAGAGCGCUCGCCUCACACGCGAGAGGUCGUAGGUUCAAGUCCUACGCCGCCCACCA'
    dbn = '(((((((..((((....[..)))).(((((.......))))).....(((((..]....))))))))))))....'
    
    #seq = "GGUUCCAC"

    score, pred = Nussinov(seq)

    print(seq)
    print(dbn)
    print(pred)
    print(score)
