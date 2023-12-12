import numpy as np

def Nussinov(seq):

    seq = seq.upper()

    bps = {"GU":1,"UG":1,
           "AU":2,"UA":2,
           "GC":3,"CG":3}


if __name__ == "__main__":

    seq = "ACGUACGCUAGCUGCUCGAUCGUCGAUCGAUCGACGCUAGCGCGUCGGGU"

    dbn = Nussinov(seq)

    print(seq)
    print(dbn)
