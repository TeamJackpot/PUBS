import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def reversecomp(some_seq):
    out = ''
    pairs = {'A':'T','T':'A','C':'G','G':'C'}
    for i in some_seq:
        out = out + pairs[i]
    return(out[::-1])

def getbarcodes(fastq):
    count = 0
    barcodes = []
    print('Starting')
    t = min(len(fastq),10000000)
    for i in range(t):
        if (i % 4 == 1):
            if 'N' not in fastq[i][0:18]:
                barcodes.append(reversecomp(fastq[i][0:18]))
    return(barcodes)

alleles = pickle.load(open("allele_dic_with_WT.pkl",'rb'))
amino = pickle.load(open("aminotonumber.pkl",'rb'))
translate = pickle.load(open("translate.pkl",'rb'))

def make_grid(path):

    df = pd.DataFrame()

    f = open(path, 'r')
    fastq = f.readlines()
    df['Barcode'] = getbarcodes(fastq)
    f.close()

    df[["Positions","Codons"]] = df['Barcode'].map(alleles).apply(pd.Series)
    df["Codons"] = df['Codons'].str.replace('T','U')
    df["Amino Acid"] = df["Codons"].map(translate)
    df["Amino Acid Numeric"] = df["Amino Acid"].map(amino)

    grid = pd.crosstab(df["Amino Acid Numeric"],df["Positions"])
    wt_count = df['Codons'].value_counts()["WU"]

    return(grid, wt_count)

ts1 = make_grid('OutputCGTGAT.fastq') #CobaltT0R1
ts2 = make_grid('OutputACATCG.fastq') #CobaltT1R1
ts3 = make_grid('OutputGCCTAA.fastq') #CobaltT2R1
ts4 = make_grid('OutputTGGTCA.fastq') #CobaltT0R2
ts5 = make_grid('OutputCACTGT.fastq') #CobaltT1R2
ts6 = make_grid('OutputATTGGC.fastq') #CobaltT2R2
ts7 = make_grid('OutputATTCCG.fastq') #ControlT0R1
ts8 = make_grid('OutputAGCTAG.fastq') #ControlT1R1
ts9 = make_grid('OutputGTATAG.fastq') #ControlT2R1

grids = {
        'CobaltT0R1':ts1,
        'CobaltT1R1':ts2,
        'CobaltT2R1':ts3,
        'CobaltT0R2':ts4,
        'CobaltT1R2':ts5,
        'CobaltT2R2':ts6,
        'ControlT0R1':ts7,
        'ControlT1R1':ts8,
        'ControlT2R1':ts9,
        }

pickle.dump(grids,open( "/Users/student/Desktop/grids.p", "wb" ))
