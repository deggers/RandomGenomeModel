from random import randint
import numpy as np 
import re
import seaborn as sns

global prom_len, cod_len, genome_len, lower_int, upper_int
prom_len = 5
cod_len = 7 
genome_len = 1e5

lower_int = 0
upper_int = 3 


def genRandomGenome():
    genom = ""
    for _ in range(int(genome_len)):
        genom += str(randint(lower_int,upper_int))
    return genom

def genRandProm_VariableLength():
    prom_seq = ""
    length = randint(lower_prom_len, upper_prom_len)

    for _ in range(length):
        prom_seq += str(randint(lower_int, upper_int))
    return prom_seq

def genRandProm_FixedLength(prom_len):
    prom_seq = ""
    for _ in range(prom_len):
        prom_seq += str(randint(lower_int, upper_int))
    return prom_seq

def Lreg_VariablePromLen(genome, fixedProm_Seq=False):
    regulatoryRegion = []
    prom_seq = genRandProm_VariableLength()
    if fixedProm_Seq:
        regulatoryRegion = re.compile(prom_seq + "\d{"+ str(cod_len) +"}(.*?)" + prom_seq).split(genome)
    else: 
        while len(genome.split(prom_seq,1)) == 2:
            prom_seq = genRandProm_VariableLength()
            splitResult = genome.split(prom_seq,1)
            if len(splitResult) == 2:
                genome = splitResult[1]
                regulatoryRegion.append(splitResult[0])
    return regulatoryRegion
