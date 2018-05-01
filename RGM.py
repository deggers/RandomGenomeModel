from random import randint
import numpy as np
import re
import seaborn as sns

global prom_len, cod_len, genome_len, lower_int, upper_int, lower_prom_len, upper_prom_len
prom_len_fixed = 5
cod_len = 7
genome_len = 1e5

lower_int, upper_int = 0, 3
lower_prom_len, upper_prom_len = 5, 8


def genRandomGenome():
    genome = ""
    for _ in range(int(genome_len)):
        genome += str(randint(lower_int, upper_int))
    return genome


def genRandProm_VariableLength():
    prom_seq = ""
    length = randint(lower_prom_len, upper_prom_len)
    for _ in range(length):
        prom_seq += str(randint(lower_int, upper_int))
    return prom_seq


def genRandProm_FixedLength():
    prom_seq = ""
    for _ in range(prom_len_fixed):
        prom_seq += str(randint(lower_int, upper_int))
    return prom_seq


def Lreg_VariablePromLen(genome, fixedProm_Seq=False):
    regulatoryRegion = []
    prom_seq = genRandProm_VariableLength()
    if fixedProm_Seq:
        regulatoryRegion = re.compile(prom_seq + "\d{" + str(cod_len) + "}(.*?)" + prom_seq).split(genome)
    else:
        while len(genome.split(prom_seq, 1)) == 2:
            prom_seq = genRandProm_VariableLength()
            splitResult = genome.split(prom_seq, 1)
            if len(splitResult) == 2:
                genome = splitResult[1]
                regulatoryRegion.append(splitResult[0])
    return regulatoryRegion


tmp = []
for _ in range(100):
    genome = genRandomGenome()
    regulatory_seq = Lreg_VariablePromLen(genome)
    tmp.append(len(regulatory_seq))
ax = sns.boxplot(x=tmp)
ax.set_title('#Lreg with variable [PromSeq')

for _ in range(100):
    genome = genRandomGenome()
    regulatory_seq = Lreg_VariablePromLen(genome,fixedProm_Seq=True)
    tmp.append(len(regulatory_seq))
ax = sns.boxplot(x=tmp)
ax.set_title('#Regulatory Regions with fixed Promoter Sequence and length (5)')
sns.plotting_context()

