from Bio.Seq import Seq
import operator
from collections import deque

codon_table = {'AAA':'K', 'AAC':'N', 'AAG':'K', 'AAU':'N', 'ACA':'T', 'ACC':'T',
               'ACG':'T', 'ACU':'T', 'AGA':'R', 'AGC':'S', 'AGG':'R', 'AGU':'S',
               'AUA':'I', 'AUC':'I', 'AUG':'M', 'AUU':'I', 'CAA':'Q', 'CAC':'H',
               'CAG':'Q', 'CAU':'H', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCU':'P',
               'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGU':'R', 'CUA':'L', 'CUC':'L',
               'CUG':'L', 'CUU':'L', 'GAA':'E', 'GAC':'D', 'GAG':'E', 'GAU':'D',
               'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCU':'A', 'GGA':'G', 'GGC':'G',
               'GGG':'G', 'GGU':'G', 'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V',
               'UAA':'*', 'UAC':'Y', 'UAG':'*', 'UAU':'Y', 'UCA':'S', 'UCC':'S',
               'UCG':'S', 'UCU':'S', 'UGA':'*', 'UGC':'C', 'UGG':'W', 'UGU':'C',
               'UUA':'L', 'UUC':'F', 'UUG':'L', 'UUU':'F'}

IntegerMass = {'G':57, 'A':71, 'S':87, 'P':97, 'V':99, 'T':101, 'C':103,
               'I':113, 'L':113, 'N':114, 'D':115, 'K':128, 'Q':128, 'E':129,
               'M':131, 'H':137, 'F':147, 'R':156, 'Y':163, 'W':186}

AminoAcid = ['G', 'A', 'S', 'P', 'V', 'T', 'C', 'I', 'L', 'N', 'D', 'K', 'Q', 'E', 'M', 'H', 'F', 'R', 'Y', 'W']
AminoAcidMass = [57, 71, 87, 97, 99, 101, 103, 113, 113, 114, 115, 128, 128, 129, 131, 137, 147, 156, 163, 186]

AminoAcidNoRedundant = ['G', 'A', 'S', 'P', 'V', 'T', 'C', 'I', 'N', 'D', 'K', 'E', 'M', 'H', 'F', 'R', 'Y', 'W']
MassNoRedundant = [57, 71, 87, 97, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186]

MassList = list(xrange(57,201))

def ProteinTranslation(Text):
    Protein  = ''
    for i in xrange(0, len(Text)/3):
        Protein = Protein + codon_table[Text[3*i:3*i+3]]
    return Protein

def RevCom(Text):
    seq = Seq(Text).reverse_complement()
    return str(seq)
    
def PeptideEncoding(Text, Peptide):
    DNA_length = len(Peptide) * 3
    DNA_array = []
    for i in xrange(0, len(Text) - DNA_length + 1):
        seq = Text[i:i+DNA_length]
        seq = seq.replace('T','U')
        revcom_seq = RevCom(Text[i:i+DNA_length])
        revcom_seq = revcom_seq.replace('T','U')
        if ProteinTranslation(seq) == Peptide or ProteinTranslation(revcom_seq) == Peptide:
            DNA_array.append(Text[i:i+DNA_length])
    DNA_array = list(set(DNA_array))  #remove duplicates
    return DNA_array

def SubpeptidesCount(n):
    return n * (n - 1)

def PeptideMass(Text):
    mass = 0
    for i in xrange(0, len(Text)):
        mass = mass + IntegerMass[Text[i]]
    return mass

def CyclospectrumOld(Text):
    SpectrumArray = [0]
    SpectrumArray.append(PeptideMass(Text))
    for i in xrange(0, len(Text)):
        for j in xrange(1,len(Text)):
            if i+j <= len(Text):
                subpeptide = Text[i:i+j]
            else:
                subpeptide = Text[i:] + Text[0:j-len(Text[i:])]
            SpectrumArray.append(PeptideMass(subpeptide))
    SpectrumArray.sort()
    return SpectrumArray

def LinearspectrumOld(Text):
    SpectrumArray = [0]
    SpectrumArray.append(PeptideMass(Text))
    for i in xrange(0, len(Text)):
        for j in xrange(1,len(Text)):
            if i+j <= len(Text):
                subpeptide = Text[i:i+j]
            else:
                continue
            SpectrumArray.append(PeptideMass(subpeptide))
    SpectrumArray.sort()
    return SpectrumArray

def Cyclospectrum(Peptide, AminoAcidRange):
    PrefixMass = {}
    PrefixMass[0] = 0
    Peptide_list = Peptide.split('-')
    for i in xrange(1, len(Peptide_list) + 1):
        for j in AminoAcidRange:
            if str(j) == Peptide_list[i-1]:
                PrefixMass[i] = PrefixMass[i - 1] + j
    CyclicSpectrum = [0]
    if Peptide_list[0]:
        peptideMass = PrefixMass[len(Peptide_list)]
        for i in xrange(0, len(Peptide_list)):
            for j in xrange(i + 1, len(Peptide_list) + 1):
                CyclicSpectrum.append(PrefixMass[j] - PrefixMass[i])
                if i > 0 and j < len(Peptide_list):
                    CyclicSpectrum.append(peptideMass - (PrefixMass[j] - PrefixMass[i]))
        CyclicSpectrum.sort()
    return CyclicSpectrum

def Linearspectrum(Peptide, AminoAcidRange):
    PrefixMass = {}
    PrefixMass[0] = 0
    Peptide_list = Peptide.split('-')
    for i in xrange(1, len(Peptide_list) + 1):
        for j in AminoAcidRange:
            if str(j) == Peptide_list[i-1]:
                PrefixMass[i] = PrefixMass[i - 1] + j
    LinearSpectrum = [0]
    for i in xrange(0, len(Peptide_list)):
        for j in xrange(i + 1, len(Peptide_list) + 1):
            LinearSpectrum.append(PrefixMass[j] - PrefixMass[i])
    LinearSpectrum.sort()
    return LinearSpectrum

def ExpandPeptide(Peptides, AminoAcidRange):
    NewPeptides = []
    for Peptide in Peptides:
        for amino_acid in AminoAcidRange:
            if Peptide:
                NewPeptides.append(Peptide + '-' + str(amino_acid))
            else:
                NewPeptides.append(str(amino_acid))
    return NewPeptides

def contained(candidate, container):
    temp = container[:]
    try:
        for v in candidate:
            temp.remove(v)
        return True
    except ValueError:
        return False
    
def CYCLOPEPTIDESEQUENCING(Spectrum, AminoAcidRange):
    Peptides = ['']
    candidates = []
    candiate_mass = []
    while len(Peptides) > 0:
        Peptides = ExpandPeptide(Peptides, AminoAcidRange)
        for Peptide in Peptides[:]:
            peptide_list = str(Peptide).split('-')
            mass_peptide = sum([int(x) for x in peptide_list])
            if mass_peptide == max(Spectrum):
                if Cyclospectrum(Peptide, AminoAcidRange) == Spectrum:
                    candidates.append(Peptide)
                Peptides.remove(Peptide)
            elif contained(Linearspectrum(Peptide, AminoAcidRange),Spectrum):
                continue
            else:
                Peptides.remove(Peptide)
    return candidates

#print CYCLOPEPTIDESEQUENCING([0, 71, 101, 113, 131, 184, 202, 214, 232, 285, 303, 315, 345, 416], MassNoRedundant)
print CYCLOPEPTIDESEQUENCING([0, 71, 99, 101, 103, 128, 129, 199, 200, 204, 227, 230, 231, 298, 303, 328, 330, 332, 333], MassNoRedundant)

def CyclopeptideScoring(Peptide, Spectrum, AminoAcidRange):
    Pep_spec = Cyclospectrum(Peptide, AminoAcidRange)
    score = 0
    for peptide_mass in Spectrum:
        if int(peptide_mass) in Pep_spec:
            score += 1
            Pep_spec.remove(int(peptide_mass))
    return score
print CyclopeptideScoring('131-71-131-71', [0, 71, 98, 99, 131, 202, 202, 202, 202, 202, 299, 333, 333, 333, 503], MassNoRedundant)

def LinearpeptideScoring(Peptide, Spectrum, AminoAcidRange):
    Pep_spec = Linearspectrum(Peptide, AminoAcidRange)
    score = 0
    for peptide_mass in Spectrum:
        if int(peptide_mass) in Pep_spec:
            score += 1
            Pep_spec.remove(int(peptide_mass))
    return score
print LinearpeptideScoring('97-129-129-97', [0, 97, 97, 129, 194, 196, 226, 226, 244, 258, 323, 323, 452], MassNoRedundant)


def Trim(Leaderboard, Spectrum, N, AminoAcidRange):
    pep_score = {}
    for Peptide in Leaderboard:
        pep_score[Peptide] = LinearpeptideScoring(Peptide, Spectrum, AminoAcidRange)
    sorted_list = sorted(pep_score.items(), key=operator.itemgetter(1), reverse=True)
    list_len = len(sorted_list)
    if list_len > N:
        for j in xrange(N, list_len):
            if sorted_list[j][1] < sorted_list[N-1][1]:
                sorted_list = sorted_list[:j]
                break
    list_len = len(sorted_list)
    neibour = [sorted_list[i][0] for i in xrange(0, list_len)]
    return neibour

def LEADERBOARDCYCLOPEPTIDESEQUENCING(Spectrum, N, AminoAcidRange):
    Leaderboard = [0]
    LeaderPeptide = ''
    pep_score = {}
    while len(Leaderboard) > 0:
        Leaderboard = ExpandPeptide(Leaderboard, AminoAcidRange)
        for Peptide in Leaderboard[:]:
            peptide_list = str(Peptide).split('-')
            mass_peptide = sum([int(x) for x in peptide_list])
            if mass_peptide == max(Spectrum):
                if CyclopeptideScoring(Peptide, Spectrum, AminoAcidRange) > CyclopeptideScoring(LeaderPeptide, Spectrum, AminoAcidRange):
                    LeaderPeptide = Peptide
                    pep_score[Peptide] = CyclopeptideScoring(Peptide, Spectrum, AminoAcidRange)
                elif CyclopeptideScoring(Peptide, Spectrum, AminoAcidRange) == CyclopeptideScoring(LeaderPeptide, Spectrum, AminoAcidRange):
                    pep_score[Peptide] = CyclopeptideScoring(Peptide, Spectrum, AminoAcidRange)
            elif mass_peptide > max(Spectrum):
                Leaderboard.remove(Peptide)
        Leaderboard = Trim(Leaderboard, Spectrum, N, AminoAcidRange)
    return pep_score

def SpectralConvolution(Spectrum):
    ConvoSpec = []
    for i in xrange(0, len(Spectrum) - 1):
        for j in xrange(i + 1, len(Spectrum)):
            if Spectrum[j] - Spectrum[i] == 0:
                continue
            ConvoSpec.append(Spectrum[j] - Spectrum[i])
    ConvoSpec.sort()
    return ConvoSpec
print SpectralConvolution([0, 57, 118, 179, 236, 240, 301])
def ConvolutionCyclopeptideSequencing(M, N, Spectrum):
    convolution_list = SpectralConvolution(Spectrum)
    freq_list = {}
    for mass in convolution_list:
        if mass >= 57 and mass <= 200:
            freq_list[mass] = 0
    for mass in convolution_list:
        if mass >= 57 and mass <= 200:
            freq_list[mass] += 1
    MFreqList = sorted(freq_list.items(), key=operator.itemgetter(1), reverse=True)
    for i in xrange(M, len(MFreqList)):
        if MFreqList[i][1] < MFreqList[M-1][1]:
            MFreqList = MFreqList[:i]
            break
    AminoAcidRange = [MFreqList[i][0] for i in xrange(0, len(MFreqList))]
    AminoAcidRange = list(set(AminoAcidRange))
    return LEADERBOARDCYCLOPEPTIDESEQUENCING(Spectrum, N, AminoAcidRange)
    
        
input_file = 'input.txt'
with open(input_file, 'r') as f:
    for line in f:
        mytext = line.rstrip('\n')
        mylist = mytext.split(' ')
        mylist = [ int(x) for x in mylist ]    # change to integer
        mylist.sort()
        print ConvolutionCyclopeptideSequencing(20, 1000, mylist)
        