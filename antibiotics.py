RNACodon = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L","UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S", "UAU":"Y", "UAC":"Y", "UAA":"*", "UAG":"*", "UGU":"C", "UGC":"C", "UGA":"*", "UGG":"W",
    "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L","CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P","CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q","CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
    "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
    "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G"}
def ProteinTranslationProblem(RNA):
    """
Input: An RNA string Pattern and the array GeneticCode.
Output: The translation of Pattern into an amino acid string Peptide.
    """
    keys = RNACodon.keys()
    peptide = ""
    for i in range(0,len(RNA)-2, 3):
        amino= RNA[i:i+3]
        if RNACodon[amino] == "*":
            return peptide
        peptide += RNACodon[amino]
    return peptide
def Reverse(Pattern):
    """
     Input:  A string Pattern
     Output: The reverse of Pattern
    """
    revPattern =  ""
    for i in range(len(Pattern)-1 , -1, -1):
        revPattern += Pattern[i]
    return revPattern

def Complement(Pattern):
    """
    Input:  A DNA string Pattern
    Output: The complementary string of Pattern (with every nucleotide replaced
    by its complement).
    """
    comPattern = ""
    for i in range(len(Pattern)):
        if Pattern[i] == "T":
            comPattern += "A"
        elif Pattern[i] == "A":
            comPattern += "T"
        elif Pattern[i] == "C":
            comPattern += "G"
        else:
            comPattern += "C"
    return comPattern

def ReverseComplement(Pattern):
    """
    Input: A DNA string Pattern.
    Output: The reverse complement of Pattern
    """
    Pattern = Reverse(Pattern) # reverse all letters in a string
    Pattern = Complement(Pattern) # complement each letter in a string
    return Pattern


def PeptideEncodingProblem(DNA, AA, RNACodon):
    """
    Input: A DNA string Text, an amino acid string Peptide, and the array GeneticCode.
    Output: All substrings of Text encoding Peptide (if any such substrings exist).
    """
    seq = []
    DNA2= ReverseComplement(DNA)
    for i in range(len(DNA)-len(AA)*3+1):
        sub = DNA[i:i+len(AA)*3]
        RNA = sub.replace("T", "U")
        if ProteinTranslationProblem(RNA) == AA:
            seq.append(sub)
    for i in range(len(DNA2)-len(AA)*3+1):
        sub = DNA2[i:i+len(AA)*3]
        RNA = sub.replace("T", "U")
        if ProteinTranslationProblem(RNA) == AA:
            seq.append(ReverseComplement(sub))
    return seq
def PeptideEncodingProblem2(DNA, AA, RNACodon):
    """
    Input: A DNA string Text, an amino acid string Peptide, and the array GeneticCode.
    Output: All substrings of Text encoding Peptide (if any such substrings exist).
    """
    seq = []
    for i in range(len(DNA)-len(AA)*3+1):
        sub = DNA[i:i+len(AA)*3]
        rev = ReverseComplement(sub)
        RNA = sub.replace("T", "U")
        RNA2 = rev.replace("T", "U")
        if ProteinTranslationProblem(RNA) == AA:
            seq.append(sub)
        if ProteinTranslationProblem(RNA2) == AA:
            seq.append(sub)
    return seq
mass = {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114, 'D': 115, 'K': 128, 'Q': 128, 'E': 129, 'M': 131, 'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186}
mass2  = {57: ['G'], 71: ['A'], 87: ['S'], 97: ['P'], 99: ['V'], 101: ['T'], 103: ['C'], 113: ['I', 'L'], 114: ['N'], 115: ['D'], 128: ['K', 'Q'], 129: ['E'], 131: ['M'], 137: ['H'], 147: ['F'], 156: ['R'], 163: ['Y'], 186: ['W']}
aminoAcid = ['G', 'A', 'S', 'P', 'V', 'T', 'C',"I", 'L', 'N', 'D',"K", 'Q', 'E', 'M', 'H', 'F', 'R', 'Y', 'W']
mass4  = {57: 'G', 71: 'A', 87: 'S', 97: 'P', 99: 'V', 101: 'T', 103: 'C', 113: 'L', 114: 'N', 115: 'D', 128: 'Q', 129: 'E', 131: 'M', 137: 'H', 147: 'F', 156: 'R', 163: 'Y', 186: 'W'}

def LinearSpectrum(Peptide, Alphabet, AminoAcidMass):
    """
    Input: An amino acid string Peptide.
    Output: The linear spectrum of Peptide.
    """
    PrefixMass = {}
    PrefixMass[0]= 0
    for i in range(len(Peptide)):
        for s in Alphabet:
            if s == Peptide[i] :
                PrefixMass[i+1] = PrefixMass[i] + AminoAcidMass[s]
    LinearSpectrum = [0]
    for i in range(len(PrefixMass)):
        for j in range(i+1, len(PrefixMass)):
            a  =PrefixMass[j] - PrefixMass[i]
            LinearSpectrum.append(a)
    LinearSpectrum.sort()
    return LinearSpectrum


def CyclicSpectrum(Peptide, Alphabet, AminoAcidMass):
    """
    Input: An amino acid string Peptide.
    Output: The cyclic spectrum of Peptide.
    """
    PrefixMass = {}
    PrefixMass[0]= 0
    for i in range(len(Peptide)):
        for s in Alphabet:
            if s == Peptide[i] :
                PrefixMass[i+1] = PrefixMass[i] + AminoAcidMass[s]
    peptideMass = PrefixMass[len(Peptide)]
    CyclicSpectrum = [0]
    for i in range(len(PrefixMass)):
        for j in range(i+1, len(PrefixMass)):
            a  =PrefixMass[j] - PrefixMass[i]
            CyclicSpectrum.append(a)
            if i > 0 and j < (len(PrefixMass)-1):
                b = peptideMass - PrefixMass[j] + PrefixMass[i]
                CyclicSpectrum.append(b)
    CyclicSpectrum.sort()
    return CyclicSpectrum


#print(*data, sep=' ')
def CountingPeptidesGivenMass(m):
    """
    Input: An integer m.
    Output: The number of linear peptides having integer mass m.
    """
    table = {}
    count = 0
    a = list(mass2.keys())
import copy

def Mass(pep):
    mass3 = 0
    for p in pep:
        mass3 += mass[p]
    return mass3
def Expand(Peptid):
    pept = []
    for pep in Peptid:
        for am in aminoAcid:
            ne = pep + am
            pept.append(ne)
    return pept
def consistent(peptide, spectrum):
    b = copy.deepcopy(spectrum)
    nk = LinearSpectrum(peptide, aminoAcid, mass)
    for nnk in nk:
        if nnk not in b:
            return "no"
        else:
            b.remove(nnk)
    return "yes"

def CyclopeptideSequencing(Spectrum):
    CandidatePeptides =[""]
    FinalPeptides =[]
    while CandidatePeptides != []:
        CandidatePeptides = Expand(CandidatePeptides)
        a = copy.deepcopy(CandidatePeptides)
        for Peptide in a:
            if Mass(Peptide) == Spectrum[-1]:
                if CyclicSpectrum(Peptide, aminoAcid, mass) == Spectrum and Peptide not in FinalPeptides:
                    FinalPeptides.append(Peptide)
                CandidatePeptides.remove(Peptide)
            elif consistent(Peptide, Spectrum) == "no":
                CandidatePeptides.remove(Peptide)
    finalspec = []
    for final in FinalPeptides:
        go = []
        for fina in final:
            go.append(str(mass[fina]))
        go2 = "-".join(go)
        if go2 not in finalspec:
            finalspec.append(go2)
    return finalspec
#data = CyclicSpectrum("RQAAGKGWAQLWLH", aminoAcid, mass)
#print(*data, sep=' ')
def CyclopeptideScoring(Peptide, Spectrum):
    """
    Input: An amino acid string Peptide and a collection of integers Spectrum.
    Output: The score of Peptide against Spectrum, Score(Peptide, Spectrum).
    """
    score = 0
    spec = CyclicSpectrum(Peptide, aminoAcid, mass)
    b = copy.deepcopy(Spectrum)
    for s in spec:
        if s in b:
            b.remove(s)
            score+= 1
    return score
def LinearScore(Peptide, Spectrum):
    """
    """
    score = 0
    spec = LinearSpectrum(Peptide, aminoAcid, mass)
    b = copy.deepcopy(Spectrum)
    for s in spec:
        if s in b:
            b.remove(s)
            score+= 1
    return score


def Trim(Leaderboard, Spectrum, N):
    """
    """
    if Leaderboard == {}:
        return {}
    scores = []
    leaderscore = {}
    finalleader = {}
    newpepts = []
    for leader in list(Leaderboard.keys()):
        a = LinearScore(leader, Spectrum)
        if a not in leaderscore:
            leaderscore[a] = [leader]
        else:
            leaderscore[a].append(leader)
        scores.append(a)
    if len(scores) < N:
        return Leaderboard
    scores.sort(reverse = True)
    a= scores[0]
    newpepts.extend(leaderscore[a])
    for i in range(1, N):
        a = scores[i]
        olda = scores[i-1]
        if olda == a :
            continue
        newpepts.extend(leaderscore[a])
    for pep in newpepts:
        finalleader[pep] = Leaderboard[pep]
    return finalleader

    # else:
    #     scores.sort(reverse = True)
    # n = scores[N]
    # for e in scores:
    #     if e > N:
    #         for lea in leaderscore[e]:
    #             finalleader[lea]=Leaderboard[lea]
    # return finalleader
aminoAcid2 = ['G', 'A', 'S', 'P', 'V', 'T', 'C', 'L', 'N', 'D', 'Q', 'E', 'M', 'H', 'F', 'R', 'Y', 'W']
def Expand3(Peptid):
    pept = {}
    for pep in list(Peptid.keys()):
        for am in aminoAcid2:
            ne = pep + am
            m= Peptid[pep]+ mass[am]
            pept[ne] = m
    return pept
def LeaderboardCyclopeptideSequencing(N, Spectrum):
    Leaderboard = {"":0}
    LeaderPeptide = ""
    leaderpeptides = []
    while Leaderboard != {}:
        Leaderboard = Expand3(Leaderboard)
        Leaderboard2 = copy.deepcopy(Leaderboard)
        for Peptide in list(Leaderboard2.keys()):
            mmmm = Leaderboard2[Peptide]
            if mmmm == Spectrum[-1]:
                if CyclopeptideScoring(Peptide, Spectrum) > CyclopeptideScoring(LeaderPeptide, Spectrum):
                    LeaderPeptide = Peptide
                    leaderpeptides = [LeaderPeptide]
                    Leaderboard.pop(Peptide)
                elif CyclopeptideScoring(Peptide, Spectrum) == CyclopeptideScoring(LeaderPeptide, Spectrum):
                    leaderpeptides.append(Peptide)
                    Leaderboard.pop(Peptide)
            elif mmmm > Spectrum[-1] or (Spectrum[-1] -mmmm)< 57 :
                # exclude (Spectrum[-1] -mmmm)< 57 for Spectrum25 in interactive book
                Leaderboard.pop(Peptide)
        Leaderboard = Trim(Leaderboard, Spectrum, N)
    return LeaderPeptide, leaderpeptides

# file2 =  open("hi2.txt", encoding="UTF-16").read()
# pa = file2.split(" ")
# pa3  = []
# for p in pa:
#     pa3.append(int(p))
# re = LeaderboardCyclopeptideSequencing(1000, pa3)
# re2 = []
# re22 = []
# print(len(re[1]))
# for r1 in re[1]:
#     re2 = []
#     for r in r1:
#         re2.append(str(mass[r]))
#     re2 = "-".join(re2)
#     if re2 not in re22:
#         re22.append(re2)
# print(*re22, sep=' ')



def LeaderboardCyclopeptideSequencing2(N, Spectrum):
    Leaderboard = {"":0}
    LeaderPeptide = ""
    while Leaderboard != {}:
        Leaderboard = Expand3(Leaderboard)
        Leaderboard2 = copy.deepcopy(Leaderboard)
        for Peptide in list(Leaderboard2.keys()):
            mmmm = Leaderboard2[Peptide]
            if mmmm == Spectrum[-1]:
                if LinearScore(Peptide, Spectrum) > LinearScore(LeaderPeptide, Spectrum):
                    LeaderPeptide = Peptide
                Leaderboard.pop(Peptide)
            elif mmmm > Spectrum[-1] or (Spectrum[-1] -mmmm)< 57 :
                Leaderboard.pop(Peptide)
        Leaderboard = Trim(Leaderboard, Spectrum, N)
    return LeaderPeptide

#print(pa2)
#res = CyclopeptideSequencing(pa2)
#print(*res, sep=' ')
#print(res)
#file2 =  open("hi2.txt", encoding="UTF-16").read()

#pa = file2.split(" ")
#pa3  = []
#for p in pa:
#    pa3.append(int(p))
#re = LeaderboardCyclopeptideSequencing2(210, pa3)
#re2 = []
#for r in re:
#    re2.append(str(mass[r]))
#re2 = "-".join(re2)
#print(re2)
# file2 =  open("hi2.txt", encoding="UTF-16").read()
# pa = file2.split(" ")
# pa3  = []
# for p in pa:
#     pa3.append(int(p))
# re = LeaderboardCyclopeptideSequencing2(50, pa3)
# re2 = []
# for r in re:
#     re2.append(str(mass[r]))
# re2 = "-".join(re2)
# print(re2)

# THE FOLLOWING CODE IS FOR ALL AMINACID MASSES BETWEEN 57 AND 200
nonprot = [57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200]
def LinearSpectrumnonp(Peptide, Alphabet):
    """
    Input: An amino acid string Peptide.
    Output: The linear spectrum of Peptide.
    """
    PrefixMass = {}
    PrefixMass[0]= 0
    for i in range(len(Peptide)):
        for s in Alphabet:
            if s == Peptide[i] :
                PrefixMass[i+1] = PrefixMass[i] + int(Peptide[i])
    LinearSpectrum = [0]
    for i in range(len(PrefixMass)):
        for j in range(i+1, len(PrefixMass)):
            a  =PrefixMass[j] - PrefixMass[i]
            LinearSpectrum.append(a)
    LinearSpectrum.sort()
    return LinearSpectrum


def CyclicSpectrumnonp(Peptide, Alphabet):
    """
    Input: An amino acid string Peptide.
    Output: The cyclic spectrum of Peptide.
    """
    PrefixMass = {}
    PrefixMass[0]= 0
    for i in range(len(Peptide)):
        for s in Alphabet:
            if s == Peptide[i] :
                PrefixMass[i+1] = PrefixMass[i] + int(Peptide[i])
    peptideMass = PrefixMass[len(Peptide)]
    CyclicSpectrum = [0]
    for i in range(len(PrefixMass)):
        for j in range(i+1, len(PrefixMass)):
            a  =PrefixMass[j] - PrefixMass[i]
            CyclicSpectrum.append(a)
            if i > 0 and j < (len(PrefixMass)-1):
                b = peptideMass - PrefixMass[j] + PrefixMass[i]
                CyclicSpectrum.append(b)
    CyclicSpectrum.sort()
    return CyclicSpectrum
def CyclopeptideScoringnonp(Peptide, Spectrum):
    """
    Input: An amino acid string Peptide and a collection of integers Spectrum.
    Output: The score of Peptide against Spectrum, Score(Peptide, Spectrum).
    """
    score = 0
    spec = CyclicSpectrumnonp(Peptide, nonprot)
    b = copy.deepcopy(Spectrum)
    for s in spec:
        if s in b:
            b.remove(s)
            score+= 1
    return score
def LinearScorenonp(Peptide, Spectrum):
    """
    """
    score = 0
    spec = LinearSpectrumnonp(Peptide, nonprot)
    b = copy.deepcopy(Spectrum)
    for s in spec:
        if s in b:
            b.remove(s)
            score+= 1
    return score


def Trimnonp(Leaderboard, Spectrum, N):
    """
    """
    if Leaderboard == {}:
        return {}
    scores = []
    leaderscore = {}
    finalleader = {}
    newpepts = []
    for leader in list(Leaderboard.keys()):
        a = LinearScorenonp(leader, Spectrum)
        if a not in leaderscore:
            leaderscore[a] = [leader]
        else:
            leaderscore[a].append(leader)
        scores.append(a)
    if len(scores) < N:
        return Leaderboard
    scores.sort(reverse = True)
    a= scores[0]
    newpepts.extend(leaderscore[a])
    for i in range(1, N):
        a = scores[i]
        olda = scores[i-1]
        if olda == a :
            continue
        newpepts.extend(leaderscore[a])
    for pep in newpepts:
        finalleader[pep] = Leaderboard[pep]
    return finalleader

    # else:
    #     scores.sort(reverse = True)
    # n = scores[N]
    # for e in scores:
    #     if e > N:
    #         for lea in leaderscore[e]:
    #             finalleader[lea]=Leaderboard[lea]
    # return finalleader
def Expandnonp(Peptid, alphapet):
    pept = {}
    for pep in list(Peptid.keys()):
        for am in alphapet:
            ne = (*pep, am)
            m= Peptid[pep]+ am
            pept[ne] = m
    return pept
def LeaderboardCyclopeptideSequencingnonp(N, Spectrum, alphapet):
    Leaderboard = {():0}
    LeaderPeptide = ()
    leaderpeptides = []
    while Leaderboard != {}:
        Leaderboard = Expandnonp(Leaderboard, alphapet)
        Leaderboard2 = copy.deepcopy(Leaderboard)
        for Peptide in list(Leaderboard2.keys()):
            mmmm = Leaderboard2[Peptide]
            if mmmm == Spectrum[-1]:
                if CyclopeptideScoringnonp(Peptide, Spectrum) > CyclopeptideScoringnonp(LeaderPeptide, Spectrum):
                    LeaderPeptide = Peptide
                    leaderpeptides = [LeaderPeptide]
                    Leaderboard.pop(Peptide)
                elif CyclopeptideScoringnonp(Peptide, Spectrum) == CyclopeptideScoringnonp(LeaderPeptide, Spectrum):
                    leaderpeptides.append(Peptide)
                    Leaderboard.pop(Peptide)
            elif mmmm > Spectrum[-1]:
                # exclude (Spectrum[-1] -mmmm)< 57 for Spectrum25 in interactive book
                Leaderboard.pop(Peptide)
        Leaderboard = Trimnonp(Leaderboard, Spectrum, N)
    return LeaderPeptide, leaderpeptides
# file2 =  open("hi2.txt", encoding="UTF-16").read()
# pa = file2.split(" ")
# pa3  = []
# for p in pa:
#     pa3.append(int(p))
# re = LeaderboardCyclopeptideSequencingnonp(1000, pa3)
# print(len(re), re)
# re2 = []
# re22 = []
# for r1 in re:
#     r3 = []
#     for r2 in r1:
#         r3.append(str(r2))
#     re2 = "-".join(r3)
#     if re2 not in re22:
#         re22.append(re2)
# print(*re22, sep=' ')
def SpectralConvolution(Spectrum):
    """
    Input: A collection of integers Spectrum.
    Output: The list of elements in the convolution of Spectrum. If an element has multiplicity k, it should appear exactly k times; you may return the elements in any order.
    """
    finalspec = []
    spec = Spectrum[:]
    for e in spec:
        for i in Spectrum:
            a = e - i
            if a > 0:
                finalspec.append(a)
    finalspec.sort(reverse = True)
    return finalspec
# file2 =  open("hi2.txt", encoding="UTF-16").read()
#
# pa = file2.split(" ")
# pa3  = []
# for p in pa:
#     pa3.append(int(p))
# res = SpectralConvolution(pa3)
# print(*res, sep=' ')
from collections import Counter
def ConvolutionCyclopeptideSequencing(M, N, Spectrum):
    """
    Input: An integer M, an integer N, and a collection of (possibly repeated) integers Spectrum.
    Output: A cyclic peptide LeaderPeptide with amino acids taken only from the
    top M elements (and ties) of the convolution of Spectrum that fall between
     57 and 200, and where the size of Leaderboard is restricted to the top N (and ties).
    """
    conv = SpectralConvolution(Spectrum)
    occurence_count = Counter(conv)
    frequency = list(occurence_count.values())
    frequency = set(frequency)
    frequency = list(frequency)
    frequency.sort(reverse = True)
    bestM = []
    for freq in frequency:
        if len(bestM) < M:
            for key in occurence_count:
                if occurence_count[key] == freq and 56<key<201:
                    bestM.append(key)
        else:
            break
    return LeaderboardCyclopeptideSequencingnonp(N, Spectrum, bestM)
file2 =  open("hi2.txt", encoding="UTF-16").read()

pa = file2.split(" ")
pa3  = []
for p in pa:
    pa3.append(int(p))
re = ConvolutionCyclopeptideSequencing(20, 1000, pa3)[1]
print(len(re))
re2 = []
re22 = []
for r1 in re:
    r3 = []
    for r2 in r1:
        r3.append(str(r2))
    re2 = "-".join(r3)
    if re2 not in re22:
        re22.append(re2)
print(*re22, sep=' ')
