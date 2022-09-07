def HammingDistance(p, q):
    """
    Input: Two strings of equal length.
    Output: An integer value representing the Hamming Distance between p and q.
    """
    hamdis= 0
    for i in range(len(p)):
        if p[i] != q[i]:
            hamdis += 1
    return hamdis

def Neighbors(Pattern, d):
    """
    Input: pattern of dna, and d an integer of most accepted mismatches.
    Output: all strings with d mismatched to pattern at most.
    """
    if d == 0:
        return [Pattern]
    if len(Pattern) == 1:
        return ["A", "C", "G", "T"]
    Neighborhood = []
    SuffixNeighbors = Neighbors(Pattern[1:], d)
    for string in SuffixNeighbors:
        if HammingDistance(Pattern[1:], string) < d:
            for symbol in "ACGT":
                pat = symbol + string
                Neighborhood.append(pat)
        else:
            pat =  Pattern[0] + string
            Neighborhood.append(pat)
    return Neighborhood

def NumberToPattern(index, k):
    """
    input: index whoch is number in numerical base, k length of pattern.
    output: pattern corresponding to the value given.
    we change base from decimal to base 4
    """
    num = index
    pat = ""
    lis = []
    for i in range(k):
        remainder = num % 4
        num = num // 4
        lis.insert(0, remainder)
    for i in lis:
        if i == 0:
            pat += "A"
        if i == 1:
            pat += "C"
        if i == 2:
            pat += "G"
        if i == 3:
            pat += "T"
    return pat
def MotifEnumeration(Dna, k, d):
    """
    Input: Integers k and d, followed by a collection of strings Dna.
    Output: All (k, d)-motifs in Dna using brute force algorithm.
    """
    Patterns = []
    frequents = {}
    count = {}
    for i in range(4**k):
        Pattern = NumberToPattern(i, k)
        if d >0:
            frequents[Pattern] = Neighbors(Pattern, d)
    for frequent in frequents.keys():
        values = frequents[frequent]
        count[frequent] = [0]*len(Dna)
        for value in values:
            for a in range(len(Dna)):
                if value in Dna[a]:
                    count[frequent][a] += 1
    for frequent in count.keys():
            if 0 not in count[frequent]:
                Patterns.append(frequent)
    return list(dict.fromkeys(Patterns).keys())

#a = MotifEnumeration(["AAAAA", "AAAAA", "AACAA"], 3, 0)
#print(*a, sep = " ")
def Count(Motifs):
    """
    Input:  A set of kmers Motifs (list of strings)
    Output: Count(Motifs) (a dictionary of lists)
    see bioinformatis for beginners(dropbox) for count
    """
    count = {}
    j = len(Motifs[0])
    for sympol in "ACTG":
        count[sympol] = [0] * j
    for i in Motifs:
        for e in range(j):
            if i[e] == "A":
                count["A"][e] += 1
            elif i[e] == "T":
                count["T"][e] += 1
            elif i[e] == "C":
                count["C"][e] += 1
            else:
                count["G"][e] += 1
    return count

def Profile(Motifs):
    """
    Input:  A set of kmers Motifs (list of strings)
    Output: the profile matrix of Motifs, as a dictionary of lists.
    """
    profile = Count(Motifs)
    k = len(Motifs[0])
    t = len(Motifs)
    for ele in profile:
        for i in range(k):
            profile[ele][i] = profile[ele][i]/t
    return profile

def Consensus(Motifs):
    """
    Input:  A set of kmers Motifs
    Output: A consensus STRING of Motifs.
    """
    count = Count(Motifs)
    k = len(Motifs[0])
    t = len(Motifs)
    con = ""
    for i in range(k):
        most = 0
        msympol = ""
        for key in count.keys():
            if count[key][i] > most:
                most = count[key][i]
                msympol += key
        con += msympol[-1]
    return con

def Score(Motifs):
    """
    Input:  A set of k-mers Motifs
    Output: The score of these k-mers.
    """
    con = Consensus(Motifs)
    score = 0
    for i in range(len(con)):
        sum = 0
        for kmer in Motifs:
            if kmer[i] != con[i]:
                sum += 1
        score += sum
    return score

def Pr(Text, Profile):
    """
    Input:  String Text and profile matrix Profile
    Output: Pr(Text, Profile)
    """
    prob = 1
    for i in range(len(Text)):
        if Text[i] == "A":
            prob *= Profile["A"][i]
        if Text[i] == "C":
            prob *= Profile["C"][i]
        if Text[i] == "G":
            prob *= Profile["G"][i]
        if Text[i] == "T":
            prob *= Profile["T"][i]
    return prob


def ProfileMostProbableKmer(Text, k, Profile):
    """
    Input: A string Text, an integer k, and a 4 x k matrix Profile.
    Output: A Profile-most probable k-mer in Text.
    """
    Mpr = -1
    # if probability is zero we don't get qn error
    for i in range(len(Text)-k+1):
        if Pr(Text[i:i+k], Profile) > Mpr:
            Mpr = Pr(Text[i:i+k], Profile)
            mk = Text[i:i+k]
    return mk

def GreedyMotifSearch(Dna, k, t):
    """
    Input:  A list of kmers Dna, and integers k and t (where t is the number of kmers in Dna)
    Output: GreedyMotifSearch(Dna, k, t)
    """
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = Profile(Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
                BestMotifs = Motifs
    return BestMotifs


#file1 = open("dors.txt","r")
#Dna = file1.readlines()
#a = GreedyMotifSearch(Dna, 15, len(Dna))
#print(a, Score(a))
import math
def entropy(Profile):
    """
    Input : profile (dictionary of four neucleotides and there probability for each column)
    output: entropy (the sum of entropy of each column) and a list of entropies of every column
    """
    entropy = 0
    entropyc = []
    keys = list(Profile.keys())
    column = len(list(Profile.values())[0])
    for i in range(column):
        ent = 0
        for key in keys:
            if Profile[key][i] > 0:
                ent += Profile[key][i]* math.log(Profile[key][i], 2)
        entropy += - ent
        entropyc.append(-ent)
    return entropy, entropyc

def CountWithPseudocounts(Motifs):
    """
    Input:  A set of kmers Motifs (list of strings)
    Output:CountWithPseudocounts(Motifs) (a dictionary of lists)
    """
    count = {}
    j = len(Motifs[0])
    for sympol in "ACTG":
        count[sympol] = [1] * j
    for i in Motifs:
        for e in range(j):
            if i[e] == "A":
                count["A"][e] += 1
            elif i[e] == "T":
                count["T"][e] += 1
            elif i[e] == "C":
                count["C"][e] += 1
            else:
                count["G"][e] += 1
    return count

def ProfileWithPseudocounts(Motifs):
    """
    Input:  A set of kmers Motifs (list of strings)
    Output: ProfileWithPseudocounts(Motifs), as a dictionary of lists.
    """
    profile = CountWithPseudocounts(Motifs)
    k = len(Motifs[0])
    t = len(Motifs)
    for ele in profile:
        for i in range(k):
            profile[ele][i] = profile[ele][i]/(t+4)
    return profile

def ScoreWithPseudocounts(Motifs):
    """
    Input:  A set of k-mers Motifs
    Output: The score of these k-mers.
    """
    con = ConsensusWithPseudocounts(Motifs)
    score = 0
    for i in range(len(con)):
        sum = 0
        for kmer in Motifs:
            if kmer[i] != con[i]:
                sum += 1
        score += sum
    return score

def ConsensusWithPseudocounts(Motifs):
    """
    Input:  A set of kmers Motifs
    Output: A consensus STRING of Motifs.
    """
    count = CountWithPseudocounts(Motifs)
    k = len(Motifs[0])
    t = len(Motifs)
    con = ""
    for i in range(k):
        most = 0
        msympol = ""
        for key in count.keys():
            if count[key][i] > most:
                most = count[key][i]
                msympol += key
        con += msympol[-1]
    return con

def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    """
    Input:  A list of kmers Dna, and integers k and t (where t is the number of kmers in Dna)
    Output: A collection of strings BestMotifs resulting from applying
    GreedyMotifSearch(Dna, k, t) with pseudocounts.
    """
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = ProfileWithPseudocounts(Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
        if ScoreWithPseudocounts(Motifs) < ScoreWithPseudocounts(BestMotifs):
                BestMotifs = Motifs
    return BestMotifs

#x = open("hii.txt").read()
#x1 = x.split("\n")
#a =GreedyMotifSearchWithPseudocounts(x1, 12, 25)
#print(*a, sep = " ")
def Motifss(Profile, Dna):
    """
    Input:  profile matrix Profile corresponding to a list of strings Dna.
    Output: a list of the Profile-most probable k-mers in each string from Dna.
    """
    list = []
    length = len(Profile["A"])
    for e in range(len(Dna)):
        score = 0
        for i in range(len(Dna[0])- length+1):
            if Pr(Dna[e][i:i+length], Profile) > score:
                score = Pr(Dna[e][i:i+length], Profile)
                mkmer = Dna[e][i:i+length]
        list.append(mkmer)
    return list

import random
def RandomMotifs(Dna, k, t):
    """
    Input: k:length of k-mers required, t: number of different strings Dna.
    Output: random list of motifs.
    """
    motifs = []
    j = len(Dna[0])
    for i in range(t):
        ran = random.randint(0, j-k)
        sub = Dna[i][ran: ran+k]
        motifs.append(sub)
    return(motifs)

def RandomizedMotifSearch(Dna, k, t):
    """
    Input:  Positive integers k and t, followed by a list of strings Dna.
    Output: best random motifs
    """
    M = RandomMotifs(Dna, k, t)
    BestMotifs = M
    while True:
        Profile = ProfileWithPseudocounts(M)
        M = Motifss(Profile, Dna)
        if ScoreWithPseudocounts(M) < ScoreWithPseudocounts(BestMotifs):
            BestMotifs = M
        else:
            return BestMotifs

x = open("hii.txt").read()
x1 = x.split("\n")
#a = RandomizedMotifSearch(x1, 8, 5)
#print(*a, sep = " ")

# if N is too low and i want to make sure random run enough times i can set aloop
#N = 100
#BestMotifs= []
#score = 10000
#for n in range(N):
#    M = RandomizedMotifSearch(x1, 15, 20)
#    sm = ScoreWithPseudocounts(M)
#    if sm < score:
#        score = sm
#        BestMotifs= M
#a = BestMotifs
#for ans in a:
#    print(ans)
#print(*a, sep = " ")
#profile = {"A" : [0.8, 0.0, 0.0, 0.2], "C" : [0.0, 0.6, 0.2, 0.0], "G" : [0.2, 0.2, 0.8, 0.0], "T" : [0.0, 0.2, 0.0, 0.8]}
#dna = ["GGCGTTCAGGCA", "AAGAATCAGTCA", "CAAGGAGTTCGC", "CACGTCAATCAC", "CAATAATATTCG"]
#print(RandomizedMotifSearch(dna, 3, 5))

def Normalize(Probabilities):
    """
    Input: A dictionary Probabilities, where keys are k-mers and values are the
     probabilities of these k-mers (which do not necessarily sum up to 1)
    Output: A normalized dictionary where the probability of each k-mer was
     divided by the sum of all k-mers' probabilities
    """
    sum = 0
    pr = {}
    for key in Probabilities.keys():
        sum += Probabilities[key]
    for key in Probabilities.keys():
        pr[key] = Probabilities[key]/sum
    return pr


def WeightedDie(Probabilities):
    """
    Input:  A dictionary Probabilities whose keys are k-mers and whose values
     are the probabilities of these kmers
    Output: A randomly chosen k-mer with respect to the values in Probabilities
    """
    kmer = ''
    lis = []
    keys = []
    su = 0
    for key in Probabilities.keys():
        lis.append(Probabilities[key])
        keys.append(key)
    p = random.uniform(0, 1)
    for i in range(len(lis)):
        if su < p < lis[i]+su:
            return keys[i]
        su += lis[i]

def ProfileGeneratedString(Text, profile, k):
    """
    Input:  A string Text, a profile matrix Profile, and an integer k
    Output: a randomly generated k-mer
     from Text whose probabilities are generated from profile,
    """
    n = len(Text)
    probabilities = {}
    for i in range(0, n-k+1):
        probabilities[Text[i:i+k]] = Pr(Text[i:i+k], profile)
    Normprobabilities = Normalize(probabilities)
    return WeightedDie(Normprobabilities)

def GibbsSampler(Dna, k, t, N):
    """
    Input:  Integers k, t, and N, followed by a collection of strings Dna
    # Output: Best Motifs
    """
    motifs= RandomMotifs(Dna, k, t)
    bestmotifs = motifs.copy()
    for j in range(1, N):
        i = random.randint(1, t-1)
        bestmotifs = bestmotifs[:i]+bestmotifs[i+1:]
        profile = ProfileWithPseudocounts(bestmotifs)
        bestmotifs.insert(i, ProfileGeneratedString(Dna[i], profile, k))
    if ScoreWithPseudocounts(motifs) < ScoreWithPseudocounts(bestmotifs):
            bestmotifs = motifs
    return bestmotifs


#x = open("hii.txt").read()
#x1 = x.split("\n")
#bestmotifs = GibbsSampler(x1, 3, 4, 100)
#score = ScoreWithPseudocounts(bestmotifs)
#for i in range(20):
#    a = GibbsSampler(x1, 3, 4, 100)
#    if ScoreWithPseudocounts(a)<= ScoreWithPseudocounts(bestmotifs):
#        bestmotifs = a
#        score = ScoreWithPseudocounts(a)
#for e in bestmotifs:
#    print(e)
def dis(Pattern, Text):
    """
    Input: a string pattern and string text.
    Output: the minimum d with each kmer in text.
    """
    Patterns = []
    k = len(Pattern)
    ham = k + 1
    kmer = []
    for i in range(len(Text)- k +1):
        Patterns.append(Text[i: i +k])
    for Pat in Patterns:
        if HammingDistance(Pattern, Pat) < ham:
            ham = HammingDistance(Pattern, Pat)
    return ham
def motifs(Pattern, Text):
        Patterns = []
        k = len(Pattern)
        ham = k + 1
        kmer = []
        for i in range(len(Text)- k +1):
            Patterns.append(Text[i: i +k])
        for Pat in Patterns:
            if HammingDistance(Pattern, Pat) < ham:
                ham = HammingDistance(Pattern, Pat)
                kmer = [Pat]
            if HammingDistance(Pattern, Pat) == ham:
                kmer.append(Pat)
        return kmer[0]
def distance(Pattern, Dna):
    count = 0
    for DNA in Dna:
        count +=dis(Pattern, DNA)
    return count


#x = open("hii.txt").readlines()
#x1 = x[0].split(" ")
def medianstring(Dna, k):
    """
    Input: A collection of strings Dna and an integer k.
    Output: All k-mers Pattern that minimizes d(Pattern, Dna) among all possible choices of k-mers.
    """
    dis = k + 1
    patterns = {}
    for i in range(4**k):
        pattern = NumberToPattern(i, k)
        patterns[pattern] = distance(pattern, Dna)
    FrequentNumbers = [key for m in [min(patterns.values())] for key,val in patterns.items() if val == m]
    return FrequentNumbers

#x = open("hii.txt").read()
#x1 = x.split("\n")
#a = medianstring(x1, 6)
#print(*a, sep = " ")
