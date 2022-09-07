from __future__ import print_function
def PatternCount(Pattern, Text):
    """
    input:
    Text is a string of DNA
    Pattern is a string (k mer we are looking for)
    output: count of Pattern occurences in Text
    """
    count = 0
    for i in range(len(Text) - len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count + 1
    return count


def FrequencyMap(Text, k):
    """
    input:
    Text is a string of DNA
    k is the length of substring
    output: a dictionary of k-mers with length k and their frequencies
     as values
    """
    freq = {}
    n = len(Text)
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        if Pattern not in freq:
            freq[Pattern] = 1
        else:
            freq[Pattern]+= 1
    return freq
#text1 ="ACGCGGCTCTGAAA"
#k1 = 2
#print(FrequencyMap(text1, 3))
def FrequentWords(Text, k):
    """
    input:
    Text is a string of DNA
    k is the length of substring
    output: the maximum frequency, a list of the most frquent k-mers
    """
    freq = []
    map = FrequencyMap(Text, k)
    maxi = max(map.values())
    for key in map:
        if map[key] == maxi:
            freq.append(key)
    return maxi, freq

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
def PatternMatching(Pattern, Genome):
    """
    Input: Strings Pattern and Genome.
    Output: All starting positions in Genome where Pattern appears as
    a substring.
    """
    positions= []
    for i in range(len(Genome)):
        if Genome[i: i+len(Pattern)] == Pattern:
            positions += [i]
    return positions
#myfile = open('Vibrio_cholerae.txt', 'r')
#data = myfile.read()
def ClumpFinding(Genome, k, L, t):
    """
    Input: A string Genome, and integers k, L, and t.
    Output: All distinct k-mers forming (L, t)-clumps in Genome
    """
    li = []
    for i in range(len(Genome)):
        window = Genome[i:i+L]
        a = FrequencyMap(window, k)
        for key in a.keys():
            if a[key] == t:
                if key not in li:
                    li.append(key)
    return li

def SymbolArray(Genome, symbol):
    """
    Input:  Strings Text and Pattern
    Output: The number of times Pattern appears in Text
    """
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    for i in range(n):
        array[i] = PatternCount(symbol, ExtendedGenome[i:i+(n//2)])
    return array

def FasterSymbolArray(Genome, symbol):
    """
    Input:  Strings Text and Pattern
    Output: The number of times Pattern appears in Text
    """
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    array[0] = PatternCount(symbol, ExtendedGenome[0:0+(n//2)])
    str = len(symbol)
    count = array[0]
    for e in range(1, n):
        if ExtendedGenome[e-1: e+str-1] == symbol:
            count = count - 1
        if ExtendedGenome[e+(n//2)-str:e+(n//2)] == symbol:
            count += 1
        array[e] = count
    return array
def FasterSymbolArrayIns(Genome, symbol):
    """
    form instructors
    mine is more correct because it considers differnet lengths of sympol
    """
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]

    # look at the first half of Genome to compute first array value
    array[0] = PatternCount(symbol, Genome[0:n//2])

    for i in range(1, n):
        # start by setting the current array value equal to the previous array value
        array[i] = array[i-1]

        # the current array value can differ from the previous array value by at most 1
        if ExtendedGenome[i-1] == symbol:
            array[i] = array[i]-1
        if ExtendedGenome[i+(n//2)-1] == symbol:
            array[i] = array[i]+1
    return array
def SkewArray(Genome):
    """
    Input:  A String Genome
    Output: The skew array of Genome as a list.
    """
    sarray = {}
    sarray[0] = 0
    for i in range(len(Genome)):
        if Genome[i] == "C":
            sarray[i+1] = sarray[i] - 1
        elif Genome[i] == "G":
            sarray[i+1] = sarray[i] + 1
        else:
            sarray[i+1] = sarray[i]
    return sarray

def MinimumSkew(Genome):
    """
    Input: A DNA string Genome.
    Output: A list of positions where difference between G and C is min
    (All integer(s) i minimizing Skewi (Genome) among all values of i
    (from 0 to |Genome|).)
    """
    skew = SkewArray(Genome)
    a = min(list(skew.values()))
    li = []
    for i in skew.keys():
        if skew[i] == a:
            li.append(i)
    return li

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


def ApproximatePatternMatching(Text, Pattern, d):
    """
    Input:  Strings Pattern and Text along with an integer d
    Output: A list containing all starting positions where Pattern appears
    as a substring of Text with at most d mismatches
    """
    lis = []
    for i in range(len(Text)):
        if len(Pattern) == len(Text[i: i+len(Pattern)]):
            if HammingDistance(Pattern, Text[i: i+len(Pattern)]) <= d:
                lis.append(i)
    return lis

def ApproximatePatternCount(Pattern, Text, d):
    """
    Input:  Strings Pattern and Text, and an integer d
    Output: The number of times Pattern appears in Text with at most d mismatches
    """
    count = 0
    for i in range(len(Text) - len(Pattern)+1):
        if HammingDistance(Pattern, Text[i: i+len(Pattern)]) <= d:
            count = count + 1
    return count
def PatternToNumber(Pattern):
    """
    input: a Dna Pattern
    output: converts the string to an integer in a 4 numeric base ("ACGT")
    more on int fun :https://stackoverflow.com/questions/23190060/what-does-base-value-do-in-int-function
    """
    pattern2 = ''
    for i in Pattern:
        if i == "A":
            pattern2 += "0"
        if i == "C":
            pattern2 += "1"
        if i == "G":
            pattern2 += "2"
        if i == "T":
            pattern2 += "3"
    return int(pattern2, 4)

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

def ComputingFrequencies(Text, k):
    """
    input: text is a string of Dna, k is the length of k-mers
    output: frequencies of evere possible arrangment of "ACGT" with that length
    in the text
    """
    FrequencyArray = {}
    for i in range(4**k):
        FrequencyArray[i] = 0
    for i in range(len(Text) - k+1):
        Pattern = Text[i : i+ k]
        j  = PatternToNumber(Pattern)
        FrequencyArray[j] += 1
    return FrequencyArray


def FrequentWordswithMismatchesProblem(Text, k, d):
    """
    Input: A string Text as well as integers k and d. (You may assume k ≤ 12 and d ≤ 3.)
    Output: All most frequent k-mers with up to d mismatches in Text.
    """
    kmers = []
    counts = []
    count = 0
    for i in range(4**k):
        kmers.append(i)
    for kmer in kmers:
        pattern = NumberToPattern(kmer, k)
        co = ApproximatePatternCount(pattern, Text, d)
        if co > count:
            counts = [pattern]
            count = co
        elif co == count:
            counts.append(pattern)
    return counts


def Neighbors(Pattern, d):
    """
    Input: pattern of dna, and d an integer of most accepted mismatches.
    Output: all strings with d mismatched to pattern at most.
    """
    if d == 0:
        return Pattern
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


def ModifiedNeighbors(Pattern, d):
    """
    Input: pattern of dna, and d an integer of most accepted mismatches.
    Output: all strings with exactly d mismatched to pattern.
    """
    if d == 0:
        return Pattern
    if len(Pattern) == 1:
        return ["A", "C", "G", "T"]
    Neighborhood = []
    firstsymbol = Pattern[0]
    SuffixNeighbors = Neighbors(Pattern[1:], d)
    for string in SuffixNeighbors:
        if HammingDistance(Pattern[1:], string) == d - 1:
            for symbol in "ACGT":
                if symbol != firstsymbol:
                    pat = symbol + string
                    Neighborhood.append(pat)
        else:
            pat =  Pattern[0] + string
            Neighborhood.append(pat)
    return Neighborhood

def ComputingFrequenciesWithMismatches(Text, k, d):
    """
    Input: text, astring of Dna, k is an integer for length of pattern, d is an
    integer of most mismatches.
    Output: frequency array of every pattern with up to d mismatches.
    """
    frequencyarray = {}
    for i in range(4**k):
        frequencyarray[i] = 0
    for i in range(len(Text)- k + 1):
        Pattern = Text[i: i+k]
        Neighborhood = Neighbors(Pattern, d)
        for string in Neighborhood:
            j = PatternToNumber(string)
            frequencyarray[j] = frequencyarray[j] + 1
    return frequencyarray


def FasterFrequentWordsWithMismatches(Text, k, d):
    """
    Input: A string Text as well as integers k and d.
    Output: All most frequent k-mers with up to d mismatches in Text.
    """
    freq = []
    map = ComputingFrequenciesWithMismatches(Text, k, d)
    maxi = max(map.values())
    for key in map:
        if map[key] == maxi:
            freq.append(NumberToPattern(key, k))
    return freq

def FrequentWordswithMismatchesandReverseComplementsProblem(Text, k, d):
    """
    Input: A string Text as well as integers k and d. (You may assume k ≤ 12 and d ≤ 3.)
    Output: All most frequent k-mers and reverse complement  with up to d mismatches in Text.
    """
    kmers = []
    counts = []
    count = 0
    for i in range(4**k):
        kmers.append(i)
    for kmer in kmers:
        pattern = NumberToPattern(kmer, k)
        repattern = ReverseComplement(pattern)
        co1 = ApproximatePatternCount(pattern, Text, d)
        co2 = ApproximatePatternCount(repattern, Text, d)
        co = co1+co2
        if pattern == repattern:
            co = co/2
        if co > count:
            counts = [pattern, repattern]
            count = co
        elif co == count:
            counts.append(pattern)
    return list(dict.fromkeys(counts))

def ComputingFrequenciesWithMismatchesandReverseComplementsProblem(Text, k, d):
    """
    Input: text, astring of Dna, k is an integer for length of pattern, d is an
    integer of most mismatches.
    Output: frequency array of every pattern and its reverse complement
    with up to d mismatches.
    """
    frequencyarray = {}
    for i in range(4**k):
        frequencyarray[i] = 0
    for i in range(len(Text)- k + 1):
        Pattern = Text[i: i+k]
        repattern = ReverseComplement(Pattern)
        Neighborhood = Neighbors(Pattern, d)
        Neighborhood2 = Neighbors(repattern, d)
        for string in Neighborhood:
            j = PatternToNumber(string)
            frequencyarray[j] = frequencyarray[j] + 1
        for string in Neighborhood2:
            j = PatternToNumber(string)
            frequencyarray[j] = frequencyarray[j] + 1
    return frequencyarray

def FasterFrequentWordsWithMismatchesandReverseComplementsProblem(Text, k, d):
    """
    Input: A string Text as well as integers k and d.
    Output: All most frequent k-mers and their reverse complemetn
    with up to d mismatches in Text.
    """
    freq = []
    map = ComputingFrequenciesWithMismatchesandReverseComplementsProblem(Text, k, d)
    maxi = max(map.values())
    for key in map:
        if map[key] == maxi:
            freq.append(NumberToPattern(key, k))
    return list(dict.fromkeys(freq).keys())


def ImmediateNeighbors(Pattern):
    Neighborhood = [Pattern]
    for i in range(len(Pattern)):
        symbol = Pattern[i]
        for x in "ACGT":
            if x != symbol:
                Neighbor= Pattern [:i] + x +Pattern[i+1:]
                Neighborhood.append(Neighbor)
    return Neighborhood
def ModifiedNeighbors(Pattern, d):
    """
    Input: pattern of dna, and d an integer of most accepted mismatches.
    Output: all strings with exactly d mismatched to pattern.
    """
    if d == 0:
        return Pattern
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

def IterativeNeighbors(Pattern, d):
    Neighborhood = [Pattern]
    for j in range(d):
        for string in Neighborhood:
            neighbor = ImmediateNeighbors(Pattern)
        Neighborhood.extend(neighbor)
    return list(dict.fromkeys(Neighborhood))

def BetterClumpFinding(Genome, k, L, t):
    FrequentPatterns = []
    Clump = {}
    for i in range(4**k):
            Clump[i] = 0
    Text = Genome[0: L]
    FrequencyArray = list(ComputingFrequencies(Text, k))
    for i in range(4**k):
        if FrequencyArray[i] >= t:
            Clump[i] = 1
    for i in range(1 ,len(Genome) - L +1):
        FirstPattern = Genome[i - 1: i-1+k]
        index = PatternToNumber(FirstPattern)
        FrequencyArray[index] = FrequencyArray[index] - 1
        LastPattern = Genome[i + L - k: i+L]
        index = PatternToNumber(LastPattern)
        FrequencyArray[index] = FrequencyArray[index] + 1
        if FrequencyArray[index] >= t:
            Clump[index] = 1
    for i in range(4**k):
        if Clump[i] == 1:
            Pattern = NumberToPattern(i, k)
            if Pattern not in FrequentPatterns:
                FrequentPatterns.append(Pattern)
    return FrequentPatterns
#myfile = open('E_coli.txt', 'r')
#data = myfile.read()
#a = BetterClimpFinding(data ,9 ,500 ,3)
#print(len(a))
#print(*a, sep = " ")


# to print list without brackets or commas with the import at beggining
#print(*a, sep = " ")
def FasterFrequentWords(Text, k):
    """
    input: text is a string of Dna, k is the length of k-mers
    output: a list of the most frequent patterns or pattern.
    """
    FrequentPatterns = []
    FrequencyArray = list(ComputingFrequencies(Text, k))
    maxCount = max(FrequencyArray)
    for i in range(4**k):
        if FrequencyArray[i] == maxCount:
            Pattern = NumberToPattern(i, k)
            FrequentPatterns.append(Pattern)
    return FrequentPatterns
import operator
def FindingFrequentWordsBySorting(Text , k):
    """
    input: text is a string of Dna, k is the length of k-mers
    output: list of most frequent k-mers.
    """
    patt = []
    indexd = {}
    for i in range(len(Text)- k +1):
        pat = Text[i:i+k]
        indexd[PatternToNumber(pat)] = indexd.get(PatternToNumber(pat), 0) +1
    tup = list(indexd.items())
    tup.sort(key = operator.itemgetter(1))
    count = tup[-1][1]
    for item in tup:
        if item[1] == count:
            patt.append(NumberToPattern(item[0], k))
    return patt


def FrequentWordsWithMismatchesSorting(Text, k, d):

    Neighborhoods =[]
    NeighborhoodArray = []
    countarray ={}
    for i in range(len(Text) - k +1):
        Neighborhoods.extend(Neighbors(Text[i: i+k], d))
    for i in Neighborhoods:
        NeighborhoodArray.append(PatternToNumber(i))
    NeighborhoodArray.sort()
    pat = NeighborhoodArray[0]
    countarray[pat] = 1
    for i in NeighborhoodArray[1:]:
        if i == pat:
            countarray[pat] = countarray.get(pat, 0) + 1
        else:
            pat = i
    FrequentNumbers = [key for m in [max(countarray.values())] for key,val in countarray.items() if val == m]
    FrequentPatterns = []
    for num in FrequentNumbers:
        FrequentPatterns.append(NumberToPattern(int(num), k))
    return FrequentPatterns

#getting a list of keys with highest value
#FrequentNumbers = [key for m in [max(countarray.values())] for key,val in countarray.items() if val == m]
#getting a list with distinct elements
#list(dict.fromkeys(freq).keys())
