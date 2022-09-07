def StringComposition(Text, k):
    """
    Input: An integer k and a string Text.
    Output: a list of all k-mers in lexocaphical order
    """
    composition = []
    for i in range(len(Text)-k +1):
        composition.append(Text[i:i+k])
    composition.sort()
    return composition
#myfile = open('hi.txt', 'r')
#data = myfile.read()
#a = StringCompositionProblem(data, 100)
#print(*a, sep = "\n")
def PathToGenome(path):
    """
    Input: A sequence path of k-mers Pattern1, … ,Patternn such that the last
    k - 1 symbols of Patterni are equal to the first k-1 symbols of Patterni+1 for 1 ≤ i ≤ n-1.
    Output: A string Text of length k+n-1 such that the i-th k-mer in Text is equal to Patterni (for 1 ≤ i ≤ n).
    """
    genome = path[0]
    for pat in path[1:]:
        genome += pat[-1]
    return genome

def Overlap(Patterns):
    """
    Input: A collection Patterns of k-mers.
    Output: The overlap graph Overlap(Patterns)), in the form of an adjacency list.
    """
    Patterns.sort()
    overlap = {}
    k = len(Patterns[0])
    Pattern2 = Patterns[:]
    for pattern in Patterns:
        suffix = pattern[1:]
        lt = []
        Pattern2.remove(pattern)
        for pattern2 in Pattern2:
            prefix = pattern2[:k-1]
            if suffix == prefix:
                lt.append(pattern2)
        Pattern2.append(pattern)
        if len(lt)> 0:
            overlap[pattern] = lt
    return overlap
#file1 =  open("hi.txt").read()
#pa = file1.split("\n")
#a = Overlap(pa)
#for pattern, adjacencies in a.items():
#        print(pattern, '->', ','.join(adjacencies))
def DeBruijnGraph(Text, k):
    """
     Input: An integer k and a string Text.
     Output: DeBruijnk(Text).
    """
    kmer = []
    deb = {}
    for i in range(len(Text)- k + 1):
        kmer.append(Text[i:i+k])
    for i in range(len(kmer)):
        n1 = kmer[i]
        if n1[:k-1] not in deb:
            deb[n1[:k-1]] =[]
        lenh = len(kmer)-1
        if i == lenh:
            n2 = kmer[i]
            deb[n1[:k-1]].append(n2[1:k])
        if i != lenh:
            n2 = kmer[i+1]
            deb[n1[:k-1]].append(n2[:k-1])
    return deb

#myfile = open('hi.txt', 'r')
#data = myfile.read()
#a = DeBruijnGraph(data, 12)
#for pattern, adjacencies in a.items():
#        print(pattern, '->', ','.join(adjacencies))
def DeBruijnGraphkmers(Patterns):
    """
     Input: A collection of k-mers Patterns.
     Output: The adjacency list of the de Bruijn graph DeBruijn(Patterns).
    """
    deb = {}
    k = len(Patterns[0])
    for km in Patterns:
        prefix1 = km[:k-1]
        suffix1 = km[1:]
        if prefix1 in deb:
            deb[prefix1].append(suffix1)
        else:
            deb[prefix1] = [suffix1]
    return deb
from collections import Counter
def MaximalNonBranchingPaths(Graph):
    """
    Input: The adjacency list of a graph whose nodes are integers.
     Output: The collection of all maximal nonbranching paths in this graph.
    """
    Paths = []
    wlists = list(Graph.values())
    wlist = []
    for lis in wlists:
        wlist.extend(lis)
    b= Counter(wlist)
    a = list(Graph.keys())
    for v in a:
        if len(Graph[v])!=1 or b[v]!= 1:
            if len(Graph[v]) > 0:
                for i in range(len(Graph[v])):
                    w =  Graph[v][i]
                    Path = [v, w]
                    while w in Graph and len(Graph[w]) == 1 and b[w] == 1:
                        u = Graph[w][0]
                        a.remove(w)
                        Path.append(u)
                        w = u
                    Paths.append(Path)
    a.reverse()
    for v in a:
        if len(Graph[v]) ==1 and  b[v] == 1:
            cycle= [v]
            w = Graph[v][0]
            cycle.append(w)
            while w in Graph and len(Graph[w]) == 1 and b[w] == 1:
                u = Graph[w][0]
                cycle.append(u)
                if u == v:
                    Paths.append(cycle)
                    for cy in cycle[1:]:
                        a.remove(cy)
                    break
                w = u
    return Paths

#file1 = open("hi.txt").read()
#content = file1.split("\n")
#dict1 = {}
#for line in content: #creates dict, a dictionary, dict[nodes] = connectednodes
#    toAdd = line.split(' -> ')[0]
#    toAdd2 = line.split(' -> ')[1]
#    toAdd2 = toAdd2.split(',')
#    toAdd3 = []
#    for ele in toAdd2:
#        toAdd3.append(int(ele))
#    dict1[int(toAdd)] = toAdd2
#pr = MaximalNonBranchingPaths(dict1)
#for element in pr:
#    element2= []
#    for ele in element:
#        element2.append(str(ele))
#    element3="->".join(element2)
#    print(element3)




def ContigGeneration(Patterns):
    """
     Input: A collection of k-mers Patterns.
     Output: All contigs in DeBruijn(Patterns).
    """
    graph = DeBruijnGraphkmers(Patterns)
    paths = MaximalNonBranchingPaths(graph)
    contings = []
    for path in paths:
        contings.append(PathToGenome(path))
    return contings



#a = ContigGeneration(pa)
#a.sort()
#res = ""
#for a2 in a:
#    res+= a2 +" "
#print(res)
#for pattern, adjacencies in a.items():
#    print(pattern, '->', ','.join(adjacencies))

def EulerianCycle(Graph):
    """
    Input: The adjacency list of an Eulerian directed graph.
    Output: An Eulerian cycle in this graph.
    """
    stack = []
    location = list(Graph.keys())[0]
    circuit = []
    while Graph != {}:
        while location in Graph:
            stack.append(location)
            location = Graph[location][0]
            if len(Graph[stack[-1]]) == 1:
                Graph.pop(stack[-1])
            else:
                Graph[stack[-1]].remove(location)
        circuit.append(location)
        location = stack[-1]
        stack = stack[:-1]
    stack.append(location)
    stack.reverse()
    circuit.extend(stack)
    circuit.reverse()
    return circuit

from collections import Counter
def EulerianPath(Graph):
    """
    Input: The adjacency list of a directed graph that has an Eulerian path.
    Output: An Eulerian path in this graph.
    """
    wlists = list(Graph.values())
    wlist = []
    for lis in wlists:
        wlist.extend(lis)
    a = {}
    for ele in list(Graph.keys()):
        a[ele] = len(Graph[ele])
    b= Counter(wlist)
    for ele in a:
        if ele not in b or a[ele]!= b[ele]:
            v = ele
    for ele in b:
        if ele not in a:
            w = ele
            Graph[w]= [v]
            break
        if a[ele]!= b[ele]:
            w = ele
            Graph[w].append(v)
    stack = []
    location = v
    circuit = []
    while Graph != {}:
        while location in Graph:
            stack.append(location)
            location = Graph[location][0]
            if len(Graph[stack[-1]]) == 1:
                Graph.pop(stack[-1])
            else:
                Graph[stack[-1]].remove(location)
        circuit.append(location)
        location = stack[-1]
        stack = stack[:-1]
    stack.append(location)
    stack.reverse()
    circuit.extend(stack)
    circuit.reverse()
    return circuit[:-1]
def StringReconstruction(k, Patterns):
        dB = DeBruijnGraphkmers(Patterns)
        path = EulerianPath(dB)
        tex = PathToGenome(path)
        return tex
file1 =  open("hi.txt").read()
pa = file1.split("\n")
print(StringReconstruction(4, pa))
#dict = {}
#for i in content: #creates dict, a dictionary, dict[nodes] = connectednodes
#    list2 = i.split() #split by whitespace
#    x = list2[2] #takes the last section
#    dict[int(list2[0])] = [int(s) for s in x.split(',') if s.isdigit()]
#pr = EulerianPath(dict)
#out = ""
#for i in range(0, len(pr)-1):
#    out = out+ str(pr[i]) + "->"
#out = out + str(pr[-1])
#print(out)
import itertools
def kUniversalCircularString(k):
    """
     Input: An integer k.
     Output: A k-universal circular string.
    """
    a = ["".join(seq) for seq in itertools.product("01", repeat=k)]
    dB = DeBruijnGraphkmers(a)
    path = EulerianCycle(dB)
    tex = PathToGenome(path)
    return tex[:len(tex)-k+1], len(tex[:len(tex)-k+1])
#file1 =  open("hi.txt").read()
#pa = file1.split("\n")
#a = kUniversalCircularString(8)
#print(a)
def PairedComposition(k,d,Text):
    """
    Input:
    Output:
    """
    reads = []
    for i in range(len(Text)-2*k-d+1):
        pat1 = Text[i:i+k]
        pat2 = Text[i+k+d:i+k+d+k]
        read = [pat1+"|"+ pat2]
        reads.append(read)
    reads.sort()
    return reads
def deBruijnReadsgraph(k, d, Reads):
    """
    Integers k and d followed by a collection of paired k-mers PairedReads.
    The adjacency list of the de Bruijn graph DeBruijn(Reads)
    """
    Graph = {}
    for read in Reads:
        read = read[1:-1]
        r = read.split("|")
        prefix1 = r[0][:k-1]
        prefix2 = r[1][:k-1]
        suffix1 = r[0][1:]
        suffix2 = r[1][1:]
        prefix = prefix1+"|"+prefix2
        suffix = suffix1+"|"+suffix2
        if prefix in Graph:
            Graph[prefix].append(suffix)
        else:
            Graph[prefix] = [suffix]
    return Graph

def StringSpelledbyaGappedGenomePath(path, k, d):
    """
    Input: A sequence of (k, d)-mers (a1|b1), ..., (an|bn) such that Suffix((ai|bi)) = Prefix((ai+1|bi+1)) for 1 ≤ i ≤ n - 1.
    Output: A string Text of length k + d + k + n - 1 such that the i-th (k, d)-mer of Text is equal to (ai|bi) for 1 ≤ i ≤ n (if such a string
    exists).
    """
    p1 = path[0].split("|")
    prefixstr = p1[0]
    sufixstr = p1[1]
    for p in path[1:]:
        p1 = p.split("|")
        prefixstr += p1[0][-1]
        sufixstr+= p1[1][-1]
    for i in range(k + d + 1, len(prefixstr)):
        if prefixstr[i] != sufixstr[i - k - d]:
            return "there is no string spelled by the gapped patterns"
    return prefixstr+ sufixstr[-k-d:]

def StringReconstructionReadPairs(k, d, Reads):
    """
     Input: Integers k and d followed by a collection of paired k-mers PairedReads.
     Output: A string Text with (k,d)-mer composition equal to PairedReads (if such a string exists).
    """
    graph = deBruijnReadsgraph(k, d, Reads)
    path = EulerianPath(graph)
    gen = StringSpelledbyaGappedGenomePath(path, k, d)
    return gen

#a =  open("hi.txt").read()
#b = a.split("\n")
#c = b[0].split(" ")
#k,d = int(c[0]), int(c[1])
#reads = b[1:]
#print(StringReconstructionReadPairs(k, d, reads))
