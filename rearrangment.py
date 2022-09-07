import copy
def GreedySorting(p):
    """
    Input: A permutation P.
    Output: The sequence of permutations corresponding to applying GreedySorting to P, ending with the identity permutation.
    """
    steps= []
    for i in range(len(p)):
        if p[i] != (i+1):
            if p[i]!= (-i-1):
                if i+1 in p:
                    a = p.index(i+1)
                else:
                    a = p.index(-i-1)
                b = p[i:a+1]
                c= b[::-1]
                d = [j * -1 for j in c]
                p = p[:i] + d +p[a+1:]
                steps.append(p)
            if p[i] != (i+1):
                p = p[:i]+ [i+1]+ p[i+1:]
                steps.append(p)
    return steps
# file2 =  open("hi.txt").read()
# pa = file2.split(" ")
# pa3  = []
#
# for p1 in pa:
#     pa3.append(int(p1))
# re =GreedySorting(pa3)
# print(len(re))
# for r in re:
#     r2 = []
#     for element in r:
#         if element>0:
#             r2.append("+"+str(element))
#         else:
#             r2.append(str(element))
#     print(*r2, sep=' ')
def NumberofBreakpoints(permutation):
    """
    Input: A permutation.
    Output: The number of breakpoints in this permutation.
    """
    adjacencies = 0
    permutation= [0]+permutation+[len(permutation)+1]
    for i in range(len(permutation)-1):
        if permutation[i+1]-permutation[i] == 1:
            adjacencies +=1
    breakpoints = len(permutation) - 1-adjacencies
    return breakpoints
# file2 =  open("hi.txt").read()
# pa = file2.split(" ")
# pa3  = []
# for p1 in pa:
#     pa3.append(int(p1))
# print(NumberofBreakpoints(pa3))
def ChromosomeToCycle(chromosome):
    """
    Input: A chromosome Chromosome containing n synteny blocks.
    Output: The sequence Nodes of integers between 1 and 2n resulting from applying ChromosomeToCycle to Chromosome.
    """
    nodes = []
    for c in chromosome:
        if c > 0:
            nodes.extend([c*2-1, c*2])
        elif c <0:
            nodes.extend([c*-2, c*-2-1])
    return nodes
print(ChromosomeToCycle((+11, +3, +13, +7, -2, +4, +12, +5, +9, +8, -14, +10, -1, +6)))
print(ChromosomeToCycle((-3, +13, -4, +6, -11, +1,+7, +14, -9, -12, -10, +8, +5,+2)))
def CycleToChromosome(nodes):
    """
    Input: A sequence Nodes of integers between 1 and 2n.
    Output: The chromosome Chromosome containing n synteny blocks resulting from applying CycleToChromosome to Nodes.
    """
    print(nodes)
    Chromosome = []
    for j in range(1, len(nodes)//2+1):
        if nodes[2*j-2] < nodes[2*j-1]:
            Chromosome.append(nodes[2*j-2]//2)
        else:
            Chromosome.append(nodes[2*j-2]//-2)
    return Chromosome
def CycleToChromosome2(nodes):
    """
    Input: A sequence Nodes of integers between 1 and 2n.
    Output: The chromosome Chromosome containing n synteny blocks resulting from applying CycleToChromosome to Nodes.
    """
    print(nodes)
    chromosome= []
    for node in range(0, len(nodes), 2):
        if nodes[node+1] < nodes[node]:
            chromosome.append((nodes[node])//-2)
        else:
            chromosome.append((nodes[node]+1)//2)
    return chromosome
# r2 = []
# for r in pa4:
#     if r>0:
#         r2.append("+"+str(r))
#     else:
#         r2.append(str(r))
# print(*r2, sep=' ')
def ColoredEdges(p):
    """
    Input: A genome P.
    Output: The collection of colored edges in the genome graph of P in the form (x, y).
    """
    chromosomes = []
    colorededges = []
    for ps in p:
        chromosome = ChromosomeToCycle(ps)
        edge = []
        for c in range(1, len(chromosome)-1, 2):
            edge.append([chromosome[c], chromosome[c+1]])
        edge.append([chromosome[len(chromosome)-1], chromosome[0]])
        colorededges.extend(edge)
    return colorededges

# file2 =  open("hi.txt").read()
# pa = file2.split(")")
# pa3 = []
# for p in pa:
#     p = p[1:]
#     p2 = p.split(" ")
#     pa2 = []
#     for p1 in p2:
#       pa2.append(int(p1))
#     pa3.append(pa2)
# pa4  = ColoredEdges(pa3)
# print(*pa4, sep = ", ")

def GraphToGenome(genomegraph):
    """
    Input: The colored edges ColoredEdges of a genome graph.
    Output: The genome P corresponding to this genome graph.
    """
    cycles = []
    while genomegraph != []:
        start = genomegraph[0][0]
        if start %2 == 0:
            end = start-1
        else:
            end = start +1
        cycle = genomegraph[0]
        oldgenome = genomegraph
        for gen in range(1, len(genomegraph)):
            # if genomegraph[gen][1] == start+1 or genomegraph[gen][1] == start-1:
            if end in genomegraph[gen]:
                if end == genomegraph[gen][1]:
                    cycle.extend(genomegraph[gen])
                else:
                    cycle.extend([genomegraph[gen][1], genomegraph[gen][0]])
                genomegraph= genomegraph[gen+1:]
                break
            # elif genomegraph[gen][0] == start+1 or genomegraph[gen][0] == start-1:
            #     cycle.extend([genomegraph[gen][1], genomegraph[gen][0]])
            #     genomegraph= genomegraph[gen+1:]
            #     break
            else:
                cycle.extend(genomegraph[gen])
        cycles.append(cycle)
        newgenome = genomegraph
        if oldgenome == newgenome:
            genomegraph = []
    p = []
    for cycle in cycles:
        cycle = [cycle[len(cycle)-1]] + cycle[:len(cycle)-1]
        p.append(CycleToChromosome2(cycle))
    return p

def GraphToGenome3(genomegraph):
    """
    Input: The colored edges ColoredEdges of a genome graph.
    Output: The genome P corresponding to this genome graph.
    """
    cycles = []
    while genomegraph != []:
        start = genomegraph[0][0]
        end = 0
        cycle = [genomegraph[0]]
        for gen in range(len(genomegraph)-1, 1, -1):
            if genomegraph[gen][1] == start+1 or genomegraph[gen][1] == start-1:
                if gen == len(genomegraph)-1:
                    cycle.extend(genomegraph[1:])
                    genomegraph= []
                else:
                    cycle.extend(genomegraph[1:gen+1])
                    cycle = [cycle[len(cycle)-1]] + cycle[:len(cycle)-1]
                    genomegraph= genomegraph[gen+1:]
                break
        cycles.append(cycle)
    p = []
    for cycle in cycles:
        cycle2 = []
        for node in cycle:
            cycle2.extend(node)
        p.append(CycleToChromosome(cycle2))
    return p
def GraphToGenome2(genomegraph):
    """
    Input: The colored edges ColoredEdges of a genome graph.
    Output: The genome P corresponding to this genome graph.
    """
    oldcycles = 0
    cycles = 0
    while genomegraph != []:
        start = genomegraph[0][0]
        end = 0
        cycle = genomegraph[0]
        for gen in range(len(genomegraph)-1, 1, -1):
            if genomegraph[gen][1] == start+1 or genomegraph[gen][1] == start-1 or genomegraph[gen][1] == start:
                genomegraph= genomegraph[:gen]
                cycles+=1
                break
        print(cycles)
        if oldcycles == cycles:
            cycles +=1
            return cycles
        else:
            oldcycles = cycles

# file2 =  open("hi2.txt").read()
# pa = file2.split("), ")
# pa3 = []
# for p in pa:
#     p = p[1:]
#     p2 = p.split(", ")
#     pa2 = []
#     for p1 in p2:
#       pa2.append(int(p1))
#     pa3.append(pa2)
# pa4  = GraphToGenome(pa3)
# for p4 in pa4:
#     r2 = []
#     for element in p4:
#         if element>0:
#             r2.append("+"+str(element))
#         else:
#             r2.append(str(element))
#     print(*r2, sep = " ")
def BreakDistance(p, q):
    """
    Input: Genomes P and Q.
    Output: The 2-break distance d(P, Q).
    """
    p.extend(q)
    maxi = 0
    for a in p:
        b = max(a)
        if b > maxi:
            maxi = b
    edge =ColoredEdges(p)
    cycles = GraphToGenome2(edge)
    return maxi-cycles
print(BreakDistance([+11, +3, +13, +7, -2, +4, +12, +5, +9, +8, -14, +10, -1, +6], [-3, +13, -4, +6, -11, +1, +7, +14, -9, -12, -10, +8, +5, +2]))
# file2 =  open("hi.txt").read()
# file = file2.split('\n')
# pa = file[0].split(")")
# pa4 = file[1].split(")")
# p = []
# q = []
# for p3 in pa:
#     p3 = p3[1:]
#     p2 = p3.split(" ")
#     pa2 = []
#     for p1 in p2:
#       pa2.append(int(p1))
#     p.append(pa2)
# for p3 in pa4:
#     p3 = p3[1:]
#     p2 = p3.split(" ")
#     pa2 = []
#     for p1 in p2:
#       pa2.append(int(p1))
#     q.append(pa2)
# print(BreakDistance(p,q))
def BreakOnGenomeGraph(GenomeGraph, i1 , i2 , i3 , i4):
    """
    Input: The colored edges of a genome graph GenomeGraph, followed by indices i1 , i2 , i3 , and i4 .
    Output: The colored edges of the genome graph resulting from applying the 2-break operation 2-BreakOnGenomeGraph(GenomeGraph, i1 , i2 , i3 , i4 ).
    """
    print(GenomeGraph)
    for gen in range(len(GenomeGraph)):
        if i1 in GenomeGraph[gen]:
            i1num = gen
            if i1 - GenomeGraph[gen-1][1] == 1 or i1 - GenomeGraph[gen-1][1] == -1:
                node = [i1, i3]
            else:
                node = [i3, i1]
        if i3 in GenomeGraph[gen]:
            i2num = gen
            if i2 - GenomeGraph[gen-1][1] == 1 or i2 - GenomeGraph[gen-1][1] == -1:
                node2 = [i2, i4]
            else:
                node2 = [i4, i2]
            # if GenomeGraph[gen][0]== i2:
            #     node2 = [i4, i2]
            # else:
            #     node2 = [i2, i4]
        # if GenomeGraph[gen][1] == back:
        #     node = [i1, i3]
        # if GenomeGraph[gen][1] == back2:
        #     node2 = [i2, i4]
    if i1num < i2num:
        # GenomeGraph2 = GenomeGraph[:i1num] + GenomeGraph[i1num+1:i2num]+ GenomeGraph[i2num+1:]+[[i1, i3]] + [[i2, i4]]
        #GenomeGraph2 = GenomeGraph[:i1num] + GenomeGraph[i2num+1:] + [node] + GenomeGraph[i1num+1:i2num] + [node2]
        GenomeGraph2 = GenomeGraph[:i1num] +[node] +GenomeGraph[i2num+1:] + GenomeGraph[i1num+1:i2num] + [node2]
    else:
        # GenomeGraph2 = GenomeGraph[:i2num]  + GenomeGraph[i2num+1:i1num] + GenomeGraph[i1num+1:] + [[i1, i3]] +[[i2, i4]]
        #GenomeGraph2 = GenomeGraph[:i2num] + GenomeGraph[i1num+1:]+ [node2] + GenomeGraph[i2num+1:i1num] +[node]
        GenomeGraph2 = GenomeGraph[:i2num] + [node] +GenomeGraph[i1num+1:] + GenomeGraph[i2num+1:i1num] +[node2]
    return GenomeGraph2
# file2 =  open("hi.txt").read()
# file = file2.split('\n')
# pa = file[0].split("), ")
# p = []
# for p3 in pa:
#     p3 = p3[1:]
#     p2 = p3.split(", ")
#     pa2 = []
#     for p1 in p2:
#       pa2.append(int(p1))
#     p.append(pa2)
# pa4 = file[1].split(", ")
# print(p, pa4)
# print(BreakOnGenomeGraph(p, int(pa4[0]), int(pa4[1]), int(pa4[2]), int(pa4[3])))
def BreakOnGenome(p, i1 , i2 , i3 , i4):
    """
    Input: A genome P, followed by indices i1 , i2 , i3 , and i4 .
    Output: The genome P' resulting from applying the 2-break operation 2-BreakOnGenome(GenomeGraph i1 , i2 , i3 , i4 ).
    """
    chromosome = ChromosomeToCycle(p)
    edge = []
    for c in range(1, len(chromosome)-1, 2):
        edge.append([chromosome[c], chromosome[c+1]])
    edge.append([chromosome[len(chromosome)-1], chromosome[0]])
    edges2 = BreakOnGenomeGraph(edge, i1 , i2 , i3 , i4)
    geaph = GraphToGenome(edges2)
    return geaph
# file2 =  open("hi2.txt").read()
# file = file2.split('\n')
# pa = file[0].split(" ")
# p = [int(pa[0][1:])]
# for a in pa[1:]:
#     p.append(int(a))
# result = BreakOnGenome(p, 36, 82, 76, 9)
# for res in result:
#     a = "("
#     for r  in res:
#         if r > 0:
#             a+="+"
#             a+= str(r)
#             a+= " "
#         else:
#             a+= str(r)
#             a+= " "
#     a += ")"
#     print(a)
    # print(*res, sep = " ")

# def ShortestRearrangementScenario(P, Q):
#     BlueEdges ← ColoredEdges(Q)
#      while BreakpointGraph has a non-trivial cycle Cycle
#           RedEdges = ColoredEdges(P)
#           BreakpointGraph = the graph formed by RedEdges and BlueEdges
#           (i, j) ← an arbitrary edge from BlueEdges in a nontrivial red-blue cycle from BreakpointGraph
#           (i, i') ← an edge from RedEdges originating at node i
#          (j, j') ← an edge from RedEdges originating at node j
#           BreakpointGraph ← the graph formed by RedEdges and BlueEdges
#           P ← 2-BreakOnGenome(P, i, i′, j, j′)
#           output P


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
def Sharedkmers(k, kmer1, kmer2):
    """
    Input: An integer k and two strings.
    Output: All k-mers shared by these strings, in the form of ordered pairs (x, y) corresponding to starting positions of these k-mers in the respective strings.
    """
    result = []
    pattern = {}
    patterns = {}
    for i in range(len(kmer1)-k+1):
        s = kmer1[i:i+k]
        if s not in pattern:
            pattern[s] = [[i]]
        else:
            pattern[s][0].append(i)
        a = ReverseComplement(kmer1[i:i+k])
        if a not in patterns:
            patterns[a] = [[i]]
        else:
            patterns[a][0].append(i)
    for j in range(len(kmer2)-k+1):
        ks1= kmer2[j:j+k]
        if ks1 in list(pattern.keys()):
            pattern[ks1].append(j)
        elif ks1 in list(patterns.keys()):
            patterns[ks1].append(j)
    for pat in list(pattern.keys()):
        if len(pattern[pat]) > 1:
            for l in range(1, len(pattern[pat])):
                for m in pattern[pat][0]:
                    result.append((m, pattern[pat][l]))
    for pat in list(patterns.keys()):
        if len(patterns[pat]) > 1:
            for l in range(1, len(patterns[pat])):
                for m in patterns[pat][0]:
                    result.append((m, patterns[pat][l]))
    result.sort()
    return result
#
#


# def two_break_sorting(P,Q):
#     '''
#     CODE CHALLENGE: Solve the 2-Break Sorting Problem.
#     2-Break Sorting Problem: Find a shortest transformation
#     of one genome into another via 2-breaks.
#     Input: Two genomes with circular chromosomes on the same
#     set of synteny blocks.
#     Output: The collection of genomes resulting from applying
#     a shortest sequence of 2-breaks transforming one genome into the other.
#     '''
#     red = colored_edges(Q)
#     path = [P]
#     while two_break_distance(P,Q) > 0:
#         cycles = colored_edges_cycles(colored_edges(P),red)
#         for c in cycles:
#             if len(c) >= 4:
#                 P = two_break_on_genome(P,c[0],c[1],c[3],c[2])
#                 path.append(P)
#                 break
#     return path
