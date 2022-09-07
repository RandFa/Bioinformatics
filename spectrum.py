mass = {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114, 'D': 115, 'K': 128, 'Q': 128, 'E': 129, 'M': 131, 'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186}
mass2  = {57: ['G'], 71: ['A'], 87: ['S'], 97: ['P'], 99: ['V'], 101: ['T'], 103: ['C'], 113: ['L',"I"], 114: ['N'], 115: ['D'], 128: ['Q',"K"], 129: ['E'], 131: ['M'], 137: ['H'], 147: ['F'], 156: ['R'], 163: ['Y'], 186: ['W']}
aminoAcid = ['G', 'A', 'S', 'P', 'V', 'T', 'C',"I", 'L', 'N', 'D',"K", 'Q', 'E', 'M', 'H', 'F', 'R', 'Y', 'W']
mass4  = {57: 'G', 71: 'A', 87: 'S', 97: 'P', 99: 'V', 101: 'T', 103: 'C', 113: 'L', 114: 'N', 115: 'D', 128: 'Q', 129: 'E', 131: 'M', 137: 'H', 147: 'F', 156: 'R', 163: 'Y', 186: 'W'}

def spectgrapg(spec):
    graph = []
    for ele in range(len(spec)):
        for ele2 in spec[ele:]:
            if ele2 - spec[ele] in mass2.keys():
                graph.append([spec[ele],ele2, mass2[ele2 - spec[ele]][0]])
    return sorted(graph)

# file2 =  open("hi2.txt").read()
# file = file2.split(' ')
# spec = [0]
# for f in file:
#     spec.append(int(f))
# result = spectgrapg(spec)
# for r in result:
#     print(str(r[0])+"->"+str(r[1])+":"+r[2])

def nodes(graph):
    nodes = {}
    for edge in graph:
        if edge[1] in nodes.keys():
            nodes[edge[1]].append(edge[0])
        else:
            nodes[edge[1]] = [edge[0]]
    return nodes
def pathes(nodes, graph):
    path = [graph[-1][1]]
    source = graph[0][0]
    sink = graph[-1][1]
    while source != sink and sink in nodes.keys():
        if len(nodes[sink])>1:
            sink2 = nodes[sink][0]
            nodes[sink] = nodes[sink][1:]
        else:
            sink2 = nodes[sink][0]
            nodes.pop(sink)
        path.append(sink)
        sink= sink2
    path.append(sink)
    if path[-1] == 0:
        return path
    else:
        return pathes(nodes, graph)


def DecodingIdealSpectrum(spec, graph):
     mynodes = nodes(graph)
     paths = pathes(mynodes, graph)
     paths = paths[::-1]
     str = ""
     for i in range(len(paths)-2):
         str += mass2[paths[i+1]-paths[i]][0]
     return str
# file2 =  open("hi2.txt").read()
# file = file2.split(' ')
# spec = [0]
# for f in file:
#     spec.append(int(f))
# graphbi = spectgrapg(spec)
# print(DecodingIdealSpectrum(spec,graphbi))

def peptidetovetor(peptide):
    vector = []
    for i in range(len(peptide)):
        le = mass[peptide[i]]
        for i in range(le-1):
            vector.append(0)
        vector.append(1)
    return vector
# a= peptidetovetor("EIWYKKELSMKMEPSFDGDEFVYMHRYYITLRL")
# print(*a)
def vectortopeptide(vector):
    peptide = ""
    score = 0
    for i in vector:
        if i==0:
            score +=1
        else:
            peptide += mass2[score+1][0]
            score=0
    return peptide



# file2 =  open("hi2.txt").read()
# file = file2.split(' ')
# vector = []
# for f in file:
#     vector.append(int(f))
# print(vectortopeptide(vector))

mass5 = {4:'X', 5:"Z"}
def PeptideSequencing(spectrum):
    """
     Input: A spectral vector Spectrum'.
     Output: An amino acid string Peptide that maximizes Score(Peptide', Spectrum') among all possible amino acid strings.
    """
    backward = {}
    visited = []
    paths = [[spectrum[-1],len(spectrum)-1]]
    for i in mass4.keys():
        if i < len(spectrum):
            backward[i-1] = [0]
    while len(spectrum)-1 not in visited:
        backwardtemp = {}
        for back in backward.keys():
            if back not in visited:
                visited.append(back)
                for i in mass4.keys():
                    if back+i < len(spectrum):
                        if back+i not in backwardtemp:
                            backwardtemp[back+i] = [back]
                        elif back+i in backwardtemp and back not in backwardtemp[back+i]:
                            backwardtemp[back+i].append(back)
        for tem in backwardtemp.keys():
            if tem in backward.keys():
                backward[tem].extend(backwardtemp[tem])
            else:
                backward[tem] = backwardtemp[tem]
    finals = [10000000]
    while max(finals)!= 0:
        finals= []
        pathtemp = []
        for path in paths:
            if path[-1] != 0:
                for i in backward[path[-1]]:
                    path2 = [path[0]+spectrum[i]]+ path[1:]+[i]
                    pathtemp.append(path2)
                    finals.append(i)
            else:
                pathtemp.append(path)
        paths = pathtemp
    paths.sort()
    print(paths)
    win = paths[-1][1:]
    win.reverse()
    pep = ""
    for i in range(len(win)-1):
        if i != 0:
            pep+= mass4[win[i+1]-win[i]]
        else:
            pep+= mass4[win[i+1]-win[i]+1]
    return pep
    # nodesv = {}
    # nodes = []
    # finals = []
    # edges = []
    # backtrack
    # for s in range(len(spectrum)):
    #     nodesv[s] = spectrum[s]
    #     nodes.append(s)
    # paths=[[nodesv[0], 0]]
    # arb = nodes[-1]
    # for i in range(arb):
    #     temp= []
    #     for path in paths:
    #         source = path[-1]
    #         if source != len(nodes)-1:
    #             for n in mass4.keys():
    #                 if n+source in nodes:
    #                     temp.append([path[1]+nodesv[source+n]]+path[1:]+[source+n])
    #         else:
    #             finals.append(path)
    #     paths = temp
    # finals.sort()
    # finalfinal= ""
    # final = [finals[-1][0]]+finals[-1][1:]
    # for f in range(len(final)-2):
    #     finalfinal += mass4[final[f+2]-final[f+1]]
    # return finalfinal


file2 =  open("hi2.txt").read()
file = file2.split(' ')
spectrum = [0]
for f in file:
    spectrum.append(int(f))
print(PeptideSequencing(spectrum))
