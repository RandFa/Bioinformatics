def ProbabilityHiddenPath(path,states, transitionm):
    """
    Input: A hidden path π followed by the states States and transition matrix Transition of an HMM (Σ, States, Transition, Emission).
    Output: The probability of this path, Pr(π).
    """
    prob = 1/len(states)
    for i in range(len(path)-1):
        prob *= transitionm[(path[i],path[i+1])]
    return prob

# file2 =  open("hi2.txt").read()
# file = file2.split('\n')
# path = file[0]
# states = file[2].split(" ")
# transitionm = {}
# first = file[4].split("\t")[1:]
# rest = file[5:]
# for i in range(len(first)):
#     temp = rest[i].split("\t")
#     for j in range(len(first)):
#         transitionm[(temp[0],first[j])] = float(temp[j+1])
# print(path, states, transitionm)
# print(ProbabilityHiddenPath(path, states, transitionm))

def ProbabilityofanOutcomeGivenaHiddenPath(string, alphabet, path, states, emissionm):
    """
    Input: A string x, followed by the alphabet from which x was constructed, followed by a hidden path π, followed by the states States and emission matrix Emission of an HMM (Σ, States, Transition, Emission).
    Output: The conditional probability Pr(x|π) that x will be emitted given that the HMM follows the hidden path π.
    """
    prob = 1
    for i in range(len(string)):
        prob*= emissionm[(path[i],string[i])]
    return prob

# file2 =  open("hi2.txt").read()
# file = file2.split('\n')
# string = file[0]
# alphabet = file[2].split(' ')
# path = file[4]
# states = file[6].split(" ")
# emissionm = {}
# first = file[8].split("\t")[1:]
# rest = file[9:]
# for i in range(len(states)):
#     temp = rest[i].split("\t")
#     for j in range(len(alphabet)):
#         emissionm[(temp[0],alphabet[j])] = float(temp[j+1])
# print(string, alphabet, path, states, emissionm)
# print(ProbabilityofanOutcomeGivenaHiddenPath(string, alphabet, path, states, emissionm))


def ViterbialgorithmsolvingtheDecodingProblem(string, alphabet, states,transitionm, emissionm):
    """
    Input: A string x, followed by the alphabet from which x was constructed, followed by the states States, transition matrix Transition, and emission matrix Emission of an HMM (Σ, States, Transition, Emission).
    Output: A path that maximizes the (unconditional) probability Pr(x, π) over all possible paths π.
    """
    path = ""
    prob = 1
    probs = {}
    for i in range(len(string)):
        probs = {}
        for j in states:
            if path == "":
                probs[j] =  emissionm[(j, string[i])]*1/len(states)
            else:
                probs[j]= prob* emissionm[(j, string[i])]*transitionm[path[-1],j]
        mx = max(probs.values())
        for key in probs.keys():
            if probs[key] == mx:
                path+=key
                prob = mx
                break
    return path


# def ViterbialgorithmsolvingtheDecodingProblem(string, alphabet, states,transitionm, emissionm):
#     """
#     Input: A string x, followed by the alphabet from which x was constructed, followed by the states States, transition matrix Transition, and emission matrix Emission of an HMM (Σ, States, Transition, Emission).
#     Output: A path that maximizes the (unconditional) probability Pr(x, π) over all possible paths π.
#     """
#     prob = 1
#     probs = {}
#     for s in states:
#         probs[s] = [[1,1]]
#     for i in range(len(string)):
#         for j in states:
#             temp = []
#             for jj in range(len(states)):
#                 temp1 = probs[j][-1][jj] * emissionm[(states[jj], string[i])] * transitionm[j, states[jj]]
#                 temp.append(temp1)
#     path = ""
#     final = []
#     first = ""
#     max = -1000
#     for s in states:
#         if max(probs[s][-1]) > max:
#             max = max(probs[s][-1])
#
#     mx = final.index(max(final))
#     return path


# file2 =  open("hi2.txt").read()
# file = file2.split('\n')
# string = file[0]
# alphabet = file[2].split(' ')
# states = file[4].split(" ")
# emissionm = {}
# transitionm = {}
# first = file[6].split("\t")[1:]
# rest1 = file[7:7+len(states)+1]
# rest2 = file[7+len(states)+2:]
# for i in range(len(states)):
#     temp = rest1[i].split("\t")
#     for j in range(len(states)):
#         transitionm[(states[j],temp[0])] = float(temp[j+1])
# for i in range(len(states)):
#     temp = rest2[i].split("\t")
#     for j in range(len(alphabet)):
#         emissionm[(temp[0],alphabet[j])] = float(temp[j+1])
# print(string, alphabet, states, transitionm, emissionm)
# print(ViterbialgorithmsolvingtheDecodingProblem(string, alphabet, states,transitionm, emissionm))


def seedAlignment(alignment, threshold):
    """
    """
    k =len(alignment[0])
    seeda = []
    seednotseed = []
    for a in alignment:
        if a.count("-")/k < threshold:
            seeda.append(a)
            seednotseed.append(1)
        else:
            seednotseed.append(0)
    return seeda, seednotseed

# def profile(seedalignment, alphabet):
#     """
#     """
#     profile = []
#     for a in seedalignment:
#         c = []
#         for b in alphabet:
#             c.append(a.count(b)/len(a))
#         profile.append(c)
#     return profile
#
# def match(salignment):
#     """
#     """
#     match = {}
#     for i in range(len(salignment)):
#         if i != len(salignment)-1:
#             match["M"+str(i+1)] = ["M"+str(i+2),"I"+str(i+1), "D"+str(i+2)]
#             match["I"+str(i+1)] = ["M"+str(i+2), "I"+str(i+1), "D"+str(i+2)]
#             match["D"+str(i+1)] = ["M"+str(i+2), "D"+str(i+2),"I"+str(i+1)]
#         else:
#             match["M"+str(i+1)] = ["I"+str(i+1),"E"]
#             match["I"+str(i+1)] = ["I"+str(i+1),"E"]
#             match["D"+str(i+1)] = ["I"+str(i+1), "E"]
#     match["I"+str(0)] = ["M"+str(1), "I"+str(0)]
#     match["S"] =  ["M"+str(1), "I"+str(0), "D"+str(1)]
#     return match
def tprobability(alignment, salignment, seednotseed, alphabet):
    """
    """
    paths = []
    countt = {}
    counte = {}
    transmission = []
    emission = []
    transmissioncount={}
    emissioncount = {}
    nodes  = ["S","I0"]
    for i in range(3*len(salignment)+3):
        transmission.append([0]*(3*len(salignment)+3))
        emission.append([0]*len(alphabet))
    for i in range(len(salignment)):
        nodes.append("M"+str(i+1))
        nodes.append("D"+str(i+1))
        nodes.append("I"+str(i+1))
    nodes.append("E")
    for i in range(len(alignment[0])):
        m=1
        paths.append(["S"])
        for a in range(len(alignment)):
            if seednotseed[a] == 1:
                if alignment[a][i] != "-":
                    paths[i].append("M"+str(m))
                    if (paths[i][-1], alignment[a][i]) in counte.keys():
                        counte[(paths[i][-1], alignment[a][i])] = counte[(paths[i][-1], alignment[a][i])] +1
                    else:
                        counte[(paths[i][-1], alignment[a][i])] = 1
                    if paths[i][-1] in emissioncount:
                        emissioncount[paths[i][-1]]= emissioncount[paths[i][-1]]+1
                    else:
                        emissioncount[paths[i][-1]]=1
                else:
                    paths[i].append("D"+str(m))
                m+=1
            else:
                if alignment[a][i] != "-":
                    if paths[i][-1][0] == "I":
                        paths[i].append(paths[i][-1])
                    else:
                        paths[i].append("I"+str(m-1))
                    if (paths[i][-1], alignment[a][i]) in counte.keys():
                        counte[(paths[i][-1], alignment[a][i])] = counte[(paths[i][-1], alignment[a][i])] +1
                    else:
                        counte[(paths[i][-1], alignment[a][i])] = 1
                    if paths[i][-1] in emissioncount:
                        emissioncount[paths[i][-1]]= emissioncount[paths[i][-1]]+1
                    else:
                        emissioncount[paths[i][-1]]=1
            if (paths[i][-2],paths[i][-1]) in countt.keys():
                countt[(paths[i][-2],paths[i][-1])] =countt[(paths[i][-2],paths[i][-1])]+1
            else:
                countt[(paths[i][-2],paths[i][-1])] =1
            if paths[i][-2] in transmissioncount:
                transmissioncount[paths[i][-2]]= transmissioncount[paths[i][-2]]+1
            else:
                transmissioncount[paths[i][-2]]=1
        paths[i].append("E")
        if (paths[i][-2],paths[i][-1]) in countt.keys():
            countt[(paths[i][-2],paths[i][-1])] =countt[(paths[i][-2],paths[i][-1])]+1
        else:
            countt[(paths[i][-2],paths[i][-1])] =1
        if paths[i][-2] in transmissioncount:
            transmissioncount[paths[i][-2]]= transmissioncount[paths[i][-2]]+1
        else:
            transmissioncount[paths[i][-2]]=1
    for key,value in countt.items():
        a1 = transmission[:nodes.index(key[0])]
        a3= transmission[nodes.index(key[0])+1:]
        a2 =[transmission[nodes.index(key[0])][:nodes.index(key[1])]+[value/transmissioncount[key[0]],3]+transmission[nodes.index(key[0])][nodes.index(key[1])+1:]]
        transmission = a1+a2+a3
    for key,value in counte.items():
        a1= emission[:nodes.index(key[0])]
        a3= emission[nodes.index(key[0])+1:]
        a2= [emission[nodes.index(key[0])][:alphabet.index(key[1])]+[value/emissioncount[key[0]],3]+emission[nodes.index(key[0])][alphabet.index(key[1])+1:]]
        emission =a1+a2+a3
    return transmission, emission, nodes


#
# file2 =  open("hi2.txt").read()
# file = file2.split('\n')
# threshold = float(file[0])
# alphabet = file[2].split(" ")
# alignment = []
# for i in range(len(file[4])):
#     alignment.append([])
#     for f in file[4:]:
#         alignment[i].append(f[i])
# print(threshold, alphabet, alignment)
# print(seedAlignment(alignment,threshold))
# seedAlignment, seednotseed = seedAlignment(alignment,threshold)
# a,b,c =tprobability(alignment,seedAlignment,seednotseed, alphabet)
# f= open("guru99.txt","w+")
# c = [" "]+c
# alphabet = [" "]+alphabet
# for word in c:
#     f.write(str(word) + ' ')
# f.write('\n')
# for i in range(len(a)):
#     d = [c[i+1]]+a[i]
#     for word in d:
#         f.write(str(word) + ' ')
#     f.write('\n')
# f.write("--------")
# f.write('\n')
# for word in alphabet:
#     f.write(str(word) + ' ')
# f.write('\n')
# for i in range(len(b)):
#     d = [c[i+1]]+b[i]
#     for word in d:
#         f.write(str(word) + ' ')
#     f.write('\n')


# print(profile(seedAlignment, alphabet))
# print(match(seedAlignment))


def HMMParameterEstimation(string, alphabet, path, states):
    """
    Input: A string x of symbols emitted from an HMM, followed by the HMM's alphabet Σ, followed by a path π, followed by the collection of states of the HMM.
    Output: A transition matrix Transition followed by an emission matrix Emission that maximize Pr(x, π) over all possible transition and emission matrices.
    """
    countt = {}
    counte = {}
    totalt = []
    totale = []
    transmission= []
    emission = []
    for i in range(len(states)):
        totalt.append(0)
        totale.append(0)
        transmission.append([0]*len(states))
        emission.append([0]*len(alphabet))
    for i in range(len(string)):
        if (path[i],string[i]) in countt.keys():
            countt[(path[i],string[i])] = countt[(path[i],string[i])]+1
        else:
            countt[(path[i],string[i])] = 1
        totalt = totalt[:states.index(path[i])] + [totalt[states.index(path[i])]+1] + totalt[states.index(path[i])+1:]
        if i != len(string)-1:
            if (path[i],path[i+1]) in counte.keys():
                counte[path[i],path[i+1]] = counte[path[i],path[i+1]]+1
            else:
                counte[path[i],path[i+1]] = 1
            totale = totale[:states.index(path[i])] + [totale[states.index(path[i])]+1] + totale[states.index(path[i])+1:]
    for key,value in counte.items():
        a1 = transmission[:states.index(key[0])]
        a3= transmission[states.index(key[0])+1:]
        a2 =[transmission[states.index(key[0])][:states.index(key[1])]+[round(value/totale[states.index(key[0])],3)]+transmission[states.index(key[0])][states.index(key[1])+1:]]
        transmission = a1+a2+a3
    for key,value in countt.items():
        a1= emission[:states.index(key[0])]
        a3= emission[states.index(key[0])+1:]
        a2= [emission[states.index(key[0])][:alphabet.index(key[1])]+[round(value/totalt[states.index(key[0])],3)]+emission[states.index(key[0])][alphabet.index(key[1])+1:]]
        emission =a1+a2+a3
    for i in range(len(transmission)):
        if transmission[i] == [0]*len(states):
            a1 = transmission[:i]
            a3= transmission[i+1:]
            a2 = [[round(1/len(transmission[i]),3)]* len(states)]
            transmission = a1+a2+a3
        if emission[i] == [0]*len(alphabet):
            a1 = emission[:i]
            a3= emission[i+1:]
            a2 = [[round(1/len(emission[i]),3)]* len(alphabet)]
            emission = a1+a2+a3
    return transmission, emission

# file2 =  open("hi2.txt").read()
# file = file2.split('\n')
# string = file[0]
# alphabet = file[2].split(' ')
# path = file[4]
# states = file[6].split(" ")
transmission, emission =HMMParameterEstimation("HTH",["H","T"],"BBB",["B","F"])
print(transmission, emission)
# states= [" "]+states
# alphabet = [" "]+alphabet
# print(*states, sep= "\t")
# for i in range(len(transmission)):
#     a = [states[i+1]]+transmission[i]
#     print(*a, sep= "\t")
# print("--------")
# print(*alphabet, sep= "\t")
# for i in range(len(emission)):
#     a = [states[i+1]]+emission[i]
#     print(*a, sep= "\t")



def ViterbilearningforestimatingtheparametersofanHMM(n, string, alphabet, states,transitionm, emissionm):
    """
    Input: A number of iterations j, followed by a string x of symbols emitted by an HMM, followed by the HMM's alphabet Σ, followed by the HMM's states, followed by initial transition and emission matrices for the HMM.
    Output: Emission and transition matrices resulting from applying Viterbi learning for j iterations.
    """
    for nn in range(n):
        path = ViterbialgorithmsolvingtheDecodingProblem(string, alphabet, states,transitionm, emissionm)
        print(path)
        transitiontemp, emissiontemp = HMMParameterEstimation(string, alphabet, path, states)
        print(transitiontemp, emissiontemp)
        transitionm = {}
        emissionm = {}
        for iii in range(len(states)):
            for jjj in range(len(states)):
                transitionm[(states[iii], states[jjj])] = transitiontemp[iii][jjj]
        for iii in range(len(states)):
            for jjj in range(len(alphabet)):
                emissionm[(states[iii], alphabet[jjj])] = emissiontemp[iii][jjj]
        print(transitionm, emissionm)
    return transitiontemp, emissiontemp

file2 =  open("hi2.txt").read()
file = file2.split('\n')
emissionm = {}
transitionm = {}
n = int(file[0])
string = file[2]
alphabet = file[4].split(' ')
states = file[6].split(" ")
rest1 = file[9:9+len(states)]
rest2 = file[9+len(states)+2:]
for i in range(len(states)):
    temp = rest1[i].split("\t")
    for j in range(len(states)):
        transitionm[(states[j],temp[0])] = float(temp[j+1])
for i in range(len(states)):
    temp = rest2[i].split("\t")
    for j in range(len(alphabet)):
        emissionm[(temp[0],alphabet[j])] = float(temp[j+1])
# transmission, emission =ViterbilearningforestimatingtheparametersofanHMM(n, string, alphabet, states,transitionm, emissionm)
# states= [" "]+states
# alphabet = [" "]+alphabet
# print(*states, sep= "\t")
# for i in range(len(transmission)):
#     a = [states[i+1]]+transmission[i]
#     print(*a, sep= "\t")
# print("--------")
# print(*alphabet, sep= "\t")
# for i in range(len(emission)):
#     a = [states[i+1]]+emission[i]
#     print(*a, sep= "\t")
