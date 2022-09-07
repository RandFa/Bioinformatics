import copy
def leafinter(matrix):
    edges  = list(matrix.keys())
    matrix1= sorted(edges)
    count= 1
    leaves = []
    internal = []
    for m in range(1, len(matrix1)):
        if matrix1[m][0] == matrix1[m-1][0]:
            count+=1
            a = matrix1[m][0]
            if a not in internal:
                internal.append(a)
        else:
            if count == 1:
                leaves.append(matrix1[m-1][0])
            count = 1
    return leaves, internal
def interroute(internal, matrix):
    '''
    not working use interroute2
    '''
    # routes = [internal[0]]
    routes = []
    route = [internal[0]]
    internal = internal[1:]
    while internal != []:
        before = len(route)
        for int in internal:
            if int[0] == route[len(route)-1]:
                route.append(int[1])
                internal.remove(int)
                break
            elif int[1] == route[len(route)-1]:
                route.append(int[0])
                internal.remove(int)
                break
            elif int[1] == route[0]:
                route= [int[0]]+ route
                internal.remove(int)
                break
            elif int[0] == route[0]:
                route= [int[1]]+ route
                internal.remove(int)
                break
        after = len(route)
        if after== before:
            routes.pop(0)
            routes.append(route)
            if internal != []:
                routes = [internal[0]] +routes
                route = routes[0]
                internal = internal[1:]
    return routes
def interroute2(internal, matrix):
    # routes = [internal[0]]
    routes = []
    route = [internal[0]]
    internal = internal[1:]
    hi = list(matrix.keys())
    while internal != []:
        before = len(route)
        for int in internal:
            if (int,route[len(route)-1]) in hi:
                route.append(int)
                internal.remove(int)
                break
            elif (route[len(route)-1], int) in hi:
                route.append(int)
                internal.remove(int)
                break
            elif  (route[0], int) in hi:
                route= [int]+ route
                internal.remove(int)
                break
            elif (int,route[0]) in hi:
                route= [int]+ route
                internal.remove(int)
                break
        after = len(route)
        if internal == []:
            routes.append(route)
        if after == before:
            routes.append(route)
            if internal != []:
                route = [internal[0]]
                internal = internal[1:]
    return routes
def DistancesBetweenLeavesProblem(n, matrix):
    '''
    Input:  An integer n followed by the adjacency list of a weighted tree with n leaves.
    Output: An n x n matrix (di,j), where di,j is the length of the path between leaves i and j.
    '''
    leaves, internal = leafinter(matrix)
    print(internal)
    hi = list(matrix.keys())
    routes = interroute2(internal, matrix)
    print(leaves, routes)
    distmat = []
    for i in range(n):
        lis = []
        for j in range(n):
            score = 0
            if i != j:
                if i in leaves:
                    edge = [item for item in hi if item[0] == i]
                    score += matrix[edge[0]]
                    ii = edge[0][1]
                else:
                    ii = i
                if j in leaves:
                    edge = [item for item in hi if item[0] == j]
                    score += matrix[edge[0]]
                    jj = edge[0][1]
                else:
                    jj = j
                if ii != jj:
                    for r in routes:
                        if ii in r and jj in r:
                            f = r.index(ii)
                            l = r.index(jj)
                            for nn in range(min([f,l]), max([f, l])):
                                score += matrix[(r[nn], r[nn+1])]
            lis.append(score)
        distmat.append(lis)
    return distmat

# file2 =  open("hi2.txt").read()
# file = file2.split('\n')
# n = int(file[0])
# mat = {}
# mat1= file[1:]
# for edge in mat1:
#     edge2 = edge.split(':')
#     edge3 = edge2[0].split('->')
#     edge3 = [int(edge3[0]), int(edge3[1])]
#     edge3 = tuple(edge3)
#     mat[edge3] = int(edge2[1])
# solution = DistancesBetweenLeavesProblem(n, mat)
# for s in solution:
#     print(*s, sep = " ")
def LimbLengthProblem(n, j, matrix):
    '''
    Input: An integer n, followed by an integer j between 0 and n - 1, followed by a space-separated additive distance matrix D (whose elements are integers).
    Output: The limb length of the leaf in Tree(D) corresponding to row j of this distance matrix (use 0-based indexing).
    '''
    lengths = []
    for i in range(n):
        for k in range(n):
            if k!=j and i!=j:
                leng = (matrix[i][j]+ matrix[j][k]- matrix[i][k])//2
                lengths.append(leng)
    return min(lengths)
# file2 =  open("hi2.txt").read()
# file = file2.split('\n')
# n = int(file[0])
# j = int(file[1])
# mat1= file[2:]
# mat = []
# for edge in mat1:
#     edge2 = edge.split(' ')
#     edge3 = []
#     for e in edge2:
#         edge3.append(int(e))
#     mat.append(edge3)
# print(LimbLengthProblem(n, j, mat))

def path(i, j, matrix):
    if (i, j) in list(matrix.keys()):
        return [i, j]
    hi = list(matrix.keys())
    leaves, internal = leafinter(matrix)
    if len(internal)>2:
        routes = interroute2(internal, matrix)
    else:
        routes = [internal]
    path = []
    path2 = []
    if i in leaves:
        edge = [item for item in matrix if item[0] == i]
        path.append(i)
        ii = edge[0][1]
    else:
        ii = i
    if j in leaves:
        edge = [item for item in matrix if item[0] == j]
        path2.append(j)
        jj = edge[0][1]
    else:
        jj = j
    if ii != jj :
        for r in routes:
            if ii in r and jj in r:
                f = r.index(ii)
                l = r.index(jj)
                m, n = min([f,l]), max([f, l])
                if m == f:
                    path.extend(r[m:n+1])
                else:
                    rr = r[m:n+1]
                    rr.reverse()
                    path.extend(rr)
    else:
        path.append(ii)
    if path2 != []:
        path.extend(path2)
    print('path', path)
    return path
def findpath(D, x,i, k):
    path1 = path(i, k, D)
    start = 0
    node = path1[start]
    sum = 0
    while sum < x:
        sum += D[(path1[start], path1[start+1])]
        start+= 1
        node = path1[start]
    if sum == x:
        return path1[start], D
    else:
        global globaln
        globaln += 1
        D[(globaln, path1[start])] = sum - x
        D[( path1[start], globaln)] = sum - x
        sum -= D[(path1[start], path1[start-1])]
        D[(globaln, path1[start-1])] = x - sum
        D[( path1[start-1], globaln)] = x -sum
        D.pop((path1[start-1], path1[start]))
        D.pop((path1[start], path1[start-1]))
        return globaln,D
def AdditivePhylogeny(matrix):
    '''
    Input: An integer n followed by a space-separated n x n distance matrix.
    Output: A weighted adjacency list for the simple tree fitting this matrix.
    '''
    n = len(matrix)
    if n == 2:
        return {(0, 1): matrix[0][1], (1,0): matrix[0][1]}
    limb = LimbLengthProblem(n, n-1, matrix)
    for m in range(n-1):
        matrix[m][n-1] = matrix[m][n-1]-limb
        matrix[n-1][m] = matrix[n-1][m]- limb
    for i in range(n-1):
        for k in range(n-1):
            if matrix[i][k] == matrix[i][n-1]+matrix[n-1][k]:
                node = (i,k)
                break
        if matrix[i][k] == matrix[i][n-1]+matrix[n-1][k]:
            break
    x = matrix[i][n-1]
    matrix.pop()
    for ma in matrix:
        ma.pop()
    T = AdditivePhylogeny(matrix)
    v, T = findpath(T, x, i, k)
    T[(n-1,v)] = limb
    T[(v,n-1)] = limb
    return T
# file2 =  open("hi2.txt").read()
# file = file2.split('\n')
# globaln = int(file[0]) -1
# mat1= file[1:]
# mat = []
# for edge in mat1:
#     edge2 = edge.split(' ')
#     edge3 = []
#     for e in edge2:
#         edge3.append(int(e))
#     mat.append(edge3)
# matrix1 = mat[:]
# solution = AdditivePhylogeny(mat)
# solutions = []
# for k, v in solution.items():
#     s= list(k)
#     s.append(v)
#     solutions.append(s)
# solutions.sort()
# for s in solutions:
#     print(str(s[0]) +'->' +str(s[1])+ ':' +str(s[2]))
def substituematrix(matrix1, node, new0, internode, memory):
    row = copy.deepcopy(matrix1)
    for k,v in list(matrix1.items()):
        if k != node and k!= (node[1], node[0]):
            if node[0] == k[0] and k[1] != node[0]:
                val = (matrix1[k[1],node[0]] + matrix1[k[1], node[1]])/2
                row[(new0, k[1])] = val
                row[(k[1], new0)] = val
                memory[k] = row[k]
                row.pop(k)
                memory[(k[1],k[0])] = row[(k[1],k[0])]
                row.pop((k[1],k[0]))
                memory[(k[1], node[1])] = row[(k[1], node[1])]
                row.pop((k[1], node[1]))
                if k[1] != node[1]:
                    memory[(node[1], k[1])] = row[(node[1], k[1])]
                    row.pop((node[1], k[1]))
    memory[node] = row[node]
    row.pop(node)
    if (node[1], node[0]) in list(row.keys()):
        memory[(node[1], node[0])] = row[(node[1], node[0])]
        row.pop((node[1], node[0]))
    return row, memory


# file2 =  open("hi2.txt").read()
# file = file2.split('\n')
# n = int(file[0])
# mat1= file[1:]
# mat = []
# print(n)
# for edge in mat1:
#     edge2 = edge.split('\t')
#     print(edge2)
#     edge3 = []
#     for e in edge2:
#         edge3.append(int(e))
#     mat.append(edge3)
# result = UPGMA(n, mat)
# result.sort()
# for r in result:
#     print(str(r[0])+'->'+str(r[1])+':'+str(round(r[2], 3)))



def UPGMA(D, n):
    '''
    Input: An integer n followed by a space separated n x n distance matrix.
    Output: An adjacency list for the ultrametric tree returned by UPGMA.
    Edge weights should be accurate to two decimal places
    (answers in the sample dataset below are provided to three decimal places).
    '''
    clusters = {}
    T= []
    T2={}
    age = {}
    for i in range(n):
        clusters[i]= []
        T.append(i)
        age[str(i)] = 0.000
    while len(clusters)> 1:
        pre = 10000000
        for blah in list(D.keys()):
            if D[blah]!= 0 and D[blah]<pre:
                pre= D[blah]
                aa= blah[0]
                bb= blah[1]
        age[str(n)] = round(pre/2,3)
        clusters[n]= clusters[aa]+ clusters[bb]+[aa]+[bb]
        T2[(n,aa)] = age[str(n)]- age[str(aa)]
        T2[(aa,n)] = age[str(n)]- age[str(aa)]
        T2[(n,bb)] = age[str(n)]- age[str(bb)]
        T2[(bb,n)] = age[str(n)]- age[str(bb)]
        T.remove(aa)
        T.remove(bb)
        for i in T:
            if len(clusters[aa]) == 0 and len(clusters[bb]) == 0:
                bla = round((D[(i,aa)]+D[(i,bb)])/2,3)
            elif len(clusters[aa]) != 0 and len(clusters[bb]) == 0:
                bla = round((D[(i,aa)]*len(clusters[aa])+ D[(i,bb)])/(len(clusters[aa])+1),3)
            elif len(clusters[aa]) == 0 and len(clusters[bb]) != 0:
                bla = round((D[(i,bb)]*len(clusters[bb])+ D[(i,aa)])/(len(clusters[bb])+1),3)
            else:
                bla = round((D[(i,bb)]*len(clusters[bb])+ D[(i,aa)]*len(clusters[aa]))/(len(clusters[bb])+len(clusters[aa])),3)
            D[(i,n)] = bla
            D[(n,i)]= bla
            D.pop((i,aa))
            D.pop((i,bb))
            D.pop((aa,i))
            D.pop((bb,i))
        clusters.pop(aa)
        clusters.pop(bb)
        T.append(n)
        n+=1
        D.pop((aa,bb))
        D.pop((bb,aa))
    return T2

# file2 =  open("hi2.txt").read()
# file = file2.split('\n')
# n = int(file[0])
# mat1= file[1:]
# mat = {}
# for edge in range(len(mat1)):
#     edge2 = mat1[edge].split(' ')
#     # edge2=edge2[:len(edge2)-1]
#     edge3 = []
#     for e in range(len(edge2)):
#         mat[(edge,e)] = int(edge2[e])
# result = UPGMA(mat,n)
# result2 = []
# for key, value in result.items():
#     result2.append([key,value])
# result2.sort()
# for r in result2:
#     print(str(r[0][0])+'->'+str(r[0][1])+':'+str(round(r[1], 3)))

def TotalDistance(D):
    TD = {}
    nli= []
    for ele in D.keys():
        if ele[0] not in nli:
            nli.append(ele[0])
    for nl in nli:
        dd = 0
        for key,value in D.items():
            if key[0]==nl:
                dd+= value
        TD[nl] = dd
    return TD
def NeighborMatrix(D, TD,n):
    D2 = {}
    for key, value in list(D.items()):
        i = key[0]
        j = key[1]
        if i == j:
            D2[key] =0
        else:
            D2[(i,j)]= round((n-2)*D[(i,j)]-TD[i]-TD[j],2)
            D2[(j,i)]= round((n-2)*D[(i,j)]-TD[i]-TD[j],2)
    return D2


def NeighborJoining(D,n,m):
    T ={}
    if n == 2:
        for key, value in D.items():
            if key[0]!=key[1]:
                T[key]=value
        return T
    TD = TotalDistance(D)
    D2 = NeighborMatrix(D,TD,m)
    pre = 1000000000
    for key,value in D2.items():
        if value < pre and key[0]!= key[1]:
            pre = value
            i = key[0]
            j=key[1]
    mtlt = (TD[i] - TD[j]) /(n - 2)
    limbLengthi = round((D[(i,j)] + mtlt)/2,2)
    limbLengthj = round((D[(i,j)] - mtlt)/2,2)
    temp = []
    tempd = copy.deepcopy(D)
    for key,value in tempd.items():
        if key[0] != i and key[0]!=j:
            D[(key[0],m)]=round((D[(key[0],i)] + D[(key[0],j)] - D[(i,j)])/2,2)
            D[(m,key[0])]=round((D[(key[0],i)] + D[(key[0],j)] - D[(i,j)])/2,2)
    for key,value in tempd.items():
        if key[0] == i or key[0]==j or key[1]==i or key[1]==j:
            D.pop(key)
    T = NeighborJoining(D,n-1, m+1)
    T[(i,m)]= limbLengthi
    T[(m,i)]=limbLengthi
    T[(j,m)]= limbLengthj
    T[(m,j)]= limbLengthj
    return T

def NeighborJoining2(D, n,m):
    '''
    Input: An integer n followed by a space separated n x n distance matrix.
    Output: An adjacency list for the ultrametric tree returned by UPGMA.
    Edge weights should be accurate to two decimal places
    (answers in the sample dataset below are provided to three decimal places).
    '''
    clusters = {}
    T= []
    T2={}
    age = {}

    for i in range(n):
        clusters[i]= []
        T.append(i)
        age[str(i)] = 0.000
    while len(clusters)> 2:
        TD = TotalDistance(D)
        D2 = NeighborMatrix(D,TD,m)
        pre = 10000000
        for blah in list(D2.keys()):
            if D2[blah]!= 0 and D2[blah]<=pre:
                pre= D2[blah]
                i= blah[0]
                j= blah[1]
        clusters[m]= clusters[i]+ clusters[j]+[i]+[j]
        mtlt = (TD[i] - TD[j]) /(n - 2)
        limbLengthi = round((D[(i,j)] + mtlt)/2,2)
        limbLengthj = round((D[(i,j)] - mtlt)/2,2)
        T2[(i,m)]= limbLengthi
        T2[(m,i)]=limbLengthi
        T2[(j,m)]= limbLengthj
        T2[(m,j)]= limbLengthj
        T.remove(i)
        T.remove(j)
        for k in T:
            D[(k,m)] = round((D[(k,i)] + D[(k,j)] - D[(i,j)])/2,2)
            D[(m,k)]= round((D[(k,i)] + D[(k,j)] - D[(i,j)])/2,2)
            D.pop((i,k))
            D.pop((k,j))
            D.pop((k,i))
            D.pop((j,k))
        clusters.pop(i)
        clusters.pop(j)
        T.append(m)
        n-=1
        m+=1
        D.pop((i,j))
        D.pop((j,i))
    i,j= list(clusters.keys())[0],list(clusters.keys())[1]
    T2[(i,j)]=D[(i,j)]
    T2[(j,i)]=D[(j,i)]
    return T2


# file2 =  open("hi2.txt").read()
# file = file2.split('\n')
# n = int(file[0])
# m=n
# mat1= file[1:]
# mat = {}
# for edge in range(len(mat1)):
#     edge2 = mat1[edge].split(' ')
#     # edge2=edge2[:len(edge2)-1]
#     edge3 = []
#     for e in range(len(edge2)):
#         mat[(edge,e)] = int(edge2[e])

# result = NeighborJoining2(mat,n,m)
# result2 = []
# for key, value in result.items():
#     result2.append([key,value])
# result2.sort()
# for r in result2:
#     print(str(r[0][0])+'->'+str(r[0][1])+':'+str(round(r[1],2)))


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


def SmallParsimony(t, tree, leaves, leavesv, nodes):
    '''
    Input: An integer n followed by an adjacency list for a rooted binary tree with n leaves labeled by DNA strings.
    Output: The minimum parsimony score of this tree, followed by the adjacency list of a tree corresponding to labeling
    internal nodes by DNA strings in order to minimize the parsimony score of the tree.  You may break ties however you like.
    A C G T
    '''
    finals = {}
    scoref = 0
    for l in range(nodes[-1]+1):
        finals[l]=""
    for i in range(len(leavesv[0])):
        m=t
        tags= {}
        values ={}
        score = {}
        for n in leaves:
            tags[n]=1
            values[n] = leavesv[n][i]
            if values[n] == "A":
                score[n] = [0,1,1,1]
            if values[n] == "C":
                score[n] = [1,0,1,1]
            if values[n] == "G":
                score[n] = [1,1,0,1]
            if values[n] == "T":
                score[n] = [1,1,1,0]
        for n in nodes:
            tags[n]=0
        while m in nodes:
            son =tree[m][0]
            daughter =tree[m][1]
            score[m]=[]
            for ba in range(4):
                add = 0
                if score[son][ba]== min(score[son]) and score[daughter][ba]== min(score[daughter]):
                    add+= min(score[son])+min(score[daughter])
                elif score[son][ba]!= min(score[son]) and score[daughter][ba]!= min(score[daughter]):
                    add+= min(score[son])+min(score[daughter])+2
                else:
                    add+= min(score[son])+min(score[daughter])+1
                score[m].append(add)
            m+=1
        scoref += min(score[m-1])
        if score[m-1][0] == min(score[m-1]):
            finals[m-1] += "A"
        elif score[m-1][1] == min(score[m-1]):
            finals[m-1] += "C"
        elif score[m-1][2] == min(score[m-1]):
            finals[m-1] += "G"
        elif score[m-1][3] == min(score[m-1]):
            finals[m-1] += "T"
        for l in range(m-2,0,-2):
            if score[l][0] == min(score[l]) and finals[m-1][-1]=="A":
                finals[l] += "A"
            elif score[l][1] == min(score[l]) and finals[m-1][-1]=="C":
                finals[l] += "C"
            elif score[l][3] == min(score[l])and finals[m-1][-1]=="T":
                finals[l] += "T"
            elif score[l][2] == min(score[l]) and finals[m-1][-1]=="G":
                finals[l] += "G"
            elif score[l][0] == min(score[l]):
                finals[l] += "A"
            elif score[l][1] == min(score[l]):
                finals[l] += "C"
            elif score[l][3] == min(score[l]):
                finals[l] += "T"
            elif score[l][2] == min(score[l]):
                finals[l] += "G"
            if score[l-1][0] == min(score[l-1]) and finals[m-1][-1]=="A":
                finals[l-1] += "A"
            elif score[l-1][1] == min(score[l-1]) and finals[m-1][-1]=="C":
                finals[l-1] += "C"
            elif score[l-1][3] == min(score[l-1]) and finals[m-1][-1]=="T":
                finals[l-1] += "T"
            elif score[l-1][2] == min(score[l-1]) and finals[m-1][-1]=="G":
                finals[l-1] += "G"
            elif score[l-1][0] == min(score[l-1]):
                finals[l-1] += "A"
            elif score[l-1][1] == min(score[l-1]):
                finals[l-1] += "C"
            elif score[l-1][3] == min(score[l-1]):
                finals[l-1] += "T"
            elif score[l-1][2] == min(score[l-1]):
                finals[l-1] += "G"
            m-=1
    return(scoref, finals)


#for rooted
file2 =  open("hi2.txt").read()
file = file2.split('\n')
n = int(file[0])
mat1= file[1:]
hi = 0
mat={}
leaves = []
leavesv = {}
nodes = []
for edge in range(0,len(mat1),2):
    edge2 = mat1[edge].split('->')
    mat[int(edge2[0])]=[]
    leavesv[int(edge2[0])]= ""
    nodes.append(int(edge2[0]))
    if len(edge2[1])>3:
        leaves.append(hi)
        leavesv[hi]= edge2[1]
        mat[int(edge2[0])].append(hi)
        hi+=1
    else:
        mat[int(edge2[0])].append(int(edge2[1]))
    edge2 = mat1[edge+1].split('->')
    if len(edge2[1])>3:
        leaves.append(hi)
        leavesv[hi]= edge2[1]
        mat[int(edge2[0])].append(hi)
        hi+=1
    else:
        mat[int(edge2[0])].append(int(edge2[1]))
print(mat)
res, result = SmallParsimony(n, mat, leaves, leavesv, nodes)
print(res)
for key in mat.keys():
    for val in mat[key]:
        ham= HammingDistance(result[val], result[key])
        print(result[val]+"->"+result[key]+":"+str(ham))
        print(result[key]+"->"+result[val]+":"+str(ham))

#for unrooted
# file2 =  open("hi2.txt").read()
# file = file2.split('\n')
# n = int(file[0])
# mat1= file[1:]
# hi = 0
# mat={}
# leaves = []
# leavesv = {}
# nodes = []
# for edge in range(0,len(mat1),2):
#     edge2 = mat1[edge].split('->')
#     mat[int(edge2[0])]=[]
#     leavesv[int(edge2[0])]= ""
#     nodes.append(int(edge2[0]))
#     if len(edge2[1])>3:
#         leaves.append(hi)
#         leavesv[hi]= edge2[1]
#         mat[int(edge2[0])].append(hi)
#         hi+=1
#     else:
#         mat[int(edge2[0])].append(int(edge2[1]))
#     edge2 = mat1[edge+1].split('->')
#     if len(edge2[1])>3:
#         leaves.append(hi)
#         leavesv[hi]= edge2[1]
#         mat[int(edge2[0])].append(hi)
#         hi+=1
#     else:
#         mat[int(edge2[0])].append(int(edge2[1]))
# print(mat)
# res, result = SmallParsimony(n, mat, leaves, leavesv, nodes)
# print(res)
# print(result)
# for key in mat.keys():
#     for val in mat[key]:
#         ham= HammingDistance(result[val], result[key])
#         print(result[val]+"->"+result[key]+":"+str(ham))
#         print(result[key]+"->"+result[val]+":"+str(ham))
