import math
def FarthestFirstTraversal(Data, k,m):
    '''
    Input: Integers k and m followed by a set of points Data in m-dimensional space.
    Output: A set Centers consisting of k points (centers) resulting from applying FarthestFirstTraversal(Data, k),
     where the first point from Data is chosen as the first center to initialize the algorithm.
    '''
    Centers = [Data[0]]
    while len(Centers) < k:
        bigdiss = []
        for d in Data:
            diss=[]
            for c in Centers:
                dis = 0
                for mm in range(m):
                    dis+=(d[mm]-c[mm])**2
                dis = math.sqrt(dis)
                diss.append(dis)
            bigdiss.append(min(diss))
        DataPoint = Data[bigdiss.index(max(bigdiss))]
        Centers.append(DataPoint)
    return Centers





# file2 =  open("hi2.txt").read()
# file = file2.split('\n')
# n = file[0].split(' ')
# k= int(n[0])
# m= int(n[1])
# mat = []
# for edge in file[1:]:
#     edge3 = []
#     edge2 = edge.split(' ')
#     for e in edge2:
#         edge3.append(float(e))
#     mat.append(edge3)
# solution =FarthestFirstTraversal(mat,k,m)
# for s in solution:
#     print(*s, sep = " ")
def SquaredErrorDistortionProblem(centers, data, k, m):
    """
    Input: A set of points Data and a set of centers Centers. 
    Output: The squared error distortion Distortion(Data, Centers).
    """
    bigdiss = []
    ssum = 0
    for d in data:
        diss=[]
        for c in centers:
            dis = 0
            for mm in range(m):
                dis+=(d[mm]-c[mm])**2
            dis = math.sqrt(dis)
            diss.append(dis)
        ssum+= min(diss)**2
    return ssum/len(data)


# file2 =  open("hi2.txt").read()
# file = file2.split('\n')
# split = file.index('--------')
# n= file[0].split(' ')
# k= int(n[0])
# m= int(n[1])
# mat2 = []
# for f in file[1:split]:
#     edge3 = []
#     edge2 = f.split(' ')
#     for e in edge2:
#         edge3.append(float(e))
#     mat2.append(edge3)
# mat = []
# for edge in file[split+1:]:
#     edge3 = []
#     edge2 = edge.split(' ')
#     for e in edge2:
#         edge3.append(float(e))
#     mat.append(edge3)
# print(SquaredErrorDistortionProblem(mat2,mat,k,m))


def Lloydalgorithm(data, k, m):
    """
    Input: Integers k and m followed by a set of points Data in m-dimensional space.
Output: A set Centers consisting of k points (centers) resulting from applying the
Lloyd algorithm to Data and Centers, where the first k points from Data are selected as the first k centers.
  You should report your answers to at least three decimal points.
    """
    clusters = []
    for i in range(k):
        clusters.append([])
    centers = data[:k]
    for d in data:
        diss=[]
        for c in centers:
            dis = 0
            for mm in range(m):
                dis+=(d[mm]-c[mm])**2
            dis = math.sqrt(dis)
            diss.append(dis)
        clusters[diss.index(min(diss))].append(d)
    newcenters = []
    for cluster in clusters:
        center = []
        for mm in range(m):
            s = 0
            for c in cluster:
                s+=c[mm]
            center.append(s/len(cluster))
        newcenters.append(center)
    while sorted(newcenters) != sorted(centers):
        centers = newcenters
        clusters = []
        for i in range(k):
            clusters.append([])
        for d in data:
            diss=[]
            for c in centers:
                dis = 0
                for mm in range(m):
                    dis+=(d[mm]-c[mm])**2
                dis = math.sqrt(dis)
                diss.append(dis)
            clusters[diss.index(min(diss))].append(d)
        newcenters = []
        for cluster in clusters:
            center = []
            for mm in range(m):
                s = 0
                for c in cluster:
                    s+=c[mm]
                center.append(s/len(cluster))
            newcenters.append(center)
    return centers






# file2 =  open("hi2.txt").read()
# file = file2.split('\n')
# n = file[0].split(' ')
# k= int(n[0])
# m= int(n[1])
# mat = []
# for edge in file[1:]:
#     edge3 = []
#     edge2 = edge.split(' ')
#     for e in edge2:
#         edge3.append(float(e))
#     mat.append(edge3)
# solution =Lloydalgorithm(mat,k,m)
# for s in solution:
#     print(*s, sep = " ")


def HierarchicalClustering(D, n):
    """
    Input: An integer n, followed by an n x n distance matrix.
Output: The result of applying HierarchicalClustering to this distance matrix (using Davg),
 with each newly created cluster listed on each line.
    """
    Clusters = {}
    remaining = []
    for i in range(n):
        Clusters[i]=[i]
        remaining.append(i)
    T = []
    m = n
    while len(Clusters.keys())>1:
        min = 1000000000000
        for key,value in D.items():
            if value<= min:
                min = value
                minin = key
        T.append(Clusters[minin[0]]+Clusters[minin[1]])
        Clusters[m] = Clusters[minin[0]]+Clusters[minin[1]]
        remaining.remove(minin[0])
        remaining.remove(minin[1])
        for nn in remaining:
            D[nn,m]= (D[minin[0],nn]*len(Clusters[minin[0]])+D[minin[1],nn]*len(Clusters[minin[1]]))/(len(Clusters[minin[0]])+len(Clusters[minin[1]]))
            D[m,nn]= (D[minin[0],nn]*len(Clusters[minin[0]])+D[minin[1],nn]*len(Clusters[minin[1]]))/(len(Clusters[minin[0]])+len(Clusters[minin[1]]))
        blah= []
        for key in D.keys():
            if key[0]== minin[0] or key[0]== minin[1] or key[1] == minin[0] or key[1]== minin[1]:
                blah.append(key)
        for bblah in blah:
            del D[bblah]
        del Clusters[minin[0]]
        del Clusters[minin[1]]
        remaining.append(m)
        m +=1
    return T


file2 =  open("hi2.txt").read()
file = file2.split('\n')
n = int(file[0])
mat = {}
file = file[1:]
for edge in range(len(file)):
    edge2 = file[edge].split(' ')
    for wedge in range(len(edge2)):
        if edge != wedge:
            mat[(edge,wedge)] = float(edge2[wedge])
solution =HierarchicalClustering(mat, n)
print(solution)
for s in solution:
    s2 = []
    for ss in s:
        s2.append(ss+1)
    print(*s2, sep = " ")
