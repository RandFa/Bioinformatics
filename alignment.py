
def GreedyChange(money, coins):
    """
    input: a number representing money and a list of coins that can be used to construct the number
    output: a list "change" that make up money using greedy algorithm
    """
    change = []
    coins.sort(reverse = True)
    while money > 0:
        for coin in coins:
            if coin <= money:
                change.append(coin)
                money = money - coin
                break
    return change


def RecursiveChange(money, Coins):
    if money == 0:
        return 0
    MinNumCoins = 100000000
    for i in range(len(Coins)):
        if money >= Coins[i]:
            NumCoins = RecursiveChange(money - Coins[i], Coins)
            if (NumCoins + 1) < MinNumCoins:
                MinNumCoins = NumCoins + 1
    return MinNumCoins



def DPChange(money, Coins):
    """
    Input: An integer money and an array Coins = (coin1, ..., coind).
    Output: The minimum number of coins with denominations Coins that changes money.
    """
    MinNumCoins = {0:[0, []]}
    for m in range(1, money+1):
        MinNumCoins[m] = [99999999999999999, []]
        for coin in Coins:
                if m >= coin:
                    if MinNumCoins[m- coin][0] + 1 < MinNumCoins[m][0]:
                        MinNumCoins[m][0] = MinNumCoins[m- coin][0] + 1
                        MinNumCoins[m][1]= MinNumCoins[m- coin][1] + [coin]
        if len(MinNumCoins.keys())> Coins[0]:
            newmin = {}
            for a in range(m-Coins[0]+1, m+1):
                newmin[a] = MinNumCoins[a]
            MinNumCoins= newmin
    return MinNumCoins[money]



def ManhattanTourist(n, m, Down, Right):
    """
    Input: Integers n and m, followed by an n × (m + 1) matrix Down and an (n + 1) × m matrix Right. The two matrices are separated by the "-" symbol.
    Output: The length of a longest path from source (0, 0) to sink (n, m) in the rectangular grid whose edges are defined by the matrices Down and Right.
    """
    s = {(0,0): 0}
    maps = {}
    for i in range(1, n+1):
        s[(i, 0)] = s[(i-1, 0)] + int(Down[i-1][0])
    for i in range(1, m+1):
        s[(0, i)] = s[(0, i-1)] + int(Right[0][i-1])
    for i in range(1, n+1):
        for j in range(1, m+1):
            south = s[(i-1,j)] + int(Down[i-1][j])
            east = s[(i, j-1)] + int(Right[i][j-1])
            maxi = max([east, south])
            if maxi == east:
                s[(i,j)] = east
                maps[(i,j)] = (i, j-1)
            if maxi == south:
                s[(i,j)] = south
                maps[(i,j)] = (i-1,j)
    return s[(n,m)], maps



def ModifiedManhattanTourist(n, m, Down, Right, Horizantal):
    """
    Input: Integers n and m, followed by an n × (m + 1) matrix Down and an (n + 1) × m matrix Right. The two matrices are separated by the "-" symbol.
    Output: The length of a longest path from source (0, 0) to sink (n, m) in the rectangular grid whose edges are defined by the matrices Down and Right.
    """
    s = {(0,0): 0}
    maps = {}
    for i in range(1, n+1):
        s[(i, 0)] = s[(i-1, 0)] + int(Down[i-1][0])
    for i in range(1, m+1):
        s[(0, i)] = s[(0, i-1)] + int(Right[0][i-1])
    for i in range(1, n+1):
        for j in range(1, m+1):
            south = s[(i-1,j)] + int(Down[i-1][j])
            east = s[(i, j-1)] + int(Right[i][j-1])
            horz = s[(i-1, j-1)] + int(Horizantal[i-1][j-1])
            maxi = max([east, south, horz])
            if maxi == east:
                s[(i,j)] = east
                maps[(i,j)] = (i, j-1)
            if maxi == south:
                s[(i,j)] = south
                maps[(i,j)] = (i-1,j)
            if maxi == east:
                s[(i,j)] = south
                maps[(i,j)] = (i-1, j-1)
    return s[(n,m)]




import sys
sys.setrecursionlimit(3000)
def LCSBackTrack(v, w):
    s = {}
    Backtrack = {}
    for i in range(len(v)+1):
        s[(i, 0)] = 0
    for j in range(len(w)+1):
        s[(0, j)] = 0
    for i in range(1, len(v)+1):
        for j in range(1, len(w)+1):
            match = 0
            if v[i-1] == w[j-1]:
                match = 1
            s[(i, j)] = max([s[(i-1, j)] , s[(i,j-1)] , (s[(i-1, j-1)] + match)])
            if s[(i,j)] == s[(i-1,j)]:
                Backtrack[(i, j)] = "↓"
            elif s[(i, j)] == s[(i, j-1)]:
                Backtrack[(i, j)] = "→"
            elif s[(i, j)] == (s[(i-1, j-1)] + match):
                    Backtrack[(i, j)] = "↘"
    return Backtrack
def OutputLCS(backtrack, v, i, j):
    """
    Input: Two strings s and t.
    Output: A longest common subsequence of s and t. (Note: more than one solution may exist, in which case you may output any one.)
    """
    while i > 0 and j > 0:
        if backtrack[(i,j)] == "↘":
            LCS= LCS + v[i-1]
            i = i-1
            j = j-1
        elif backtrack[(i,j)] == "↓":
            i = i-1
        else:
            j = j-1
    return LCS[::-1]
# v = "PLEASANTLY"
# w = "MEANLY"
# backtrack=LCSBackTrack(v,w)
# print(OutputLCS(backtrack, v, len(v), len(w)))


def LongestPathDAG(st, en , edges):
    """
    Input: An integer representing the starting node to consider in a graph, followed by an integer representing the ending node to consider, followed by a list of edges in the graph. The edge notation "0->1:7" indicates that an edge connects node 0 to node 1 with weight 7.  You may assume a given topological order corresponding to nodes in increasing order.
    Output: The length of a longest path in the graph, followed by a longest path. (If multiple longest paths exist, you may return any one.)
    """
    longpath = [en]
    score = 0
    while st != en:
        edges2 = {}
        for edge in edges:
            if edge[1] == en:
                edges2[edge] = edges[edge]
        maxi = max(list(edges2.values()))
        score = score + maxi
        for key in edges2:
            if edges2[key] == maxi:
                longpath = [key[0]]+ longpath
                break
        en = longpath[0]
    return score, longpath



scorematrix = {'A': {'A': 4, 'C': 0, 'E': -1, 'D': -2, 'G': 0, 'F': -2, 'I': -1, 'H': -2, 'K': -1, 'M': -1, 'L': -1, 'N': -2, 'Q': -1, 'P': -1, 'S': 1, 'R': -1, 'T': 0, 'W': -3, 'V': 0, 'Y': -2}, 'C': {'A': 0, 'C': 9, 'E': -4, 'D': -3, 'G': -3, 'F': -2, 'I': -1, 'H': -3, 'K': -3, 'M': -1, 'L': -1, 'N': -3, 'Q': -3, 'P': -3, 'S': -1, 'R': -3, 'T': -1, 'W': -2, 'V': -1, 'Y': -2}, 'E': {'A': -1, 'C': -4, 'E': 5, 'D': 2, 'G': -2, 'F': -3, 'I': -3, 'H': 0, 'K': 1, 'M': -2, 'L': -3, 'N': 0, 'Q': 2, 'P': -1, 'S': 0, 'R': 0, 'T': -1, 'W': -3, 'V': -2, 'Y': -2}, 'D': {'A': -2, 'C': -3, 'E': 2, 'D': 6, 'G': -1, 'F': -3, 'I': -3, 'H': -1, 'K': -1, 'M': -3, 'L': -4, 'N': 1, 'Q': 0, 'P': -1, 'S': 0, 'R': -2, 'T': -1, 'W': -4, 'V': -3, 'Y': -3}, 'G': {'A': 0, 'C': -3, 'E': -2, 'D': -1, 'G': 6, 'F': -3, 'I': -4, 'H': -2, 'K': -2, 'M': -3, 'L': -4, 'N': 0, 'Q': -2, 'P': -2, 'S': 0, 'R': -2, 'T': -2, 'W': -2, 'V': -3, 'Y': -3}, 'F': {'A': -2, 'C': -2, 'E': -3, 'D': -3, 'G': -3, 'F': 6, 'I': 0, 'H': -1, 'K': -3, 'M': 0, 'L': 0, 'N': -3, 'Q': -3, 'P': -4, 'S': -2, 'R': -3, 'T': -2, 'W': 1, 'V': -1, 'Y': 3}, 'I': {'A': -1, 'C': -1, 'E': -3, 'D': -3, 'G': -4, 'F': 0, 'I': 4, 'H': -3, 'K': -3, 'M': 1, 'L': 2, 'N': -3, 'Q': -3, 'P': -3, 'S': -2, 'R': -3, 'T': -1, 'W': -3, 'V': 3, 'Y': -1}, 'H': {'A': -2, 'C': -3, 'E': 0, 'D': -1, 'G': -2, 'F': -1, 'I': -3, 'H': 8, 'K': -1, 'M': -2, 'L': -3, 'N': 1, 'Q': 0, 'P': -2, 'S': -1, 'R': 0, 'T': -2, 'W': -2, 'V': -3, 'Y': 2}, 'K': {'A': -1, 'C': -3, 'E': 1, 'D': -1, 'G': -2, 'F': -3, 'I': -3, 'H': -1, 'K': 5, 'M': -1, 'L': -2, 'N': 0, 'Q': 1, 'P': -1, 'S': 0, 'R': 2, 'T': -1, 'W': -3, 'V': -2, 'Y': -2}, 'M': {'A': -1, 'C': -1, 'E': -2, 'D': -3, 'G': -3, 'F': 0, 'I': 1, 'H': -2, 'K': -1, 'M': 5, 'L': 2, 'N': -2, 'Q': 0, 'P': -2, 'S': -1, 'R': -1, 'T': -1, 'W': -1, 'V': 1, 'Y': -1}, 'L': {'A': -1, 'C': -1, 'E': -3, 'D': -4, 'G': -4, 'F': 0, 'I': 2, 'H': -3, 'K': -2, 'M': 2, 'L': 4, 'N': -3, 'Q': -2, 'P': -3, 'S': -2, 'R': -2, 'T': -1, 'W': -2, 'V': 1, 'Y': -1}, 'N': {'A': -2, 'C': -3, 'E': 0, 'D': 1, 'G': 0, 'F': -3, 'I': -3, 'H': 1, 'K': 0, 'M': -2, 'L': -3, 'N': 6, 'Q': 0, 'P': -2, 'S': 1, 'R': 0, 'T': 0, 'W': -4, 'V': -3, 'Y': -2}, 'Q': {'A': -1, 'C': -3, 'E': 2, 'D': 0, 'G': -2, 'F': -3, 'I': -3, 'H': 0, 'K': 1, 'M': 0, 'L': -2, 'N': 0, 'Q': 5, 'P': -1, 'S': 0, 'R': 1, 'T': -1, 'W': -2, 'V': -2, 'Y': -1}, 'P': {'A': -1, 'C': -3, 'E': -1, 'D': -1, 'G': -2, 'F': -4, 'I': -3, 'H': -2, 'K': -1, 'M': -2, 'L': -3, 'N': -2, 'Q': -1, 'P': 7, 'S': -1, 'R': -2, 'T': -1, 'W': -4, 'V': -2, 'Y': -3}, 'S': {'A': 1, 'C': -1, 'E': 0, 'D': 0, 'G': 0, 'F': -2, 'I': -2, 'H': -1, 'K': 0, 'M': -1, 'L': -2, 'N': 1, 'Q': 0, 'P': -1, 'S': 4, 'R': -1, 'T': 1, 'W': -3, 'V': -2, 'Y': -2}, 'R': {'A': -1, 'C': -3, 'E': 0, 'D': -2, 'G': -2, 'F': -3, 'I': -3, 'H': 0, 'K': 2, 'M': -1, 'L': -2, 'N': 0, 'Q': 1, 'P': -2, 'S': -1, 'R': 5, 'T': -1, 'W': -3, 'V': -3, 'Y': -2}, 'T': {'A': 0, 'C': -1, 'E': -1, 'D': -1, 'G': -2, 'F': -2, 'I': -1, 'H': -2, 'K': -1, 'M': -1, 'L': -1, 'N': 0, 'Q': -1, 'P': -1, 'S': 1, 'R': -1, 'T': 5, 'W': -2, 'V': 0, 'Y': -2}, 'W': {'A': -3, 'C': -2, 'E': -3, 'D': -4, 'G': -2, 'F': 1, 'I': -3, 'H': -2, 'K': -3, 'M': -1, 'L': -2, 'N': -4, 'Q': -2, 'P': -4, 'S': -3, 'R': -3, 'T': -2, 'W': 11, 'V': -3, 'Y': 2}, 'V': {'A': 0, 'C': -1, 'E': -2, 'D': -3, 'G': -3, 'F': -1, 'I': 3, 'H': -3, 'K': -2, 'M': 1, 'L': 1, 'N': -3, 'Q': -2, 'P': -2, 'S': -2, 'R': -3, 'T': 0, 'W': -3, 'V': 4, 'Y': -1}, 'Y': {'A': -2, 'C': -2, 'E': -2, 'D': -3, 'G': -3, 'F': 3, 'I': -1, 'H': 2, 'K': -2, 'M': -1, 'L': -1, 'N': -2, 'Q': -1, 'P': -3, 'S': -2, 'R': -2, 'T': -2, 'W': 2, 'V': -1, 'Y': 7}}
def LCSBackTrackGlobal(v, w):
    s = {(0,0): 0}
    Backtrack = {}
    for i in range(1, len(v)+1):
        s[(i, 0)] = s[(i-1, 0)] -5
    for j in range(1, len(w)+1):
        s[(0, j)] = s[(0, j-1)] -5
    for i in range(1, len(v)+1):
        for j in range(1, len(w)+1):
            a =s[(i-1, j)] -5
            b = s[(i,j-1)] - 5
            c = s[(i-1, j-1)] + scorematrix[v[i-1]][w[j-1]]
            s[(i, j)] = max([ a, b ,c ])
            if s[(i,j)] == a:
                Backtrack[(i, j)] = "↓"
            elif s[(i, j)] == b:
                Backtrack[(i, j)] = "→"
            elif s[(i, j)] == c:
                Backtrack[(i, j)] = "↘"
    i = len(v)
    j = len(w)
    score2 = s[(i, j)]
    return Backtrack, score2

def GlobalAlignmentProblem(v, w):
    """
    Input: Two protein strings written in the single-letter amino acid alphabet.
    Output: The maximum alignment score of these strings followed by an alignment achieving this maximum score. Use the BLOSUM62 scoring matrix for matches and mismatches as well as the indel penalty σ = 5.
    """
    backtrack, score = LCSBackTrackGlobal(v, w)
    i = len(v)
    j= len(w)
    newV = ""
    newW= ""
    while i > 0 and j > 0:
        if backtrack[(i,j)] == "↘":
            newV= newV + v[i-1]
            newW= newW + w[j-1]
            i = i-1
            j = j-1
        elif backtrack[(i,j)] == "↓":
            newV= newV + v[i-1]
            newW= newW + "-"
            i = i-1
        else:
            newV= newV + "-"
            newW= newW + w[j-1]
            j = j-1
    if newV[::-1][0] != v[0]:
        newV= newV + v[0]
        newW= newW + "-"
    elif newW[::-1][0] != w[0]:
        newW= newW + w[0]
        newV= newV + "-"
    return newV[::-1], newW[::-1], score



pam250 = {'A': {'A': 2, 'C': -2, 'E': 0, 'D': 0, 'G': 1, 'F': -3, 'I': -1, 'H': -1, 'K': -1, 'M': -1, 'L': -2, 'N': 0, 'Q': 0, 'P': 1, 'S': 1, 'R': -2, 'T': 1, 'W': -6, 'V': 0, 'Y': -3}, 'C': {'A': -2, 'C': 12, 'E': -5, 'D': -5, 'G': -3, 'F': -4, 'I': -2, 'H': -3, 'K': -5, 'M': -5, 'L': -6, 'N': -4, 'Q': -5, 'P': -3, 'S': 0, 'R': -4, 'T': -2, 'W': -8, 'V': -2, 'Y': 0}, 'E': {'A': 0, 'C': -5, 'E': 4, 'D': 3, 'G': 0, 'F': -5, 'I': -2, 'H': 1, 'K': 0, 'M': -2, 'L': -3, 'N': 1, 'Q': 2, 'P': -1, 'S': 0, 'R': -1, 'T': 0, 'W': -7, 'V': -2, 'Y': -4}, 'D': {'A': 0, 'C': -5, 'E': 3, 'D': 4, 'G': 1, 'F': -6, 'I': -2, 'H': 1, 'K': 0, 'M': -3, 'L': -4, 'N': 2, 'Q': 2, 'P': -1, 'S': 0, 'R': -1, 'T': 0, 'W': -7, 'V': -2, 'Y': -4}, 'G': {'A': 1, 'C': -3, 'E': 0, 'D': 1, 'G': 5, 'F': -5, 'I': -3, 'H': -2, 'K': -2, 'M': -3, 'L': -4, 'N': 0, 'Q': -1, 'P': 0, 'S': 1, 'R': -3, 'T': 0, 'W': -7, 'V': -1, 'Y': -5}, 'F': {'A': -3, 'C': -4, 'E': -5, 'D': -6, 'G': -5, 'F': 9, 'I': 1, 'H': -2, 'K': -5, 'M': 0, 'L': 2, 'N': -3, 'Q': -5, 'P': -5, 'S': -3, 'R': -4, 'T': -3, 'W': 0, 'V': -1, 'Y': 7}, 'I': {'A': -1, 'C': -2, 'E': -2, 'D': -2, 'G': -3, 'F': 1, 'I': 5, 'H': -2, 'K': -2, 'M': 2, 'L': 2, 'N': -2, 'Q': -2, 'P': -2, 'S': -1, 'R': -2, 'T': 0, 'W': -5, 'V': 4, 'Y': -1}, 'H': {'A': -1, 'C': -3, 'E': 1, 'D': 1, 'G': -2, 'F': -2, 'I': -2, 'H': 6, 'K': 0, 'M': -2, 'L': -2, 'N': 2, 'Q': 3, 'P': 0, 'S': -1, 'R': 2, 'T': -1, 'W': -3, 'V': -2, 'Y': 0}, 'K': {'A': -1, 'C': -5, 'E': 0, 'D': 0, 'G': -2, 'F': -5, 'I': -2, 'H': 0, 'K': 5, 'M': 0, 'L': -3, 'N': 1, 'Q': 1, 'P': -1, 'S': 0, 'R': 3, 'T': 0, 'W': -3, 'V': -2, 'Y': -4}, 'M': {'A': -1, 'C': -5, 'E': -2, 'D': -3, 'G': -3, 'F': 0, 'I': 2, 'H': -2, 'K': 0, 'M': 6, 'L': 4, 'N': -2, 'Q': -1, 'P': -2, 'S': -2, 'R': 0, 'T': -1, 'W': -4, 'V': 2, 'Y': -2}, 'L': {'A': -2, 'C': -6, 'E': -3, 'D': -4, 'G': -4, 'F': 2, 'I': 2, 'H': -2, 'K': -3, 'M': 4, 'L': 6, 'N': -3, 'Q': -2, 'P': -3, 'S': -3, 'R': -3, 'T': -2, 'W': -2, 'V': 2, 'Y': -1}, 'N': {'A': 0, 'C': -4, 'E': 1, 'D': 2, 'G': 0, 'F': -3, 'I': -2, 'H': 2, 'K': 1, 'M': -2, 'L': -3, 'N': 2, 'Q': 1, 'P': 0, 'S': 1, 'R': 0, 'T': 0, 'W': -4, 'V': -2, 'Y': -2}, 'Q': {'A': 0, 'C': -5, 'E': 2, 'D': 2, 'G': -1, 'F': -5, 'I': -2, 'H': 3, 'K': 1, 'M': -1, 'L': -2, 'N': 1, 'Q': 4, 'P': 0, 'S': -1, 'R': 1, 'T': -1, 'W': -5, 'V': -2, 'Y': -4}, 'P': {'A': 1, 'C': -3, 'E': -1, 'D': -1, 'G': 0, 'F': -5, 'I': -2, 'H': 0, 'K': -1, 'M': -2, 'L': -3, 'N': 0, 'Q': 0, 'P': 6, 'S': 1, 'R': 0, 'T': 0, 'W': -6, 'V': -1, 'Y': -5}, 'S': {'A': 1, 'C': 0, 'E': 0, 'D': 0, 'G': 1, 'F': -3, 'I': -1, 'H': -1, 'K': 0, 'M': -2, 'L': -3, 'N': 1, 'Q': -1, 'P': 1, 'S': 2, 'R': 0, 'T': 1, 'W': -2, 'V': -1, 'Y': -3}, 'R': {'A': -2, 'C': -4, 'E': -1, 'D': -1, 'G': -3, 'F': -4, 'I': -2, 'H': 2, 'K': 3, 'M': 0, 'L': -3, 'N': 0, 'Q': 1, 'P': 0, 'S': 0, 'R': 6, 'T': -1, 'W': 2, 'V': -2, 'Y': -4}, 'T': {'A': 1, 'C': -2, 'E': 0, 'D': 0, 'G': 0, 'F': -3, 'I': 0, 'H': -1, 'K': 0, 'M': -1, 'L': -2, 'N': 0, 'Q': -1, 'P': 0, 'S': 1, 'R': -1, 'T': 3, 'W': -5, 'V': 0, 'Y': -3}, 'W': {'A': -6, 'C': -8, 'E': -7, 'D': -7, 'G': -7, 'F': 0, 'I': -5, 'H': -3, 'K': -3, 'M': -4, 'L': -2, 'N': -4, 'Q': -5, 'P': -6, 'S': -2, 'R': 2, 'T': -5, 'W': 17, 'V': -6, 'Y': 0}, 'V': {'A': 0, 'C': -2, 'E': -2, 'D': -2, 'G': -1, 'F': -1, 'I': 4, 'H': -2, 'K': -2, 'M': 2, 'L': 2, 'N': -2, 'Q': -2, 'P': -1, 'S': -1, 'R': -2, 'T': 0, 'W': -6, 'V': 4, 'Y': -2}, 'Y': {'A': -3, 'C': 0, 'E': -4, 'D': -4, 'G': -5, 'F': 7, 'I': -1, 'H': 0, 'K': -4, 'M': -2, 'L': -1, 'N': -2, 'Q': -4, 'P': -5, 'S': -3, 'R': -4, 'T': -3, 'W': 0, 'V': -2, 'Y': 10}}
def LCSBackTrackLocal(v, w):
    s = {(0,0): 0}
    Backtrack = {}
    score = 0
    for i in range(1, len(v)+1):
        s[(i, 0)] = 0
    for j in range(1, len(w)+1):
        s[(0, j)] = 0
    for i in range(1, len(v)+1):
        for j in range(1, len(w)+1):
            a =s[(i-1, j)] -5
            b = s[(i,j-1)] -5
            c = s[(i-1, j-1)] + pam250[v[i-1]][w[j-1]]
            s[(i, j)] = max([ a, b ,c, 0])
            if s[(i,j)] == a:
                Backtrack[(i, j)] = "↓"
            elif s[(i, j)] == b:
                Backtrack[(i, j)] = "→"
            elif s[(i, j)] == c:
                Backtrack[(i, j)] = "↘"
            elif s[(i, j)] == 0:
                Backtrack[(i, j)] = 0
            # if s[(i, j)] > score:
            #     score = s[(i, j)]
    return Backtrack, s

def LocalAlignmentProblem(v, w):
    """
    Input: Two protein strings written in the single-letter amino acid alphabet.
    Output: The maximum score of a local alignment of the strings, followed by a local alignment of these strings achieving the maximum score. Use the PAM250 scoring matrix for matches and mismatches as well as the indel penalty σ = 5.
    """
    backtrack, s = LCSBackTrackLocal(v, w)
    newV = ""
    newW= ""
    score = max(list(s.values()))
    for track in s:
        if s[track] == score:
            i, j = track[0], track[1]
    weight = s[(i,j)]
    while weight !=0:
        if backtrack[(i,j)] == "↘":
            newV= newV + v[i-1]
            newW= newW + w[j-1]
            i = i-1
            j = j-1
        elif backtrack[((i, j))] == "↓":
            newV= newV + v[i-1]
            newW= newW + "-"
            i = i-1
        elif backtrack[(i, j)] == "→":
            newV= newV + "-"
            newW= newW + w[j-1]
            j = j-1
        if i == 1 and j == 1:
            break
        weight = backtrack[(i,j)]
    return newV[::-1], newW[::-1], score
# print(LocalAlignmentProblem("MEANLY", "PENALTY"))



def EditDistanceProblemBacktrack(v, w):
    """
    Input: Two strings.
    Output: The edit distance between these strings.
    """
    s = {(0,0): 0}
    Backtrack = {}
    distace = 0
    for i in range(1, len(v)+1):
        s[(i, 0)] = s[(i-1, 0)] -1
    for j in range(1, len(w)+1):
        s[(0, j)] = s[(0, j-1)]-1
    for i in range(1, len(v)+1):
        for j in range(1, len(w)+1):
            if v[i-1] == w[j-1]:
                c = s[(i-1, j-1)]
            else:
                c = s[(i-1, j-1)] -1
            a =s[(i-1, j)] -1
            b = s[(i,j-1)] -1
            s[(i, j)] = max([ a, b ,c ])
            if s[(i, j)] == c:
                Backtrack[(i, j)] = "↘"
            elif s[(i,j)] == a:
                Backtrack[(i, j)] = "↓"
            elif s[(i, j)] == b:
                Backtrack[(i, j)] = "→"
    return Backtrack



def EditDistanceProblem(v, w):
    """
    Input: Two strings.
    Output: The edit distance between these strings.
    """
    backtrack= EditDistanceProblemBacktrack(v, w)
    i = len(v)
    j= len(w)
    newV = ""
    newW= ""
    distance = 0
    while i > 0 and j > 0:
        if backtrack[(i,j)] == "↘":
            newV= newV + v[i-1]
            newW= newW + w[j-1]
            if v[i-1] != w[j-1]:
                distance += 1
            i = i-1
            j = j-1
        elif backtrack[(i,j)] == "↓":
            newV= newV + v[i-1]
            newW= newW + "-"
            i = i-1
            distance+=1
        else:
            newV= newV + "-"
            newW= newW + w[j-1]
            j = j-1
            distance+=1
    if newV[::-1][0] != v[0]:
        newV= newV + v[0]
        newW= newW + "-"
        distance += 1
    elif newW[::-1][0] != w[0]:
        newW= newW + w[0]
        newV= newV + "-"
        distance +=1
    return distance, newV[::-1], newW[::-1]
# print(EditDistanceProblem("AG","T"))

def LCSBackTrackFitting(v, w):
    s = {(0,0): 0}
    Backtrack = {}
    score = 0
    for i in range(1, len(v)+1):
        s[(i, 0)] = 0
    for j in range(1, len(w)+1):
        s[(0, j)] = s[(0,j-1)] -1
    for i in range(1, len(v)+1):
        for j in range(1, len(w)+1):
            a =s[(i-1, j)] -1
            b = s[(i,j-1)] -1
            if v[i-1] == w[j-1]:
                c = s[(i-1, j-1)] + 1
            else:
                c = s[(i-1, j-1)] - 1
            s[(i, j)] = max([ a, b ,c, 0])
            if s[(i,j)] == a:
                Backtrack[(i, j)] = "↓"
            elif s[(i, j)] == b:
                Backtrack[(i, j)] = "→"
            elif s[(i, j)] == c:
                Backtrack[(i, j)] = "↘"
            elif s[(i, j)] == 0:
                Backtrack[(i, j)] = 0
            # if s[(i, j)] > score:
            #     score = s[(i, j)]
    return Backtrack, s

def FittingAlignmentProblem(v, w):
    """
    Input: Two nucleotide strings v and w, where v has length at most 1000 and w has length at most 100.
    Output: A highest-scoring fitting alignment between v and w. Use the simple scoring method in which matches count +1 and both the mismatch and indel penalties are 1.
    """
    backtrack, s = LCSBackTrackFitting(v, w)
    newV = ""
    newW= ""
    score = 0
    maxsc = []
    for column in s:
        a, b = column[0], column[1]
        if a>= len(w) and b== len(w):
            if s[column] > score:
                score = s[column]
                i,j = a,b
                maxsc = [(i , j)]
            elif s[column] == score:
                i,j = a,b
                maxsc.append((i , j))
    best = {score:[], score-1: []}
    print(maxsc)
    for maxc in maxsc:
        meas = len(w)
        i, j = maxc[0], maxc[1]
        weight = s[(i,j)]
        newV = ""
        newW = ""
        while meas != 0:
            if backtrack[(i,j)] == "↘" or backtrack[(i,j)] ==0:
                newV= newV + v[i-1]
                newW= newW + w[j-1]
                i = i-1
                j = j-1
                meas = meas-1
            elif backtrack[((i, j))] == "↓":
                newV= newV + v[i-1]
                newW= newW + "-"
                i = i-1
            elif backtrack[(i, j)] == "→":
                newV= newV + "-"
                newW= newW + w[j-1]
                j = j-1
                meas = meas-1
            if i <= 1 and j <= 1:
                newV= newV + v[0]
                newW= newW + w[0]
                break
        # newV= newV + v[i-1]
        # newW= newW + w[j-1]
        if newV[-1] != newW[-1]:
            best[score-1].append([newV[::-1], newW[::-1]])
        else:
            best[score].append([newV[::-1], newW[::-1]])

    print(best)
    if best[score] !=[]:
        bestv, bestw = best[score][-1][0], best[score][-1][1]
    else:
        bestv, bestw = best[score-1][0][0], best[score-1][0][1]
        score = score-1
    return bestv, bestw, score


def LCSBackTrackAlignmentOverlap(v, w):
    s = {(0,0): 0}
    Backtrack = {}
    score = 0
    for i in range(1, len(v)+1):
        s[(i, 0)] = 0
    for j in range(1, len(w)+1):
        s[(0, j)] = s[(0,j-1)] -2
    for i in range(1, len(v)+1):
        for j in range(1, len(w)+1):
            a =s[(i-1, j)] -2
            b = s[(i,j-1)] -2
            if v[i-1] == w[j-1]:
                c = s[(i-1, j-1)] + 1
            else:
                c = s[(i-1, j-1)] - 2
            s[(i, j)] = max([ a, b ,c, 0])
            if s[(i,j)] == a:
                Backtrack[(i, j)] = "↓"
            elif s[(i, j)] == b:
                Backtrack[(i, j)] = "→"
            elif s[(i, j)] == c:
                Backtrack[(i, j)] = "↘"
            elif s[(i, j)] == 0:
                Backtrack[(i, j)] = 0
            # if s[(i, j)] > score:
            #     score = s[(i, j)]
    return Backtrack, s

def OverlapAlignmentProblem(v, w):
    """
    Input: Two nucleotide strings v and w, where v has length at most 1000 and w has length at most 100.
    Output: A highest-scoring fitting alignment between v and w. Use the simple scoring method in which matches count +1 and both the mismatch and indel penalties are 1.
    """
    backtrack, s = LCSBackTrackAlignmentOverlap(v, w)
    newV = ""
    newW= ""
    maxsc = []
    score =0
    for node in s:
        if node[0] == len(v) and node[1]!= 0:
            if s[node]> score:
                maxsc = [node]
                score = s[node]
            elif s[node] == score:
                maxsc.append(node)
    print(maxsc)
    best = {score:[], score-1: []}
    for maxc in maxsc:
        i, j = maxc[0], maxc[1]
        meas = min([i, len(w)])

        weight = s[(i,j)]
        newV = ""
        newW = ""
        while meas != 0:
            if backtrack[(i,j)] == "↘" or backtrack[(i,j)] ==0:
                newV= newV + v[i-1]
                newW= newW + w[j-1]
                i = i-1
                j = j-1
                meas = meas-1
            elif backtrack[((i, j))] == "↓":
                newV= newV + v[i-1]
                newW= newW + "-"
                i = i-1
            elif backtrack[(i, j)] == "→":
                newV= newV + "-"
                newW= newW + w[j-1]
                j = j-1
                meas = meas-1
            if i <= 1 or j <= 1:
                newV= newV + v[i-1]
                newW= newW + w[0]
                break
        # newV= newV + v[i-1]
        # newW= newW + w[j-1]
        if newV[-1] != newW[-1]:
            best[score-1].append([newV[::-1], newW[::-1]])
        else:
            best[score].append([newV[::-1], newW[::-1]])

    if best[score] !=[]:
        bestv, bestw = best[score][0][0], best[score][0][1]
    else:
        bestv, bestw = best[score-1][0][0], best[score-1][0][1]
        score = score-2
    # if newV[::-1][0] != v[0]:
    #     newV= newV + v[0]
    #     newW= newW + "-"
    # elif newW[::-1][0] != w[0]:
    #     newW= newW + w[0]
    #     newV= newV + "-"
    return bestv, bestw, score

BLOSUM62={'A': {'A': 4, 'C': 0, 'E': -1, 'D': -2, 'G': 0, 'F': -2, 'I': -1, 'H': -2, 'K': -1, 'M': -1, 'L': -1, 'N': -2, 'Q': -1, 'P': -1, 'S': 1, 'R': -1, 'T': 0,  'W': -3, 'V': 0, 'Y': -2},

          'C': {'A': 0, 'C': 9, 'E': -4, 'D': -3, 'G': -3, 'F': -2, 'I': -1, 'H': -3, 'K': -3, 'M': -1, 'L': -1, 'N': -3, 'Q': -3, 'P': -3, 'S': -1, 'R': -3, 'T': -1, 'W': -2, 'V': -1, 'Y': -2},

          'E': {'A': -1, 'C': -4, 'E': 5, 'D': 2, 'G': -2, 'F': -3, 'I': -3, 'H': 0, 'K': 1, 'M': -2, 'L': -3, 'N': 0, 'Q': 2, 'P': -1, 'S': 0, 'R': 0, 'T': -1, 'W': -3, 'V': -2, 'Y': -2},

          'D': {'A': -2, 'C': -3, 'E': 2, 'D': 6, 'G': -1, 'F': -3, 'I': -3, 'H': -1, 'K': -1, 'M': -3, 'L': -4, 'N': 1, 'Q': 0, 'P': -1, 'S': 0, 'R': -2, 'T': -1,  'W': -4, 'V': -3, 'Y': -3},

          'G': {'A': 0, 'C': -3, 'E': -2, 'D': -1, 'G': 6, 'F': -3, 'I': -4, 'H': -2, 'K': -2, 'M': -3, 'L': -4, 'N': 0, 'Q': -2, 'P': -2, 'S': 0, 'R': -2, 'T': -2, 'W': -2, 'V': -3, 'Y': -3},

          'F': {'A': -2, 'C': -2, 'E': -3, 'D': -3, 'G': -3, 'F': 6, 'I': 0, 'H': -1, 'K': -3, 'M': 0, 'L': 0, 'N': -3, 'Q': -3, 'P': -4, 'S': -2, 'R': -3, 'T': -2, 'W': 1, 'V': -1, 'Y': 3},

          'I': {'A': -1, 'C': -1, 'E': -3, 'D': -3, 'G': -4, 'F': 0, 'I': 4, 'H': -3, 'K': -3, 'M': 1, 'L': 2, 'N': -3, 'Q': -3, 'P': -3, 'S': -2, 'R': -3, 'T': -1, 'W': -3, 'V': 3, 'Y': -1},

          'H': {'A': -2, 'C': -3, 'E': 0, 'D': -1, 'G': -2, 'F': -1, 'I': -3, 'H': 8, 'K': -1, 'M': -2, 'L': -3, 'N': 1, 'Q': 0, 'P': -2, 'S': -1, 'R': 0, 'T': -2, 'W': -2, 'V': -3, 'Y': 2},

          'K': {'A': -1, 'C': -3, 'E': 1, 'D': -1, 'G': -2, 'F': -3, 'I': -3, 'H': -1, 'K': 5, 'M': -1, 'L': -2, 'N': 0, 'Q': 1, 'P': -1, 'S': 0, 'R': 2, 'T': -1,  'W': -3, 'V': -2, 'Y': -2},

          'M': {'A': -1, 'C': -1, 'E': -2, 'D': -3, 'G': -3, 'F': 0, 'I': 1, 'H': -2, 'K': -1, 'M': 5, 'L': 2, 'N': -2, 'Q': 0, 'P': -2, 'S': -1, 'R': -1, 'T': -1,  'W': -1, 'V': 1, 'Y': -1},

          'L': {'A': -1, 'C': -1, 'E': -3, 'D': -4, 'G': -4, 'F': 0, 'I': 2, 'H': -3, 'K': -2, 'M': 2, 'L': 4, 'N': -3, 'Q': -2, 'P': -3, 'S': -2, 'R': -2, 'T': -1,  'W': -2, 'V': 1, 'Y': -1},

          'N': {'A': -2, 'C': -3, 'E': 0, 'D': 1, 'G': 0, 'F': -3, 'I': -3, 'H': 1, 'K': 0, 'M': -2, 'L': -3, 'N': 6, 'Q': 0, 'P': -2, 'S': 1, 'R': 0, 'T': 0, 'W': -4, 'V': -3, 'Y': -2},

          'Q': {'A': -1, 'C': -3, 'E': 2, 'D': 0, 'G': -2, 'F': -3, 'I': -3, 'H': 0, 'K': 1, 'M': 0, 'L': -2, 'N': 0, 'Q': 5, 'P': -1, 'S': 0, 'R': 1, 'T': -1, 'W': -2,'V': -2, 'Y': -1},

          'P': {'A': -1, 'C': -3, 'E': -1, 'D': -1, 'G': -2, 'F': -4, 'I': -3, 'H': -2, 'K': -1, 'M': -2, 'L': -3, 'N': -2, 'Q': -1, 'P': 7, 'S': -1, 'R': -2, 'T':  -1, 'W': -4, 'V': -2, 'Y': -3},

          'S': {'A': 1, 'C': -1, 'E': 0, 'D': 0, 'G': 0, 'F': -2, 'I': -2, 'H': -1, 'K': 0, 'M': -1, 'L': -2, 'N': 1, 'Q': 0, 'P': -1, 'S': 4, 'R': -1, 'T': 1, 'W': -3, 'V': -2, 'Y': -2},

          'R': {'A': -1, 'C': -3, 'E': 0, 'D': -2, 'G': -2, 'F': -3, 'I': -3, 'H': 0, 'K': 2, 'M': -1, 'L': -2, 'N': 0, 'Q': 1, 'P': -2, 'S': -1, 'R': 5, 'T': -1,  'W': -3, 'V': -3, 'Y': -2},

          'T': {'A': 0, 'C': -1, 'E': -1, 'D': -1, 'G': -2, 'F': -2, 'I': -1, 'H': -2, 'K': -1, 'M': -1, 'L': -1, 'N': 0, 'Q': -1, 'P': -1, 'S': 1, 'R': -1, 'T': 5,  'W': -2, 'V': 0, 'Y': -2},

          'W': {'A': -3, 'C': -2, 'E': -3, 'D': -4, 'G': -2, 'F': 1, 'I': -3, 'H': -2, 'K': -3, 'M': -1, 'L': -2, 'N': -4, 'Q': -2, 'P': -4, 'S': -3, 'R': -3, 'T': -2, 'W': 11, 'V': -3, 'Y': 2},

          'V': {'A': 0, 'C': -1, 'E': -2, 'D': -3, 'G': -3, 'F': -1, 'I': 3, 'H': -3, 'K': -2, 'M': 1, 'L': 1, 'N': -3, 'Q': -2, 'P': -2, 'S': -2, 'R': -3, 'T': 0, 'W': -3, 'V': 4, 'Y': -1},

          'Y': {'A': -2, 'C': -2, 'E': -2, 'D': -3, 'G': -3, 'F': 3, 'I': -1, 'H': 2, 'K': -2, 'M': -1, 'L': -1, 'N': -2, 'Q': -1, 'P': -3, 'S': -2, 'R': -2, 'T':  -2, 'W': 2, 'V': -1, 'Y': 7}}
def AffineGapPenaltiesBacktrack(v, w):

    supper ={(0,0): 0}
    slower={(0,0): 0}
    smiddle ={(0,0): 0}
    for i in range(1, len(v)+1):
        supper[i, 0]= -11 + (i-1)*-1
        slower[i, 0] = -1000
        smiddle[i, 0] = -11 + (i-1)*-1
    for j in range(1, len(w)+1):
        slower[0, j]= -11 + (j-1)*-1
        supper[0, j] = -1000
        smiddle[0, j] = -11 + (j-1)*-1
    backtrack = {}
    backtracklow = {}
    backtrackup = {}
    for j in range(1, len(w)+1):
        for i in range(1, len(v)+1):
            slower[i, j] = max([slower[(i-1,j)]-1, smiddle[(i-1,j)]-11])
            if slower[(i,j)] == slower[(i-1,j)]-1:
                backtracklow[i,j] = "↓"
            elif slower[(i,j)] == smiddle[(i-1,j)]-11:
                backtracklow[i,j] = "↘"
            supper[i,j] = max([supper[(i,j-1)]-1, smiddle[(i,j-1)]-11])
            if supper[(i,j)] == supper[(i,j-1)]-1:
                backtrackup[i,j] = "→"
            elif supper[(i,j)] == smiddle[(i,j-1)]-11:
                backtrackup[i,j] = "↘"
            smiddle[i, j]= max([supper[(i,j)], slower[(i,j)], smiddle[i-1, j-1]+BLOSUM62[v[i-1]][w[j-1]]])
            if smiddle[i, j] == smiddle[i-1, j-1]+BLOSUM62[v[i-1]][w[j-1]]:
                backtrack[i,j] = "↘"
            elif smiddle[i, j] == slower[(i,j)]:
                backtrack[i,j] = "↓"
            elif smiddle[i, j] == supper[(i,j)]:
                backtrack[i,j] = "→"
    return backtrack, backtracklow, backtrackup, smiddle[len(v), len(w)]
def AffineGapPenaltiesProblem(v, w):
    """
    Input: Two amino acid strings v and w (each of length at most 100).
    Output: The maximum alignment score between v and w, followed by an alignment of v and w achieving this maximum score. Use the BLOSUM62 scoring matrix, a gap opening penalty of 11, and a gap extension penalty of 1.
    """
    backtrack, backtracklow, backtrackup, score = AffineGapPenaltiesBacktrack(v, w)
    i = len(v)
    j = len(w)
    newV, newW = "", ""
    node = (i,j)
    while node != (0, 0):
        if backtrack[i,j] == "↘":
            newV = newV+ v[i-1]
            newW = newW + w[j-1]
            i = i-1
            j = j-1
        elif backtrack[i,j] == "→":
            newV= newV + "-"
            newW= newW + w[j-1]
            i1, j1 = i, j
            while backtrackup[i1, j1] != "↘":
                j1 = j1-1
                newV= newV + "-"
                newW= newW + w[j1-1]
            j = j1 -1
        elif backtrack[i,j] == "↓":
            newV= newV + v[i-1]
            newW= newW + "-"
            i2, j2 = i, j
            while backtracklow[i2, j2] != "↘":
                i2 = i2- 1
                newV= newV + v[i2-1]
                newW= newW + "-"
            i = i2 - 1
        node = (i,j)
    return newV[::-1], newW[::-1], score



def FromSource(v, w):
    s = {(0,0): 0}
    maxColumn = {}
    best = -10000
    bestn = {}
    if w == "":
        for i in range(1, len(v)+1):
            s[(i, 0)] = s[(i-1, 0)] -5
            bestn = s
    else:
        for i in range(1, len(v)+1):
            s[(i, 0)] = s[(i-1, 0)] -5
        for j in range(1, len(w)+1):
            s[(0, j)] = s[(0, j-1)] -5
        for i in range(1, len(v)+1):
            for j in range(1, len(w)+1):
                a =s[(i-1, j)] -5
                b = s[(i,j-1)] - 5
                c = s[(i-1, j-1)] + BLOSUM62[v[i-1]][w[j-1]]
                s[(i, j)] = max([ a, b ,c ])
                if j == len(w):
                    bestn[(i,j)]= s[(i,j)]

    return bestn
def ToSink(v, w):
    v= v[::-1]
    w = w[::-1]
    s = {(0,0): 0}
    maxColumn = {}
    best = -10000
    bestn = {}
    for i in range(1, len(v)+1):
        s[(i, 0)] = s[(i-1, 0)] -5
    for j in range(1, len(w)+1):
        s[(0, j)] = s[(0, j-1)] -5
    for i in range(1, len(v)+1):
        for j in range(1, len(w)+1):
            a =s[(i-1, j)] -5
            b = s[(i,j-1)] - 5
            c = s[(i-1, j-1)] + BLOSUM62[v[i-1]][w[j-1]]
            s[(i, j)] = max([ a, b ,c ])
    return bestn


def MiddleEdgeinLinearSpaceProblem(v, w):
    """
    Input: Two amino acid strings.
    Output: A middle edge in the alignment graph in the form "(i, j) (k, l)", where (i, j) connects to (k, l). To compute scores, use the BLOSUM62 scoring matrix and a (linear) indel penalty equal to 5.
    """
    if len(w) == 1:
        best1 = FromSource(v, "")
        best2 = ToSink(v,w, b)
    else:
        if len(w)%2 == 0:
            best1 = FromSource(v, w[0:len(w)//2])
            best2 = ToSink(v,w[len(w)//2:])
        else:
            best1 = FromSource(v, w[0:len(w)//2])
            best2 = ToSink(v,w[len(w)//2:])
    return best1, best2
# print(MiddleEdgeinLinearSpaceProblem("PLEASANTLY", "MEASNLY"))

def LinearSpaceAlignment(v, w):
    """
    Input: Two long (10000 amino acid) protein strings written in the single-letter amino acid alphabet.
    Output: The maximum alignment score of these strings, followed by an alignment achieving this maximum score. Use the BLOSUM62 scoring matrix and indel penalty σ = 5.
    """

    backtrack = ""
    top , bottom = 1, len(w)
    right, left = 1 , len(v)
    if right == left:
        v = ""
        return (bottom - top)*"-"
    if top == bottom:
        w = ""
        return (right-left)*"-"
    a, b, edge, score = MiddleEdgeinLinearSpaceProblem(v, w)
    v1= v[:b[1]]
    v2 = v[a[1]:]
    w1 = w[:b[0]]
    w2 = w[a[0]:]
    backtrack+= edge
    return LinearSpaceAlignment(v1, w1) + backtrack + LinearSpaceAlignment(v2, w2)
