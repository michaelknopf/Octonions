### This module contains functions I used often while testing hypthoses
### regarding integral octonions.


# Instantiate units
o.e = list()
for i in range(8):
    o.e.append(o())
for i in range(8):
    o.e[i][i] = 1
e = deepcopy(o.e)



def o_listnorms(n, order = 2):
    """
    Generates all elements of the desired order.

    order:
    Integral Octonions: 2
    Integral Quaternions: 1
    """

    onorms = list()
    for i in range(n+1):
        onorms.append(list())

    def norm(x):
        m = 0
        for i in x:
            m = m+i^2
        return m

    half4 = o(0,1/2,1/2,1/2,0,1/2,0,0)
    half8 = o(1/2,1/2,1/2,1/2,1/2,1/2,1/2,1/2)

    for a1 in range(0, floor(sqrt(n)) + 1):
        for a2 in range(a1, floor(sqrt(n - norm([a1]))) + 1):
            for a3 in range(a2, floor(sqrt(n - norm([a1, a2]))) + 1):
                for a4 in range(a3, floor(sqrt(n - norm([a1, a2, a3]))) + 1):
                    for a5 in range(a4, floor(sqrt(n - norm([a1, a2, a3, a4]))) + 1):
                        for a6 in range(a5, floor(sqrt(n - norm([a1, a2, a3, a4,a5]))) + 1):
                            for a7 in range(a6, floor(sqrt(n - norm([a1, a2, a3, a4, a5, a6]))) + 1):
                                for a8 in range(a7, floor(sqrt(n - norm([a1, a2, a3, a4, a5, a6, a7]))) + 1):
                                    x = o(a1,a2,a3,a4,a5,a6,a7,a8)
                                    onorms[norm(x)].append(x)
                                    if order >= 1:
                                        y = x + half8
                                        m = norm(y)
                                        if m <= n:
                                            onorms[m].append(y)
    # Halving sets
    if order == 2:
        for a1 in range(0, floor(sqrt(n)) + 1):
            for a2 in range(a1, floor(sqrt(n - norm([a1]))) + 1):
                for a3 in range(a2, floor(sqrt(n - norm([a1, a2]))) + 1):
                    for a4 in range(a3, floor(sqrt(n - norm([a1, a2, a3]))) + 1):
                        for a5 in range(0, floor(sqrt(max(0, n - norm([a1, a2, a3, a4]))) - 1/2) + 1):
                            for a6 in range(a5, floor(sqrt(max(0, n - norm([a1, a2, a3, a4, a5+1/2]))) - 1/2) + 1):
                                for a7 in range(a6, floor(sqrt(max(0, n - norm([a1, a2, a3, a4, a5+1/2, a6+1/2]))) - 1/2) + 1):
                                    for a8 in range(a7, floor(sqrt(max(0, n - norm([a1, a2, a3, a4, a5+1/2, a6+1/2, a7+1/2]))) - 1/2) + 1):
                                        x = o(a1,a5,a6,a7,a2,a8,a3,a4) + half4
                                        onorms[norm(x)].append(x)
    return onorms


def all_unique(elems):
    """
    Utility function.  Checks whether the list contains duplicates.
    """
    return len(set([tuple(x) for x in elems])) == len(elems)

def sets_equal(l):
    """
    Utility function.  Checks whether the sets within l are all equal.
    """
    if len(l) <= 1:
        return True
    if len(l) == 2:
        s1 = deepcopy(l[0])
        s2 = deepcopy(l[1])
        for i in range(len(s1)):
            s1[i] = tuple(s1[i])
        for i in range(len(s2)):
            s2[i] = tuple(s2[i])
        return(set(s1) == set(s2))
    else:
        return(sets_equal([l[0],l[len(l) - 1]]) and sets_equal(l[:-1]))


def findlast(onorms):
    i = 0
    for j in range(len(onorms)):
        if len(onorms[j]) > 0:
            i = j
    return i

def get_pos_twins(x, order = 3):
    positive_elems = list()
    elems = list()

    if order == 3:
        hsets = o.halving_sets
    if order == 2:
        hsets = [[1,2,3,5], [0,4,6,7]]

    half = x.halving_set()
    if len(half) == 4:
        halves = Permutations([x[i] for i in half])
        ints = Permutations([x[i] for  i in list(set(range(8)) - set(half))])
        for halving_set in hsets:
            inting_set = list(set(range(8)) - set(halving_set))
            for int_perm in ints:
                for half_perm in halves:
                    elem = o()
                    for i in range(4):
                        elem[inting_set[i]] = int_perm[i]
                        elem[halving_set[i]] = half_perm[i]
                    elems.append(o(list(elem)))
    else:
        for elem in Permutations(x):
            elems.append(o(list(elem)))

#     elems2 = list()

#     for elem in elems:
#         print(elem)
#         elems2.extend(pm(elem))

    return elems

def pm(beta):
    normlist = [beta]

    def pmi(elems,i):
        newelems = list()
        for elem in elems:
            newelem = deepcopy(elem)
            newelem[i] = -newelem[i]
            newelems.append(newelem)
        return elems + newelems

    for i in range(len(beta)):
        if beta[i] != 0:
            normlist = pmi(normlist,i)
    return normlist

def get_twins(x, order = 3):
    elems = list()
    for elem in get_pos_twins(x, order):
        elems.extend(pm(elem))
    return elems

o.normlist = o_listnorms(20)
o.units = list()
for elem in o.normlist[1]:
    o.units.extend(get_twins(elem))
o.units[0:16] = list(reversed(o.units[0:16]))
for i in range(8):
    u = o.units[2*i]
    o.units[2*i] = o.units[2*i+1]
    o.units[2*i + 1] = u


def all_norm(n):
    elems = list()
    normlist = o.normlist
    print(str(len(normlist[n])) + " elements")

    count = 1

    for elem in normlist[n]:
        print(count)
        count += 1
        elems.extend(get_twins(elem))
    return elems

def count_norms(n, proof = False):

    if proof == False:

        count = 0

        for elem in o.normlist[n]:
            half_set = elem.halving_set()

            if len(half_set) == 4:
                int_set = list(set(range(8)) - set(half_set))

                ints = [elem[i] for i in int_set]
                halves = [elem[i] for i in half_set]

                ints = [ints.count(i) for i in list(set(ints))]
                halves = [halves.count(i) for i in list(set(halves))]

                m = 1
                for i in ints:
                    m = m*factorial(i)
                for i in halves:
                    m = m*factorial(i)

                m = 14 * factorial(4)^2 / m
                m = m* 2^(8 - elem.count(0))

                count = count + m

            else:
                components = [elem.count(i) for i in list(set(elem))]
                m = 1
                for i in components:
                    m = m*factorial(i)
                m = factorial(8) / m
                m = m* 2^(8 - elem.count(0))

                count = count + m
    else:

        count = 0

        for i in divisors(n):
            count += i^3

        count = count * 240

    return count


def o_embed(tri, order = 3):
    normlist = o.normlist
    a = tri.sides[0]^2
    b = tri.sides[1]^2
    c = tri.sides[2]^2
    counter = 0
    for alpha in normlist[a]:
        counter  = counter + 1
        if counter % 25 == 0:
            print(counter)
        for elem in normlist[b]:
            for beta in get_twins(elem, order):
                norm = 0
                for i in range(len(alpha)):
                    norm = norm + (alpha[i] - beta[i])^2
                if norm == c:
                    return([alpha,beta])
        for elem in normlist[c]:
            for beta in get_twins(elem, order):
                norm = 0
                for i in range(len(alpha)):
                    norm = norm + (alpha[i] - beta[i])^2
                if norm == b:
                    return([alpha,beta])
    return 0


def embed_tris(char, start = 0, m = None, order = 3):
    if len(tris[char]) == 0:
        return None
    print(len(tris[char]))
    if m == None:
        m = len(tris[char])
    for i in range(start, min(m+1, len(tris[char]))):
        tri = tris[char][i]
        tri = Triangle(tri.sides[0],tri.sides[1],tri.sides[2])
        emb = tri.lattice_embed()
        if emb == None or emb == 'No rational embedding':
            print(str(i) + ":   " + str(o_embed(tri, order)))
        else:
            print(str(i) + ":   " + str("embeds"))


def init_vars(vectors = True):
    global x,y,z,w

    x1 = var('x1')
    x2 = var('x2')
    x3 = var('x3')
    x4 = var('x4')
    x5 = var('x5')
    x6 = var('x6')
    x7 = var('x7')
    x8 = var('x8')
    y1 = var('y1')
    y2 = var('y2')
    y3 = var('y3')
    y4 = var('y4')
    y5 = var('y5')
    y6 = var('y6')
    y7 = var('y7')
    y8 = var('y8')
    z1 = var('z1')
    z2 = var('z2')
    z3 = var('z3')
    z4 = var('z4')
    z5 = var('z5')
    z6 = var('z6')
    z7 = var('z7')
    z8 = var('z8')
    w1 = var('w1')
    w2 = var('w2')
    w3 = var('w3')
    w4 = var('w4')
    w5 = var('w5')
    w6 = var('w6')
    w7 = var('w7')
    w8 = var('w8')

    if vectors:
        x = o(x1,x2,x3,x4,x5,x6,x7,x8)
        y = o(y1,y2,y3,y4,y5,y6,y7,y8)
        z = o(z1,z2,z3,z4,z5,z6,z7,z8)
        w = o(w1,w2,w3,w4,w5,w6,w7,w8)

def test(char):
    for tri in tris[char]:
        s = 0
        for side in tri.sides:
            s += side^2
        if s % 2 == 0:
            print(True)

        else:
            z = 16*tri.area^2 / char
            if z % 2 == 0:
                print(True)
            else:
                print(tri.sides)

def getPureMults(t,m):
    l = list()
    x = list()
    for a in range(m):
        print(a)
        for b in range(m):
            for c in range(m):
                for d in range(m):
                    f = (o(a,b,c,d)*t) % m
                    if f[0] == 0:
                        l.append(o(a,b,c,d))
                        x.append(tuple(f))
    n = list(set(x))
    k = list()
    for i in range(len(n)):
        k.append([])
    for i in range(len(n)):
        for j in range(len(l)):
            if tuple(x[j]) == n[i]:
                a = l[j]
                if a.count(0) == 6 and a[0] != 0:
                    k[i].append(a)
    return (l,x,n,k)


def save(self, filename = "tris"):
        """
        Saves self to default directory.  Default filename is "tris".
        """
        save_path = path + filename
        save(self, save_path)

def normdata(max, data = [[0]], saving = True):

    def input_thread(L):
        raw_input()
        L.append(None)

    try:
        o.normlist[max]
    except:
        o.normlist = o_listnorms(max)

    L = []
    thread.start_new_thread(input_thread, (L,))

    for n in range(len(data),max+1):

        if L:
            break

        print("Current iteration: " + str(n))
        normn = all_norm(n)
        for i in range(len(normn)):
            normn[i] = tuple(normn[i])
        data.append(normn)
    return data

def tuple_data(data):
    for i in range(len(data)):
        i = int(i)
        for j in range(len(data[i])):
            j = int(j)
            print(data[i][j])
            data[i][j] = tuple(data[i][j])




# a = [-3, -3/2, -3, 3/2, -3, 3/2, -3, 3/2]
# b = [0, 0, 0, 0, 0, 0, 3, 3]