import time
import thread

#path = "/Users/kida/Workspace/Embeddings/"
path = "/Users/kida/Desktop/Workspace/Embeddings/"

class ListOfTrianglesByCharacteristic(list):
    """
    This class is a list specialized for storing and handling Triangle objects.
    These Triangles are organized by characteristic.  Methods are provided for
    embedding triangles in a systematic way, then storing and outputting the
    results.
    """
    def __init__(self, side_upper_bound = 0):
        list.__init__(self)
        self.side_upper_bound = side_upper_bound
        self.searching_time = 0
        self.embedding_time = 0
        self.num_tris = []
        self.ratemb = []
        self.latemb = []
        self.ratprop = []
        self.latprop = []
        proof.number_field(False)
        self.invariant_factors = load("/Users/kida/Workspace/Embeddings/invariantFactors.sobj")
        ## Store class numbers and the invariant factors in the
        ## decomposition of the class groups corresponding to
        ## each characteristic.

    def __str__(self, tris = False, emb = True, class_num = True, inv = True):
        """
        The string method.

        Keyword arguments:

        tris -- determines whether lists of triangles are to be printed
        emb  -- determines whether the embedding proportion is to be printed
        """

        s = ""
        for i in range(len(self)):
            if i !=0 and squarefree_part(i) == i:
                s += self.char_str(i, tris, emb, class_num, inv)
        s += "Searching time: " + str(self.searching_time) + "\n" + "Embedding time: " + str(self.embedding_time)
        return s

    def char_str(self, char, tris = True, emb = True, class_num = False, inv = True):
        """
        Returns a string containing information about the stored
        triangles of this characteristic.

        Keyword arguments:
        char -- the characteristic of the triangles
        tris -- if true, then the triangles are printed
        emb  -- if true, then the embeddings are printed
        """

        s = ""

        ## Check that this characteristic list has triangles.
        ## Check if the embedding proportion is supposed to be printed and,
        ## if so, then check that triangles have been embedded - else,
        ## embedding proportion would be trivial.

        if self.num_tris[char] != 0:
            s += "\n"
            s += "Characteristic: " + str(char) + "\n"
            s += "# of Triangles: " + str(self.num_tris[char]) + "\n"
            if self.latprop[char] != None and emb != False:
                s += "% Field   Embedded: " + str(n(self.ratprop[char])) + "\n"
                s += "% Lattice Embedded: " + str(n(self.latprop[char])) + "\n"
            if class_num:
                s += "Class Number: " + str(self.class_num(char)) + "\n"
            if inv:
                s += "Invariants: " + str(self.invariant_factors[char]) + "\n"
            s += "\n"
            if tris:
                for j in self[char]:
                    s += str(j) + "\n"
            s += "\n"
        return s

#     def embed_char(self, char, proof=True):
#         """
#         Attempts to lattice embed all triangles of characteristic char.
#         Also stores the proportion of triangles that embedded.  This
#         probability is safely updated if new triangles are later added to
#         this list and then embedded.  Also stores the time taken by process.
#         """
#
#         start_time = time.time()
#
#         if char < len(self) and len(self[char]) != 0:
#             if self.probs[char] == None:
#                 self.probs[char] = 0
#             num_new = 0
#             num_embedded = 0
#             startTime = time.time()
#             for j in self[char]:
#                 if j.embedding == 0:
#                     if j.lattice_embed(proof = proof) != None:
#                         num_embedded += 1
#                     if j.rational_embedding != None:
#                         num_new += 1
#                         self.ratemb[char] += 1
#                     print("Characteristic: " + str(char) + "               " + str(j))
#             self.probs[char] = (self.probs[char]*(self.ratemb[char] - num_new) + num_embedded)/self.ratemb[char]
#             self.embedding_time += time.time() - start_time

    def embed_char(self, char, proof=True):
        """
        Attempts to lattice embed all triangles of characteristic char.
        Also stores the proportion of triangles that embedded.  This
        probability is safely updated if new triangles are later added to
        this list and then embedded.  Also stores the time taken by process.
        """

        if char < len(self) and len(self[char]) != 0:
            start_time = time.time()
            num_new = 0
            num_newrat = 0
            num_newlat = 0
            for j in self[char]:
                if j.embedding == 0:
                    if j.rational_embed() != None:
                        num_newrat += 1
                        if j.lattice_embed(proof = proof) != None:
                            num_newlat += 1
                    print("Characteristic: " + str(char) + "               " + str(j))

            self.ratemb[char] += num_newrat
            self.latemb[char] += num_newlat

            self.ratprop[char] = self.ratemb[char] / self.num_tris[char]

            if self.ratemb[char] != 0:
                self.latprop[char] = self.latemb[char] / self.ratemb[char]

            print(self.num_tris[char], self.ratemb[char], self.latemb[char])
            print(self.ratprop[char], self.latprop[char])

            self.embedding_time += time.time() - start_time

    def embed_range(self, start, end, proof = True):
        """
        Embeds all characteristics from start to end (both inclusive).
        """
        self.embed_set(range(start, end+1), proof = proof)

    def embed_set(self, A, proof = True):
        """
        Embeds all characteristics listed in A.  Will abort on press
        of "return" without corrupting data.
        """
        def input_thread(L):
            raw_input()
            L.append(None)

        def run_loop():
            L = []
            thread.start_new_thread(input_thread, (L,))
            for char in A:
                if L: break
                self.embed_char(char, proof = proof)
        run_loop()

    def class_num(self, n):
        """
        Returns the class number of characteristic n.
        """
        class_num = 1
        for factor in self.invariant_factors[n]:
            class_num = class_num * factor
        return class_num
#         return QQ[sqrt(-n)].class_group(proof=False).order()

    def inv(self, n):
        """
        Returns the invariant factors of the class group for
        characteristic n.
        """
        return self.invariant_factors[n]
#         return QQ[sqrt(-n)].class_group(proof=False).invariants()

    def cyclic_twos(self, char):
        """
        Returns True if the class group of char is a direct product
        of cyclic two groups.
        """
        cyclic_twos = True
        for factor in self.inv(char):
            if factor != 2 and factor !='':
                cyclic_twos=False
        return cyclic_twos

    def save(self, filename = "tris"):
        """
        Saves self to default directory.  Default filename is "tris".
        """
        save_path = path + filename
        save(self, save_path)




###############################################################################
###############################################################################
###############################################################################
###############################################################################

class Triangle:

    """This class contains a triangle and its attributes."""
    def __init__(self, a, b, c):
        """ Initializes dictionary of side lengths.
            Sorts side lengths by size.
            Calls update to set area and characteristic.
        """
        self.sides = sorted([a^2,b^2,c^2])
        self.sides = [sqrt(self.sides[0]), sqrt(self.sides[1]), sqrt(self.sides[2])]
        self.base = None
        self.update()
        self.embedding = 0
        self.embeds = None

    def __str__(self):
        """
        Returns a string of this triangle's list of sides.
        """
        sides = str(self.sides)
        if self.embedding != 0:
            embedding = str(self.embedding)
        else:
            embedding = "-"
        return '{0:50}{1:100}'.format(sides, embedding)

    def update(self):
        """
        Updates the area and characteristic.
        """
        a = self.sides[0]
        b = self.sides[1]
        c = self.sides[2]
        radicand = Integer(-a^4 - b^4 - c^4 + 2*(a^2*b^2 + a^2*c^2 + b^2*c^2))
        self.characteristic = radicand.squarefree_part()
        self.area = sqrt(radicand)/4
        self.rational_embedding = self.rational_embed()

    def class_group(self):
        """Returns class group of this triangle's characteristic."""
        return QQ[sqrt(-self.characteristic)].class_group()

    def class_number(self):
        """Returns class number of this triangle's characteristic."""
        return self.class_group().order()

    def rational_embed(self):
        """
        Returns a list of this triangle's coordinates when
        embedded in QQ[sqrt(-D)]

        Keyword arguments:
        """
        a = self.sides[2]
        b = self.sides[1]
        c = self.sides[0]
        Q.<x> = QQ[sqrt(-self.characteristic)]
        Z = Q.maximal_order()

        alpha = gen_norm_reps(a^2, self.characteristic)

        if len(alpha) == 0:
            return None
        else:
            alpha = alpha[0]

        h = simplify(2*self.area/a)
        w = simplify(sqrt(b^2 - h^2))

        beta_imagcoef = QQ(simplify(h/simplify(sqrt(self.characteristic)*sqrt(alpha.norm()))))
        beta_real = QQ(sqrt(simplify(w/sqrt(alpha.norm()))^2))
        beta = Q(alpha) * (beta_real + beta_imagcoef*x)
#         beta = Q(alpha) * (QQ(simplify(w/sqrt(alpha.norm()))) + QQ(simplify(h/simplify(sqrt(self.characteristic)*sqrt(alpha.norm()))))*x)
        return [0, alpha, beta]



    def lattice_embed(self, base_ind = None, proof = True):
        """
        Attempts to store and return a lattice embedding for this triangle.
        If the triangle cannot even embed in the rational field extension,
        then 'No rational embedding' is returned.
        If the triangle embeds rationally, but not in the maximal order, then
        'None' is returned.
        Else, the lattice embedding is returned.
        """

        if self.rational_embedding == None:
            self.embedding = "No rational embedding"
            return "No rational embedding"

        if proof == True:
            self.embedding = self.proof_embed()
            return self.embedding

        Q.<x> = QQ[sqrt(-self.characteristic)]
        Z = Q.maximal_order()
        new_emb = [None, None, None]
        A = gen_norm_reps(self.sides[self.base]**2, self.characteristic)
        ## Main algorithm
        try:
            if base_ind == None:
                rat_emb = self.rational_embedding
            else:
                rat_emb = self.rational_embed(base_ind)
            b = QQ(rat_emb[2][0]) + QQ(rat_emb[2][1]/sqrt(self.characteristic))*x
            ## Loop through all elements of magnitude a
            for a in A:
                transformation = a/sqrt(norm(a))
                b = QQ(rat_emb[2][0]) + QQ(rat_emb[2][1]/sqrt(self.characteristic))*x
                new_b = transformation * b
                if new_b.is_integral():
                    new_emb[0] = 0 + 0*x
                    new_emb[1] = a
                    new_emb[2] = new_b
                    self.embedding = new_emb
                    return new_emb

                transformation = Z(a.conjugate())/sqrt(norm(a))
                new_b = transformation * b
                if new_b.is_integral():
                    new_emb[0] = 0 + 0*x
                    new_emb[1] = Z(a.conjugate())
                    new_emb[2] = new_b
                    self.embedding = new_emb
                    return new_emb

        except TypeError:
            pass
        self.embedding = None
        return None



###############################################################################

#                                FUNCTIONS                                    #

###############################################################################

    def IJK(self):
        """
        Returns (I,J,K) or (I,K) where:

        I is the factor I identified in the proof of the main lemma.

        If triangle has a lattice embedding, then returns the factor
        J identified in the proof of the main theorem.

        K is the GCD of K_1, ..., K_n as identified in the main theorem.
        """
        Q.<x> = QQ[sqrt(-self.characteristic)]
        Z = Q.maximal_order()

        rational = self.rational_embedding

        if rational == None:
            return None

        point = rational[2]
        r = point.denominator()

        ## Check D=3(mod4) case
        if r % 2 == 0 and self.characteristic % 4 == 3:
            r = r/2

        alpha = Z.ideal(Z(point*r))
        base = Z.ideal(rational[1] * r)
        r = Z.ideal(r)
        I = Z.ideal(1)
        factors = factor(r)


        for pair in factors:
            fact = pair[0]
            mult = pair[1]
            if (fact^2).divides(alpha):
                I = I*fact^mult
                alpha = alpha/fact^(2*mult)
        alpha_K = alpha
        base_K = base / I^2
        K = Z.ideal(alpha_K.gens() + base_K.gens())
        divisors = get_divisors(K)
        for J in divisors:
            prod = (I*J)^2
            if prod.is_principal():
                return (I,J,K)
        return (I, None, K)

    def check_embeds(self):
        if self.embeds == True:
            return True
        if self.embedding == None:
            self.embeds = False
            return False
        IJK = self.IJK()
        J = IJK[1]
        if J == None:
            self.embeds = False
            return False
        else:
            self.embeds = True
            return True

    def proof_embed(self):
        Q.<x> = QQ[sqrt(-self.characteristic)]
        Z = Q.maximal_order()

        IJK = self.IJK()
        I = IJK[0]
        J = IJK[1]

        if J == None:
            return None

        IJconj_gens = []
        for gen in (I*J).gens_reduced():
            IJconj_gens.append(gen.conjugate())

        IJconj = Z.ideal(IJconj_gens)
        rot = (IJconj / (I*J)).gens_reduced()[0]

        return [0, rot*self.rational_embedding[1], rot*self.rational_embedding[2]]





def generate(max_side, max_char, tris_by_char=None, integers=False, primitive = True):
    """
    Generates all primitive, incongruent integer triangles
    with side lengths not exceeding max_side and
    characteristics not exceeding max_char.  Organizes these
    triangles into a ListOfTrianglesByCharacteristic

    Keyword arguments:

    max_side -- the maximum side length
    max_char -- the maximum characteristic
    tris_by_char -- a list of triangles organized by characteristic.
                  If a previously generated list is given as an argument,
                  that list will be extended by this method.  Otherwise,
                  a new list is initialized with a max_side attribute of 0.
    """

    if tris_by_char == None:
        tris_by_char = ListOfTrianglesByCharacteristic()

    start_time = time.time()

    ## Setup functions based on whether integer triangles or integer-normed
    ## triangles are to be generated.
    if integers:
        f = lambda x: x
        c_start = lambda a,b: floor(a-b) + 1
        triangle = lambda a,b,c: Triangle(c,b,a)
#         def good_tri(a,b,c):
#             return integers and gcd((a,b,c)) == 1
    else:
        f = lambda x: sqrt(x)
        c_start = lambda a,b: floor(a-2*sqrt(a*b)+b) + 1
        triangle = lambda a,b,c: Triangle(sqrt(c),sqrt(b),sqrt(a))
#         def good_tri(a,b,c):
#             return gcd([a/ZZ(a).squarefree_part(),b/ZZ(b).squarefree_part(),c/ZZ(c).squarefree_part()]) == 1 \
#                 and (ZZ(a).is_square() or ZZ(b).is_square() or ZZ(c).is_square())

    def good_tri(a,b,c):
        if primitive:
            return gcd([a,b,c]) == 1
        else:
            return True

    def generateLayer(a, tris_by_char):
        """
        Generates all primitive, incongruent "layer"
        of triangles with greatest side length a.
        """
        ## Generating Algorithm
        for b in range(1,a+1):
            for c in range(c_start(a,b), b+1):
                ## Check primitive and contains at least one integer side
                if good_tri(a,b,c):
                    print(f(a),f(b),f(c))
                    tri = triangle(a,b,c)
                    char = tri.characteristic
                    if char <= max_char:
                        tris_by_char[char].append(tri)

        return tris_by_char

    ## Beginning of outer function
    start = tris_by_char.side_upper_bound + 1   ## The starting value for a
    tris_by_char.side_upper_bound = max_side    ## In case this list is later extended
    ## Populate lists with empty lists and zeros
    while len(tris_by_char) <= max_char:
        tris_by_char.append([])
        tris_by_char.num_tris.append(0)
        tris_by_char.ratemb.append(0)
        tris_by_char.latemb.append(0)
        tris_by_char.ratprop.append(None)
        tris_by_char.latprop.append(None)

    ## Outer loop of triangle generation
    def input_thread(L):
        raw_input()
        L.append(None)

    def run_loop(tris):
        L = []
        thread.start_new_thread(input_thread, (L,))
        for a in range(start, max_side + 1):
            if L:
                tris_by_char.side_upper_bound = a-1    ## In case this list is later extended
                break
            tris = generateLayer(a, tris)
        return tris

    tris = run_loop(tris_by_char)

    tris_by_char.searching_time = tris_by_char.searching_time + time.time() - start_time

    for char in range(len(tris_by_char)):
        tris_by_char.num_tris[char] = len(tris_by_char[char])

    return tris_by_char





def gen_norm_reps(n, d):
    """
    This function returns a list of all elements in QQ[sqrt(-d)]
    that have norm n, but are neither associates nor conjugates of
    associates of each other.
    The algorithm only checks points which are within a predetermined
    arc ranging counter-clockwise from the x-axis.  Within this arc,
    only one representative of each "conjugate-associate class"
    exists.
    The algorithm checks each half-integer within the appropriate
    range of x-values, then evaluates their corresponding y-coordinate
    on the circle of radius n centered at the origin.  If the number
    x + y*sqrt(-d) that is produced is an element of the integer ring,
    then this element is stored in a list.  This list is returned.

    Keyword arguments:

    n -- all representative elements of norm n will be returned
    d -- the discriminant (the function can take both D or 4D;
         the 4 will be divided out if present)

    returns
    """

    ## Initialize rings and list of elements
    Q.<x> = QQ[sqrt(-d)]
    Z = Q.maximal_order()
    elements = []

    ## Trivial case
    if n == 0:
        return [Z(0)]

    ## Initialize the angle t for Gaussian case, Eisenstein case,
    ## and all other cases
    if   d == 1:
        t = cos(pi/4)
    elif d == 3:
        t = cos(pi/6)
    else:
        t = cos(pi/2)

    ## Main algorithm.
#    for A in range(floor(2*sqrt(n)*t + 1), floor(2*sqrt(n))+1):
    if d%4 == 3:
        for A in range(ceil(2*sqrt(n)*t), floor(2*sqrt(n))+1):
            a = A/2
            b = sqrt(n-a**2)
            ## Try to coerce a+b*x into rational ring
            try:
                r = QQ(a) + QQ(b/sqrt(d))*x
                if r.is_integral():
                    elements.append(r)
            except TypeError:
                pass
        return elements

    else:
        for a in range(ceil(sqrt(n)*t), floor(sqrt(n))+1):
            b = sqrt(n-a**2)
            ## Try to coerce a+b*x into rational ring
            try:
                r = QQ(a) + QQ(b/sqrt(d))*x
                if r.is_integral():
                    elements.append(r)
            except TypeError:
                pass
        return elements


def get_divisors(prod):
    """
    Returns a list of all divisors of prod.

    Keyword Arguments:
    prod -- the product being factored, assumed to be an ideal
    """

    Z = prod.ring()
    f = factor(prod)
    factors = []            ## List of prime factors of prod
    for fac in f:
        for mult in range(fac[1]):
            factors.append(fac[0])
    power = list(powerset(factors))     ## Multiset of factors
#     print(power)
    divisors = []
    for ls in power:
        prod = Z.ideal(1)
        for fac in ls:
            prod = prod*fac
        divisors.append(prod)
    divisors = list(set(divisors))
    return divisors

def cyclic_twos(D):
        """
        Returns True if the class group of char is a direct product
        of cyclic two groups.
        """
        inv = load("/Users/kida/Workspace/Embeddings/invariantFactors.sobj")
        is_cyclic_twos = True
        for fact in inv[D]:
            if fact != 2 and fact !='':
                is_cyclic_twos=False
        return is_cyclic_twos


# tris = generate(75,199, primitive = True)
# tris.embed_range(1,199)