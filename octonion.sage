### This module defines the Octonion class.


import itertools
import thread
import cPickle

class o(list):
    """
    Octonion class.
    Representation of an element of the 8-dimensional octonion algebra.
    The first 4 coordinates correspond to 1, i, j, and k of the standard quaternion subalgebra.

    This class overloads many standard operators to provide intuitive methods for performing
    computations within the octonion algebra.
    """

    halving_sets = [[1,2,4,5],[1,4,3,6], [1,3,5,7], [0,5,6,1], [1,6,7,2], [0,7,1,4], [0,1,2,3],
                    [0,3,6,7], [0,2,5,7], [0,2,4,6], [2,4,3,7], [0,4,3,5], [2,3,5,6], [4,5,6,7]]
    for i in range(len(halving_sets)):
        halving_sets[i] = list.sort(halving_sets[i])
    """
    This defines a collection of halving sets used to define the (nonassociative) ring of octonion integers.
    """

    def __init__(self,a0=0,a1=0,a2=0,a3=0,a4=0,a5=0,a6=0,a7=0):
        if type(a0) == type(list()) or type(a0) == type(o([0,0,0,0,0,0,0,0])):
            if len(a0) != 8:
                raise TypeError
            else:
                list.__init__(self)
                self.extend(a0)
        else:
            list.__init__(self)
            self.extend([a0,a1,a2,a3,a4,a5,a6,a7])

        self.divs = dict()
        self.pdivs = dict()
        self.eucs = dict()

    def norm(self):
        n = 0
        for i in self:
            n = n + i^2
        return n

    def halving_set(self):
        """
        Returns the halving set of self if it is an octonion integer.
        Otherwise, returns None.
        """
        ints = list()
        for i in range(8):
            if self[i] in ZZ:
                ints.append(1)
            elif 2*self[i] not in ZZ:
                return None
            else:
                ints.append(0)
        return [i for i in range(8) if ints[i] == 0]

    def __neg__(self):
        a = deepcopy(self)
        for i in range(8):
            a[i] = -a[i]
        return a

    def __add__(self, b):
        b = o(b)
        c = o()
        for i in range(8):
            c[i] = self[i] + b[i]
        return c

    def __radd__(self,b):
        return self + o(b)

    def __sub__(self,b):
        return self + -o(b)

    def __rsub__(self,b):
        return o(b) - self

    def __mul__(self,b):
        if b in RR or type(b) == sage.symbolic.expression.Expression:
            c = o()
            for i in range(8):
                c[i] = self[i]*b
            return c
        elif self in o.e and b in o.e:
            a = o.e.index(self)
            d = o.e.index(b)
            if a == 0:
                return o.e[d]
            elif d == 0:
                return o.e[a]
            elif a == d:
                return -o.e[0]
            elif a == 1:
                if d == 2:
                    return o.e[3]
                elif d == 3:
                    return -o.e[2]
                elif d == 4:
                    return o.e[5]
                elif d == 5:
                    return -o.e[4]
                elif d == 6:
                    return -o.e[7]
                elif d == 7:
                    return o.e[6]
            elif a == 2:
                if d == 1:
                    return -o.e[3]
                elif d == 3:
                    return o.e[1]
                elif d == 4:
                    return o.e[6]
                elif d == 5:
                    return o.e[7]
                elif d == 6:
                    return -o.e[4]
                elif d == 7:
                    return -o.e[5]
            elif a == 3:
                if d == 1:
                    return o.e[2]
                elif d == 2:
                    return -o.e[1]
                elif d == 4:
                    return o.e[7]
                elif d == 5:
                    return -o.e[6]
                elif d == 6:
                    return o.e[5]
                elif d == 7:
                    return -o.e[4]
            elif a == 4:
                if d == 1:
                    return -o.e[5]
                elif d == 2:
                    return -o.e[6]
                elif d == 3:
                    return -o.e[7]
                elif d == 5:
                    return o.e[1]
                elif d == 6:
                    return o.e[2]
                elif d == 7:
                    return o.e[3]
            elif a == 5:
                if d == 1:
                    return o.e[4]
                elif d == 2:
                    return -o.e[7]
                elif d == 3:
                    return o.e[6]
                elif d == 4:
                    return -o.e[1]
                elif d == 6:
                    return -o.e[3]
                elif d == 7:
                    return o.e[2]
            elif a == 6:
                if d == 1:
                    return o.e[7]
                elif d == 2:
                    return o.e[4]
                elif d == 3:
                    return -o.e[5]
                elif d == 4:
                    return -o.e[2]
                elif d == 5:
                    return o.e[3]
                elif d == 7:
                    return -o.e[1]
            elif a == 7:
                if d == 1:
                    return -o.e[6]
                elif d == 2:
                    return o.e[5]
                elif d == 3:
                    return o.e[4]
                elif d == 4:
                    return -o.e[3]
                elif d == 5:
                    return -o.e[2]
                elif d == 6:
                    return o.e[1]
        else:
            c = o()
            for i in range(8):
                for j in range(8):
                    c = c + (o.e[i]*o.e[j])*(simplify(self[i]*b[j]))
            for i in range(8):
                if c[i] not in QQ:
                    c[i] = c[i].maxima_methods().rootscontract().simplify()
            return c

    def __rmul__(self,b):
        return o(b)*self

    def real(self):
        return o(self[0])

    def imag(self):
        return self - o(self[0])

    def conj(self):
        return self - 2*self.imag()

    def __abs__(self):
        return sqrt(self.norm())

    def __div__(self,b):
        return self * o(b).conj() * (1/o(b).norm())

    def __rdiv__(self,b):
        return o(b) / self

    def __pow__(self,b):
        if b == -1:
            return o(1)/self
        elif b in ZZ:
            if b == 0:
                return o(1)
            if b > 0:
                return self * (self^(b-1))
            if b < 0:
                return (o(1)/self) * self^(b+1)

    def dot(self, b):
        b = o(b)
        c = 0
        for i in range(8):
            c += self[i]*b[i]
        return c

    def angle(self, b):
        return arccos(self.dot(b) / (abs(self) * abs(b)))

    def __mod__(self, m):
        a = o()
        for i in range(8):
            a[i] = self[i] % m
        return a

    def is_integer(self):
        """
        Return true if and only if self is an octonion integer.
        """
        half = []
        for i in range(8):
            try:
                ZZ(self[i])
            except:
                try:
                    ZZ(2*self[i])
                    half.append(i)
                except:
                    return False
        return len(half) == 0 or len(half) == 8 or half in o.halving_sets

    def ldivides(self,b):
        return ((1/self)*b).is_integer()

    def rdivides(self,b):
        return (b*(1/self)).is_integer()

    def denominator(self):
        """
        Returns the smallest possible rational integer denominator
        of self when expressed as self/d.
        """
        if self.is_integer():
            return o(1)
        dens = []
        for i in self:
            dens.append(denominator(i))
        r = lcm(dens)
        if r % 2 == 0 and (self*(r/2)).is_integer():
            return o(r/2)
        else:
            return o(r)

    def numerator(self):
        """
        Returns the numerator when self is expressed as above.
        """
        return self*self.denominator()

    # def get_divs(self, n):
    #     if ZZ(n).divides(self.norm()) == False:
    #         return None
    #     ldivs = []
    #     rdivs = []
    #     for el in all_norm(n):
    #         if el.ldivides(self):
    #             ldivs.append(el)
    #         if el.rdivides(self):
    #             rdivs.append(el)
    #     return [ldivs,rdivs]
    # def princ_div(self, norm, right = False):
    #     divs = self.get_divs(norm)[int(right)]

    #     i = 0
    #     found = False

    #     while found == False:
    #         wrong = False
    #         div = divs[i]
    #         j = 0
    #         while wrong == False:
    #             u = o.units[j]
    #             if (u*div).ldivides(self) == False:
    #                 wrong = True
    #             j += 1
    #         if wrong == False:
    #             found = True
    #         i += 1
    #     return div

    def is_princ_div(self, div, norm, right = False):
        """
        Returns true if and only if div is a principal divisor
        (by default, left) of self, of the given norm.
        """
        if norm in self.pdivs:
            return(div in self.princ_div(norm))

        divs = self.get_divs(norm)[int(right)]
        j = 0
        wrong = False
        while wrong == False and j < len(divs):
            div2 = divs[j]
            if not div.ldivides(div2):
                wrong = True
            j += 1
        return(not wrong)

    def princ_div(self, norm, right = False):
        """
        Returns a principle divisor (by default, left) of self
        of the given norm.
        """
        princ = []
        divs = self.get_divs(norm)[int(right)]

        if norm in self.pdivs:
            for i in self.pdivs[norm]:
                princ.append(divs[i])
            return princ
        else:
            for i in range(len(divs)):
                print(str(i) + "   " + str(len(princ)))
                div = divs[i]
                if self.is_princ_div(div,norm, right = right):
                    princ.append(i)
        self.pdivs[norm] = princ
        return princ

    def round(self):
        """
        Returns the octonion integer that lies closest to self,
        by standard Euclidean distance
        """
        choice = []
        for i in range(8):
            x = floor(2*self[i])
            choice.append([x/2,(x+1)/2])
        for y in itertools.product(*choice):
            y = o(list(y))
            if y.is_integer() and (self - y).norm() <= 1/2:
                return y

    def floor_div(self, div, right = False):
        """
        Divides self by div and returns the integral
        part of the quotient.
        """
    	if right:
    		return (self/div).round()
    	else:
        	return ((div)^(-1) * self).round()

    def euc_remainders(self, M):
        """
        Returns (gamma, rho, m) triples.
        """
        if M in self.eucs:
        	return self.eucs[M]
        else:
	        rho = self
	        m = M
	        trips = []
	        while rho != o(0):
	            gamma = rho.floor_div(m)
	            trips.append([gamma,rho, m])
	            rho = (rho - m*gamma).conj()
	            m = rho.norm() / m
	    	self.eucs[M] = trips
        return trips


    def get_divs(self, M, right = False):
        """
        Returns the list of all octonion divisors of self
        of norm M.
        """
        if M not in self.divs:
            if right == False:
                M = self.norm() / M
            trips = list(reversed((self.euc_remainders(M))))
            trip = trips[0]
            m = trip[2]

            norms = []

            if m == 1:
                norms = o.units

            ldivs = []
            rdivs = []

            for mu_N in norms:
                mus = [mu_N]

                gamma = trip[0]
                rho = trip[1]
                m = trip[2]

                mus.append(gamma*mu_N)

                for tri in trips[1:]:
                    gamma = tri[0]
                    rho = tri[1]
                    m = tri[2]
                    k = len(mus)

                    mus.append(gamma * mus[k - 1] + mus[k-2])

                ldivs.append(mus[len(mus) - 1])
                rdivs.append((mus[len(mus) - 2]).conj())

            if right == False:
                M = self.norm() / M

            self.divs[M] = [ldivs, rdivs]
        
        return self.divs[M]

    def common_divs(self, x, m, right = False):
        """
        Returns a list of all common divisors of
        self and x having norm m.
        """
        common = []
        sdivs = self.get_divs(m)[int(right)]
        xdivs = x.get_divs(m)[int(right)]
        for a in sdivs:
            if a in xdivs:
                common.append(a)
        return common

    def euc_divide(self, divisor, right = True):
        """
        Returns the quotient and remainder obtained
        when dividing self by divisor.
        """
    	if right == False:
    		gamma = self.floor_div(divisor)
    		rho = divisor*(divisor^(-1)*self - gamma)
    	else:
    		gamma = self.floor_div(divisor, True)
    		rho = (self*divisor^(-1) - gamma)*divisor
    	return [gamma,rho]

