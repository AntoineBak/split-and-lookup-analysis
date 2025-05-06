
from collections import Counter, defaultdict
from sage.modules.free_module_integer import IntegerLattice

class SplitAndLookup():
    """
    Class implementing Split-and-lookups.

    Constructor:
    * SplitAndLookup(Sboxes, buckets, p)

    Getters:
    * get_p(): int, return the modulus p.
    * get_n(): int, return the number of lookups n.
    * get_small_S(k: int): list of int, return the small S-box Sk.
    * eval_small_S(k: int, x: int): int, evaluate Sk(x).
    * get_buckets(): list of int, return the list of bucket sizes.
    * get_bucket(k): int, return the bucket size sk.

    Methods:
    * eval(x: int): int, evaluate S(x).
    * is_permutation(): Boolean, check if S is a permutation according to the criteria on small S-boxes.
    * inverse(): SplitAndLookup, if S is a permutation, return a SplitAndLookup Object representing the functionnal inverse of S.
    * linear_approximation(a: int, b: int, c: int): int, return the number of 0<=x<p such that a*S(x) + b*x = c.
    * linear_correlation(a: int, b: int): CC, return the (a, b) Fourier coefficient of S.
    """

    def _fill_sbox(S, s):
        """
        Make the S-box length match the size of the bucket.
        * if len(S) > s, truncate S to its first s values
        * if len(S) < s, pad S with S(x) = x for x >= s
        """
        l = len(S)

        if l < s:
            return S + list(range(l, s))
        else:
            return copy(S[:s])

    ## Class constructor.
    def __init__(self, Sboxes, buckets, p):
        """
        Input:
            Sboxes: (list of integers) or (list of list of integers). S-boxes of the S&L function.
            buckets: list of integers. Size of the buckets.
            p: integer. Prime field over which the object is defined.
        """
        # If buckets is an integer, set each bucket size to the same value.
        if buckets in ZZ:
            self._n = ceil(log(p, buckets))
            self._buckets = [buckets for i in range(self._n - 1)]
            self._buckets.append(ceil(p / (buckets**(self._n - 1))))
        else:
            self._n = len(buckets)
            self._buckets = copy(buckets)

        self._p = p

        # If Sboxes is a list of integers, we let it be the S-box for all buckets.
        if Sboxes[0] in ZZ:
            self._Sboxes = [SplitAndLookup._fill_sbox(Sboxes, self._buckets[i]) for i in range(self._n)]
        else:
            self._Sboxes = [SplitAndLookup._fill_sbox(Sboxes[i], self._buckets[i]) for i in range(self._n)]

    def __str__(self):
        """
        Print the prime modulus and the bucket sizes.
        """
        return "SplitAndLookup object defined modulo {} decomposed into buckets {}".format(self.get_p(), tuple(self.get_buckets()))

    ## Getters.

    def get_p(self):
        """
        Get the modulus p over which the S-box is defined.
        """
        return self._p

    def get_n(self):
        """
        Get the number n of small S-boxes in parallel.
        """
        return self._n

    def get_small_S(self, k):
        """
        Get the k-th small S-box, 
        that is the one applied to the k-th least significant bucket.
        """
        return self._Sboxes[k]

    def eval_small_S(self, k, x):
        """
        Evaluatie the k-th small S-box on the value x, 
        that is the one applied to the k-th least significant bucket.
        """
        return self._Sboxes[k][x]
    
    def get_buckets(self):
        """
        Get the size of the buckets.
        """
        return self._buckets

    def get_bucket(self, k):
        """
        Get the size of the k-th least significant bucket.
        """
        return self._buckets[k]

    ## Evaluation of the split-and-lookup.

    def _decomp(self, x):
        """
        Decompose x using the buckets.
        """
        x_list = []

        for si in self.get_buckets():
            x_list.append(x % si)
            x //= si

        return x_list


    def _recomp(self, x_list):
        """
        Recompose x using the buckets.
        """
        n = self.get_n()
        x = 0

        for i in range(n):
            x *= self.get_bucket(n - 1 - i)
            x += x_list[n - 1 - i]

        return x

    def _apply_Sboxes(self, x_list):
        """
        Apply the S-boxes in parallel.
        """
        n = self.get_n()
        y_list = []

        for i in range(n):
            xi = x_list[i]
            yi = self.eval_small_S(i, xi)

            y_list.append(yi)

        return y_list

    def eval(self, x):
        """
        Compute the value S(x).
        """
        p = self.get_p()

        x_red = x % p

        x_list = self._decomp(x_red)
        y_list = self._apply_Sboxes(x_list)
        y = self._recomp(y_list)

        y_red = y % p
        return y_red

    ## Check if the function is a permutation.

    def _small_permutations(self):
        """
        Check if all small S-boxes are permutations.
        """
        n = self.get_n()
        small_perm = True

        for i in range(n):
            si = self.get_bucket(i)
            Si = self.get_small_S(i)

            small_perm &= (Counter(Si) == Counter(range(si)))
        return small_perm
    
    def _partial_ordering(self, epsilon):
        """
        Check if all S-boxes are partially ordered according to epsilon, that is:
        - if p - epsilon = Sum_(k<n) Prod_(i<k) si * vk .
        - for each k:
            S(xk) < vk when xk < vk .
            S(vk) = vk .
            S(xk) > vk when xk > vk .
        """
        partially_ordered = True

        v_list = self._decomp(self.get_p() - epsilon)

        for i in range(self.get_n()):
            vi = v_list[i]
            Si = self.get_small_S(i)
            si = self.get_bucket(i)

            for xi in range(vi):
                partially_ordered &= (Si[xi] < vi)

            partially_ordered &= (Si[vi] == vi)

            for xi in range(vi + 1, si):
                partially_ordered &= (Si[xi] > vi)

        return partially_ordered


    def is_permutation(self):
        """
        Check if the split-and-lookup is a permutation using the "partially-ordered" criteria. 
        """
        small_perm = self._small_permutations()
        partially_ordered = self._partial_ordering(0) or self._partial_ordering(1)

        return (small_perm and partially_ordered)

    ## Compute the inverse of the function.

    def _inverse_small(self, k):
        """
        Compute the inverse of a small permutation.
        """
        sk = self.get_bucket(k)
        Sk = self.get_small_S(k)
        
        Sk_inv = [0 for x in range(sk)]

        for x in range(sk):
            Sk_inv[Sk[x]] = x

        return Sk_inv

    def inverse(self):
        """
        Compute the inverse of the split-and-lookup
        (returns an error if it is not a permutation).
        """
        if not(self.is_permutation()):
            print("[Error] this SplitAndLookup is not a permutation.")
            return None
        Sboxes_inv = [self._inverse_small(k) for k in range(self._n)]
        return SplitAndLookup(Sboxes_inv, self.get_buckets(), self.get_p())


    ## Linear approximation.

    def _gen_lis(self, i, a, b, vals):
        """
        Generate a table:
        * keys are the values achieved by the map a*Sk(x)+b*x where x in vals
        * values are lists of x's achieving this value
        """
        lis = defaultdict(int)
        Si = self.get_small_S(i)
        
        for x in vals:
            y = a*Si[x] + b*x
            lis[y] += 1
            
        return lis

    def _const_to_mat(self, lis, a, b, c, k):
        """
        Compute the carry matrix whose coefficients are
            m_ij = #{x | a*Sk(x)+b*x = c+i*s - j} .
        """
        s = self.get_bucket(k)
        low = -(max(0, -a) + max(0,  -b))
        up  =   max(0,  a) + max(0, b) -1
        
        mat = matrix(ZZ, [[lis[c-i+j*s] for i in range(low, up+1)]
                        for j in range(low, up+1)])
        
        return mat

    def _list_to_mat(self, lis, a, b, c_list):
        """
        Compute the carry matrix over a list of small S-boxes, 
        according to a list of constants.
        """
        mat = self._const_to_mat(lis, a, b, c_list[0], 0)
        for k in range(1, self.get_n()):
            mat = mat * self._const_to_mat(lis, a, b, c_list[k], k)
        return mat

    def _mat(self, a, b, c):
        """
        Generate the carry matrix of S.
        """
        p = self.get_p()
        n = self.get_n()

        c_list = self._decomp(c)
        v_list = self._decomp(p-1)

        # Contains the lists for [sj], cj .
        lis_list_full_bucket = []
        # Contains the lists for [vj], cj .
        lis_list_less_v = []
        # Contains the lists for {vj}, cj .
        lis_list_only_v = []
        
        for k in range(n):
            lis_list_full_bucket.append(self._gen_lis(k, a, b, range(self.get_bucket(k))))
            lis_list_less_v.append(self._gen_lis(k, a, b, range(v_list[k])))
            lis_list_only_v.append(self._gen_lis(k, a, b, [v_list[k]]))
        
        # Contains the M^(Sn,{vn})_(a,b,cn) * ... * M^(Sj,{vj})_(a,b,cj) .
        mat_list_only_v_prod = [self._const_to_mat(lis_list_only_v[n-1], a, b, c_list[n-1], n-1)]
        # Contains the M^(Sj,[vj])_(a,b,cj) .
        mat_list_less_v = []
        # Contains the M^(Sj,[sj])_(a,b,cj) * ... * M^(S1,[s1])_(a,b,c1) .
        mat_list_full_bucket_prod = [self._const_to_mat(lis_list_full_bucket[0], a, b, c_list[0], 0)]
        
        for k in range(n):
            mat_list_less_v.append(self._const_to_mat(lis_list_less_v[k], a, b, c_list[k], k))
        
        for k in range(1, n):
            mat_list_only_v_prod.append(mat_list_only_v_prod[-1]*self._const_to_mat(lis_list_only_v[n-k-1], a, b, c_list[n-k-1], n-k-1))
            mat_list_full_bucket_prod.append(self._const_to_mat(lis_list_full_bucket[k], a, b, c_list[k], k)*mat_list_full_bucket_prod[-1])
        
        mat_list_only_v_prod.reverse()
        
        # Sum over j the
        #     M^(Sn,{vn})_(a,b,cn) * ... * M^(S(j+1),{v(j+1)})_(a,b,c(j+1)) * M^(Sj,[vj])_(a,b,cj)
        #   * M^(S(j-1),[s(j-1)])_(a,b,c(j-1)) * ... * M^(S1,[s1])_(a,b,c1) .
        
        mat = mat_list_less_v[n-1]*mat_list_full_bucket_prod[n-2]
        
        for k in range(n-2):
            mat += mat_list_only_v_prod[k+2]*mat_list_less_v[k+1]*mat_list_full_bucket_prod[k]
            
        mat += mat_list_only_v_prod[1]*mat_list_less_v[0]
        mat += mat_list_only_v_prod[0]
        
        return mat

    def _get_small_ab(self, a, b, c):
        """
        Find lambda such that |a'|+|b'| is small where:
            a' = lambda*a mod p, b' = lambda*b mod p .
        """
        p = self.get_p()

        vector_mod_matrix = matrix(ZZ, [[a, b], [p, 0], [0, p]])
        vector_mod_lattice = IntegerLattice(vector_mod_matrix)

        small_ab = vector_mod_lattice.shortest_vector()

        a_prime, b_prime = small_ab[0], small_ab[1]
        if a_prime < 0:
            a_prime *= -1
            b_prime *= -1

        c_prime = (c*a_prime / a) % p

        return a_prime, b_prime, c_prime

    def _linear_approx(self, a, b, c):
        """
        Return the number of solutions of the equation:
            a*S(x) + b*x = c .
        """
        p = self.get_p()
        c %= p

        low = -(max(0, -a) + max(0,  -b))
        up  =   max(0,  a) + max(0, b) -1

        bucket_prod = prod(self.get_buckets())
        
        score = 0
        
        for k in range(low, up+1):
            c_prime  = (c + k*p)
            k_prime  = c_prime // bucket_prod
            c_prime %= bucket_prod
            
            mat_k = self._mat(a, b, c_prime)

            score += mat_k[k_prime -low, -low]
            
        return score
    
    def linear_approximation(self, a, b, c):
        """
        Normalize a, b, c to get a, b to be small integers and
        return the number of solutions of the equation:
            a*S(x) + b*x = c .
        """
        a_prime, b_prime, c_prime = self._get_small_ab(a, b, c)
        return self._linear_approx(a_prime, b_prime, c_prime)

    ## Upper bounds for the linear approximations.

    def _lower_bound_ab_k(self, a, b, k):
        """
        Given a, b compute the maximal coefficient of
            M^Sk (a, b, ck) ,
        over 0 <= ck < sk.
        """
        sk = self.get_bucket(k)
        lis = self._gen_lis(k, a, b, range(sk))

        return max(lis.values())

    def _lower_bound_ab(self, a, b):
        """
        Given a, b provides a lower bound for the number of solutions to
            a*S(x) + b*x = c .
        
        This lower bound is computed using the trivial coordinate-wise bound.
        """
        n = self.get_n()
        bound = 1

        for k in range(n):
            maxi_k = 0
            maxi_k = self._lower_bound_ab_k(a, b, k)
            bound *= maxi_k

        return bound

    def lower_bound(self, max_ab):
        """
        Provide a lower bound for the number of solutions to
            a*S(x) + b*x = c , |a|, |b| < max_ab.
        """
        maxi_bound = 0

        for a in range(1, max_ab):
            for b in range(max_ab):
                norm_ab  = self._lower_bound_ab(a,  b)
                norm_amb = self._lower_bound_ab(a, -b)

                maxi_bound = max([maxi_bound, norm_ab, norm_amb])
        
        return maxi_bound

    ## Linear correlation.

    def _exp_summation(self, i, a, b, vals):
        """
        Sum over a single small S-box:
            Sum_(x in vals) e^(2*pi*i (a*Si(xi) + b*xi)*Prod_(j<i) sj /p) .
        """
        p = self.get_p()
        exp_sum = 0
        Si = self.get_small_S(i)
        # Compute Prod_(k<i) sk .
        s_cur = prod(self.get_buckets()[:i])
        
        zeta_p = CC.zeta(p)
        
        for x in vals:
            exp_sum += zeta_p**((a*Si[x]+b*x)*s_cur % p)
        return exp_sum


    def linear_correlation(self, a, b):
        """
        Compute the linear correlation of the function:
            Sum_x e^(2*pi*i (a*S(x) + b*x)/p) .
        """
        p = self.get_p()
        n = self.get_n()
        v_list = self._decomp(p-1)
        exp_sum_p = 0
        
        # Contains the Sum_(x<sj) e^(2*pi*i (a*Sj(x) + b*x)*Prod_(k<j) sk/p) .
        exp_sum_list_full_bucket = []
        # Contains the Sum_(x<vj) e^(2*pi*i (a*Sj(x) + b*x)*Prod_(k<j) sk/p) .
        exp_sum_list_less_v = []
        # Contains the          e^(2*pi*i (a*Sj(vj) + b*vj)*Prod_(k<j) sk/p) .
        exp_sum_list_only_v = []
        
        for k in range(n):
            exp_sum_list_full_bucket.append(self._exp_summation(k, a, b, range(self.get_bucket(k))))
            exp_sum_list_less_v.append(self._exp_summation(k, a, b, range(v_list[k])))
            exp_sum_list_only_v.append(self._exp_summation(k, a, b, [v_list[k]]))
        
        # Contains the products Prod_(k<=j) exp_sum_list_full_bucket[k] .
        exp_sum_list_full_bucket_prod = [exp_sum_list_full_bucket[0]]
        # Contains the products Prod_(k>=j) exp_sum_list_only_v_prod[k] .
        exp_sum_list_only_v_prod = [exp_sum_list_only_v[n-1]]
        
        for k in range(1,n):
            exp_sum_list_only_v_prod.append(exp_sum_list_only_v[n-k-1]*exp_sum_list_only_v_prod[-1])
            exp_sum_list_full_bucket_prod.append(exp_sum_list_full_bucket_prod[-1]*exp_sum_list_full_bucket[k])
        
        exp_sum_list_only_v_prod.reverse()
        
        exp_sum_p += exp_sum_list_only_v_prod[0]
        exp_sum_p += exp_sum_list_only_v_prod[1]*exp_sum_list_less_v[0]
        
        for k in range(n-2):
            exp_sum_p += exp_sum_list_only_v_prod[k+2]*exp_sum_list_less_v[k+1]*exp_sum_list_full_bucket_prod[k]
            
        exp_sum_p += exp_sum_list_less_v[n-1]*exp_sum_list_full_bucket_prod[n-2]
        
        return exp_sum_p


## Tests on the S-boxes of Tip5, Monolith, Reinforced Concrete.

# Fixed parameters.

p_goldi  = 2**64 - 2**32 + 1
p_mers   = 2**31 - 1
p_BLS381 = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
p_BN254  = 0x30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000001
p_ST     =  0x3fa000000000000000000000000000000000000000000000000000000000001

bytesize = 256


# Case of Monolith.

def circular_shift(x, k, n):
    """
    Shift circularly the bits of an n-bit word of k positions towards the most significant bit.
    """
    return ((x  << k) | (x >> (n-k))) & (2**n-1)

def bitwise_not(x, n):
    """
    Compute the bitwise not of an n-bit word.
    """
    return x ^^ (2**n -1)

def chi8(x):
    """
    Compute the 8-bit Daemen chi function.
    """
    return x ^^ (circular_shift(bitwise_not(x, 8), 1, 8) & circular_shift(x, 2, 8) & circular_shift(x, 3, 8))

def chi7(x):
    """
    Compute the 7-bit Daemen chi function.
    """
    return x ^^ (circular_shift(bitwise_not(x, 7), 1, 7) & circular_shift(x, 2, 7))

def test_monolith():
    """
    Tests on the 64-bits and 31-bits Monolith S-boxes.
    """
    small_S_monolith8 = [circular_shift(chi8(x), 1, 8) for x in range(bytesize)]
    small_S_monolith7 = [circular_shift(chi7(x), 1, 7) for x in range(128)]

    S_monolith_64 = SplitAndLookup(small_S_monolith8, bytesize, p_goldi)
    S_monolith_31 = SplitAndLookup([small_S_monolith8, small_S_monolith8, small_S_monolith8, small_S_monolith7], bytesize, p_mers)

    print("Case of the 64-bit permutation for Monolith:")
    print("Is a permutation: {}".format(S_monolith_64.is_permutation()))
    print("Has {} fixed points.".format(S_monolith_64.linear_approximation(1, -1, 0)))
    print("Best linear approximation found is S(x) = 2*x: {} solutions.".format(S_monolith_64.linear_approximation(1, -2, 0)))
    print("Best linear correlation found is for S(x) - 2*x: {}.".format(abs(S_monolith_64.linear_correlation(1, -2)/p_goldi)))


    print("Case of the 31-bit permutation for Monolith:")
    print("Is a permutation: {}".format(S_monolith_31.is_permutation()))
    print("Has {} fixed points.".format(S_monolith_31.linear_approximation(1, -1, 0)))
    print("Best linear approximation found is S(x) = 2*x + 2^24: {} solutions.".format(S_monolith_31.linear_approximation(1, -2, 2**24)))
    print("Best linear correlation found is for 128*S(x) - 256*x: {}.".format(abs(S_monolith_31.linear_correlation(128, -256)/p_mers)))


# Case of Tip5.

def test_tip_five():
    """
    Tests on the S-box of Tip5.
    """
    small_S_tip5 = [((x+1)**3 -1) % 257 for x in range(bytesize)]

    S_tip5 = SplitAndLookup(small_S_tip5, bytesize, p_goldi)

    print("Case of Tip5:")
    print("Is a permutation: {}".format(S_tip5.is_permutation()))
    print("Has {} fixed points.".format(S_tip5.linear_approximation(1, -1, 0)))
    print("Best linear approximation found is S(x) = -x + 0x086a0896f76a0894: {} solutions.".format(S_tip5.linear_approximation(1, 1, 0x086a0896f76a0894)))
    print("Best linear correlation found is for 22*S(x) - 14*x: {}.".format(abs(S_tip5.linear_correlation(22, -14)/p_goldi)))


# Case of Reinforced Concrete.

s_BLS381 = [693, 696, 694, 668, 679, 695, 691, 693, 700, 
            688, 700, 694, 701, 694, 699, 701, 701, 701, 
            695, 698, 697, 703, 702, 691, 688, 703, 679]

small_S_RC_BLS381 = [
        171, 178, 483, 527, 653, 408, 197, 599, 300, 607, 403, 511, 579, 520, 591, 412, 261, 559,
        551, 154, 180, 138, 596, 150, 276, 271, 48, 168, 362, 637, 467, 164, 536, 554, 287, 530,
        431, 92, 654, 518, 323, 572, 624, 4, 258, 439, 430, 495, 534, 222, 545, 31, 44, 18, 80, 55,
        399, 328, 505, 313, 441, 586, 501, 598, 566, 568, 77, 496, 106, 563, 537, 78, 50, 450, 445,
        166, 237, 617, 185, 404, 621, 578, 133, 517, 646, 98, 86, 492, 267, 193, 33, 476, 207, 17,
        487, 643, 52, 384, 74, 148, 121, 657, 633, 528, 269, 611, 567, 601, 391, 231, 226, 658,
        331, 191, 354, 23, 474, 277, 390, 341, 279, 442, 422, 638, 15, 196, 329, 377, 36, 433, 398,
        72, 256, 352, 253, 550, 635, 142, 343, 176, 500, 588, 413, 569, 266, 42, 283, 535, 410,
        538, 647, 85, 27, 423, 558, 61, 356, 348, 43, 19, 625, 291, 238, 274, 432, 448, 100, 642,
        260, 587, 622, 608, 366, 420, 477, 316, 605, 254, 130, 407, 471, 174, 631, 34, 652, 628,
        175, 134, 122, 192, 531, 217, 32, 257, 145, 307, 262, 83, 509, 440, 600, 589, 359, 522,
        268, 143, 498, 512, 333, 651, 151, 183, 126, 351, 39, 246, 242, 630, 543, 574, 610, 655,
        25, 494, 456, 612, 123, 315, 340, 296, 580, 503, 281, 428, 62, 10, 76, 203, 288, 91, 426,
        128, 629, 29, 218, 292, 447, 161, 117, 388, 540, 364, 245, 541, 224, 502, 370, 229, 90,
        466, 636, 208, 51, 562, 259, 344, 334, 111, 235, 488, 632, 577, 54, 386, 75, 181, 463, 421,
        24, 96, 406, 156, 158, 265, 5, 310, 37, 124, 88, 155, 480, 593, 202, 451, 1, 497, 645, 457,
        187, 56, 206, 179, 640, 249, 99, 240, 460, 490, 163, 369, 293, 186, 553, 46, 449, 41, 219,
        308, 7, 234, 336, 373, 372, 347, 215, 481, 542, 146, 357, 656, 136, 330, 595, 516, 592,
        273, 365, 8, 47, 641, 81, 484, 573, 614, 437, 533, 0, 282, 184, 400, 49, 114, 374, 280,
        499, 418, 139, 382, 613, 233, 345, 393, 575, 508, 299, 101, 582, 360, 285, 2, 376, 548,
        189, 648, 214, 618, 385, 371, 425, 552, 204, 286, 443, 210, 294, 211, 241, 461, 275, 165,
        350, 59, 583, 159, 434, 252, 71, 436, 529, 236, 475, 339, 367, 147, 170, 110, 22, 298, 506,
        172, 247, 513, 73, 230, 314, 239, 157, 116, 65, 11, 570, 40, 620, 205, 251, 594, 468, 69,
        489, 109, 452, 465, 312, 383, 129, 379, 335, 353, 602, 546, 243, 57, 473, 486, 320, 162,
        526, 115, 26, 560, 107, 458, 519, 169, 97, 358, 504, 414, 13, 459, 132, 167, 402, 14, 491,
        571, 105, 112, 363, 581, 194, 84, 349, 201, 462, 289, 53, 603, 209, 396, 303, 317, 102, 82,
        131, 639, 3, 435, 378, 415, 539, 223, 30, 510, 199, 479, 397, 45, 248, 561, 67, 213, 438,
        20, 405, 557, 120, 89, 584, 555, 264, 419, 525, 429, 392, 311, 68, 446, 270, 585, 113, 627,
        472, 38, 375, 327, 127, 417, 547, 12, 108, 368, 95, 250, 322, 198, 380, 149, 104, 87, 332,
        135, 28, 318, 482, 221, 188, 58, 544, 521, 93, 324, 64, 272, 297, 644, 453, 225, 606, 295,
        216, 152, 411, 361, 444, 469, 427, 507, 395, 609, 153, 381, 464, 424, 94, 9, 564, 321, 615,
        21, 227, 137, 70, 326, 549, 556, 565, 416, 470, 255, 60, 604, 590, 305, 35, 278, 6, 125,
        387, 220, 597, 63, 454, 401, 119, 302, 309, 342, 16, 619, 493, 290, 616, 173, 304, 195,
        524, 263, 212, 649, 626, 409, 338, 306, 389, 79, 160, 66, 177, 232, 478, 514, 650, 455,
        103, 144, 355, 182, 346, 284, 200, 634, 244, 140, 337, 325, 319, 532, 394, 118, 485, 301,
        623, 190, 523, 515, 576, 141, 228
    ]

s_BN254  = [651, 658, 656, 666, 663, 654, 668, 677, 681, 
            683, 669, 681, 680, 677, 675, 668, 675, 683, 
            681, 683, 683, 655, 680, 683, 667, 678, 673]

small_S_RC_BN254 = [
        377, 222, 243, 537, 518, 373, 152, 435, 526, 352, 2, 410, 513, 545, 567, 354, 405, 80, 233,
        261, 49, 240, 568, 74, 131, 349, 146, 278, 330, 372, 43, 432, 247, 583, 105, 203, 637, 307,
        29, 597, 633, 198, 519, 95, 148, 62, 68, 312, 616, 357, 234, 433, 154, 90, 163, 249, 101,
        573, 447, 587, 494, 103, 608, 394, 409, 73, 317, 305, 346, 562, 262, 313, 303, 550, 64,
        102, 259, 400, 495, 572, 238, 40, 612, 236, 586, 15, 361, 386, 138, 136, 107, 33, 190, 423,
        176, 161, 460, 35, 202, 589, 32, 160, 444, 517, 490, 515, 144, 195, 269, 332, 25, 308, 192,
        276, 623, 180, 626, 217, 329, 66, 392, 431, 12, 478, 67, 232, 258, 355, 94, 191, 632, 181,
        298, 1, 301, 79, 618, 523, 627, 484, 306, 610, 635, 619, 544, 420, 408, 158, 328, 61, 406,
        299, 442, 178, 625, 621, 497, 465, 574, 143, 54, 57, 89, 322, 135, 96, 605, 599, 473, 97,
        85, 133, 200, 93, 291, 525, 529, 206, 614, 319, 196, 482, 17, 168, 70, 104, 441, 159, 364,
        603, 78, 150, 230, 116, 31, 630, 132, 69, 499, 532, 218, 492, 112, 505, 437, 333, 457, 456,
        439, 639, 398, 16, 436, 264, 450, 211, 241, 524, 294, 235, 126, 165, 527, 452, 212, 157,
        272, 208, 469, 611, 338, 83, 326, 151, 139, 607, 285, 585, 58, 14, 193, 71, 440, 511, 542,
        390, 470, 155, 413, 606, 142, 367, 371, 174, 5, 60, 289, 297, 336, 370, 76, 209, 622, 453,
        257, 555, 44, 430, 345, 335, 548, 459, 47, 426, 591, 559, 417, 284, 552, 137, 277, 281,
        463, 631, 350, 265, 323, 108, 290, 169, 634, 609, 414, 130, 6, 166, 316, 207, 592, 280,
        391, 274, 20, 300, 593, 549, 3, 602, 418, 472, 419, 296, 41, 46, 615, 638, 388, 553, 282,
        356, 327, 462, 115, 325, 121, 399, 273, 334, 383, 488, 292, 55, 628, 9, 19, 601, 496, 228,
        201, 576, 374, 558, 153, 162, 341, 353, 84, 220, 461, 221, 547, 344, 507, 577, 140, 485,
        471, 11, 175, 13, 53, 543, 270, 120, 30, 584, 384, 368, 397, 239, 4, 483, 620, 189, 522,
        540, 510, 149, 245, 533, 283, 256, 369, 302, 571, 128, 253, 448, 446, 183, 99, 438, 468,
        42, 594, 487, 403, 23, 172, 340, 106, 481, 251, 363, 295, 489, 474, 337, 87, 86, 246, 215,
        376, 315, 415, 117, 286, 600, 56, 145, 91, 358, 429, 411, 516, 310, 213, 598, 10, 395, 111,
        506, 237, 170, 512, 82, 147, 579, 402, 501, 343, 38, 434, 214, 314, 360, 77, 565, 320, 385,
        404, 199, 331, 351, 466, 596, 365, 231, 477, 604, 254, 268, 539, 424, 167, 378, 491, 535,
        141, 267, 177, 27, 546, 219, 556, 216, 451, 387, 28, 50, 569, 255, 288, 156, 449, 379, 508,
        528, 531, 624, 581, 554, 59, 171, 252, 0, 595, 185, 51, 520, 575, 475, 113, 187, 194, 428,
        500, 617, 188, 321, 179, 263, 110, 467, 18, 401, 22, 164, 342, 21, 382, 381, 127, 52, 570,
        45, 445, 36, 534, 339, 98, 293, 244, 266, 629, 229, 122, 123, 48, 88, 225, 173, 100, 114,
        536, 636, 205, 34, 425, 502, 514, 304, 613, 530, 118, 75, 561, 582, 81, 480, 92, 498, 464,
        224, 479, 563, 223, 640, 521, 427, 503, 250, 375, 186, 72, 242, 125, 380, 271, 204, 407,
        366, 197, 119, 7, 493, 26, 109, 65, 359, 396, 311, 309, 458, 134, 393, 557, 476, 324, 421,
        275, 37, 39, 580, 184, 560, 8, 455, 509, 422, 24, 287, 590, 182, 416, 318, 260, 578, 454,
        389, 129, 566, 63, 486, 541, 362, 210, 551, 348, 279, 538, 347, 504, 124, 564, 443, 412,
        226, 227, 248, 588
    ]

s_ST     = [1023] + [1024 for i in range(24)]

small_S_RC_ST = [
        849, 68, 27, 909, 988, 687, 828, 507, 847, 380, 656, 379, 340, 296, 974, 3, 338, 355, 263,
        968, 754, 119, 442, 231, 629, 634, 938, 484, 73, 954, 704, 20, 1006, 447, 977, 591, 528,
        593, 103, 69, 236, 45, 843, 461, 762, 158, 908, 661, 751, 874, 545, 96, 35, 802, 738, 495,
        597, 560, 956, 518, 262, 991, 54, 156, 821, 646, 620, 581, 454, 470, 753, 617, 550, 91,
        647, 481, 475, 992, 287, 141, 523, 7, 233, 51, 89, 614, 336, 126, 857, 882, 194, 806, 55,
        793, 443, 584, 213, 967, 110, 673, 645, 979, 446, 116, 621, 795, 760, 1, 473, 543, 185,
        1008, 399, 105, 344, 205, 914, 830, 851, 927, 393, 290, 716, 17, 906, 170, 918, 895, 638,
        26, 327, 409, 161, 371, 559, 363, 513, 67, 61, 121, 549, 886, 62, 822, 925, 747, 357, 618,
        201, 624, 464, 665, 892, 317, 302, 943, 235, 64, 642, 416, 104, 799, 521, 839, 875, 220,
        623, 921, 361, 522, 234, 625, 562, 128, 800, 117, 275, 81, 90, 313, 834, 176, 554, 82, 168,
        928, 504, 637, 764, 721, 532, 193, 100, 911, 434, 890, 28, 160, 565, 541, 397, 901, 996,
        23, 922, 146, 301, 844, 303, 697, 107, 136, 768, 869, 494, 347, 428, 798, 949, 957, 735,
        929, 0, 294, 619, 677, 299, 548, 8, 43, 284, 202, 232, 260, 109, 745, 982, 976, 78, 695,
        845, 790, 826, 375, 755, 524, 823, 450, 512, 1000, 540, 948, 856, 568, 269, 712, 771, 873,
        816, 999, 195, 811, 199, 708, 348, 539, 765, 133, 774, 500, 49, 492, 24, 162, 891, 211,
        258, 582, 729, 346, 18, 788, 360, 217, 820, 448, 249, 1010, 405, 316, 430, 228, 410, 803,
        692, 4, 852, 224, 777, 752, 22, 950, 455, 883, 97, 557, 488, 221, 585, 124, 879, 342, 458,
        981, 670, 827, 387, 219, 120, 858, 930, 414, 932, 411, 207, 558, 123, 696, 47, 369, 920,
        813, 351, 503, 809, 343, 268, 664, 505, 118, 70, 970, 30, 324, 325, 863, 570, 987, 789, 76,
        936, 903, 190, 218, 401, 706, 2, 276, 514, 632, 247, 705, 805, 586, 794, 993, 32, 34, 15,
        502, 84, 672, 214, 733, 984, 417, 724, 72, 866, 66, 842, 685, 717, 297, 469, 668, 636, 192,
        12, 145, 1003, 627, 700, 756, 281, 635, 385, 783, 893, 298, 11, 251, 131, 819, 931, 31,
        641, 285, 429, 178, 19, 868, 186, 792, 530, 689, 106, 366, 730, 169, 739, 10, 538, 872,
        796, 786, 39, 833, 273, 563, 271, 200, 453, 283, 825, 462, 1007, 657, 727, 139, 419, 280,
        740, 720, 898, 889, 510, 832, 423, 383, 256, 942, 33, 841, 613, 319, 471, 48, 779, 406,
        198, 564, 924, 465, 770, 650, 535, 413, 330, 590, 933, 1001, 734, 651, 432, 534, 436, 486,
        876, 111, 596, 345, 531, 177, 41, 95, 245, 552, 606, 653, 743, 667, 837, 767, 138, 744,
        203, 659, 307, 648, 723, 726, 569, 997, 980, 44, 154, 227, 797, 499, 881, 153, 609, 382,
        511, 812, 763, 439, 216, 125, 323, 566, 900, 517, 818, 305, 814, 293, 400, 728, 829, 815,
        241, 854, 592, 304, 913, 261, 407, 370, 533, 703, 403, 372, 761, 229, 75, 587, 150, 669,
        575, 252, 985, 463, 164, 438, 542, 181, 516, 526, 288, 306, 680, 595, 248, 556, 425, 959,
        426, 386, 519, 311, 561, 38, 995, 690, 746, 191, 182, 758, 679, 962, 896, 772, 130, 877,
        870, 951, 536, 701, 732, 134, 553, 114, 490, 989, 174, 242, 860, 940, 941, 939, 574, 736,
        601, 420, 389, 71, 888, 183, 359, 986, 332, 681, 451, 567, 412, 267, 188, 675, 58, 576,
        850, 817, 722, 894, 77, 152, 395, 112, 944, 750, 551, 277, 135, 254, 240, 459, 189, 196,
        684, 259, 643, 87, 333, 132, 741, 749, 74, 376, 115, 907, 599, 364, 92, 171, 6, 508, 244,
        600, 79, 952, 824, 209, 477, 631, 958, 710, 501, 80, 694, 660, 1011, 14, 289, 52, 537, 808,
        339, 801, 836, 445, 766, 963, 384, 546, 250, 615, 965, 785, 487, 686, 640, 468, 264, 698,
        525, 910, 945, 489, 926, 418, 579, 961, 731, 328, 835, 611, 318, 13, 616, 255, 855, 865,
        341, 444, 312, 737, 257, 491, 529, 639, 279, 326, 662, 682, 496, 804, 482, 983, 583, 960,
        775, 572, 424, 916, 5, 603, 711, 757, 335, 912, 396, 278, 295, 791, 848, 947, 781, 172,
        972, 282, 50, 699, 971, 719, 184, 633, 449, 966, 25, 147, 381, 368, 197, 274, 742, 655,
        649, 955, 748, 1004, 605, 588, 905, 969, 457, 493, 978, 309, 666, 479, 885, 21, 140, 973,
        917, 644, 923, 466, 871, 678, 476, 769, 626, 472, 573, 934, 437, 300, 878, 374, 98, 688,
        149, 485, 838, 143, 594, 57, 1005, 1012, 919, 179, 456, 377, 773, 349, 85, 880, 718, 243,
        433, 166, 520, 93, 408, 807, 715, 83, 390, 334, 467, 864, 398, 337, 155, 602, 352, 388,
        173, 440, 208, 391, 422, 709, 612, 331, 266, 497, 707, 270, 431, 460, 478, 598, 810, 483,
        392, 498, 350, 129, 365, 713, 953, 246, 610, 904, 787, 663, 435, 29, 226, 127, 902, 452,
        759, 652, 215, 506, 782, 362, 676, 37, 674, 265, 9, 831, 163, 862, 180, 671, 358, 167, 42,
        238, 60, 780, 272, 59, 315, 314, 946, 658, 94, 367, 322, 884, 113, 175, 210, 122, 212, 607,
        555, 310, 994, 329, 778, 204, 63, 714, 1002, 101, 40, 1009, 630, 474, 142, 225, 230, 165,
        292, 404, 702, 846, 693, 421, 157, 102, 108, 137, 515, 86, 222, 580, 402, 577, 36, 353,
        148, 46, 320, 975, 887, 776, 990, 628, 683, 151, 56, 544, 937, 480, 604, 571, 159, 622,
        356, 861, 53, 99, 608, 589, 65, 784, 691, 239, 867, 853, 441, 859, 998, 237, 527, 899, 725,
        394, 373, 16, 223, 253, 354, 509, 378, 578, 187, 291, 308, 415, 964, 427, 915, 547, 144,
        897, 935, 88, 840, 286, 206, 321, 654
    ]

def test_RC_BLS():
    """
    Tests on the S-box of the BLS381 version of Reinforced Concrete.
    """
    S_RC_BLS381 = SplitAndLookup(small_S_RC_BLS381, s_BLS381, p_BLS381)

    print("Case of Reinforced Concrete-BLS381:")
    print("Is a permutation: {}".format(S_RC_BLS381.is_permutation()))
    print("Has {} fixed points.".format(S_RC_BLS381.linear_approximation(1, -1, 0)))
    print("Best linear correlation found is for 4*S(x) - 6*x: {}.".format(abs(S_RC_BLS381.linear_correlation(4, -6)/p_BLS381)))


def test_RC_BN():
    """
    Tests on the S-box of the BN254 version of Reinforced Concrete.
    """
    S_RC_BN254  = SplitAndLookup(small_S_RC_BN254, s_BN254, p_BN254)

    print("Case of Reinforced Concrete-BN254:")
    print("Is a permutation: {}".format(S_RC_BN254.is_permutation()))
    print("Has {} fixed points.".format(S_RC_BN254.linear_approximation(1, -1, 0)))
    print("Best linear correlation found is for 5*S(x) - 7*x: {}.".format(abs(S_RC_BN254.linear_correlation(5, -7)/p_BN254)))

def test_RC_ST():
    """
    Tests on the S-box of the ST version of Reinforced Concrete.
    """
    S_RC_ST     = SplitAndLookup(small_S_RC_ST, s_ST, p_ST)

    print("Case of Reinforced Concrete-ST:")
    print("Is a permutation: {}".format(S_RC_ST.is_permutation()))
    print("Has {} fixed points.".format(S_RC_ST.linear_approximation(1, -1, 0)))
    print("Best linear correlation found is for 4*S(x) + x: {}.".format(abs(S_RC_ST.linear_correlation(4, 1)/p_ST)))


## Our bit-shuffle proposals for constant-time split-and-lookups.

def int_to_vec(x, n):
    """
    Convert a n-bit integer into its bit representation.
    """
    x_vec = []
    for i in range(n):
        x_vec.append(x % 2)
        x //= 2
    return x_vec

def vec_to_int(x_vec):
    """
    Convert a bit vector into an integer.
    """
    x = 0
    n = len(x_vec)
    
    for i in range(n):
        x += 2**i * x_vec[i]
    return x

def shuffle(x, sigma):
    """
    Shuffle the bits of x using the sigma permutation.
    """
    n = len(sigma)
    
    x_vec = int_to_vec(x, n)
    y_vec = [x_vec[sigma_i] for sigma_i in sigma]
    
    y = vec_to_int(y_vec)
    return y


def test_shuffle():
    """
    Test the proposed 64-bit and 31-bit shuffle S-boxes.
    """
    sigma_8 = [6, 4, 3, 7, 0, 2, 5, 1]
    sigma_7 = [4, 0, 3, 1, 6, 2, 5]

    small_S_shuffle_8 = [shuffle(x, sigma_8) for x in range(256)]
    small_S_shuffle_7 = [shuffle(x, sigma_7) for x in range(128)]

    print("Value table for the 8-bit permutation: {}".format(small_S_shuffle_8))
    print("Value table for the 7-bit permutation: {}".format(small_S_shuffle_7))

    S_shuffle_64 = SplitAndLookup(small_S_shuffle_8, bytesize, p_goldi)
    S_shuffle_31 = SplitAndLookup([small_S_shuffle_8, small_S_shuffle_8, small_S_shuffle_8, small_S_shuffle_7], bytesize, p_mers)

    print("Case of the 64-bit shuffle S-box:")
    print("Is a permutation: {}".format(S_shuffle_64.is_permutation()))

    print("Case of the 31-bit shuffle S-box:")
    print("Is a permutation: {}".format(S_shuffle_31.is_permutation()))

    max_ab = 20

    bound_64 = S_shuffle_64.lower_bound(max_ab)
    bound_31 = S_shuffle_31.lower_bound(max_ab)

    print("Lower bounds for the linear approximations: number of solutions to a*S(x) + b*x = c, |a|, |b| < {}.".format(max_ab))
    print("Case of the 64-bit shuffle S-box: at least {} solutions.".format(bound_64))
    print("Case of the 31-bit shuffle S-box: at least {} solutions.".format(bound_31))

    a, b, c = 10, -17, 2399141888

    print("The carry propagation phenomena can increase those values:")

    print("For instance in the case a, b, c = {}, {}, {}.".format(a, b, c))
    print("we have {} solutions for the 64-bit S-box.".format(S_shuffle_64.linear_approximation(a, b, c)))

    c = 0

    print("For instance in the case a, b, c = {}, {}, {}.".format(a, b, c))
    print("we have {} solutions for the 31-bit S-box.".format(S_shuffle_31.linear_approximation(a, b, c)))


## Run the tests.
if __name__ == "__main__":
    test_monolith()
    test_tip_five()

    test_RC_BLS()
    test_RC_BN()
    test_RC_ST()

    test_shuffle()
