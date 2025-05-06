
p = 2**64-2**32+1
Fp = GF(p)

r = 2**8+1
Fr = GF(r)

## Define the S-boxes T and S.

def T(x):
    """
    Input:
        x:Fp
    Output:
        T(x):Fp
    """
    return x**7

def T_inv(x):
    """
    Inverse of T.
    """
    return x**10540996611094048183

R = Fp(2**64)

def sigma(x):
    """
    Input:
        x:Fp
    Output:
        x_list:Fr⁸, decomposition of x in base 256
    """
    b = 256
    x_int = ZZ(x)
    x_list = []
    for i in range(8):
        xi = Fr(x_int % b)
        x_list.append(xi)
        x_int //= 256
    return x_list

def L(xi):
    """
    Input:
        xi:Fr
    Output:
        L(xi):Fr
    """
    return (xi+1)**3-1

def L8(x_list):
    """
    Input:
        x_list:Fr⁸
    Output:
        L8(x_list):Fr⁸, L applied elementwise
    """
    x_L8 = [L(xi) for xi in x_list]
    return x_L8

def L_inv(xi):
    """
    Inverse of L.
    """
    return (xi+1)**171-1

def L8_inv(x_list):
    """
    Inverse of L8.
    """
    x_L8_inv = [L_inv(xi) for xi in x_list]
    return x_L8_inv

def rho(x_list):
    """
    Input:
        x_list:Fr^8
    Output:
        x:Fp, inverse operation of sigma
    """
    x = Fp(0)
    b = 256
    for i in range(8):
        x += Fp(x_list[i])*b**i
    return x
    

def S(x):
    """
    Input:
        x:Fp
    Output:
        S(x):Fp, output of the S-box S
    """
    y = 1/R * rho(L8(sigma(R*x)))
    return y

def S_inv(x):
    """
    Inverse of S.
    """
    y = 1/R * rho(L8_inv(sigma(R*x)))
    return y

## Linear part of Tip5.

Mcol =  vector(Fp, [
    61402, 1108, 28750, 33823, 7454, 43244, 53865, 12034, 56951, 27521, 41351, 40901, 12021, 59689, 26798, 17845
])
M = matrix.circulant(Mcol).transpose()

## Compute the f permutation

def f_round(X, round_constant):
    """
    Input:
        X:Fp^16, state of the permutation
        round_constant:Fp^16, constant added to the round
    Output:
        Y:Fp^16, state after one round
    """
    Y = zero_vector(Fp, 16)
    for i in range(4):
        Y[i] = S(X[i])
    for i in range(4,16):
        Y[i] = T(X[i])
    
    Y = M*Y + round_constant
    
    return Y

def f_round_inv(Y, round_constant):
    """
    Inverse of a round of f.
    """
    M_inv = M.inverse()
    X = M_inv*(Y - round_constant)
    
    for i in range(4):
        X[i] = S_inv(X[i])
    for i in range(4,16):
        X[i] = T_inv(X[i])
    
    return X


def f_inv(Y, round_constants):
    """
    Input:
        X:Fp^16, state of the permutation
        round_constant:Fp^16, constant added to the round
    Output:
        Y:Fp^16, state after one round
    """
    X = Y
    for round_constant in round_constants[::-1]:
        X = f_round_inv(X, round_constant)
    
    return X

def f(X, round_constants):
    """
    Input:
        X:Fp^16, state of the permutation
        round_constant:Fp^16, constant added to the round
    Output:
        Y:Fp^16, state after one round
    """
    Y = X
    for round_constant in round_constants:
        Y = f_round(Y, round_constant)
    
    return Y



round_constants = [
    vector(Fp, [
        13630775303355457758,
        16896927574093233874,
        10379449653650130495,
        1965408364413093495,
        15232538947090185111,
        15892634398091747074,
        3989134140024871768,
        2851411912127730865,
        8709136439293758776,
        3694858669662939734,
        12692440244315327141,
        10722316166358076749,
        12745429320441639448,
        17932424223723990421,
        7558102534867937463,
        15551047435855531404
            ]
        ),
    vector(Fp, [
        17532528648579384106,
        5216785850422679555,
        15418071332095031847,
        11921929762955146258,
        9738718993677019874,
        3464580399432997147,
        13408434769117164050,
        264428218649616431,
        4436247869008081381,
        4063129435850804221,
        2865073155741120117,
        5749834437609765994,
        6804196764189408435,
        17060469201292988508,
        9475383556737206708,
        12876344085611465020
            ]
        ),
    vector(Fp, [
        13835756199368269249,
        1648753455944344172,
        9836124473569258483,
        12867641597107932229,
        11254152636692960595,
        16550832737139861108,
        11861573970480733262,
        1256660473588673495,
        13879506000676455136,
        10564103842682358721,
        16142842524796397521,
        3287098591948630584,
        685911471061284805,
        5285298776918878023,
        18310953571768047354,
        3142266350630002035
            ]
        ),
    vector(Fp, [
        549990724933663297,
        4901984846118077401,
        11458643033696775769,
        8706785264119212710,
        12521758138015724072,
        11877914062416978196,
        11333318251134523752,
        3933899631278608623,
        16635128972021157924,
        10291337173108950450,
        4142107155024199350,
        16973934533787743537,
        11068111539125175221,
        17546769694830203606,
        5315217744825068993,
        4609594252909613081
            ]
        ),
    vector(Fp, [
        3350107164315270407,
        17715942834299349177,
        9600609149219873996,
        12894357635820003949,
        4597649658040514631,
        7735563950920491847,
        1663379455870887181,
        13889298103638829706,
        7375530351220884434,
        3502022433285269151,
        9231805330431056952,
        9252272755288523725,
        10014268662326746219,
        15565031632950843234,
        1209725273521819323,
        6024642864597845108
            ]
        )
]


def Tip5(input_vec):
    """
    Implementation of the 10-to-5 Tip5 hash function.
    
    Input:
        input_vec: Fp^10.
    Output:
        output: Fp^5.
    """
    X = zero_vector(Fp, 16)
    for i in range(10):
        X[i] = input_vec[i]
    for i in range(10,16):
        X[i] = Fp(1)
    
    Y = f(X, round_constants)
    output = Y[:5]
    return output


K.<x1,x2,x3,x4> = PolynomialRing(Fp, 4, order="degrevlex")

FpZ.<z> = PolynomialRing(Fp)

# The functions used for the collision attack on 4 rounds, c=d=5.

def complementary(indices, n):
    """
    Input:
        indices: set of integer i<n
        n: integer
    Output:
        comp: list of i<n such that i not in indices.
    """
    comp = []
    
    for i in range(n):
        if i in indices:
            continue
        else:
            comp.append(i)
    return comp


def orthogonal_vector(M, indices_a):
    """
    Input:
        M: 16x16 matrix, linear layer.
        indices: indices controlled in the vector a.
    Output:
        a: length 16 vector in Fp.
    """
    M_inv = M.inverse()
    M_inv_sub = M_inv[11:, indices_a]
    
    # basis of the vector space in which the 5 words in the capacity and the 4 wires before S-boxes are fixed.
    basis = matrix(K, M_inv_sub.right_kernel().basis()).transpose()
    
    vec = basis*vector(K, [1,x1,x2,x3,x4])
    a = zero_vector(K, 16)
    
    # shift the coordinates to get back to Fp^16.
    indices_b = complementary(indices_a, 16)
    index = 0
    for j in range(len(indices_a)):
        while index in indices_b:
            index += 1
        a[index] = T(vec[j])
        index += 1

    a = M*a
    return a

def fix_constants(M, round_constant1, round_constant2, indices_b, target):
    """
    Choose the constant of the affine subspace.
    
    Input:
        M:16x16 matrix, linear layer.
        round_constant1: length 16 vector in Fp, constants of the first round.
        indices: the wires we can control.
        target: the fixed values in the capacity of the hash function (in Fp^c).
    Output:
        b: length 16 vector in Fp.
    """
    # Describe the linear system defining b.
    c = len(target)
    M_inv = M.inverse()
    target_after_affine_layer = (M_inv * round_constant1)[16-c:] + target
    M_inv_sub = M_inv[16-c:, indices_b]
    constant = M_inv_sub.solve_right(target_after_affine_layer)
    
    # put the values obtained in the support of b.
    b = zero_vector(Fp, 16)
    index = 0
    indices_a = complementary(indices_b, 16)
    for j in range(len(indices_b)):
        while index in indices_a:
            index += 1
        if index < 4:
            b[index] = S(constant[j])
        else:
            b[index] = T(constant[j])
        index += 1
    
    # Get to the end of round 2.
    b = M*b + round_constant2
    
    return b

def attack_4_rounds(M, round_constants, a, indices_b, target):
    """
    Knowing an attack vector a with support indices, 
    gets an input with capacity target whose output lies in the right hyperplane.
    
    Input:
        M:16x16 matrix, linear layer.
        round_constants: 4 length 16 vectors in Fp.
        a: length 16 vector in Fp.
        indices: list of integers.
        target: the fixed values in the capacity of the hash function (in Fp^c).
    Output:
        x: length 16 vector in Fp.
    """
    # The vector we use to cancel out the last S-box layer.
    orth_vec = M[:5, :4].left_kernel().basis()[0]
    
    # We can fix the first value of the target (in case the polynomial has no root).
    x0 = 0
    while True:
        list_target = [x0] + list(target)
        # Apply S-boxes of first round.
        cur_target = vector(Fp, list_target)
        for i in range(6):
            cur_target[i] = cur_target[i]**7
        
        b = fix_constants(M, round_constants[0], round_constants[1], indices_b, cur_target)
        poly_1 = a*z + b
        poly_2 = zero_vector(FpZ, 16)
        for i in range(4):
            poly_2[i] = S(poly_1[i])
        for i in range(4, 16):
            poly_2[i] = poly_1[i]**7
        poly_2 = M*poly_2 + round_constants[2]
        # State after round 3.
        
        poly_3 = zero_vector(FpZ, 16)
        for i in range(4, 16):
            poly_3[i] = poly_2[i]**7
        poly_3 = M*poly_3 + round_constants[3]
        # State after round 4 (not taking into account the first 4 wires).
        
        # Solve the equation to lie in the hyperplane.
        equation = poly_3[:5].dot_product(orth_vec)
        try: 
            z_sol = equation.roots()[0][0]
            x_round_1 = poly_1(z_sol)
            x = f_inv(x_round_1, round_constants[:2])
            return x
        except:   
            x0 += 1



# We choose to fix the values of the wires 0,1,2,3 (before S-boxes) and 5,6 (before power maps).
indices_b = [0,1,2,3,5,6]
indices_a = complementary(indices_b, 16)

# a is the vector before the S-box layer of the third round.
a = orthogonal_vector(M, indices_a)

# generate the system of equation to solve in order to fix the wires of the S-boxes of the third round.
I = ideal(K, list(a[:4]))
print("The ideal we need to solve to find the vector a is:")
print(I)

print("Using MAGMA, we found the solution:")
print("x1, x2, x3, x4 = {}, {}, {}, {}".format(2857675761863926834, 12096320840878172146, 16939308632471610364, 506223986432544699))

# using MAGMA, we find the vector:
a = a(2857675761863926834, 12096320840878172146, 16939308632471610364, 506223986432544699)
print("State of a before round 3:")
print(a)

vec = M.inverse()*a

for i in range(4):
    vec[i] = S_inv(vec[i])
for i in range(4,16):
    vec[i] = T_inv(vec[i])

print("State of a before round 2:")
print(vec)
print("State of a in the input:")
print(M.inverse()*vec)

# We define the hyperplane in which th eoutputs will lie.

orth_vec = M[:5, :4].left_kernel().basis()[0]

print("We can generate outputs othogonal to the vector v:")
print(orth_vec)

# For all value of the capacity, we can fill the rate in order to have the output in some hyperplane.
target = [0,0,0,0,0]

x = attack_4_rounds(M, round_constants, a, indices_b, target)

print("The input x is:")
print(x)

y = f(x, round_constants[:4])

print("The output y is:")
print(y)

res = orth_vec.dot_product(y[:5])

print("The dot product y[:5]*v equals:")
print(res)

