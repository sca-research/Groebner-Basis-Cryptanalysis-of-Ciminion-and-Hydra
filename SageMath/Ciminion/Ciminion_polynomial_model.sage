load("Ciminion.sage")
load("../utilities.sage")

def generate_Ciminion_variables(rounds_C, rounds_E):
    """
    Generates variables for a Ciminion instance as strings.

    INPUT:
    "rounds_C" -- Rounds of C permutation.
    "rounds_E" -- Rounds of E permutation.

    OUTPUT:
    List of variables as strings.
    """
    variables = []

    for i in range(0, rounds_C - 1):
        for j in range(0, 3):
            variables.append("x_C_" + str(j + 1) + "__" + str(i + 1))

    for i in range(0, rounds_E):
        for j in range(0, 3):
            variables.append("x_E_" + str(j + 1) + "__" + str(i + 1))
    
    variables.append("y_1")
    variables.append("y_2")
    variables.append("x")

    return variables

def generate_Ciminion_polynomials(ciminion=Ciminion(),
                                  plain=None,
                                  cipher=None,
                                  nonce=None,
                                  termorder="degrevlex",
                                  field_equations=False,
                                  info_level=1):
    """
    Generates a polynomial model for the first key pair of Ciminion.

    INPUT:
    "ciminion" -- A Ciminion instance. If no argument is supplied,
                  then a random instance is generated.
    "plain" -- A 2 element list/vector of integers/field elements.
               If no plaintext is provided a random one is generated.
    "cipher" -- A 2 element list/vector of integers/field elements.
                If no ciphertext is provided a random one is generated.
    "nonce" -- An integer/field element.
               If no nonce is provided a random one is generated.
    "termorder" -- Termorder for the polynomial ring as string.
                   E.g., "lex", "degrevlex".
    "field_equations" -- Boolean value, if set to true field equations
                         are added to the Ciminion polynomials.
                         Default value "False".
    
    OUTPUT:
    Ciminion polynomial system.
    """
    print_key = False
    if nonce is None:
        nonce = ciminion.field.random_element()
    if plain is None:
        plain = [ciminion.field.random_element(),
                 ciminion.field.random_element()]
    if cipher is None:
        print_key = True
        key = [ciminion.field.random_element(),
               ciminion.field.random_element()]
        cipher = ciminion.encrypt(plain, key, nonce)
    
    if info_level > 0:
        print("Plaintext:", plain)
        if print_key:
            print("Key:", key)
        print("Nonce:", nonce)
        print("Ciphertext:", cipher)
        print("Term order:", termorder)

    variables = generate_Ciminion_variables(ciminion.rounds_C, ciminion.rounds_E)
    P = PolynomialRing(ciminion.field, variables, order=termorder)
    variables = list(P.gens())

    polynomials = []
    N = 3 * (ciminion.rounds_C - 1 + ciminion.rounds_E)
    key_variables = variables[N:N + 2]
    auxiliary_variable = variables[N + 2]
    variables_C = variables[0:3 * (ciminion.rounds_C - 1)]
    variables_E = variables[3 * (ciminion.rounds_C - 1):N]

    current_state = vector(P, 3 * [0])
    current_state[0] = nonce
    current_state[1] = key_variables[0]
    current_state[2] = key_variables[1]
    next_state = vector(P, variables_C[0:3])
    polys = ciminion.round_function(current_state, ciminion.constants_C[0]) - next_state
    polynomials = polynomials + list(polys)

    for i in range(1, ciminion.rounds_C - 1):
        current_state = next_state
        next_state = vector(P, variables_C[3 * i:3 * (i + 1)])
        polys = ciminion.round_function(current_state, ciminion.constants_C[i]) - next_state
        polynomials = polynomials + list(polys)
    current_state = next_state
    next_state = vector(P, variables_E[0:3])
    polys = ciminion.round_function(current_state, ciminion.constants_C[ciminion.rounds_C - 1]) - next_state
    polynomials = polynomials + list(polys)
    
    current_state = next_state
    next_state = vector(P, variables_E[3:6])
    polys = ciminion.round_function(current_state, ciminion.constants_E[0]) - next_state
    polynomials = polynomials + list(polys)

    for i in range(1, ciminion.rounds_E - 1):
        current_state = next_state
        next_state = vector(P, variables_E[3 * (i + 1):3 * (i + 2)])
        polys = ciminion.round_function(current_state, ciminion.constants_E[i]) - next_state
        polynomials = polynomials + list(polys)
    
    current_state = next_state
    next_state = vector(P, 3 * [0])
    next_state[0] = -plain[0] + cipher[0]
    next_state[1] = -plain[1] + cipher[1]
    next_state[2] = auxiliary_variable
    polys = ciminion.round_function(current_state, ciminion.constants_E[ciminion.rounds_E - 1]) - next_state
    polynomials = polynomials + list(polys)

    if field_equations:
        polynomials = polynomials + generate_field_equations(variables)

    return polynomials

def efficient_Ciminion_termorder(ciminion, P):
    """
    Change to a more efficient DRL order.

    INPUT:
    "ciminion" -- Ciminion instance.
    "P" -- Polynomial ring of the Ciminion polynomials.

    OUTPUT:
    A homorphism to a polynomial ring with a more efficient
    DRL term order for Ciminion.
    """
    variables = generate_Ciminion_variables(ciminion.rounds_C, ciminion.rounds_E)
    N = 3 * (ciminion.rounds_C - 1 + ciminion.rounds_E)
    variables_Q = variables[N:N + 2] + variables[0:N] + [variables[N + 2]]
    Q = PolynomialRing(ciminion.field, variables_Q, order="degrevlex")
    variables_Q = list(Q.gens())
    variables_Q = variables_Q[2:2 + N] + variables_Q[:2] + [variables_Q[N + 2]]
    return P.hom(variables_Q, Q)

def transform_Ciminion_polynomial_system(ciminion, ciminion_polys):
    """
    Inverts the matrix for every round in the Ciminion polynomial system.

    INPUT:
    "ciminion" -- Ciminion instance.
    "ciminion_polys" -- Ciminion polynomials, expects them in the same order
    as in the output of "generate_Ciminion_polynomials".

    OUTPUT:
    Ciminion polynomials with inverted matrix in every round.
    """
    transformed_polys = []
    for i in range(0, ciminion.rounds_C):
        mat = ciminion.generate_matrix(ciminion.constants_C[i][3]).inverse()
        transformed_polys = transformed_polys + (mat * vector(ciminion_polys[3 * i:3 * (i + 1)])).list()
    for i in range(0, ciminion.rounds_E):
        mat = ciminion.generate_matrix(ciminion.constants_E[i][3]).inverse()
        transformed_polys = transformed_polys + (mat * vector(ciminion_polys[3 * (ciminion.rounds_C +  i):3 * (ciminion.rounds_C + i + 1)])).list()
    return transformed_polys
