using Oscar

struct Ciminion
    field::FqField
    rounds_C::Int64
    rounds_E::Int64
    constants_C::FqMatrix
    constants_E::FqMatrix
end

function Ciminion_constructor(;field=GF(2),
                               rounds_C=3,
                               rounds_E=3,
                               constants_C=nothing,
                               constants_E=nothing,
                               info_level=0)
    """
    Initialization of Ciminion instance.

    INPUT:
    "field" -- A prime field.
    "rounds_C" -- Rounds of C permutation.
    "rounds_E" -- Rounds of E permutation.
    "constants_C" -- Constants for C permutation, expects a matrix over the
                     field of dimension rounds_C x 4. 
                     If no constants are provide, random constants are generated.
    "constants_E" -- Constants for E permutation, expects a matrix over the
                     field of dimension rounds_E x 4. 
                     and the last element must be non-zero.
                     If no constants are provide, random constants are generated.
    "info_level" -- Integer, if greater than zero, then Ciminion parameters are
                    printed in console.
    """
    if isnothing(constants_C)
        constants_C = zero_matrix(field, 4, rounds_C)
        for i in 1:rounds_C
            for j in 1:4
                constants_C[j, i] = rand(field)
            end
            while constants_C[4, i] == zero(field)
                constants_C[4, i] = rand(field)
            end
        end
    end
    if isnothing(constants_E)
        constants_E = zero_matrix(field, 4, rounds_E)
        for i in 1:rounds_E
            for j in 1:4
                constants_E[j, i] = rand(field)
            end
            while constants_E[4, i] == zero(field)
                constants_E[4, i] = rand(field)
            end
        end
    end
    if info_level > 0
        println("Ciminion parameters")
        println("Field: ", field)
        println("Rounds_C: ", rounds_C)
        println("Rounds_E: ", rounds_E)
        println("Constants_C: ", constants_C)
        println("Constants_E: ", constants_E)
    end
    return Ciminion(field, rounds_C, rounds_E, constants_C, constants_E)
end

function generate_matrix(constant::FqFieldElem)
    """
    Generates the Ciminion matrix.

    INPUT:
    "constant" -- Non-zero element in Ciminion field.

    OUTPUT:
    Invertible 3x3 matrix.
    """
    K = parent(constant)
    mat = zero_matrix(K, 3, 3)
    mat[1, 3] = one(K)
    mat[2, 1] = one(K)
    mat[2, 2] = constant
    mat[2, 3] = constant
    mat[3, 2] = one(K)
    mat[3, 3] = one(K)
    return mat
end

function round_function(v_in::Union{FqMatrix, AbstractAlgebra.Generic.MatSpaceElem{FqMPolyRingElem}},
                        constants::FqMatrix)
    """
    Ciminion round function.

    INPUT:
    "v_in" -- 3 element vector over Ciminion field/polynomial ring.
    "constants" -- round constants.

    OUTPUT:
    3 element vector over Ciminion field/polynomial ring.
    """
    R = base_ring(v_in)
    v_out = zero_matrix(R, 3, 1)
    mat = generate_matrix(constants[4, 1])

    # Toffoli gate
    v_out[1, 1] += v_in[1, 1]
    v_out[2, 1] += v_in[2, 1]
    v_out[3, 1] += v_in[3, 1] + v_in[1, 1] * v_in[2, 1]

    # Matrix multiplcation
    v_out = mat * v_out

    # Add round constants
    v_out += constants[1:3, :]

    return v_out
end

function generate_key_stream(key, nonce, ciminion::Ciminion)
    """
    Generates first Ciminion key pair.

    INPUT:
    "key" -- 2x1 matrix of integers/field elements.
    "nonce" -- An integer/field element.
    "ciminion" -- A Ciminion instance.

    OUTPUT:
    2x1 matrix of field elements.
    """
    key_stream = zero_matrix(base_ring(key), 3, 1)
    key_stream[1, 1] += nonce
    key_stream[2, 1] += key[1, 1]
    key_stream[3, 1] += key[2, 1]

    for i in 1:ciminion.rounds_C
        key_stream = round_function(key_stream, matrix(ciminion.constants_C[:, i]))
    end
    
    for i in 1:ciminion.rounds_E
        key_stream = round_function(key_stream, matrix(ciminion.constants_E[:, i]))
    end

    return key_stream
end

function encrypt(plain, key, nonce, ciminion::Ciminion)
    """
    Encryption with Ciminion.

    INPUT:
    "plain" -- 2x1 matrix of integers/field elements.
    "key" -- 2x1 matrix of integers/field elements.
    "nonce" -- An integer/field element.
    "ciminion" -- A Ciminion instance.

    OUTPUT:
    2x1 matrix of field elements.
    """
    key_stream = generate_key_stream(key, nonce, ciminion)
    
    return plain + key_stream[1:2, :]
end

function decrypt(cipher, key, nonce, ciminion::Ciminion)
    """
    Decryption with Ciminion.

    INPUT:
    "cipher" -- 2x1 matrix of integers/field elements.
    "key" -- 2x1 matrix of integers/field elements.
    "nonce" -- An integer/field element.
    "ciminion" -- A Ciminion instance.

    OUTPUT:
    2x1 matrix of field elements.
    """
    key_stream = generate_key_stream(key, nonce, ciminion)
    
    return cipher - key_stream[1:2, :]
end
