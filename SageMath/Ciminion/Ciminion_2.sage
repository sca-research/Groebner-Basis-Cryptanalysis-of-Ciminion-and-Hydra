class Ciminion_2:
    """
    Implementation of Ciminion_2.
    In the middle of the key stream it adds key additions to protect against an attack of Bariant ( https://eprint.iacr.org/2023/1283 ).
    """
    def __init__(self,
                 field=GF(2**127 + 45),
                 rounds_C=90,
                 rounds_E=14,
                 constants_C=None,
                 constants_E=None,
                 info_level=0):
        """
        Initialization of Ciminion_2 instance.

        INPUT:
        "field" -- A prime field.
        "rounds_C" -- Rounds of C permutation.
        "rounds_E" -- Rounds of E permutation.
        "constants_C" -- Constants for C permutation, expects a list of r_C elements, 
                         where each element is a list of 4 integers/field elements, 
                         and the last element must be non-zero.
                         If no constants are provide, random constants are generated.
        "constants_E" -- Constants for E permutation, expects a list of r_E elements, 
                         where each element is a list of 4 integers/field elements, 
                         and the last element must be non-zero.
                         If no constants are provide, random constants are generated.
        "info_level" -- Integer, if greater than zero, then Ciminion_2 parameters are
                        printed in console.
        """
        self.field = field
        self.rounds_C = rounds_C
        self.rounds_E = rounds_E
        
        if constants_C is None:
            self.constants_C = []
            for i in range(0, self.rounds_C):
                con = [self.field.random_element(),
                       self.field.random_element(),
                       self.field.random_element()]
                tmp = self.field.random_element()
                while tmp == self.field(0):
                    tmp = self.field.random_element()
                con.append(tmp)
                self.constants_C.append(con)
        else:
            self.constants_C = constants_C
        
        if constants_E is None:
            self.constants_E = []
            for i in range(0, self.rounds_E):
                con = [self.field.random_element(),
                       self.field.random_element(),
                       self.field.random_element()]
                tmp = self.field.random_element()
                while tmp == self.field(0):
                    tmp = self.field.random_element()
                con.append(tmp)
                self.constants_E.append(con)
        else:
            self.constants_E = constants_E
        
        if info_level > 0:
            print("Ciminion_2 parameters")
            print("Field:", self.field)
            print("Rounds_C:", self.rounds_C)
            print("Rounds_E:", self.rounds_E)
            print("Constants_C:", self.constants_C)
            print("Constants_E:", self.constants_E)
    
    def generate_matrix(self, constant):
        """
        Generates the Ciminion matrix.

        INPUT:
        "constant" -- Non-zero element in Ciminion field.

        OUTPUT:
        Invertible 3x3 matrix.
        """
        mat = matrix(self.field, 3 * [3 * [0]])
        mat[0, 2] = self.field(1)
        mat[1, 0] = self.field(1)
        mat[1, 1] = constant
        mat[1, 2] = constant
        mat[2, 1] = self.field(1)
        mat[2, 2] = self.field(1)
        return mat

    def round_function(self, v_in, constants, key=None):
        """
        Ciminion round function.

        INPUT:
        "v_in" -- 3 element vector over Ciminion field.
        "constants" -- round constants.

        OUTPUT:
        3 element vector over Ciminion field.
        """
        v_out = vector(v_in.parent().base_ring(), 3 * [0])
        mat = self.generate_matrix(constants[3])

        # Toffoli gate
        v_out[0] = v_in[0]
        v_out[1] = v_in[1]
        v_out[2] = v_in[2] + v_in[0] * v_in[1]
        # Protection agains Bariant's attack
        if key:
            v_out[2] += key[0] + key[1]

        # Matrix multiplication
        v_out = mat * v_out

        # Constant addition
        v_out += vector(self.field, constants[0:3])

        return v_out
    
    def generate_key_stream(self, key, nonce):
        """
        Generates first Ciminion_2 key pair.
        Includes protection against Bariant's attack.

        INPUT:
        "key" -- A 2 element list/vector for integers/field elements.
        "nonce" -- An integer/field element.

        OUTPUT:
        2 element list of field elements.
        """
        key_stream = vector(self.field, 3 * [0])
        key_stream[0] = self.field(nonce)
        key_stream[1] = self.field(key[0])
        key_stream[2] = self.field(key[1])

        key_stream = self.round_function(key_stream, self.constants_C[0])
        for i in range(1, self.rounds_C):
            key_stream = self.round_function(key_stream, self.constants_C[i], key=key)

        for i in range(0, self.rounds_E):
            key_stream = self.round_function(key_stream, self.constants_E[i], key=key)
        
        return list(key_stream)[0:2]
    
    def encrypt(self, plain, key, nonce):
        """
        Encryption with Ciminion_2.

        INPUT:
        "plain" -- A 2 element list/vector of integers/field elements.
        "key" -- A 2 element list/vector of integers/field elements.
        "nonce" -- An integer/field element.

        OUTPUT:
        A 2 element list of field elements.
        """
        key_stream = vector(self.field, self.generate_key_stream(key, nonce))
        return list(vector(self.field, plain) + key_stream)
    
    def decrypt(self, cipher, key, nonce):
        """
        Decryption with Ciminion_2.

        INPUT:
        "cipher" -- A 2 element list/vector of integers/field elements.
        "key" -- A 2 element list/vector of integers/field elements.
        "nonce" -- An integer/field element.

        OUTPUT:
        A 2 element list of field elements.
        """
        key_stream = vector(self.field, self.generate_key_stream(key, nonce))
        return list(vector(self.field, cipher) - key_stream)
