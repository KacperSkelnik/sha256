class SHA256():
    """
    https://perso.crans.org/besson/publis/notebooks/Benchmark_of_the_SHA256_hash_function__Python_Cython_Numba.html
    https://en.wikipedia.org/wiki/SHA-2
    https://www.movable-type.co.uk/scripts/sha256.html
    https://stackoverflow.com/questions/7321694/sha-256-implementation-in-python
    """
    def __init__(self):
        self.k = [
            0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
            0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
            0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
            0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
            0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
            0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
            0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
            0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
            0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
            0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
            0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
            0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
            0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
            0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
            0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
            0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
        ]

        self.h = [
            0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
            0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
        ]

    def leftrotate(self, x, c):
        # left rotate the number x by c bytes
        x &= 0xFFFFFFFF
        return ((x << c) | (x >> (32 - c))) & 0xFFFFFFFF

    def rightrotate(self, x, c):
        # right rotate the number x by c bytes
        x &= 0xFFFFFFFF
        return ((x >> c) | (x << (32 - c))) & 0xFFFFFFFF

    def leftshift(self, x, c):
        # left shift the number x by c bytes
        return x << c

    def rightshift(self, x, c):
        # right shift the number x by c bytes
        return x >> c

    def process_message(self, message):
        data = bytearray(message)
        dat_len_bits = (8 * len(data)) & 0xFFFFFFFFFFFFFFFF
        # Add a single '1' bit at the end
        data.append(0x80)
        # Padding with zeros as long as the input bits length ≡ 448 (mod 512)
        while len(data) % 64 != 56:
            data.append(0)
        # Append original length in bits to message
        data += dat_len_bits.to_bytes(8, byteorder='big')

        assert len(data) % 64 == 0, "Error in padding"
        return data

    def update_hash(self, data):
        # https://en.wikipedia.org/wiki/SHA-2#Pseudocode
        h0, h1, h2, h3, h4, h5, h6, h7 = self.h
        for offset in range(0, len(data), 64):
            chunks = data[offset : offset + 64]
            w = [0]*64
            
            for i in range(16):
                w[i] = int.from_bytes(chunks[4*i : 4*i + 4], byteorder='big')
            for i in range(16, 64):
                # ^ - python XOR
                s0 = (self.rightrotate(w[i-15], 7) ^ self.rightrotate(w[i-15], 18) ^ self.rightshift(w[i-15], 3)) & 0xFFFFFFFF
                s1 = (self.rightrotate(w[i-2], 17) ^ self.rightrotate(w[i-2], 19) ^ self.rightshift(w[i-2], 10)) & 0xFFFFFFFF
                w[i] = (w[i-16] + s0 + w[i-7] + s1) & 0xFFFFFFFF
            a, b, c, d, e, f, g, h = h0, h1, h2, h3, h4, h5, h6, h7

            for i in range(64):
                S1 = (self.rightrotate(e, 6) ^ self.rightrotate(e, 11) ^ self.rightrotate(e, 25)) & 0xFFFFFFFF
                # ~ - python not
                ch = ((e & f) ^ ((~e) & g)) & 0xFFFFFFFF
                temp1 = (h + S1 + ch + self.k[i] + w[i]) & 0xFFFFFFFF
                S0 = (self.rightrotate(a, 2) ^ self.rightrotate(a, 13) ^ self.rightrotate(a, 22)) & 0xFFFFFFFF
                maj = ((a & b) ^ (a & c) ^ (b & c)) & 0xFFFFFFFF
                temp2 = (S0 + maj) & 0xFFFFFFFF

                new_a = (temp1 + temp2) & 0xFFFFFFFF
                new_e = (d + temp1) & 0xFFFFFFFF
                # Rotate the 8 variables
                a, b, c, d, e, f, g, h = new_a, a, b, c, new_e, e, f, g

            # update h's values
            h0 = (h0 + a) & 0xFFFFFFFF
            h1 = (h1 + b) & 0xFFFFFFFF
            h2 = (h2 + c) & 0xFFFFFFFF
            h3 = (h3 + d) & 0xFFFFFFFF
            h4 = (h4 + e) & 0xFFFFFFFF
            h5 = (h5 + f) & 0xFFFFFFFF
            h6 = (h6 + g) & 0xFFFFFFFF
            h7 = (h7 + h) & 0xFFFFFFFF
            self.h = [h0, h1, h2, h3, h4, h5, h6, h7]

    def digest(self):
        return sum(self.leftshift(x, 32 * i) for i, x in enumerate(self.h[::-1]))

    def hexdigest(self):
        digest = self.digest()
        raw = digest.to_bytes(32, byteorder='big')
        format_str = '{:0' + str(64) + 'x}'
        return format_str.format(int.from_bytes(raw, byteorder='big'))

    def hash_message(self, message):
        msg = self.process_message(bytes(message, encoding='utf8'))
        self.update_hash(msg)
        return self.hexdigest()