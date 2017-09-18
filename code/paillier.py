"""
 File       pailler.py
 Author     Jan Henrik Ziegeldorf (ziegeldorf (at) comsys.rwth-aachen.de)
 Brief      Paillier crypto system
 
 Copyright  BLOOM: Bloom filter based outsourced oblivious matchings
            Copyright (C) 2017 Communication and Distributed Systems (COMSYS), RWTH Aachen
            
            This program is free software: you can redistribute it and/or modify
            it under the terms of the GNU Affero General Public License as published
            by the Free Software Foundation, either version 3 of the License, or
            (at your option) any later version.
            
            This program is distributed in the hope that it will be useful,
            but WITHOUT ANY WARRANTY; without even the implied warranty of
            MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
            GNU Affero General Public License for more details.
            You should have received a copy of the GNU Affero General Public License
            along with this program. If not, see <http://www.gnu.org/licenses/>.        
"""

from gmpy import mpz, gcd, invert
from random import randint, getrandbits
from math import log
import gensafeprime
import packings

def genSafePrimes(k):
    p = gensafeprime.generate(k-1)
    return (p-1)/2, p

def randomPrime(k):
    p = mpz(1)
    lower = 2**k
    upper = 2**(k+1)-1
    while not p.is_prime():
        p = mpz(randint(lower,upper))
    return p

def randomFromCyclicGroup(n):
    r = n
    while gcd(n,r) != 1:
        r = randint(1,n-1)
    return mpz(r)

def randomIntByBitlength(k):
    return mpz(getrandbits(k) + 2**k - 1)

def randomBits(n):
    return [randint(0,1) for _ in range(n)]

def crt_pow(x, e, p, q, q_inv, phi_p, phi_q):
    m1 = pow(x, e % phi_p, p)
    m2 = pow(x, e % phi_q, q)
    h = (q_inv*(m1-m2)) % p # q_inv = q^-1 mod p """
    return m2 + h*q  

class Paillier(object):     
    """
        Implementation of the Paillier crypto system.
        
        It can be used in a procedural fashion by calling to the static methods, e.g.
            * pubkey, privkey = Paillier.generateKeys(1024)
            * a = 111
            * A = Paillier.encrypt(a, pubkey, r=1234)
            * print a == Paillier.decrypt(A, privkey)
        
        Or, it can be instantiated as an object with the same funcitons but fixed to 
        a specific set of parameters and public and private key, e.g.
            * p = Paillier(20)
            * a = 111
            * A = p.Enc(a)
            * print a == p.Dec(A)
            
        Notation in source code:
            * A,Blabla,Cipher,D ... (capital letters) variables living in the encrypted domain (i.e., Z/n^2Z)
            * a,blabla,clear,d ...  (small letters) variables living in the plaintext domain (i.e., Z/nZ)
            * encrypt, decrypt ...  static methods
            * Enc, Dec ...          class methods that require instantiation
                        
        *verbose*           switch on verbosity
        *key_length*        bit length for the key to generate
        *pubkey*            instatntiate from public key (instance cannot decrypt then)        
    """
    def __init__(self, key_length=1024, pubkey=None, verbose=False):
        if pubkey:
            self.pubkey = pubkey
            self.privkey = None
        else:      
            self.pubkey, self.privkey = Paillier.generateKeys(key_length)      
        self.n = self.pubkey['n']
        self.pubkey['nsq'] = self.pubkey['n']*self.pubkey['n']
        self.nsq = self.pubkey['nsq']              

    def priv_to_pub(self):
        """
            Get a Paillier instance the holds only the public key.
        """
        return Paillier(pubkey=self.pubkey)

    @staticmethod
    def generateKeys(bit_length):
        """ 
            Generate a Paillier keypair of desired bit length.
            
            *bitlength*    the desired bitlength for the generated key
        """
        # Generate an RSA modulus
        if bit_length <= 2:
            bit_length=3
        p = randomPrime(bit_length/2)
        q = p
        n = p*q
        while p==q or gcd(n, (p-1)*(q-1)) != 1:
            q = randomPrime(bit_length/2)
            n = p*q
        
        # lm = lambda
        lm = (p-1)*(q-1)  
        # Generator g of Z*/n^2Z.
        g = n+1
        # Precompute mu for decryption
        mu = invert(lm, n)
    
        # pubkey, privkey
        psq = p**2
        qsq = q**2
        qsq_inv = invert(qsq, psq)
        return {'n':n, 'g': g, 'nsq': n*n}, {'n': n, 'nsq': n*n, 'g': g, 'lm': lm, 'mu':mu, 'p':p, 'q':q, 'psq':psq, 'phi_psq':psq-p, 'qsq':qsq, 'phi_qsq':qsq-q, 'qsq_inv':qsq_inv}
    
    def Invert(self, s):
        return invert(s, self.n)

    def Enc(self, m, r=None, rbyn=None):
        """
            Encrypt a clear text integer *m*. 
            *m* may be negative.
            
            Limitations: log(m,2) <= log(n,2)-1.
            
            *r*     randomness used for the encryption. 
                    If none, a fresh one will be generated.
            *rbyn*  randomness used for the encryption.
                    This allows for precomputation of the modular exponentiation.
                    If none, a fresh one will be generated.
        """
        m = mpz(m)        
        m %= self.n                
        return Paillier.encrypt(m, self.pubkey, r, rbyn)    

    def Dec(self, C, CRT=True):
        """
            Decrypt a Paillier encrypted ciphertext *C*.
            
            *CRT*     True, if decryption should be sped up using
                      the Chinese remainder theoreme.
        """
        if not self.privkey:
            print("ERROR: Cannot decrypt without private key.")
            return -1
        
        m = Paillier.decrypt(C, self.privkey, CRT=CRT)
        if m > self.n/2:
            m -= self.n
        return m
    
    def Add(self, A, B):
        return Paillier.add(A,B,self.nsq)
    
    def AddScalar(self,A,b,rbyn=None):
        return self.Add(A,self.Enc(b,rbyn=rbyn))
    
    def Mult(self, A, s):
        return Paillier.mult(A, s, self.n, self.nsq)
    
    def XOR(self, A, b):
        return self.AddScalar(self.Mult(A,1-2*b), b) 
    
    def Pack(self, Vals, l):
        """ 
            The ciphertexts from *Vals* are packed as *l* bit values.
            Returns a list of ciphertexts that pack *Vals*.
        """         
        return packings.pack(Vals, self, l)
    
    def Dec_unpack(self, Packings, l, max_cnt=None):
        """ 
            Decrypts each packing in *Packings* and *l* bit plain texts and returns the result as a list.
            If *max_cnt* is specified, a maximum number of *max_cnt* plain texts are extracted.
            (*max_cnt* is required if the last ciphertext is not completely filled up.)    
        """        
        tmp = [packings.unpack(P, self, l) for P in Packings]
        [t.reverse() for t in tmp]
        res = []
        bl = log(long(self.n),2)
        k = int(bl/l)
        if max_cnt:
            cnt = max_cnt
        else:
            cnt = len(Packings)*k
        for r in tmp:
            if cnt >= k:
                res += r
                cnt -= k
            else:
                res += r[-1*cnt:]
                break
        return res
    
    @staticmethod
    def encrypt(m, pubkey, r=None, rbyn=None):
        """
            Enryption uses the following optimiztaions:
                * choosing g = n+1 reduces encryption g^mr^n mod n^2 
                  to (mn + 1)r^n mod n^2 which requires only one modulare exponation
                * *rbyn*: provide a precomputed r^n
                * *r*: provide randomness for encryption 
        """
        n = pubkey['n']
        nsq = pubkey['nsq']
        if rbyn is None:
            if r is None:
                r = randomFromCyclicGroup(n)
            rbyn = pow(r,n,nsq)
        return (m*n+1)*rbyn % nsq   
    
    @staticmethod
    def decrypt(C, privkey, CRT=True):
        """
            Decryption of a Paillier encrypted ciphertext.
            
            *privkey*    the private key used for decryption
            *CRT*        True, if decryption should be sped up using
                         chinese remainder theoreme
        """
        def L(u, n):
            return (u-1)/n
        
        if CRT:
            return (L(crt_pow(C, privkey['lm'], privkey['psq'], privkey['qsq'], privkey['qsq_inv'], privkey['phi_psq'], privkey['phi_qsq']), privkey['n'])*privkey['mu']) % privkey['n']
        else: 
            return (L(pow(C, privkey['lm'], privkey['nsq']), privkey['n'])*privkey['mu']) % privkey['n']

    @staticmethod
    def add(A, B, nsq):
        return (A * B) % nsq
    
    @staticmethod
    def mult(A,s,n,nsq):
        s %= n
        return pow(A,s,nsq)          

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("a", type=int, help="first operand")    
    parser.add_argument("b", type=int, help="second operand")
    args = parser.parse_args()
    
    # generate a Paillier instance
    paillier = Paillier()
    
    # Encrypt a and b
    aEnc = paillier.Enc(args.a)
    bEnc = paillier.Enc(args.b)
    
    # Encrypted addition and scalar multiplication
    sumEnc = paillier.Add(aEnc, bEnc)
    multEnc = paillier.Mult(aEnc, args.b)
    print("%i + %i = %i" % (args.a, args.b, paillier.Dec(sumEnc)))
    print("%i * %i = %i" % (args.a, args.b, paillier.Dec(multEnc)))
        
    # Demonstrate scalar multiplication on packed ciphertexts
    inputs = [randint(0,1) for _ in range(100)]
    Inputs = map(paillier.Enc, inputs)
    packed_Inputs = paillier.Pack(Inputs, l=32)
    selectors = [1 for _ in range(len(packed_Inputs))]    
    packed_Inputs = [paillier.Mult(V, s) for (V, s) in zip(packed_Inputs, selectors)]
    unpacked_inputs = paillier.Dec_unpack(packed_Inputs, l=32, max_cnt=len(inputs))    
    print(unpacked_inputs == inputs)