"""
 File       packings.py
 Author     Jan Henrik Ziegeldorf (ziegeldorf (at) comsys.rwth-aachen.de)
 Brief      Packing methods for generic homomorphic encryption.
 
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

import math

def pack(Vals, HE, l):
    """
        Pack *l* bit cipher texts from *Vals*.
                   
        Returns a list of cipher texts, as many as
        necessary to pack all cipher texts from *Vals*
        
        The cipher texts from *Vals* will be packed in the same order,
        i.e. (x0, x1, ..., xn) is packed into
            X1 = x0     || ... || xk
            X2 = xk+1   || ... || x2k
            ...
            Xj = xn-k+1 || ... || xn            
    """    
    shift_factor = pow(2,l)  
    keylen = math.log(long(HE.pubkey['n']), 2)
    k = int(keylen/l)
    
    #print "Packing %i %i-bit ciphertexts into %i ciphertexts, compression factor = %.2f" % (len(Vals), l, math.ceil(len(Vals) / float(k)), float(len(Vals)) / math.ceil(len(Vals) / float(k)))
    
    res = []
    i = 0
    while (i < len(Vals)):
        X = Vals[i]
        for V in Vals[i+1:i+k]:
            X = HE.Mult(X,shift_factor)
            X = HE.Add(X, V)                   
        res.append(X)
        i += k
    return res 

def unpack(P, HE, l):
    """ 
        Decrypts and unpacks one packed ciphertext *P* into *l*-bit clear texts. 
        We always try to extract as many clear texts as possible from the packing,
        i.e. log(n)/l many. The caller needs to throw away those that are potentially rubbish.
    """        
    res = []
    keylen = math.log(long(HE.pubkey['n']),2)
    k = int(keylen/l)
    shift_factor = pow(2, l)
    p = HE.Dec(P)
    for _ in range(k):
        val = p % shift_factor;
        res.append(val)
        p = p / shift_factor
    return res 