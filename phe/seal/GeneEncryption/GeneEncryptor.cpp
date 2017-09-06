#include <iostream>
#include <cmath>

#include "seal.h"

using namespace std;
using namespace seal;

Ciphertext encrypt_gene(int gene, IntegerEncoder encoder, Encryptor encryptor);
int decrypt_gene(Ciphertext encryptedGene, IntegerEncoder encoder, Decryptor decryptor);
Ciphertext bit_max(Ciphertext x, Ciphertext y, Ciphertext one, Evaluator evaluator);

int main()
{
  // Create encryption parameters.
  EncryptionParameters parms;
  parms.set_poly_modulus("1x^2048 + 1");
  parms.set_coeff_modulus(ChooserEvaluator::default_parameter_options().at(2048));

  parms.set_coeff_modulus(ChooserEvaluator::default_parameter_options().at(2048));
  parms.set_plain_modulus(1 << 8);
  parms.validate();
  IntegerEncoder encoder(parms.plain_modulus());

  // Generate keys.
  cout << "Generating keys ..." << endl;
  KeyGenerator generator(parms);
  generator.generate();
  cout << "... key generation complete" << endl;
  Ciphertext public_key = generator.public_key();
  Plaintext secret_key = generator.secret_key();
  Encryptor encryptor(parms, public_key);
  Decryptor decryptor(parms, secret_key);
  Evaluator evaluator(parms);

  Ciphertext zero = encrypt_gene(0, encoder, encryptor);
  Ciphertext one = encrypt_gene(1, encoder, encryptor);

  Ciphertext m = bit_max(one, zero, one, evaluator);
  Ciphertext n = bit_max(zero, one, one, evaluator);
  Ciphertext l = bit_max(one, one, one, evaluator);
  Ciphertext p = bit_max(zero, zero, one, evaluator);

  int a = decrypt_gene(m, encoder, decryptor);
  int b = decrypt_gene(n, encoder, decryptor);
  int c = decrypt_gene(l, encoder, decryptor);
  int d = decrypt_gene(p, encoder, decryptor);

  return 0;
}

Ciphertext bit_max(Ciphertext x, Ciphertext y, Ciphertext one, Evaluator evaluator)
{
  Ciphertext a = evaluator.sub(one, x);
  Ciphertext b = evaluator.multiply(a, y);
  Ciphertext c = evaluator.add(x, b);
  return c;
}

Ciphertext encrypt_gene(int gene, IntegerEncoder encoder, Encryptor encryptor)
{

  // Encode two integers as polynomials.
  Plaintext encodedGene = encoder.encode(gene);
  cout << "Encoded " << gene << " as polynomial " << encodedGene.to_string() << endl;

  // Encrypt values.
  cout << "Encrypting values..." << endl;
  Ciphertext encryptedGene = encryptor.encrypt(encodedGene);

  //return encrypted gene
  return encryptedGene;
}

int decrypt_gene(Ciphertext encryptedGene, IntegerEncoder encoder, Decryptor decryptor)
{
  // Decrypt results.
  cout << "Decrypting results ..." << endl;
  Plaintext decryptedGene = decryptor.decrypt(encryptedGene);

  // Decode results.
  int decodedGene = encoder.decode_int32(decryptedGene);

  // Display results.
  cout << "after encryption/decryption = " << decodedGene << endl;

  // How much noise budget did we use in these operations?
  cout << "Noise budget in encryption of " << decodedGene << ": "
      << decryptor.invariant_noise_budget(encryptedGene) << " bits" << endl;

  return decodedGene;
}
