"""Public-key encryption and signature algorithms.
Public-key encryption uses two different keys, one for encryption and
one for decryption.  The encryption key can be made public, and the
decryption key is kept private.  Many public-key algorithms can also
be used to sign messages, and some can *only* be used for signatures.
========================  =============================================
Module                    Description
========================  =============================================
Crypto.PublicKey.DSA      Digital Signature Algorithm (Signature only)
Crypto.PublicKey.DH 	 (Key exchange)
========================  =============================================
"""

__all__ = ['DH', 'DSA']
