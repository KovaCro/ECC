import secrets

# Jednostavna eliptička krivulja. U pravim aplikacijama koriste se standardizirane krivulje koje su sigurne.
p = 11
a = 1
b = 6
G = (2,4) # nisam provjeravo red tocke

def is_on_curve(x, y):
    """ Checks if a point lies on the curve. """
    if x is None or y is None:
        return True
    return (y**2) % p == (x**3 + a*x + b) % p

def point_addition(P, Q):
    """ Adds two points on the curve. """
    if P is None:
        return Q
    if Q is None:
        return P
    if P == Q:
        return point_doubling(P)
    x1, y1 = P
    x2, y2 = Q

    if x1 == x2:
        return None

    lam = ((y2 - y1) * pow(x2 - x1, p - 2, p)) % p
    x3 = (lam**2 - x1 - x2) % p
    y3 = (lam * (x1 - x3) - y1) % p
    return (x3, y3)

def point_doubling(P):
    """ Doubles a point on the curve. """
    if P is None: 
        return None
    x, y = P
    if y == 0: 
        return None

    lam = ((3 * x**2 + a) * pow(2 * y, p - 2, p)) % p
    x3 = (lam**2 - 2 * x) % p
    y3 = (lam * (x - x3) - y) % p
    return (x3, y3)

def scalar_multiplication(k, P):
    """ Multiplies a point P by a scalar k. """
    result = None
    while k > 0:
        if k & 1:
            result = point_addition(result, P)
        P = point_doubling(P)
        k >>= 1
    return result

def generate_key_pair():
    """ Generates key pair. """
    sk = secrets.randbelow(p) # p treba bit red tocke G
    pk = scalar_multiplication(sk, G)
    return sk, pk

def encrypt(m, pk):
    """ Encrypts a message. """
    k = secrets.randbelow(p) # p treba bit red tocke G
    c1 = scalar_multiplication(k, G)
    c2 = point_addition(m, scalar_multiplication(k, pk))
    return c1, c2

def decrypt(c, sk):
    """ Decrypts a ciphertext using the private key. """
    c1, c2 = c
    x_tmp, y_tmp = scalar_multiplication(sk, c1)
    m = point_addition(c2, (x_tmp, -y_tmp))
    return m

private_key, public_key = generate_key_pair()
print("Private Key:", private_key)
print("Public Key:", public_key)

message = (3,5)

ciphertext = encrypt(message, public_key)
print("Ciphertext:", ciphertext)

decrypted_message = decrypt(ciphertext, private_key)
print("Decrypted Message:", decrypted_message)
