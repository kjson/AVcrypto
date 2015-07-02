    def schoof(self):
        """ Rene Schoofs point counting algorithm """
        a,b,p = self.a,self.b,self.p
        x, y = sp.symbols('x,y')

        list_of_primes=list_of_congruences=[]

        # Create list of primes
        i=product=1
        while product < 4*math.sqrt(p):
            if is_probable_prime(i):
                list_of_primes.append(i)
                product*=i
            i+=2  

        # Create list of congruences 
        if sp.gcd(x ** p - x, x ** 3 + a * x + b) != 1:
            list_of_congruences.append((2, 0))
        else:
            list_of_congruences.append((2, 1))

        for l in listOfPrimes:
            div_poly_l = E.dpoly(l).subs(y ** 2, x ** 3 - a * x + b)
            pl = p % l
            x_pl = (x - E.dpoly(pl - 1) * E.dpoly(pl + 1) / (E.dpoly(pl)) ** \
                2).subs(y ** 2, x ** 3 - a * x + b)
            y_pl = (E.dpoly(2 * pl) / E.dpoly(pl) ** 4).subs(y ** 2, x ** 3 \
                - a * x + b)
            theta_y = y_pl / y
            # Problem polynomial has massive degree
            slope = ((x ** 3 + a + x + b) * ((x ** 3 + a * x + b) ** (((p ** 2)\
                - 1) / 2) - theta_y ) ** 2) / (x ** (2 * p) - x_pl)
            x_prime = slope ** 2 - x ** (2 * p) - x_pl
            y_prime = - (x ** 3 + a + x + b)**p + slope*(x_prime - x ** (2 * p))

            for t in xrange(1, (l - 1) / 2 + 2):
                rhs = E.symbolic_scalar(t, (x, y))[0]
                if sp.div(x_prime - rhs, div_poly_l)[1] == 0:
                    if sp.div((y_prime - E.dpoly(2 * t) / 2 * E.dpoly(t) ** 4)\
                        / y, div_poly_l)[1] == 0:
                        list_of_congruences.append((l, t))
                    else:
                        list_of_congruences.append((l, -t))

        return list_of_congruences

    def symbolic_scalar(self, scalar, (x, y)):
        """ Symbolic scalar multiplication  of (x,y) """
        if scalar == 0:
            return (sp.symbols('O'), sp.symbols('O'))
        elif scalar == 1:
            return (x, y)
        elif scalar % 2 == 0:
            return self.symbolic_scalar(scalar / 2, self.symbolic_add((x, y), (x, y)))
        else:
            return self.symbolic_add((x, y), self.symbolic_scalar(scalar - 1, (x, y)))


    def symbolic_add(self, (x1, y1), (x2, y2)):
        """ Adds (x1,y1) + (x2,y2) symbolically """
        a = self.a
        if x1 != x2:
            s = (y1 - y2) / (x1 - x2)
        else:
            if y1 == -y2:
                return (sp.symbols('O'), sp.symbols('O'))
            else:
                s = (3 * x1 ** 2 + a) / (2 * y1)
        x3 = s ** 2 - x1 - x2
        y3 = - y1 + s * (x3 - x1)
        return (x3, y3)

    def dpoly(self, n):
        """ Calculate the nth division polynomial for E """
        x, y = sp.symbols('x,y')
        a,b = self.a,self.b 
        assert n >= 0

        if n == 0 or n == 1:
            return sp.poly(n, gens=x)
        elif n == 2:
            return sp.poly(2 * y)
        elif n == 3:
            return sp.poly(3 * x ** 4 + 6 * a * x ** 2 + 12 * b * x - a ** 2)
        elif n == 4:
            return sp.poly(4 * y * (x ** 6 + 5 * a * x ** 4 \
                                   + 20 * b * x ** 3 - 5 * a ** 2 * x ** 2 - \
                                   4 * a * b * x - 8 * b ** 2 - a ** 3))
        elif n % 2 == 1:
            m = n / 2
            return self.dpoly(m + 2) * self.dpoly(m) ** 3 \
                   - self.dpoly(m - 1) * self.dpoly(m + 1) ** 3
        elif n % 2 == 0:
            m = n / 2
            return (self.dpoly(m) / 2 * y) \
                   * (self.dpoly(m + 2) * self.dpoly(m - 1) ** 2 \
                      - self.dpoly(m - 2) * self.dpoly(m + 1) ** 2)


    def order(self):
        """ Order of the group of points """
        p = self.p
        if p.bit_length() < 20:
            return self.lenstra()
        elif 20 <= p.bit_length() <= 100:
            return self.schoof()
        else:
            return self.schoof_elkies_atkin()


    def lenstra(self):
        """ Lenstra's simple point counting method """
        a,b,p = self.a,self.b,self.p
        number_of_points = 1 + p
        for x in xrange(1,p):
            number_of_points += legendre((x**3 + a*x + b) % p,p)
        return number_of_points



